/*!
 * \author Daniel Sundfeld
 * \copyright MIT License
 */


#ifndef WIN32
#include <sched.h>
#endif


#include "mpi_dependencies.h"
#include "PAStar.h"
#include "backtrace.h"
#include "Coord.h"
#include "Node.h"
#include "TimeCounter.h"
#include <chrono>

#include <sstream>
// include input and output archivers
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

// include this header to serialize vectors, maps and other types
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/collection_size_type.hpp>

#define MPI_TAG_SEND_COMMON                0
#define MPI_TAG_BARRIER                  202
#define MPI_TAG_GLOBAL_CHECK             203
#define MPI_TAG_FINAL_GLOBAL_CHECK       204
#define MPI_TAG_KILL_RECEIVER            205
#define MPI_TAG_BARRIER                  206
#define MPI_TAG_BROADCAST_NODE_TO_0      207


#ifndef WIN32

#ifdef MAC
    int sched_setaffinity(pthread_t thread, size_t cpu_size, cpu_set_t *cpu_set)
    {
        thread_port_t mach_thread;
        int core = 0;

        for (core = 0; core < 8 * cpu_size; core++) {
            if (CPU_ISSET(core, cpu_set)) break;
        }
        std::cout << "binding to core " << core << " \n";
        thread_affinity_policy_data_t policy = { core };
        mach_thread = pthread_mach_thread_np(thread);
        thread_policy_set(mach_thread, THREAD_AFFINITY_POLICY,(thread_policy_t)&policy, 1);
        return 0;
    }
#endif

#endif

template < int N >
PAStar<N>::PAStar(const Node<N> &node_zero, const struct PAStarOpt &opt)
    : m_options(opt),
      nodes_count ( ),
      nodes_reopen ( )
{
    std::cout << "Running PAStar with: "
              << opt.threads_num << " threads, "
              << Coord<N>::get_hash_name() << " hash, "
              << Coord<N>::get_hash_shift() << " shift.\n";

    end_cond = false;
	end_condLocal = false;
    sync_count = 0;
    final_node.set_max();
	best_received.set_max();

	check_stop_mutex = new std::mutex;

    OpenList = new PriorityList<N>[m_options.threads_num]();
    ClosedList = new std::map< Coord<N>, Node<N> >[m_options.threads_num]();

    nodes_count = new long long int[m_options.threads_num]();
    nodes_reopen = new long long int[m_options.threads_num]();
	

	// The additional mutex and condition are used by sender thread
    queue_mutex = new std::mutex[m_options.threads_num+1];	  
    queue_condition = new std::condition_variable[m_options.threads_num+1];
    queue_nodes = new std::vector< Node<N> >[m_options.threads_num];

    OpenListFinal = NULL;
    ClosedListFinal = NULL;
    nodes_countFinal = NULL;
    nodes_reopenFinal = NULL;

	processing_mutex = new std::shared_mutex[m_options.threads_num]();
	//senderWait = new std::mutex();

    // Allocate final structures and enqueue first node if rank zero
    if (m_options.mpiRank == 0)
    {
        OpenListFinal = new PriorityList<N>[m_options.totalThreads];
        ClosedListFinal = new std::map< Coord<N>, Node<N> >[m_options.totalThreads];
        nodes_countFinal = new long long int[m_options.totalThreads]();
        nodes_reopenFinal = new long long int[m_options.totalThreads]();
		nodes_openListSizeFinal = new long long int[m_options.totalThreads]();
		nodes_closedListSizeFinal = new long long int[m_options.totalThreads]();

        OpenList[0].enqueue(node_zero);
    }
}

template < int N >
PAStar<N>::~PAStar()
{
    delete[] nodes_count;
    delete[] nodes_reopen;
    delete[] queue_mutex;
    delete[] queue_condition;
    delete[] queue_nodes;

    if (m_options.mpiRank == 0)
    {
		delete[] OpenListFinal;
		delete[] ClosedListFinal;
		delete[] nodes_closedListSizeFinal;
		delete[] nodes_openListSizeFinal;
		delete[] nodes_countFinal;
		delete[] nodes_reopenFinal;

    }
}
#ifndef WIN32
template < int N >
int PAStar<N>::set_affinity(int tid)
{
    std::cout << "setting affinity of tid " << tid << " rank " << m_options.mpiRank << std::endl;
    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(1 << tid, &mask);
    return sched_setaffinity(0, sizeof(mask), &mask);
}
#endif
/*!
 * Add a vector of nodes \a nodes to the OpenList with id \a tid. Use the
 * ClosedList information to ignore expanded nodes.
 * This function is a expensive function and should be called with no locks.
 * Parallel access should never occur on OpenList and ClosedList with
 * same tids.
 */


template < int N >
void PAStar<N>::enqueue(int tid, std::vector< Node<N> > &nodes)
{
    //std::cout << "enqueue" << std::endl;

    typename std::map< Coord<N>, Node<N> >::iterator c_search;

    for (typename std::vector< Node<N> >::iterator it = nodes.begin() ; it != nodes.end(); ++it)
    {
        if ((c_search = ClosedList[tid].find(it->pos)) != ClosedList[tid].end())
        {
            if (it->get_g() >= c_search->second.get_g())
                continue;
            ClosedList[tid].erase(it->pos);
            nodes_reopen[tid] += 1;
			//std::cout << "tid" << tid << std::endl;
        }
        //std::cout << "Adding:\t" << *it << std::endl;
        OpenList[tid].conditional_enqueue(*it);
    }
    return;
}

//! Consume the queue with id \a tid
template < int N >
void PAStar<N>::consume_queue(int tid)
{
    //std::cout << "consume queue" << std::endl;
    std::unique_lock<std::mutex> queue_lock(queue_mutex[tid]);
    std::vector< Node<N> > nodes_to_expand(queue_nodes[tid]);
	queue_nodes[tid].clear();
    queue_lock.unlock();

    //std::cout << "consuming queue" << std::endl;
    enqueue(tid, nodes_to_expand);
    return;
}

//! Wait something on the queue
template < int N >
void PAStar<N>::wait_queue(int tid)
{
    //std::cout << "wait queue:";// << std::endl;
    //std::cout << " tid " << tid << " queue size " << queue_nodes[tid].size() << std::endl;
    std::unique_lock<std::mutex> queue_lock(queue_mutex[tid]);
    if (queue_nodes[tid].size() == 0)
        queue_condition[tid].wait(queue_lock);
    return;
}

//! Wake up everyone waiting on the queue
template < int N >
void PAStar<N>::wake_all_queue()
{
	int i;
    for (i = 0; i < m_options.threads_num; ++i)
    {
        std::unique_lock<std::mutex> queue_lock(queue_mutex[i]);
        queue_condition[i].notify_one();
    }
    return;
}

//! Sync all threads, including remote ones
template < int N >
void PAStar<N>::sync_threads()
{
    //std::cout << "sync threads" << std::endl;

    std::unique_lock<std::mutex> sync_lock(sync_mutex);

    if (++sync_count < m_options.threads_num)
        sync_condition.wait(sync_lock);
    else
    {
        //When every thread is blocked, block processes in a barrier
		Node<N> nullNode;
		std::unique_lock<std::mutex> queue_lock(queue_mutex[m_options.threads_num]);
		std::unique_lock<std::mutex> lock(sync_mutex_global);

		//std::cout << m_options.mpiRank << ": send node" << std::endl;

		send_queue.push(std::make_tuple(0, MPI_TAG_BARRIER, nullNode));
		queue_condition[m_options.threads_num].notify_one();
		queue_lock.unlock();
		sync_condition_global.wait(lock);

		//Wait until sender finishes the barrier
		//std::cout << m_options.mpiRank << ": barrier passed" << std::endl;

        //After that, restart everyone
        sync_count = 0;
        sync_condition.notify_all();
    }
}

//! Sync all local threads
template < int N >
void PAStar<N>::sync_threads_local()
{
	//std::cout << "sync threads" << std::endl;

	std::unique_lock<std::mutex> sync_lock(sync_mutex);

	if (++sync_count < m_options.threads_num)
		sync_condition.wait(sync_lock);
	else
	{
		//After that, restart everyone
		sync_count = 0;
		sync_condition.notify_all();
	}
}

//! Execute the pa_star algorithm until all nodes expand the same final node
template < int N >
void PAStar<N>::worker_inner(int tid, const Coord<N> &coord_final)
{
    Node<N> current;
    std::vector< Node<N> > *neigh = new std::vector< Node<N> >[m_options.totalThreads];

    // Loop ended by process_final_node
    while (end_cond == false)
    {
        typename std::map< Coord<N>, Node<N> >::iterator c_search;

        // Start phase

		// Threads wait until every processing message is received
		while (!queue_nodes[tid].empty())
		{

			// Threads wait until every processing message is received
			std::unique_lock<std::shared_mutex> lock(processing_mutex[tid]);
			lock.unlock();
			std::this_thread::sleep_for(std::chrono::milliseconds(20));
			// Consume openlist to find if someone have a better route
			consume_queue(tid);
		}

        // Dequeue phase
        if (OpenList[tid].dequeue(current) == false)
        {
            //std::cout << "dequeueing" << std::endl;
            wait_queue(tid);
            //std::cout << "good to go" << std::endl;
            continue;
        }
        nodes_count[tid] += 1;

        // Check if better node is already found
        if ((c_search = ClosedList[tid].find(current.pos)) != ClosedList[tid].end())
        {
            //std::cout << "searching for a better node" << std::endl;
            if (current.get_g() >= c_search->second.get_g())
                //std::cout << "oops, we didn't found it, so we are good" << std::endl;
                continue;
            nodes_reopen[tid] += 1;
        }

        //std::cout << "[" << tid << "] Opening node:\t" << current << std::endl;
        ClosedList[tid][current.pos] = current;
											  
        if (current.pos == coord_final)
        {
            //std::cout << m_options.mpiRank << "-" << tid << ": achei essa porra" << std::endl;
            process_final_node(tid, current);
            continue;
        }

        // Expand phase
        current.getNeigh(neigh, m_options.totalThreads);

        int actualTid = tid + m_options.mpiMin;
        //std::cout << actualTid << std::endl;

        // Reconciliation phase
        for (int i = 0; i < m_options.totalThreads; i++)
        {
            // process local generated nodes
            if (i == actualTid)
            {
                //std::cout << "current node in action" << std::endl;
                enqueue(tid, neigh[i]);

            }
            //enqueue for other nodes
            else if (neigh[i].size() != 0)
            {
                //std::cout << "sending to different node (i:"<< i << "; actualTid:" << actualTid;
                //std::cout << "; localMin: " << m_options.mpiMin << "; localMax:"<< m_options.mpiMax << ")" << std::endl;
                //wich can be local
                
                if (i >= m_options.mpiMin && i < m_options.mpiMax)
                {
					int iLocal = i % m_options.threads_num;
                    //std::cout << "happily, a local node" << std::endl;
                    std::unique_lock<std::mutex> queue_lock(queue_mutex[iLocal]);
                    queue_nodes[iLocal].insert(queue_nodes[iLocal].end(), neigh[i].begin(), neigh[i].end());
                    queue_condition[iLocal].notify_one();
                }
                //or remote
                else
				{
					std::unique_lock<std::mutex> queue_lock(queue_mutex[m_options.threads_num]);
					for (unsigned int j = 0; j < neigh[i].size(); j++)
						send_queue.push(std::make_tuple(i, MPI_TAG_SEND_COMMON, neigh[i].at(j)));
					queue_condition[m_options.threads_num].notify_one();
                }
            }
            neigh[i].clear();
        }
    }
    delete[] neigh;
    return;
}

/*!
 * Process \a n as an possible answer. Check end phase 1.
 * When a final node is first opened, it is broadcasted in all OpenLists.
 * When all OpenList open this node, it have the lowest priority among all
 * openlists, then it must proceed to Check end phase 2.
 * This is functions does not require synchronization between the threads.
 */
template < int N >
void PAStar<N>::process_final_node(int tid, const Node<N> &n)
{
	int i = 0;
	unsigned int originRank = ((n.pos.get_id(m_options.totalThreads) - tid) / m_options.threads_num);
	//std::cout << "process final node" << std::endl;

	//Lock all threads from all nodes
	std::unique_lock<std::mutex> final_node_lock(final_node_mutex);

	// Better possible answer already found, discard n
	if (final_node.get_f() < n.get_f())
	{
		//std::cout << "newbie, we already got a better candidate" << std::endl;
		//return;
	}
	else
	{
		if (n.pos.get_id(m_options.threads_num) == ((unsigned int)tid))
		{
			//std::cout << "[proc: " << m_options.mpiRank << " - tid: " << tid << "] Possible answer found: " << n << std::endl;
			// Broadcast the node
			final_node = n;
			final_node_count = 0;

			//If the tid is shared with founder of the final node, add to other local nodes
			for (i = 0; i < m_options.threads_num; i++)
			{
				if (i != tid)
				{
					std::unique_lock<std::mutex> queue_lock(queue_mutex[i]);
					queue_nodes[i].push_back(n);
					queue_condition[i].notify_one();
				}
			}
	
				//For remote threads: send messages to be broadcasted between node threads
				std::unique_lock<std::mutex> queue_lock(queue_mutex[m_options.threads_num]);
				send_queue.push(std::make_tuple(0, MPI_TAG_BROADCAST_NODE_TO_0 + tid, n));
				queue_condition[m_options.threads_num].notify_one();
			

			final_node_lock.unlock();
		}
		// Every other worker is unlocked
		else
		{
			//std::cout << "[" << m_options.mpiRank << ":" << tid << "] Agreed with possible answer! " << n << "/" << final_node << std::endl;
			//if (n != final_node) std::cout << "BUG HERE!\n";
			final_node_lock.unlock();
		}

		// This node have the highest priority between all Openlist. Broadcast the end condition
		if (++final_node_count == m_options.threads_num)
		{
			end_cond = true;
		}
	}
	sync_threads_local();
    return;
}

/*!
 * Check end phase 2.
 * After every thread agreed that a possible answer is found, we must syncronize
 * the threads, consume the queue and check again, if the answer have the
 * lowest priority between all OpenLists
 * The queue consume and/or thread scheduling might have caused the final_node
 * to not have the lowest priority.
 * This is a very costly function, threads syncronization are called twice.
 */
template < int N >
bool PAStar<N>::check_stop(int tid)
{
	wake_all_queue();

	Node<N> n = final_node;

	while (!queue_nodes[tid].empty())
	{

		// Threads wait until every processing message is received
		std::unique_lock<std::shared_mutex> lock(processing_mutex[tid]);
		lock.unlock();
		//std::this_thread::sleep_for(std::chrono::milliseconds(20));
		// Consume openlist to find if someone have a better route
		consume_queue(tid);
	}

	//std::cout << m_options.mpiRank << "-" << tid << ": checking end with f-value = " << final_node.get_f() << std::endl;
	// If someone have a better value, we're not at the end yet	   -> we could look for minimum value inside node and then use MPI_Allreduce to find global minimum
	if (OpenList[tid].get_highest_priority() < final_node.get_f())
	{
		std::cout << "looks like i've found a better route" << std::endl;
		end_cond = false;
	}

	sync_threads(); 

	// Get if we're at global end
	if (tid == 0)
	{
		Node<N> nullNode;
		std::unique_lock<std::mutex> queue_lock(queue_mutex[m_options.threads_num]);
		send_queue.push(std::make_tuple(0, MPI_TAG_GLOBAL_CHECK, nullNode));
		queue_condition[m_options.threads_num].notify_one();		
	}

	sync_threads();

    if (end_cond == false)							  
    {
		ClosedList[tid].erase(n.pos);
		if (n.pos.get_id(m_options.totalThreads) == ((unsigned int)tid+m_options.mpiMin))
		{
			OpenList[tid].conditional_enqueue(n);
		}
        return true;
    }

    return false;
}

/*!
* Check end phase 3.
* After everyone agreed that a possible answer is found, we must exchange
* messages between mpi nodes to get the global best answer, and then 
* choose wheter to continue our finish.
* This functions allows the removal off mpi barrier of sync function
*/
template < int N >
bool PAStar<N>::check_stop_global(int tid)
{
	// Each rank 0 thread send the process final node to other rank 0
	if (tid == 0)
	{
		//std::cout << m_options.mpiRank << "-" << tid << ": global checking" << std::endl;
		Node<N> nullNode;
		std::unique_lock<std::mutex> queue_lock(queue_mutex[m_options.threads_num]);
		send_queue.push(std::make_tuple(0, MPI_TAG_FINAL_GLOBAL_CHECK, nullNode));
		queue_condition[m_options.threads_num].notify_one();

	}
	//Sync all local threads
	/*sync_threads_local();

	Node<N> n = final_node;
	if (end_condLocal == false)
	{
		ClosedList[tid].erase(n.pos);

		if (n.pos.get_id(m_options.threads_num) == ((unsigned int)tid))
		{
			OpenList[tid].conditional_enqueue(n);
		}
	}*/

	sync_threads_local();
	return !end_cond;
}

//! Execute a worker thread. This thread have id \a tid
template < int N >
int PAStar<N>::worker(int tid, const Coord<N> &coord_final)
{
#ifndef WIN32
	set_affinity(tid);
#endif
	//sync all threads from all nodes to prevent problems
	sync_threads_local();
	sync_threads();

	// worker_inner is the main inner loop
	// check_stop syncs and check if is the optimal answer
	// and check_stop_global stops if the local answer is the global optimum
    //do 
	//{
		do 
		{
			worker_inner(tid, coord_final);
		} while (check_stop(tid));
		//std::cout << m_options.mpiRank << "-" << tid << ": threads agreed in a result" << std::endl;
		sync_threads();
		//consume_queue(tid);
		//if (OpenList[tid].get_highest_priority() < final_node.get_f())
		//{
		//	continue;
		//}
	//} while (check_stop_global(tid));

	//the first local thread of each node kill sender and receiver to prevent problems with MPI exchange 
	if (tid == 0)
	{
		// awake receiver to die
		Node<N> nullNode;
		std::unique_lock<std::mutex> queue_lock(queue_mutex[m_options.threads_num]);
		send_queue.push(std::make_tuple(0, MPI_TAG_KILL_RECEIVER, nullNode));

		// awake sender to send message to receiver and after that kill itself
		queue_condition[m_options.threads_num].notify_one();
	}

	//std::cout << "waiting for sync" << std::endl;
	sync_threads();
	
    return 0;
}

enum Flags { normal, barrier, check, final_check, end_receiver};

//! Execute MPI sender thread, with tid = numThreads + 1
template < int N >
int PAStar<N>::sender()
{
	int i = 0, tag = 0;
	const char b[] = "1";
	bool goodbye = false, empty = true;
	long long int local, val, remoteValues[2], *recvBuff = nullptr;	

	
	Flags flags = normal;
	

	while (!goodbye | !(empty = send_queue.empty()))
	{
		//std::unique_lock<std::mutex> lock(*senderWait);
		//std::cout << "sender continuing? " << std::endl;
		
		// Polling queue
		if (empty == true)
		{
			switch(flags)
			{
				case barrier:
					{
						std::unique_lock<std::mutex> lock(sync_mutex_global);
						std::cout << m_options.mpiRank << ": barrier with message queue: " << send_queue.size() << std::endl;
						MPI_Barrier(MPI_COMM_WORLD);
						flags = normal;
						sync_condition_global.notify_all();
					}
					break;
				case check:
					std::cout << m_options.mpiRank << ": sender queue size at global ckeck" << send_queue.size() << std::endl;
					// If end_cond = true, then we might be in a false end where everyone agrees with its own final node;
					local = final_node.get_f();

					MPI_Allreduce(&local, &val, 1, MPI_LONG_LONG_INT, MPI_MIN, MPI_COMM_WORLD);

					// Preserve if end flag was already up
					end_condLocal = (end_cond == true) ? true : false;

					if (local != val)
					{
						end_condLocal = false;
						std::cout << m_options.mpiRank << ":local node isnt global minimum" << std::endl;
					}
					MPI_Allreduce(&end_condLocal, &end_cond, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
					flags = normal;
					break;
				case final_check:
					std::cout << m_options.mpiRank << ": sender queue size at final global ckeck" << send_queue.size() << std::endl;

					// Allgather to sync and get if we continue
					val = final_node.get_f();
					//first value is the minumum and the second the maximum
					remoteValues[0] = LLONG_MAX;
					remoteValues[1] = 0;

					if (m_options.mpiRank == 0)
						recvBuff = new long long int[m_options.mpiCommSize]();

					MPI_Gather(&val, 1, MPI_LONG_LONG_INT, recvBuff, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

					if (m_options.mpiRank == 0)
					{
						long long int min = LLONG_MAX, max = 0;
						for (int i = 0; i < m_options.mpiCommSize; i++)
						{
							if (recvBuff[i] < remoteValues[0])
								remoteValues[0] = recvBuff[i];
							if (recvBuff[i] > remoteValues[1])
								remoteValues[1] = recvBuff[i];
						}
					}
					if (m_options.mpiRank == 0)
						delete[] recvBuff;
					MPI_Bcast(&remoteValues, 2, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
					std::cout << m_options.mpiRank << " - minVal: " << remoteValues[0] << " maxVal: " << remoteValues[1] << " localVal: " << val << std::endl;

					// If all remote nodes share the same maximum and minimum f-values, the result is optimal
					if (remoteValues[0] == remoteValues[1])
					{
						end_cond = true; //end of execution
					}
					else
					{
						end_cond = false;
						if (val != remoteValues[0])
							end_condLocal = false;
					}
					flags = normal;
					break;
				case end_receiver:
					std::cout << "sender killing receiver" << std::endl;
					MPI_Send(&b, 2, MPI_CHAR, m_options.mpiRank, tag, MPI_COMM_WORLD);
					goodbye = true;
					continue;	//break to loop
					break;
				default:
					break;
			}
			//std::cout << "locking sender mutex" << std::endl;
			std::unique_lock<std::mutex> queue_lock(queue_mutex[m_options.threads_num]);
			//std::cout << "waiting for something to be queued" << std::endl;
			queue_condition[m_options.threads_num].wait(queue_lock);
			//std::cout << "awaking" << std::endl;
			continue;
		}
		
		std::unique_lock<std::mutex> queue_lock(queue_mutex[m_options.threads_num]);
		std::tuple<int, int, Node<N>> temp;
		temp = send_queue.front();
		//std::cout << "continuing tag: " << std::get<1>(temp) << std::endl;
		send_queue.pop();
		queue_lock.unlock();
		
		switch (tag = std::get<1>(temp))
		{
			case MPI_TAG_KILL_RECEIVER:
				flags = end_receiver;
				break;
			case MPI_TAG_BARRIER:
				flags = barrier;
				break;
			case MPI_TAG_FINAL_GLOBAL_CHECK:
				flags = final_check;
				break;
			case MPI_TAG_GLOBAL_CHECK:
				flags = check;
				break;
			default:
			case MPI_TAG_SEND_COMMON:
				Node<N> n = std::get<2>(temp);
				//std::cout << "current node " << n << " and current final node " << final_node << std::endl;
				if (n.get_f() <= final_node.get_f())
				{
					//std::cout << "oh, dang it, sending nodes to remote target" << std::endl;
					std::stringstream ss;
					boost::archive::text_oarchive oa{ ss };
					oa & n;

					//precisamos enviar a thread alvo junto do conteudo dos nos
					if (tag == MPI_TAG_SEND_COMMON)
					{

						//std::cout << "sending: " << ss.str() << std::endl;
						int iLocal = std::get<0>(temp) % m_options.threads_num;
						int targetNode = std::get<0>(temp) / m_options.mpiMax;
						//std::cout << "target node " << targetNode << " thread " << iLocal << std::endl;
						MPI_Send(ss.str().c_str(), ss.str().size() + 1, MPI_CHAR, targetNode, iLocal, MPI_COMM_WORLD);
					}
					if (tag >= MPI_TAG_BROADCAST_NODE_TO_0)
					{
						for (i = 0; i < m_options.mpiCommSize; i++)
							if (i != m_options.mpiRank)
								MPI_Send(ss.str().c_str(), ss.str().size() + 1, MPI_CHAR, i, tag - MPI_TAG_BROADCAST_NODE_TO_0, MPI_COMM_WORLD);
					}
			}
		}
		send++;
		//lock.unlock();
	} 
	//std::cout << m_options.mpiRank << ": finishing sender with " << send << " messages sent" << std::endl;
	return 0;
}

//! Execute a receiver thread, with tid = numThreads + 2
template < int N >
int PAStar<N>::receiver(PAStar<N> * pastar_inst)
{
    MPI_Status status;
    int n_bytes = 0;
    int sender = 0;
    int sender_tag = 0;
    int i = 0, flag = 0;
    char * buffer = NULL;
	bool goodbye = false;

	while (!goodbye)
    {
		//std::unique_lock<std::mutex> lock(*senderWait);
		// std::cout << "waiting for probe" << std::endl;
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		//MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);

		//if (!flag)
		//	continue;

		MPI_Get_count(&status, MPI_CHAR, &n_bytes);
		sender = status.MPI_SOURCE;
		sender_tag = status.MPI_TAG;

		//std::cout << m_options.mpiRank << ":probe returned node " << status.MPI_SOURCE << " sending " << n_bytes << " bytes of data" << std::endl;
		if (n_bytes == 0)
		{
			MPI_Recv(NULL, 0, MPI_CHAR, sender, sender_tag, MPI_COMM_WORLD, &status);
			continue;
		}

		buffer = new char[n_bytes]();
		//Receive thread destination plus nodes
		MPI_Recv(buffer, n_bytes, MPI_CHAR, sender, sender_tag, MPI_COMM_WORLD, &status);

		if (sender_tag == MPI_TAG_KILL_RECEIVER)
		{
			//std::cout << "we are finishing, goodbye" << std::endl;
			goodbye = true;

			//as safety measure, wake up receiver to finish
			queue_condition[m_options.threads_num].notify_one();
		}
		else
		{
			process_message(sender_tag, buffer);
			//std::thread t(&PAStar::process_message, pastar_inst, sender_tag, buffer);
			//t.detach();
			//create process_message a thread passing buffer and sender_tag
		}
		recv++;
		//lock.unlock();
	} 
	//std::cout << m_options.mpiRank << ": finishing receiver with " << recv << " messages received" << std::endl;
    return 0;
}

//! Execute a receiver thread, with tid = numThreads + 2
template < int N >
int PAStar<N>::process_message(int sender_tag, char *buffer)
{
	//TODOOOOO: use a global control variable, that is incremented when start processing a message and decremented when it finishes
	//so that the thread that is going to receive the node keep waiting with std::this_thread::yield() until its control variable is equal to 0
	//or make it sleep until process message control variable dont wake it up
	std::shared_lock<std::shared_mutex> lock (processing_mutex[sender_tag]);

	//convert buffer to something useful, like thread id and node vector
	std::stringstream ss(buffer);
	boost::archive::text_iarchive ia{ ss };
	//std::cout << "hey, rank " << m_options.mpiRank << " received " << n_bytes << " of data, saving to tid " << sender_tag << std::endl;
	Node<N> neigh;
	ia & neigh;
	//std::cout << "node being received by " << m_options.mpiRank << "-" << ": " << ss.str() << std::endl;

	
	//if (neigh.get_f() < final_node.get_f())
	//{	
	//	std::unique_lock<std::mutex> best_received_lock(best_received_mutex);
	//	best_received = neigh;
	//	//std::cout << m_options.mpiRank << ": receiver replacing best received ->" << best_received << std::endl;
	//	best_received_lock.unlock();
	//}

	// Locking queue to add received stuff
	//std::cout << "receiving: " << ss.str() << std::endl;
	std::unique_lock<std::mutex> queue_lock(queue_mutex[sender_tag]);
	queue_nodes[sender_tag].push_back(neigh);
	queue_condition[sender_tag].notify_one();
	
	delete[] buffer;

	return 0;
}

template < int N >
void PAStar<N>::print_nodes_count()
{
    long long int nodes_total = 0;
    long long int open_list_total = 0;
    long long int closed_list_total = 0;
    long long int nodes_reopen_total = 0;

    std::cout << "Total nodes count:" << std::endl;

    for (int i = 0; i < m_options.totalThreads; ++i)
    {
		std::cout << "tid " << i
			 << "\tOpenList:" << nodes_openListSizeFinal[i]
             << "\tClosedList:" << nodes_closedListSizeFinal[i]
             << "\tReopen:" << nodes_reopenFinal[i]
             << "\tTotal: " << nodes_countFinal[i] << std::endl;
        open_list_total += nodes_openListSizeFinal[i];
        closed_list_total += nodes_closedListSizeFinal[i];
        nodes_reopen_total += nodes_reopenFinal[i];
        nodes_total += nodes_countFinal[i];
    }

    std::cout << "Sum"
          << "\tOpenList:" << open_list_total
          << "\tClosedList:" << closed_list_total
          << "\tReopen:" << nodes_reopen_total
          << "\tTotal: " << nodes_total << std::endl;
}

template < int N >
void PAStar<N>::print_answer()
{
    backtrace<N>(ClosedListFinal, m_options.totalThreads);
    print_nodes_count();
	//std::cout << "best received " << best_received << std::endl;
}

/*!
 * To calculate the correct statistics and do the
 * backtrace, we need to merge work data in rank 0 node.
 */
template < int N >
void PAStar<N>::sync_pastar_data()
{
    // Before backtracing, the rank 0 receives bunch of data from every other rank
	//MPI_Status status;
	int n_bytes = 0;
	int i, k;

    // Send all node count and reopen data to rank 0
    MPI_Gather(nodes_count, m_options.threads_num, MPI_LONG_LONG_INT, nodes_countFinal, m_options.threads_num, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(nodes_reopen, m_options.threads_num, MPI_LONG_LONG_INT, nodes_reopenFinal, m_options.threads_num, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

	// Send open and closed list data from all nodes to rank 0
	long long int * nodes_openListSize, *nodes_closedListSize;
	nodes_openListSize = new long long int[m_options.threads_num];
	nodes_closedListSize = new long long int[m_options.threads_num];

	for (i = 0; i < m_options.threads_num; i++)
	{
		nodes_openListSize[i] = OpenList[i].size();
		nodes_closedListSize[i] = ClosedList[i].size();
	}

	//std::cout << "gathering list sizes" << std::endl;
	MPI_Gather(nodes_openListSize, m_options.threads_num, MPI_LONG_LONG_INT, nodes_openListSizeFinal, m_options.threads_num, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
	MPI_Gather(nodes_closedListSize, m_options.threads_num, MPI_LONG_LONG_INT, nodes_closedListSizeFinal, m_options.threads_num, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

	delete[] nodes_openListSize;
	delete[] nodes_closedListSize;

	// Send all closed lists from all nodes to rank 0

	//First we need to serialize
	std::stringstream ss;
	boost::archive::text_oarchive oa{ ss };
	for (k = 0; k < N + 1; k++)
	{
		oa << ClosedList[k];
	}

	//Then send serialized data sizes to rank 0
	int ss_size = ss.str().size()+1;
	int * ss_sizes = nullptr;
	if (m_options.mpiRank == 0)
		ss_sizes = new int[m_options.mpiCommSize]();
	MPI_Gather(&ss_size, 1, MPI_INT, ss_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//Then send everything to node 0
	char * ss_ss = nullptr;
	int * ss_disp = nullptr;
	if (m_options.mpiRank == 0)
	{
		int totalLen = 0;
		ss_disp = new int[m_options.mpiCommSize]();

		//compute displacements of data sent
		ss_disp[0] = 0;
		totalLen += ss_sizes[0];
		for (k = 1; k < m_options.mpiCommSize; k++)
		{
			ss_disp[k] = totalLen;
			totalLen += ss_sizes[k];		
		}
		ss_ss = new char[totalLen]();
	}

	//send serialized closed lists to node 0
	MPI_Gatherv(ss.str().c_str(), ss_size, MPI_CHAR, ss_ss, ss_sizes, ss_disp, MPI_CHAR, 0, MPI_COMM_WORLD);

	if (m_options.mpiRank == 0)
	{
		char * temp = nullptr;
		//for each node
		for (i = 0; i < m_options.mpiCommSize; i++)
		{
			//read the block for current node
			temp = new char[ss_sizes[i]]();
			memcpy(temp, &ss_ss[ss_disp[i]], ss_sizes[i]);

			//std::cout << temp << std::endl;
			std::stringstream ss(temp);
			boost::archive::text_iarchive ia{ ss };

			//and then unserialize each closed list into final closed list
			for (k = 0; k < m_options.threads_num; k++)
			{
				ia >> ClosedListFinal[i*m_options.threads_num+k];
				//std::cout << "closedListFinal[" << i*m_options.threads_num + k << "] is " << ClosedListFinal[i*m_options.threads_num + k].size() << std::endl;
			}
			delete[] temp;
		}
		delete[] ss_ss;
		delete[] ss_disp;
		delete[] ss_sizes;
	}
}

/*!
 * Same a_star() function usage.
 * Starting function to do a pa_star search.
 */
template < int N >
int PAStar<N>::pa_star(const Node<N> &node_zero, const Coord<N> &coord_final, const PAStarOpt &options)
{
    if (options.threads_num <= 0)
        throw std::invalid_argument("Invalid number of threads");
    Coord<N>::configure_hash(options.hash_type, options.hash_shift);

    PAStar<N> pastar_instance(node_zero, options);
    std::vector<std::thread> threads;
    TimeCounter *t = new TimeCounter("Phase 2: PA-Star running time: ");

    // Create worker threads
    for (int i = 0; i < options.threads_num; ++i)
        threads.push_back(std::thread(&PAStar::worker, &pastar_instance, i, coord_final));

	// Create sender thread
	threads.push_back(std::thread(&PAStar::sender, &pastar_instance));

	// Create receiver thread
	threads.push_back(std::thread(&PAStar::receiver, &pastar_instance, &pastar_instance));

    // Wait for the end of all threads
    for (auto& th : threads)
        th.join();
    delete t;

    //std::cout << options.mpiRank << ":olha eu aqui, mãe" << std::endl;
	
	// Don't you dare removing that barrier
    MPI_Barrier(MPI_COMM_WORLD);
    
	pastar_instance.sync_pastar_data();
	//std::cout << "preparing to print answer" << std::endl;

    if (options.mpiRank == 0)
        pastar_instance.print_answer();

	//std::cout << options.mpiRank << ": waiting for other processes to finish" << std::endl;
	
    return 0;
}

#define PASTAR_DECLARE_TEMPLATE( X ) \
template class PAStar< X >; \

MAX_NUM_SEQ_HELPER(PASTAR_DECLARE_TEMPLATE);
