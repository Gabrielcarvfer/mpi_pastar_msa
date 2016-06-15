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
#include <stdint.h>
#include <lz4.h>

#include <sstream>
// include input and output archivers
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

// include this header to serialize vectors, maps and other types
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/collection_size_type.hpp>

typedef uint32_t u4;
#define s_u4 (sizeof(u4))

typedef int64_t u8;

#define MPI_TAG_SEND_COMMON                0
#define MPI_TAG_KILL_RECEIVER            205


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
    //if (m_options.mpiRank == 0)
    //{
    std::cout << "Running PAStar with: "
              << opt.threads_num << " threads, "
              << Coord<N>::get_hash_name() << " hash, "
              << Coord<N>::get_hash_shift() << " shift.\n";
    //}

    end_cond = false;
    end_condLocal = false;
    sync_count = 0;
    recv_cnt = 0;
    final_node.set_max();

	check_stop_mutex = new std::mutex;

    OpenList = new PriorityList<N>[m_options.threads_num]();
    ClosedList = new std::map< Coord<N>, Node<N> >[m_options.threads_num]();

    nodes_count = new u8[m_options.threads_num]();
    nodes_reopen = new u8[m_options.threads_num]();
	

	// The additional mutex and condition are used by sender thread
    queue_mutex = new std::mutex[m_options.threads_num+1];	  
    queue_condition = new std::condition_variable[m_options.threads_num+1];
    queue_nodes = new std::vector< Node<N> >[m_options.threads_num];

    send_queue = new std::vector< Node<N> > [m_options.totalThreads+1]();
    end_of_transmission = false;

    OpenListFinal = NULL;
    ClosedListFinal = NULL;
    nodes_countFinal = NULL;
    nodes_reopenFinal = NULL;

	threadLookupTable = new int[m_options.totalThreads*2]();
	int i = 0, j = 0, k = 0;
	for (; k < m_options.mpiCommSize; k++)
	{
		for (j = 0; j < m_options.threads_num; j++)
		{
			threadLookupTable[i] = k;
			threadLookupTable[m_options.totalThreads + i] = j;
			//std::cout << i << "-" << threadLookupTable[i] << "-" << threadLookupTable[m_options.totalThreads + i] << std::endl;
			i++;
		}
	}
	//senderWait = new std::mutex();

    // Allocate final structures and enqueue first node if rank zero
    if (m_options.mpiRank == 0)
    {
        OpenListFinal = new PriorityList<N>[m_options.totalThreads];
        ClosedListFinal = new std::map< Coord<N>, Node<N> >[m_options.totalThreads];
        nodes_countFinal = new u8[m_options.totalThreads]();
        nodes_reopenFinal = new u8[m_options.totalThreads]();
		nodes_openListSizeFinal = new u8[m_options.totalThreads]();

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
		delete[] nodes_openListSizeFinal;
		delete[] nodes_countFinal;
		delete[] nodes_reopenFinal;

    }
}
#ifndef WIN32
template < int N >
int PAStar<N>::set_affinity(int tid)
{
    //std::cout << "setting affinity of tid " << tid << " rank " << m_options.mpiRank << std::endl;
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
void PAStar<N>::sync_threads(bool flushing)
{
    //std::cout << "sync threads" << std::endl;

    std::unique_lock<std::mutex> sync_lock(sync_mutex);

    if (++sync_count < m_options.threads_num)
        sync_condition.wait(sync_lock);
    else
    {
        //When every thread is blocked, block processes in a barrier
		if (flushing)
		{
			flush_sender();
			MPI_Barrier(MPI_COMM_WORLD);
			flush_receiver();
		}
		MPI_Barrier(MPI_COMM_WORLD);

        //After that, restart everyone
        sync_count = 0;
        sync_condition.notify_all();
    }
}

template <int N>
void PAStar<N>::flush_sender()
{
	if (!sender_empty)
	{
		std::unique_lock<std::mutex> sender_lock(sender_mutex);
		//std::cout << m_options.mpiRank << ": esperando pelo sender with queue size " << send_queue.size() << std::endl;
		queue_condition[m_options.threads_num].notify_one(); //acorda sender se estiver dormindo
		sender_condition.wait(sender_lock);
	}
	//std::cout << m_options.mpiRank << ": buffer do sender limpo" << std::endl;
	return;
}

template <int N>
void PAStar<N>::flush_receiver()
{
	std::unique_lock<std::mutex> lock(processing_mutex);   
	if (recv_cnt > 0)
	{
		lock.unlock();
		std::unique_lock<std::mutex> receiver_lock(receiver_mutex);
		//std::cout << m_options.mpiRank << ": esperando pelo receiver working with " << recv_cnt << std::endl;
		receiver_condition.wait(receiver_lock);
	}
	//std::cout << m_options.mpiRank << ": buffer do receiver limpo" << std::endl;
	return;
}

//! Sync all local threads
template < int N >
void PAStar<N>::sync_threads_local()
{
	//std::cout << "sync threads" << std::endl;

	std::unique_lock<std::mutex> sync_lock(sync_mutex);

	if (++sync_count < m_options.threads_num)
	{
		sync_condition.wait(sync_lock);
	}
	else
	{
		sync_count = 0;
		sync_condition.notify_all();
	}
}

//! Execute the pa_star algorithm until all nodes expand the same final node
template < int N >
void PAStar<N>::worker_inner(int tid, const Coord<N> &coord_final)
{
    Node<N> current;
    //std::cout << m_options.mpiRank << ": fase 6 3 worker_inner " << tid  << "fase 1" << std::endl;

    std::vector< Node<N> > *neigh = new std::vector< Node<N> >[m_options.totalThreads];

    //std::cout << m_options.mpiRank << ": fase 6 3 worker_inner " << tid  << "fase 2" << std::endl;
    // Loop ended by process_final_node
    while (end_condLocal == false)
    {
        typename std::map< Coord<N>, Node<N> >::iterator c_search;

        // Start phase
		//flush_receiver();

		// Consume openlist to find if someone have a better route
		consume_queue(tid);
//std::cout << m_options.mpiRank << ": fase 6 3 worker_inner " << tid  << "fase 3" << std::endl;
    
        // Dequeue phase
        if (OpenList[tid].dequeue(current) == false)
        {
            //std::cout << "dequeueing" << std::endl;
            wait_queue(tid);
            //std::cout << "good to go" << std::endl;
            continue;
        }
        nodes_count[tid] += 1;
//std::cout << m_options.mpiRank << ": fase 6 3 worker_inner " << tid  << "fase 4" << std::endl;
    
        // Check if better node is already found
        if ((c_search = ClosedList[tid].find(current.pos)) != ClosedList[tid].end())
        {
            //std::cout << "searching for a better node" << std::endl;
            if (current.get_g() >= c_search->second.get_g())
                //std::cout << "oops, we didn't found it, so we are good" << std::endl;
                continue;
            nodes_reopen[tid] += 1;
        }
//std::cout << m_options.mpiRank << ": fase 6 3 worker_inner " << tid  << "fase 5" << std::endl;
    
        //std::cout << "[" << tid << "] Opening node:\t" << current << std::endl;
        ClosedList[tid][current.pos] = current;
											  
        if (current.pos == coord_final)
        {
            //std::cout << m_options.mpiRank << "-" << tid << ": achei essa porra" << std::endl;
            process_final_node(tid, current);
            continue;
        }
//std::cout << m_options.mpiRank << ": fase 6 3 worker_inner " << tid  << "fase 6" << std::endl;
    
        // Expand phase
        current.getNeigh(neigh, m_options.totalThreads);

        int actualTid = tid + m_options.mpiMin;
        //std::cout << actualTid << std::endl;
//std::cout << m_options.mpiRank << ": fase 6 3 worker_inner " << tid  << "fase 7" << std::endl;
    
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
					int iLocal = threadLookupTable[i + m_options.threads_num];
                    //std::cout << "happily, a local node" << std::endl;
                    std::lock_guard<std::mutex> queue_lock(queue_mutex[iLocal]);
                    queue_nodes[iLocal].insert(queue_nodes[iLocal].end(), neigh[i].begin(), neigh[i].end());
                    queue_condition[iLocal].notify_one();
                }
                //or remote
                else
				{
					//std::cout << "adding node to sender" << std::endl;
					std::lock_guard<std::mutex> queue_lock(queue_mutex[m_options.threads_num]);
					if (neigh[i].front().get_f() <= final_node.get_f())
						send_queue[i].insert(send_queue[i].end(), neigh[i].begin(), neigh[i].end());
					queue_condition[m_options.threads_num].notify_one();
                }
            }
            neigh[i].clear();
        }
	//std::cout << m_options.mpiRank << ": fase 6 3 worker_inner " << tid  << "fase 8" << std::endl;
    
    }
//std::cout << m_options.mpiRank << ": fase 6 3 worker_inner " << tid  << "fase 9" << std::endl;
    
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
//std::cout << m_options.mpiRank << ": fase 6 3 5 process_final_node " << tid  << "fase 1" << std::endl;
    
	unsigned int originRank = threadLookupTable[n.pos.get_id(m_options.totalThreads)];
	//std::cout << "process final node" << std::endl;

	//Lock all threads from all nodes
	std::unique_lock<std::mutex> final_node_lock(final_node_mutex);
//std::cout << m_options.mpiRank << ": fase 6 3 5 process_final_node " << tid  << "fase 2" << std::endl;
   
	// Better possible answer already found, discard n
	if (final_node.get_f() < n.get_f())
	{
		//std::cout << "newbie, we already got a better candidate" << std::endl;
		//return;
	}
	else
	{
//std::cout << m_options.mpiRank << ": fase 6 3 5 process_final_node " << tid  << "fase 3" << std::endl;
   
		if (n.pos.get_id(m_options.threads_num) == ((unsigned int)tid))
		{
			std::cout << "[proc: " << m_options.mpiRank << " - tid: " << tid << "] Possible answer found: " << n << std::endl;
			// Broadcast the node
			final_node = n;
			final_node_count = 0;

			//If the tid is shared with founder of the final node, add to other local nodes
			for (i = 0; i < m_options.threads_num; i++)
			{
				if (originRank != m_options.mpiRank)
				{
					std::cout << "Wow, broadcasting remote node to fellows" << std::endl;
				}
				if (i != tid)
				{
					std::lock_guard<std::mutex> queue_lock(queue_mutex[i]);
					queue_nodes[i].push_back(n);
					queue_condition[i].notify_one();
				}
			}			

			
			//If thread found the node, share it with other processes
			if (originRank == (unsigned int) m_options.mpiRank)
			{
				std::cout << m_options.mpiRank << "- origin " << originRank << "-" << tid << ": broadcasting " << n << " to all other nodes" << std::endl;
				//For remote threads: send messages to be broadcasted between node threads
				std::lock_guard<std::mutex> queue_lock(queue_mutex[m_options.threads_num]);
				for (int i = tid; i < m_options.totalThreads; i=i+m_options.threads_num)
				{
					send_queue[i].push_back(n);
				}
				queue_condition[m_options.threads_num].notify_one();
			}

			final_node_lock.unlock();
		}
		// Every other worker is unlocked
		else
		{		  
			std::cout << "[" << m_options.mpiRank << ":" << tid << "] Agreed with possible answer! " << n << "/" << final_node << std::endl;
			//if (n != final_node) std::cout << "BUG HERE!\n";
			final_node_lock.unlock();
		}

		// This node have the highest priority between all Openlist. Broadcast the end condition
		if (++final_node_count == m_options.threads_num)
		{
			//end_cond = true;
			end_condLocal = true;
		}
	}
//std::cout << m_options.mpiRank << ": fase 6 3 5 process_final_node " << tid  << "fase 4" << std::endl;
   
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
//std::cout << m_options.mpiRank << ": fase 6 3 check_stop " << tid  << "fase 1" << std::endl;
   
	u8 local, val;
	//std::cout << m_options.mpiRank << " : entering check stop" << std::endl;
	Node<N> n = final_node;
	
	wake_all_queue();

	sync_threads(true);

	// Consume openlist to find if someone have a better route
	consume_queue(tid);
//std::cout << m_options.mpiRank << ": fase 6 3 check_stop " << tid  << "fase 2" << std::endl;
	std::cout << m_options.mpiRank << "-" << tid << ": checking end f-value = " << final_node.get_f() << " | " << OpenList[tid].get_highest_priority()<< std::endl;
	// If someone have a better value, we're not at the end yet	   -> we could look for minimum value inside node and then use MPI_Allreduce to find global minimum
	if (OpenList[tid].get_highest_priority() < final_node.get_f())
	{
		//std::cout << m_options.mpiRank << "-" << tid << ": looks like i've found a better route " << OpenList[tid].get_highest_priority() << " node: " << OpenList[tid].get_highest_priority_node() << std::endl;
		
		end_condLocal = false;
	}

	//sync_threads(false);
//std::cout << m_options.mpiRank << ": fase 6 3 check_stop " << tid  << "fase 3" << std::endl;
	if (tid == 0)
	{
		//std::cout << m_options.mpiRank << ": sender queue size at global ckeck" << send_queue.size() << std::endl;
		// If end_cond = true, then we might be in a false end where everyone agrees with its own final node;
		local = final_node.get_f();

		MPI_Allreduce(&local, &val, 1, MPI_LONG_LONG_INT, MPI_MIN, MPI_COMM_WORLD);

		if (local != val)
		{
			end_condLocal = false;
			//std::cout << m_options.mpiRank << ":local node isnt global minimum" << std::endl;
		}

		//Some MPI implementations don't like atomic booleans, so here's a possible workaround
		MPI_Allreduce(&end_condLocal, &end_cond, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);

		//bool endLocal = end_condLocal, endGlobal = end_cond;
		//MPI_Allreduce(&endLocal, &endGlobal, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
		
		std::cout << m_options.mpiRank << "end_condLocal " << end_condLocal << " while end_cond " << end_cond << std::endl;
		//std::cout << m_options.mpiRank << ": " << local << " " << val << std::endl;

	}
//std::cout << m_options.mpiRank << ": fase 6 3 check_stop " << tid  << "fase 4" << std::endl;
	sync_threads_local();
	//sync_threads(false);

	//If someone found a better node
    if (!end_cond)							  
    {
		//If it was on a local thread of that rank
		if (!end_condLocal)
		{
			//Erase previous final node
			ClosedList[tid].erase(n.pos);

			//Enqueue node
			//if (n.pos.get_id(m_options.totalThreads) == (unsigned int)tid+m_options.mpiMin)
			if (n.pos.get_id(m_options.threads_num) == (unsigned int) tid)
			{
				OpenList[tid].conditional_enqueue(n);
			}
		}
		//If from a remote rank, just finish
//std::cout << m_options.mpiRank << ": fase 6 3 check_stop " << tid  << "fase 5" << std::endl;
        return true;
    }
//std::cout << m_options.mpiRank << ": fase 6 3 check_stop " << tid  << "fase 6" << std::endl;
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
	u8 val, remoteValues[2], *recvBuff;
	sync_threads(false);
	// Each rank 0 thread send the process final node to other rank 0
	if (tid == 0)
	{
		//std::cout << m_options.mpiRank << ": sender queue size at final global ckeck " << send_queue.size() << std::endl;

		// Allgather to sync and get if we continue
		val = final_node.get_f();
		//first value is the minumum and the second the maximum
		remoteValues[0] = LLONG_MAX;
		remoteValues[1] = 0;

		if (m_options.mpiRank == 0)
			recvBuff = new u8[m_options.mpiCommSize]();

		MPI_Gather(&val, 1, MPI_LONG_LONG_INT, recvBuff, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

		if (m_options.mpiRank == 0)
		{
			u8 min = LLONG_MAX, max = 0;
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
		//std::cout << m_options.mpiRank << " - minVal: " << remoteValues[0] << " maxVal: " << remoteValues[1] << " localVal: " << val << std::endl;

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

	}
	//Sync all local threads
	sync_threads_local();

	Node<N> n = final_node;
	if (end_condLocal == false)
	{
		ClosedList[tid].erase(n.pos);

		if (n.pos.get_id(m_options.totalThreads) == ((unsigned int)tid+m_options.mpiMin))
		{
			OpenList[tid].conditional_enqueue(n);
		}
	}
													
	sync_threads_local();										
	return !end_cond;
}

//! Execute a worker thread. This thread have id \a tid
template < int N >
int PAStar<N>::worker(int tid, const Coord<N> &coord_final)
{
	//std::cout << m_options.mpiRank << ": fase 6 worker " << tid  << "fase 1" << std::endl;
#ifndef WIN32
	set_affinity(tid);
#endif
	//std::cout << m_options.mpiRank << ": fase 6 worker " << tid  << "fase 2" << std::endl;
	//sync all threads from all nodes to prevent problems
	sync_threads(false);

	//std::cout << m_options.mpiRank << ": fase 6 worker " << tid  << "fase 3" << std::endl;
	// check_stop_global stops if all the local answer is the global optimum
    //do 
	//{
		// check_stop syncs and check if local nodes agreed in local optimum
		do
		{
			// worker_inner is the main inner loop
			worker_inner(tid, coord_final);
		} while (check_stop(tid));
		
	//} while (check_stop_global(tid));

	//std::cout << "syncing global threads" << std::endl;
	//std::cout << m_options.mpiRank << ": fase 6 worker " << tid  << "fase 4" << std::endl;
	//global sync with flush to garantee that no node will be stuck because a sender/receiver thread didn't died at right time
	sync_threads(true);

	//std::cout << "preparing to exit" << std::endl;
	//the first local thread of each node kill sender and receiver to prevent problems with MPI exchange 
	if (tid == 0)
	{
		// awake receiver to die
		Node<N> nullNode;
		std::lock_guard<std::mutex> queue_lock(queue_mutex[m_options.threads_num]);
		send_queue[m_options.totalThreads].push_back(nullNode);

		// awake sender to send message to receiver and after that kill itself
		queue_condition[m_options.threads_num].notify_one();
	}

	//std::cout << m_options.mpiRank << ": fase 6 worker " << tid  << "fase 5" << std::endl;
	
	sync_threads_local();
	//if (tid == 0)
	//{	
	//std::lock_guard<std::mutex> lock(queue_mutex[m_options.threads_num]);
	//	end_of_transmission = true;	
	//	delete[] send_queue;
	//}
	//std::cout << "waiting for sync" << std::endl;	
    return 0;
}

//! Execute MPI sender thread, with tid = numThreads + 1
template < int N >
int PAStar<N>::sender()
{
	int i = 0, tag = 0;
	const char b[] = "1";
	bool goodbye = false;
	
	std::cout << m_options.mpiRank << ": sender fase 1" << std::endl;
	
	while (!goodbye)
	{
		std::unique_lock<std::mutex> queue_lock(queue_mutex[m_options.threads_num]);
		// Polling queue
		sender_empty = true;
		//std::cout << m_options.mpiRank << ": checking if queue is empty" << std::endl;
		
		//Select a thread to send all its nodes
		for (i = 0; i < m_options.totalThreads; i++)
		{
			if (!send_queue[i].empty())
			{
				sender_empty = false;

				//Copy to auxiliary buffer
				std::vector<Node<N>> temp;
				temp.insert(temp.end(), send_queue[i].begin(), send_queue[i].end());
				send_queue[i].clear();

				//Unlock send_queue for continuous operation
				queue_lock.unlock();

				//std::cout << "oh, dang it, sending nodes to remote target" << std::endl;
				std::ostringstream ss (std::ios_base::binary);
				boost::archive::binary_oarchive oa{ ss };

				oa & temp;
				temp.clear();

				//Experimental LZ4 compression
				u4 len = (u4) (ss.str().length() + 1);
				u4 lz4len = LZ4_compressBound( (u4) len );
				u4 lz4len2 = 0;
				char * lz4Data = new char[lz4len + ( s_u4*3)]();
					

				lz4len2 = LZ4_compress(ss.str().c_str(), &lz4Data[s_u4*3], (u4) len);
				((u4*)lz4Data)[0] = lz4len2+(s_u4*3);
				((u4*)lz4Data)[1] = lz4len;
				((u4*)lz4Data)[2] = len;					
				//std::cout << m_options.mpiRank << ": sended length" << len << std::endl;
				//std::cout << m_options.mpiRank << ": input sized " << len << " and maximum Lz4 sized " << lz4len2 << std::endl;
				//std::cout << "sending: " << ss.str() << std::endl;
				//myFileLz4 << new std::string(lz4Data) << std::endl;
				//precisamos enviar a thread alvo junto do conteudo dos nos

				//int iLocal = i % m_options.threads_num;
				//int targetNode = i / m_options.threads_num;
				MPI_Send(lz4Data, (int) (lz4len2 + ( s_u4 * 3 ) ), MPI_CHAR, threadLookupTable[i], threadLookupTable[i+m_options.totalThreads], MPI_COMM_WORLD);//targetNode, iLocal, MPI_COMM_WORLD);

				//Lock again for next round
				queue_lock.lock();
			}
		}

		// If empty
		if (sender_empty)
		{

			//If some node in killing queue, finish sender and receiver
			if (!send_queue[m_options.totalThreads].empty() & !goodbye)
			{
				//std::cout << m_options.mpiRank << " sending finishing message" << std::endl;
				MPI_Send((void*)&b, 2, MPI_CHAR, m_options.mpiRank, MPI_TAG_KILL_RECEIVER, MPI_COMM_WORLD);
				goodbye = true;
				send_queue[m_options.totalThreads].clear();
			}

			// If set goodbye, continue to finish
			if (goodbye)
			{
				//std::cout << m_options.mpiRank << " saying goodbye" << std::endl;
				sender_condition.notify_one();
			}
			else
			{
				// Sinaliza buffer limpo e acorda quem estiver esperando
				//std::cout << "sender without work, going to sleep" << std::endl;
				sender_condition.notify_one();
				queue_condition[m_options.threads_num].wait(queue_lock);
			}
		}
		//At the end, unlock the lock
		queue_lock.unlock();

		send++;
	} 

	delete[] threadLookupTable;

	// Last notify to unlock flusher
	sender_condition.notify_one();
	std::cout << m_options.mpiRank << ": sender fase 2" << std::endl;
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
    //std::cout << m_options.mpiRank << ": receiver fase 1" << std::endl;
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
			//std::cout << "bytes received " << n_bytes << " and marked " << ((int*)buffer)[0] << std::endl;
			if (n_bytes == ((int*)buffer)[0])
			{
			//Lock counter of messages awaiting to be processed
			std::lock_guard<std::mutex> lock(processing_mutex);
			recv_cnt++;

			//Create a processing thread passing buffer and sender_tag
			std::thread t(&PAStar::process_message, pastar_inst, sender_tag, buffer);
			t.detach();
			}
			else
			{
				delete[] buffer;
			}
		}
		recv++;
	} 
	// Last notify to unlock flushers
	std::cout << m_options.mpiRank << ": receiver fase 2" << std::endl;
	receiver_condition.notify_one();
    return 0;
}

//! Execute a receiver thread, with tid = numThreads + 2
template < int N >
int PAStar<N>::process_message(int sender_tag, char *buffer)
{
//std::cout << m_options.mpiRank << ": processer fase 1" << std::endl;
	//Experimental LZ4 decompression
	u4 lz4len = ( (u4*) buffer) [1];
	u4 len = ( (u4*) buffer) [2];

	char * tempBuff = new char[lz4len]();
	LZ4_decompress_fast(&buffer[s_u4 * 3], tempBuff, len);
	
	//Prepare to deserialize stuff into place
	std::istringstream ss(std::string(tempBuff, tempBuff + len), std::ios_base::binary);

	//convert buffer to something useful, like thread id and node vector
	boost::archive::binary_iarchive ia{ ss };
	std::vector<Node<N>> temp;
	ia & temp;
	delete[] tempBuff;

	// Locking queue to add received stuff
	std::unique_lock<std::mutex> queue_lock(queue_mutex[sender_tag]);
    queue_nodes[sender_tag].insert(queue_nodes[sender_tag].end(), temp.begin(), temp.end());
	queue_condition[sender_tag].notify_one();
	queue_lock.unlock();

	//std::cout << m_options.mpiRank << ": processer fase 4" << std::endl;

	delete[] buffer;
	temp.clear();

	// Update counter of messages awaiting to be processed
	std::lock_guard<std::mutex> processing_lock(processing_mutex);
	recv_cnt--;
	//std::cout << m_options.mpiRank << ": recv_cnt " << recv_cnt << std::endl;

	// Notify any process awaiting for receiver queue to be clean
	if (recv_cnt == 0)
	{
		//std::cout << m_options.mpiRank << ": receiver finished processing" << std::endl;
		receiver_condition.notify_one();
	}
//std::cout << m_options.mpiRank << ": processer fase 5" << std::endl;
	return 0;
}

template < int N >
void PAStar<N>::print_nodes_count()
{
    u8 nodes_total = 0;
    u8 open_list_total = 0;
    u8 closed_list_total = 0;
    u8 nodes_reopen_total = 0;

    std::cout << "Total nodes count:" << std::endl;

    for (int i = 0; i < m_options.totalThreads; ++i)
    {
		std::cout << "tid " << i
			 << "\tOpenList:" << nodes_openListSizeFinal[i]
             << "\tClosedList:" << ClosedListFinal[i].size()
             << "\tReopen:" << nodes_reopenFinal[i]
             << "\tTotal: " << nodes_countFinal[i] << std::endl;
        open_list_total += nodes_openListSizeFinal[i];
        closed_list_total += ClosedListFinal[i].size();
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
}

/*!
 * To calculate the correct statistics and do the
 * backtrace, we need to merge work data in rank 0 node.
 */
template < int N >
void PAStar<N>::sync_pastar_data()
{
	std::cout << m_options.mpiRank << ": syncing all the shit" << std::endl;
	// Before backtracing, the rank 0 receives a bunch of data from every other rank
	if (m_options.mpiRank != 0)
	{
		//Prepare structures for serialization
		std::ostringstream ss (std::ios_base::binary);
		boost::archive::binary_oarchive oa{ ss };

		//Serialize stuff to be sent
		for (int i = 0; i < m_options.threads_num; i++)
		{
			oa << ClosedList[i];
			oa << nodes_count[i];
			oa << nodes_reopen[i];
			oa << OpenList[i].size();
		}

		//Experimental lz4 compression
		u4 len = (u4) (ss.str().length() + 1);
		u4 lz4len = LZ4_compressBound((int)len);
		u4 lz4len2 = 0;
		char * lz4Data = new char[lz4len + (s_u4 * 3)]();

		lz4len2 = LZ4_compress(ss.str().c_str(), &lz4Data[s_u4 * 3], (int) len);
		((u4*)lz4Data)[0] = lz4len2 + (s_u4 * 3);
		((u4*)lz4Data)[1] = lz4len;
		((u4*)lz4Data)[2] = len;

		//Sending message
		MPI_Send(lz4Data, (int) (lz4len2 + ( s_u4 * 3 ) ), MPI_CHAR, 0, 0, MPI_COMM_WORLD);

		//Clearing buffer
		delete[] lz4Data;
	}
	else
	{
		int i, sender, n_bytes = 0;
		MPI_Status status;
		//Save stuff from node 0 into final structures
		for (i = 0; i < m_options.threads_num; i++)
		{
			nodes_countFinal[i] = nodes_count[i];
			nodes_reopenFinal[i] = nodes_reopen[i];
			ClosedListFinal[i] = ClosedList[i];
			nodes_openListSizeFinal[i] = OpenList[i].size();
		}

		//Receive remote stuff
		for (i = 1; i < m_options.mpiCommSize; i++)
		{
			//Probe and receive final message from each rank
			MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, MPI_CHAR, &n_bytes);
			sender = status.MPI_SOURCE;

			char * buffer = new char[n_bytes]();
			MPI_Recv(buffer, n_bytes, MPI_CHAR, sender, 0, MPI_COMM_WORLD, &status);

			//Experimental LZ4 decompression
			u4 lz4len = ((u4*)buffer)[1];
			u4 len    = ((u4*)buffer)[2];
			char * tempBuff = new char[lz4len]();
			LZ4_decompress_fast(&buffer[s_u4 * 3], tempBuff, len);
			delete[] buffer;

			//Unserialize data into place
			std::istringstream ss(std::string(tempBuff, tempBuff + len), std::ios_base::binary);
			boost::archive::binary_iarchive ia{ ss };

			//Calculate offset of nodes
			int offset = sender * m_options.threads_num;

			//Recover stuff in same order as saved
			for (int j = 0; j < m_options.threads_num; j++)
			{
				ia & ClosedListFinal[offset + j];
				ia & nodes_countFinal[offset + j];
				ia & nodes_reopenFinal[offset + j];
				ia & nodes_openListSizeFinal[offset + j];
			}

			//Free buffer
			delete[] tempBuff;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

/*!
 * Same a_star() function usage.
 * Starting function to do a pa_star search.
 */
template < int N >
int PAStar<N>::pa_star(const Node<N> &node_zero, const Coord<N> &coord_final, const PAStarOpt &options)
{





    //std::cout << options.mpiRank << ": fase 4" << std::endl;
    if (options.threads_num <= 0)
        throw std::invalid_argument("Invalid number of threads");
    Coord<N>::configure_hash(options.hash_type, options.hash_shift);

    //std::cout << options.mpiRank << ": fase 5" << std::endl;


    PAStar<N> pastar_instance(node_zero, options);

    //std::cout << options.mpiRank << ": fase 6" << std::endl;

    std::vector<std::thread> threads;
    TimeCounter *t = new TimeCounter("Phase 2: PA-Star running time: ");

    // Create worker threads
    for (int i = 0; i < options.threads_num; ++i)
	// Create worker threads
        threads.push_back(std::thread(&PAStar::worker, &pastar_instance, i, coord_final));

	// Create sender thread
	threads.push_back(std::thread(&PAStar::sender, &pastar_instance));

	// Create receiver thread
	threads.push_back(std::thread(&PAStar::receiver, &pastar_instance, &pastar_instance));

    // Wait for the end of all threads
    for (auto& th : threads)
        th.join();
    delete t;
	
    std::cout << options.mpiRank << ": fase 7" << std::endl;

    // Don't you dare removing that barrier
    MPI_Barrier(MPI_COMM_WORLD);

    pastar_instance.sync_pastar_data();
    //std::cout << "preparing to print answer " << rank << std::endl;

    std::cout << options.mpiRank << ": fase 8" << std::endl;


    if (options.mpiRank == 0)
        pastar_instance.print_answer();

	//std::cout << options.mpiRank << ": waiting for other processes to finish" << std::endl;
    std::cout << options.mpiRank << ": fase 9" << std::endl;



    return 0;
}

#define PASTAR_DECLARE_TEMPLATE( X ) \
template class PAStar< X >; \

MAX_NUM_SEQ_HELPER(PASTAR_DECLARE_TEMPLATE);
