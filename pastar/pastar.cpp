/*!
 * \authors {Daniel Sundfeld, Gabriel Ferreira}
 * \copyright MIT License
 */


#ifndef WIN32
#include <sched.h>
#endif

// Defining Ryzen and running with 2 MPI processes (8 thread each) vs a single 16-thread process
// the Kinase alignment of kinase file drops from 33 seconds to 22. CCX latency related.
//#define RYZEN 

#include <chrono>
#include <mpi.h>
 
#include "include/PAStar.h"
#include "include/backtrace.h"
#include "include/Coord.h"
#include "include/Node.h"
#include "include/TimeCounter.h"
#include "include/lz4sup.h"
#include "include/Sequences.h"




#define MPI_TAG_SEND_COMMON                0
#define MPI_TAG_REQ_CHECK         0x00FFFFFD
#define MPI_TAG_ACC_REQ           0x00FFFFFE
#define MPI_TAG_KILL_RECEIVER     0x00FFFFFF

#define BACKTRACE_WORKING                  0
#define BACKTRACE_LISTENING                1
#define BACKTRACE_NEXTNODE                30
#define BACKTRACE_PARTIAL_ALIGNMENT       31
#define BACKTRACE_END                     32


#ifndef WIN32

#ifdef MAC
 int sched_setaffinity(pthread_t thread, size_t cpu_size, cpu_set_t *cpu_set)
 {
    thread_port_t mach_thread;
    int core = 0;

    for (core = 0; core < 8 * cpu_size; core++) 
    {
        if (CPU_ISSET(core, cpu_set)) 
            break;
    }

    thread_affinity_policy_data_t policy = { core };
    mach_thread = pthread_mach_thread_np(thread);
    thread_policy_set(mach_thread, THREAD_AFFINITY_POLICY,(thread_policy_t)&policy, 1);
    return 0;
}
#endif

#endif

template < int N >
PAStar<N>::PAStar(const Node<N> &node_zero, const struct PAStarOpt &opt)
: m_options(opt)
{
    if (m_options.mpiRank == 0)
    {
        std::cout << "Running PAStar with: "
        << opt.totalThreads << " threads ("
        << opt.mpiCommSize << " machines with " 
        << opt.threads_num << " threads each)," 
        << Coord<N>::get_hash_name() << " hash, "
        << Coord<N>::get_hash_shift() << " shift.\n";
    }


    //Value initialization
    end_cond = false;
    end_condLocal = false;
    sync_count = 0;
    final_node.set_max();
    recv_cnt = 0;

    nodes_countFinal = NULL;
    nodes_reopenFinal = NULL;

    //Main structures and accountability
    OpenList = new PriorityList<N>[m_options.threads_num]();
    ClosedList = new std::map< Coord<N>, Node<N> >[m_options.threads_num]();

    nodes_count = new long long int[m_options.threads_num]();
    nodes_reopen = new long long int[m_options.threads_num]();


	// The additional mutex and condition are used by sender thread and local ones
    queue_mutex = new std::mutex[m_options.threads_num];
    queue_condition = new std::condition_variable[m_options.threads_num+1];
    queue_nodes = new std::vector< Node<N> >[m_options.threads_num];
    squeue_mutex = new std::mutex[m_options.totalThreads];
    send_queue = new std::vector< Node<N> >*[m_options.totalThreads+1]();

    // Lookup table to check global thread tid info

    threadLookupTable = new int[m_options.totalThreads*2]();
    int i = 0, j = 0, k = 0;
    for (; k < m_options.mpiCommSize; k++)
    {
      for (j = 0; j < m_options.threads_num; j++)
      {
         threadLookupTable[i] = k;
         threadLookupTable[m_options.totalThreads + i] = j;
         i++;
     }
 }

    // The additional mutex and condition are used by sender thread
 queue_mutex = new std::mutex[m_options.threads_num];
 queue_condition = new std::condition_variable[m_options.threads_num+1];
 queue_nodes = new std::vector< Node<N> >[m_options.threads_num];
 squeue_mutex = new std::mutex[m_options.totalThreads];

 send_queue = new std::vector< Node<N> >*[m_options.totalThreads+1]();
 for (i = 0; i < m_options.totalThreads+1; i++)
 {
    send_queue[i] = new std::vector< Node <N>>();
}

    //Preallocating pairwise_costs structure
#ifdef WIN32
pairwise_costs = new int*[m_options.threads_num];
for (i = 0; i < m_options.threads_num; i++)
    pairwise_costs[i] = new int[N * 2 * 3];

#endif

    //Initializing message processing threads and their stuff
processing_condition = new std::condition_variable[m_options.threads_num];
processing_queue = new std::vector<char*>[m_options.threads_num];
processing_mutexes = new std::mutex[m_options.threads_num];

    // Allocate final structures and enqueue first node if rank zero
if (m_options.mpiRank == 0)
{
        //ClosedListFinal = new std::map< Coord<N>, Node<N> >[m_options.totalThreads];
    nodes_countFinal = new long long int[m_options.totalThreads]();
    nodes_reopenFinal = new long long int[m_options.totalThreads]();
    nodes_closedListSizeFinal = new long long int [m_options.totalThreads]();
    nodes_openListSizeFinal = new long long int[m_options.totalThreads]();

    OpenList[0].enqueue(node_zero);
}
}

template < int N >
PAStar<N>::~PAStar()
{
    //Accountability
    delete[] nodes_count;
    delete[] nodes_reopen;

    //Local Queue stuff
    delete[] queue_mutex;
    delete[] queue_condition;
    delete[] queue_nodes;

    //Remote Qeueu stuff
    delete[] squeue_mutex;
    delete threadLookupTable;


    //Freeing sender queues
    for (int i = 0; i < m_options.totalThreads+1; i++)
    {
        delete send_queue[i];
    }
    delete[] send_queue;

    //Freeing pairwise costs structure
#ifdef WIN32
    for (int i = 0; i < m_options.threads_num; i++)
        delete[] pairwise_costs[i];
    delete[] pairwise_costs;
#endif


    delete[] ClosedList;
    if (m_options.mpiRank == 0)
    {

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
    typename std::map< Coord<N>, Node<N> >::iterator c_search;

    for (typename std::vector< Node<N> >::iterator it = nodes.begin() ; it != nodes.end(); ++it)
    {
        if ((c_search = ClosedList[tid].find(it->pos)) != ClosedList[tid].end())
        {
            if (it->get_g() >= c_search->second.get_g())
                continue;
            ClosedList[tid].erase(it->pos);
            nodes_reopen[tid] += 1;
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
    std::unique_lock<std::mutex> queue_lock(queue_mutex[tid]);
    std::vector< Node<N> > nodes_to_expand(queue_nodes[tid]);
    queue_nodes[tid].clear();
    queue_lock.unlock();

    enqueue(tid, nodes_to_expand);
    return;
}

//! Wait something on the queue
template < int N >
void PAStar<N>::wait_queue(int tid)
{
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
    std::unique_lock<std::mutex> lock (squeues_mutex);
    if (!sender_empty)
    {
        lock.unlock();
        std::unique_lock<std::mutex> sender_lock(sender_mutex);
		queue_condition[m_options.threads_num].notify_one(); //acorda sender se estiver dormindo
		sender_condition.wait(sender_lock);
	}
	else
	{
        lock.unlock();
    }
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
		receiver_condition.wait(receiver_lock);
	}
	return;
}

//! Sync all local threads
template < int N >
void PAStar<N>::sync_threads_local()
{
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

    std::vector< Node<N> > *neigh = new std::vector< Node<N> >[m_options.totalThreads];
    int actualTid = tid + m_options.mpiMin;

    // Loop ended by process_final_node
    while (end_condLocal == false)
    {
        typename std::map< Coord<N>, Node<N> >::iterator c_search;

		// Consume openlist to find if someone have a better route
        consume_queue(tid);

        // Dequeue phase
        if (OpenList[tid].dequeue(current) == false)
        {
            wait_queue(tid);
            continue;
        }
        nodes_count[tid] += 1;

        // Check if better node was already found
        if ( (c_search = ClosedList[tid].find(current.pos) ) != ClosedList[tid].end())
        {
            if (current.get_g() >= c_search->second.get_g())
                continue;
            nodes_reopen[tid] += 1;
        }

        ClosedList[tid][current.pos] = current;

        if (current.pos == coord_final)
        {
            process_final_node(tid, current);
            continue;
        }

        // Expand phase
#ifndef WIN32
        current.getNeigh(neigh, m_options.totalThreads);
#else
        current.getNeigh(neigh, m_options.totalThreads, pairwise_costs[tid]);
#endif
        // Reconciliation phase
        for (int i = 0; i < m_options.totalThreads; i++)
        {
            // process local generated nodes
            if (i == actualTid)
            {
                enqueue(tid, neigh[i]);
            }
            //enqueue for other nodes
            else if (neigh[i].size() != 0)
            {
                //wich can be local
                if (i >= m_options.mpiMin && i < m_options.mpiMax)
                {
                   int iLocal = threadLookupTable[i + m_options.totalThreads];

                   std::unique_lock<std::mutex> queue_lock(queue_mutex[iLocal]);
                   queue_nodes[iLocal].insert(queue_nodes[iLocal].end(), neigh[i].begin(), neigh[i].end());
                   queue_lock.unlock();
                   queue_condition[iLocal].notify_one();
               }
                //or remote
               else
               {
                   std::unique_lock<std::mutex> queue_lock(squeue_mutex[i]);
                   send_queue[i]->insert(send_queue[i]->end(), neigh[i].begin(), neigh[i].end());
                   queue_lock.unlock();
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
   int originRank = threadLookupTable[n.pos.get_id(m_options.totalThreads)];

	//Lock all threads from all nodes
   std::unique_lock<std::mutex> final_node_lock(final_node_mutex);

	// Better possible answer already found, discard n
   if (final_node.get_f() < n.get_f())
   {
		//return;
   }
   else
   {
      if (n.pos.get_id(m_options.threads_num) == ((unsigned int)tid))
      {
			// Broadcast the node
         final_node = n;
         final_node_count = 0;

			//If the tid is shared with founder of the final node, add to other local nodes
         for (i = 0; i < m_options.threads_num; i++)
         {
            if (i != tid)
            {
               std::lock_guard<std::mutex> queue_lock(queue_mutex[i]);
               queue_nodes[i].push_back(n);
               queue_condition[i].notify_one();
           }
       }


			//For remote threads: send messages to be broadcasted between node threads
       for (int i = tid; i < m_options.totalThreads; i=i+m_options.threads_num)
       {
        std::lock_guard<std::mutex> lock(squeue_mutex[i]);
        send_queue[i]->push_back(n);
    }
    queue_condition[m_options.threads_num].notify_one();

    final_node_lock.unlock();
}
		// Every other worker is unlocked
else
{
			//if (n != final_node) std::cout << "BUG HERE!\n";
 final_node_lock.unlock();
}

		// This node have the highest priority between all Openlist. Broadcast the end condition
if (++final_node_count == m_options.threads_num)
{
  end_condLocal = true;
}
}
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
   long long int local, val;
   Node<N> n = final_node;

   wake_all_queue();

    //Sync all local and global threads, flushing senders and receivers from all processes
   sync_threads(true);

	// Consume openlist to find if someone have a better route
   consume_queue(tid);

	// If someone have a better value, we're not at the end yet	   -> we could look for minimum value inside node and then use MPI_Allreduce to find global minimum
   if (OpenList[tid].get_highest_priority() < final_node.get_f())
   {
    end_condLocal = false;
}

    //Sync all local threads before syncing globally
sync_threads_local();

if (tid == 0)
{
  remoteLowerVal = false;
		// If end_cond = true, then we might be in a false end where everyone agrees with its own final node;
  local = final_node.get_f();

  MPI_Allreduce(&local, &val, 1, MPI_LONG_LONG_INT, MPI_MIN, MPI_COMM_WORLD);

        //Local value isnt global minimum
  if (local != val)
  {
    end_condLocal = false;
}

		//Some MPI implementations don't like atomic booleans, so here's a possible workaround
MPI_Allreduce(&end_condLocal, &end_cond, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);

}

    //Hold all threads different than tid 0, while it communicates with other processes
sync_threads_local();

	//If someone found a better node
if (!end_cond)
{
		//If it was on a local thread of that rank
  if (!end_condLocal)
  {
            //If, and only if, out local value is worse than a remote, we clean it and enqueue again
               //Erase previous final node because ours is worse
     ClosedList[tid].erase(n.pos);

    			//Enqueue node to git it a second
     if (n.pos.get_id(m_options.threads_num) == (unsigned int) tid)
     {
        OpenList[tid].conditional_enqueue(n);
    }
}


        //If from a remote rank, just finish
return true;
}
    //If everyone agreed, then exit
return false;
}

//! Execute a worker thread. This thread have id \a tid
template < int N >
int PAStar<N>::worker(int tid, const Coord<N> &coord_final)
{

#ifndef WIN32
	#ifndef RYZEN
	set_affinity(tid);
	#else
	set_affinity(tid+m_options.mpiMin);
	#endif
#endif

	//sync all threads from all nodes to prevent problems
	sync_threads(false);

	
	// check_stop syncs and check if local nodes agreed in local optimum
	do
	{
		// worker_inner is the main inner loop
		worker_inner(tid, coord_final);
	} while (check_stop(tid));


	//global sync with flush to garantee that no node will be stuck because a sender/receiver thread didn't died at right time
	sync_threads(true);


	//the first local thread of each node kill sender and receiver to prevent problems with MPI exchange
	if (tid == 0)
	{
		// awake receiver to die
		Node<N> nullNode;
		send_queue[m_options.totalThreads]->push_back(nullNode);

		// awake sender to send message to receiver and after that kill itself
		queue_condition[m_options.threads_num].notify_one();

	}

	sync_threads_local();

    return 0;
}

//! Execute MPI sender thread, with tid = numThreads + 1
template < int N >
int PAStar<N>::sender()
{
	int i = 0, tag = 0;
	const char b[] = "1";
	bool goodbye = false;

    char ** lz4Data = NULL;
    lz4Data = (char**)calloc(1, sizeof(char*));

    std::unique_lock<std::mutex> sender_lock(squeues_mutex);

    //Preallocating vector to reduce time of locked queue
    std::vector< Node<N> >* ptr = new std::vector< Node<N> >();

    while (!goodbye)
    {
      sender_empty = true;
		//Select a thread to send all its nodes
      for (i = 0; i < m_options.totalThreads; i++)
      {
         if ( !send_queue[i]->empty() )
         {

            sender_empty = false;
            if (!squeue_mutex[i].try_lock())
            {
               continue;
           }


				//Copy to auxiliary buffer
           std::vector<Node<N>>* temp = send_queue[i];

                send_queue[i] = ptr; //prealocate firstm to reduce time blocked

				//Unlock send_queue for continuous operation
                squeue_mutex[i].unlock();

                ptr = new std::vector< Node <N> >();

				//Sending nodes to remote target
                std::ostringstream ss (std::ios_base::binary);
                boost::archive::binary_oarchive oa{ ss };

                oa & *temp;

                delete temp;

                int lz4Size = pastar_lz4_en( ss.str().c_str(), lz4Data, ss.str().length() );

                MPI_Send(*lz4Data, lz4Size, MPI_CHAR, threadLookupTable[i], threadLookupTable[i+m_options.totalThreads], MPI_COMM_WORLD);

                delete[] *lz4Data;
            }
        }

		//If some node in killing queue, finish sender and receiver
        if (!send_queue[m_options.totalThreads]->empty())
        {
         if (sender_empty)
         {
         			//sending finishing message
            MPI_Send((void*)&b, 2, MPI_CHAR, m_options.mpiRank, MPI_TAG_KILL_RECEIVER, MPI_COMM_WORLD);
            goodbye = true;
            send_queue[m_options.totalThreads]->clear();
        }
        else
        {
				//received finishing node, but still have stuff to process"
            continue;
        }
    }

		// If set goodbye, continue to finish
    if (goodbye)
    {
			//Finishing
     sender_condition.notify_one();
 }
 else
 {
			// Sinaliza buffer limpo e acorda quem estiver esperando
     if (sender_empty)
     {
    			//Sender don't have working, going to sleep
         sender_condition.notify_one();
         queue_condition[m_options.threads_num].wait(sender_lock);
     }
 }
 std::this_thread::sleep_for(std::chrono::milliseconds(20));
}
delete ptr;
free(lz4Data);
	//At the end, unlock the lock
sender_lock.unlock();

	// Last notify to unlock flusher
sender_condition.notify_all();
return 0;
}

//! Execute a receiver thread, with tid = numThreads + 2
template < int N >
int PAStar<N>::receiver(PAStar<N> * pastar_inst)
{
    //Initializing MPI and working variables
    MPI_Status status;
    int n_bytes = 0;
    int sender = 0;
    int sender_tag = 0;
    int i = 0, flag = 0;
    char * buffer = NULL;
    int prob = 1;

    while (!recv_goodbye)
    {
        //Check if any message was received
      flag = 0;

      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);

        //Make receiver sleep to reduce prevent CPU usage at probing
      if (!flag)
      {
        std::this_thread::sleep_for(std::chrono::milliseconds(prob = prob << 1));
        std::this_thread::yield();
        continue;
    }
        //Reset clock "skew"
    prob = 1;

        //Get info from message received
    MPI_Get_count(&status, MPI_CHAR, &n_bytes);
    sender = status.MPI_SOURCE;
    sender_tag = status.MPI_TAG;

		//Receive message
    if (n_bytes == 0)
    {
     MPI_Recv(NULL, 0, MPI_CHAR, sender, sender_tag, MPI_COMM_WORLD, &status);
     continue;
 }

 buffer = new char[n_bytes]();

		//Receive thread destination plus nodes
 MPI_Recv(buffer, n_bytes, MPI_CHAR, sender, sender_tag, MPI_COMM_WORLD, &status);

 switch (sender_tag)
 {
    case MPI_TAG_KILL_RECEIVER:
                //end receiver
    recv_goodbye = true;
    			//as safety measure, wake up receiver to finish
    queue_condition[m_options.threads_num].notify_one();
    break;


            // the common case is to receive only nodes
    default:
    if (n_bytes == ((int*)buffer)[0])
    {
        			//Lock counter of messages awaiting to be processed
     std::unique_lock<std::mutex> lock(processing_mutex);
     recv_cnt++;
     lock.unlock();

                    //Enqueue buffer pointer to be processed later
     std::unique_lock<std::mutex> proclock (processing_mutexes[sender_tag]);
     processing_queue[sender_tag].push_back(buffer);

                    //Unlock queue and notify
     proclock.unlock();

     processing_condition[sender_tag].notify_one();
 }
 else
 {
    delete[] buffer;
}
break;
}
}

    // Wake processers then join them
for (i = 0; i < m_options.threads_num; i++)
{
    processing_condition[i].notify_one();
    proc_threads[i].join();
}

    //Clean allocated stuff
delete[] processing_queue;
delete[] processing_condition;
delete[] processing_mutexes;

	// Last notify to unlock flushers
receiver_condition.notify_one();
return 0;
}

//! Execute a receiver thread, with tid = numThreads + 2
template < int N >
//int PAStar<N>::process_message(int sender_tag, char *buffer)
int PAStar<N>::process_message(int tid)
{

#ifndef WIN32
    //set_affinity(tid);
#endif

    char ** buffer = NULL;
    buffer = (char**)calloc(1,sizeof(char*));

    while(!recv_goodbye)
    {
        //Check if there are messages to process
        std::unique_lock<std::mutex> lock (processing_mutexes[tid]);

        if (processing_queue[tid].empty())
        {
            processing_condition[tid].wait(lock);
            lock.unlock();
            continue;
        }

        //Copy pointer and unlock queue, permitting continuous operation
        *buffer = processing_queue[tid].back();
        processing_queue[tid].pop_back();
        lock.unlock();


        int len = pastar_lz4_dec(buffer);

    	//Prepare to deserialize stuff into place
        std::istringstream ss(std::string(*buffer, *buffer + len), std::ios_base::binary);

    	//convert buffer to something useful, like thread id and node vector
        boost::archive::binary_iarchive ia{ ss };
        std::vector<Node<N>> temp;
        ia & temp;
        delete[] *buffer;

    	// Locking queue to add received stuff
        std::unique_lock<std::mutex> queue_lock(queue_mutex[tid]);
        queue_nodes[tid].insert(queue_nodes[tid].end(), temp.begin(), temp.end());//temp.rbegin(), temp.rend());//temp.begin(), temp.end());
        queue_condition[tid].notify_one();
        queue_lock.unlock();


    	// Update counter of messages awaiting to be processed
        std::unique_lock<std::mutex> processing_lock(processing_mutex);
        recv_cnt--;
        processing_lock.unlock();

        // Notify any process awaiting for receiver queue to be clean
        if (recv_cnt == 0)
        {
          receiver_condition.notify_one();
      }
  }

    //We're finishing, so we need to clean up any non processed messages
  for (unsigned int i = 0; i < processing_queue[tid].size(); i++)
  {
    delete[] processing_queue[tid][i];
}
free(buffer);
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
    int i;
    
    //Distributed backtrace
    std::list<char> finalAlignments[N];

    TimeCounter *t = NULL;

    if (m_options.mpiRank == 0)
    {
        t = new TimeCounter("Phase 3 - backtrace: ");
    }
    

    Sequences *seq = Sequences::getInstance();

    int id = seq->get_final_coord<N>().get_id(m_options.totalThreads);

    Node<N> current;

    int state = BACKTRACE_LISTENING;

    //Whoever is the owner of last node start backtracing
    if ( (id >= m_options.mpiMin) && (id < m_options.mpiMax) )
    {
        current = ClosedList[id - m_options.mpiMin][seq->get_final_coord<N>()];
        state = BACKTRACE_WORKING;
        std::cout << "Final Score: " << current << std::endl;
    }


    bool finished = false;

    do
    {
        switch(state)
        {
            case BACKTRACE_WORKING:
                {

                    std::list<char> alignments[N];
                    Node<N> node;

                        //Execute partial alignment
                    node = partial_backtrace_alignment<N>(alignments, ClosedList, current, m_options.totalThreads, m_options.mpiMin, m_options.mpiMax);

                        //Transmit partial alignment
                    if (m_options.mpiRank != 0)
                    {
                     std::ostringstream ss (std::ios_base::binary);
                     boost::archive::binary_oarchive oa{ ss };

                     for (i = 0; i < N; i++)
                        oa & alignments[i];
                    MPI_Send(ss.str().c_str(), ss.str().length(), MPI_CHAR, 0, BACKTRACE_PARTIAL_ALIGNMENT, MPI_COMM_WORLD);
                }
                        //Or save to final aligment space
                else
                {
                    for (i = 0; i < N; i++)
                        finalAlignments[i].splice(finalAlignments[i].begin(), alignments[i]);
                }

                        //Check if backtrace has finished
                if (node.pos == Sequences::get_initial_coord<N>())
                {
                            //Send finishing message
                    for (i = 0; i < m_options.mpiCommSize; i++)
                        if (i != m_options.mpiRank)
                            MPI_Send(NULL, 0, MPI_CHAR, i, BACKTRACE_END, MPI_COMM_WORLD);
                        finished = true;
                    }
                    else
                    {
                            //Serialize current node
                        std::ostringstream ss (std::ios_base::binary);
                        boost::archive::binary_oarchive oa{ ss };

                        oa & node;

                            //Transmit data
                        MPI_Send(ss.str().c_str(), ss.str().length(), MPI_CHAR, threadLookupTable[node.get_parent().get_id(m_options.totalThreads)], BACKTRACE_NEXTNODE, MPI_COMM_WORLD);

                        state = BACKTRACE_LISTENING;
                    }
                }
                break;
            case BACKTRACE_LISTENING:
                {
                    MPI_Status status;
                    int n_bytes = 0;
                    int sender = 0;
                    int sender_tag = 0;

                    //Listen to messages
                    int flag = 0;
                    int tag = 0;


                    //Check if there is a partial alignment message awaiting
                    MPI_Iprobe(MPI_ANY_SOURCE, BACKTRACE_PARTIAL_ALIGNMENT, MPI_COMM_WORLD, &flag, &status);

                    if (!flag)
                    {
                        //If didnt received partial alignment, check for next_node
                        MPI_Iprobe(MPI_ANY_SOURCE, BACKTRACE_NEXTNODE, MPI_COMM_WORLD, &flag, &status);

                        if (!flag)
                        {
                            //If didnt received next_node, check for final node
                            MPI_Iprobe(MPI_ANY_SOURCE, BACKTRACE_END, MPI_COMM_WORLD, &flag, &status);

                            if (!flag)
                            {
                                //if didn't received anything, loop
                                std::this_thread::sleep_for(std::chrono::milliseconds(40));
                                break;
                            }
                        }
                    }

                    //Iprobe messages tagged with append, that should be processed earlier than aligning

                    MPI_Get_count(&status, MPI_CHAR, &n_bytes);
                    sender = status.MPI_SOURCE;
                    sender_tag = status.MPI_TAG;

                    char * buffer = (char*)calloc(sizeof(char),n_bytes);

                    MPI_Recv(buffer, n_bytes, MPI_CHAR, sender, sender_tag, MPI_COMM_WORLD, &status);

                    std::string temp = std::string(buffer, buffer + n_bytes);

                    free(buffer);

                    //Switch for message tag
                    switch(sender_tag)
                    {
                        //If received a next_node, find its previous node and then resume partial backtrace
                        case BACKTRACE_NEXTNODE:
                            {
                                std::istringstream ss(temp, std::ios_base::binary);

                                //convert buffer to something useful, like thread id and node vector
                                boost::archive::binary_iarchive ia{ ss };
                                ia & current;
                                int index = threadLookupTable[current.get_parent().get_id(m_options.totalThreads)+m_options.totalThreads];
                                current = ClosedList[index][current.get_parent()];

                                //Set state to working
                                state = BACKTRACE_WORKING;
                            }
                            break;
                        //Append partial alignment and then continue listening
                        case BACKTRACE_PARTIAL_ALIGNMENT:
                            {
                                std::list<char> partial_alignments[N];
                                std::istringstream ss(temp, std::ios_base::binary);

                                //convert buffer to something useful, like thread id and node vector
                                boost::archive::binary_iarchive ia{ ss };

                                for (i = 0; i < N; i++)
                                {
                                    ia & partial_alignments[i];
                                    finalAlignments[i].splice(finalAlignments[i].begin(), partial_alignments[i]);

                                }
                            }
                            break;
                        case BACKTRACE_END:
                            finished = true;
                            break;
                        default:
                            break;
                    }
            }
            break;
            default:
            break;
        }
    }while(!finished);

    //Every mpi rank is synchronized
    MPI_Barrier(MPI_COMM_WORLD);

    //After backtrace, rank 0 prints alignment stuff
    if (m_options.mpiRank == 0)
    {
        delete t;
        print_entire_backtrace<N>(finalAlignments);
        print_nodes_count();
    }
}

/*!
 * To calculate the correct statistics and do the
 * backtrace, we need to merge work data in rank 0 node.
 */
template < int N >
 void PAStar<N>::sync_pastar_data()
 {
    int i, j;

	// Before backtracing, the rank 0 receives a bunch of data from every other rank
    if (m_options.mpiRank != 0)
    {
		//Prepare structures for serialization
      std::ostringstream ss (std::ios_base::binary);
      boost::archive::binary_oarchive oa{ ss };

      long long int temp = 0;
		//Serialize stuff to be sent
      for (i = 0; i < m_options.threads_num; i++)
      {
        temp = ClosedList[i].size();
        oa & temp;
        oa & nodes_count[i];
        oa & nodes_reopen[i];
			//Workaround for linux
        temp = OpenList[i].size();
        oa & temp;
    }

		//Freeing memory as soon as possible
    delete[] OpenList;


    char ** lz4Data = NULL;
    lz4Data = (char**)calloc(1, sizeof(char*));

    int lz4Size = pastar_lz4_en( ss.str().c_str(), lz4Data, ss.str().length() );

		//Sending message
    MPI_Send(*lz4Data, lz4Size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);

		//Clearing buffer
    delete[] *lz4Data;
    free(lz4Data);

}
else
{
  int sender, n_bytes = 0;
  MPI_Status status;
		//Save stuff from node 0 into final structures
  for (i = 0; i < m_options.threads_num; i++)
  {
     nodes_countFinal[i] = nodes_count[i];
     nodes_reopenFinal[i] = nodes_reopen[i];
     nodes_closedListSizeFinal[i] = ClosedList[i].size();
     nodes_openListSizeFinal[i] = OpenList[i].size();
 }


		//Freeing memory as soon as possible
 delete[] OpenList;

 char ** buffer = NULL;
 buffer = (char**)calloc(1,sizeof(char*));

		//Receive remote stuff
 for (i = 1; i < m_options.mpiCommSize; i++)
 {
			//Probe and receive final message from each rank
     MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
     MPI_Get_count(&status, MPI_CHAR, &n_bytes);
     sender = status.MPI_SOURCE;

     *buffer = new char[n_bytes]();

     MPI_Recv(*buffer, n_bytes, MPI_CHAR, sender, 0, MPI_COMM_WORLD, &status);

     int len = pastar_lz4_dec(buffer);

			//Unserialize data into place
     std::istringstream ss(std::string(*buffer, *buffer + len), std::ios_base::binary);

			//Free buffer
     delete[] *buffer;
     boost::archive::binary_iarchive ia{ ss };

			//Calculate offset of nodes
     int offset = sender * m_options.threads_num;

			//Recover stuff in same order as saved
     for (j = 0; j < m_options.threads_num; j++)
     {
        ia & nodes_closedListSizeFinal[offset + j];
        ia & nodes_countFinal[offset + j];
        ia & nodes_reopenFinal[offset + j];
        ia & nodes_openListSizeFinal[offset + j];
    }
}
free(buffer);
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
    if (options.threads_num <= 0)
        throw std::invalid_argument("Invalid number of threads");

    //Configure hash
    Coord<N>::configure_hash(options.hash_type, options.hash_shift);

    //Initiate PAStar
    PAStar<N> pastar_instance(node_zero, options);

    //Allocate a vector to hold threads
    std::vector<std::thread> threads;
    TimeCounter *t = new TimeCounter("Phase 2: PA-Star running time: ");

    // Create processing threads
    for (int i = 0; i < options.threads_num; i++)
        pastar_instance.proc_threads.push_back(std::thread(&PAStar::process_message, &pastar_instance, i));

    // Create receiver thread
    threads.push_back(std::thread(&PAStar::receiver, &pastar_instance, &pastar_instance));

    // Create worker threads
    for (int i = 0; i < options.threads_num; ++i)
        threads.push_back(std::thread(&PAStar::worker, &pastar_instance, i, coord_final));

	// Create sender thread
    threads.push_back(std::thread(&PAStar::sender, &pastar_instance));



    // Wait for the end of all threads
    for (auto& th : threads)
        th.join();
    delete t;

    // Don't you dare removing that barrier
    MPI_Barrier(MPI_COMM_WORLD);

    // Sync metrics data
    pastar_instance.sync_pastar_data();

    // Distributed backtrace than print
    pastar_instance.print_answer();

    return 0;
}

#define PASTAR_DECLARE_TEMPLATE( X ) \
template class PAStar< X >; \

MAX_NUM_SEQ_HELPER(PASTAR_DECLARE_TEMPLATE);
