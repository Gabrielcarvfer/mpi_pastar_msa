/*!
 * \authors {Daniel Sundfeld, Gabriel Ferreira}
 * \copyright MIT License
 */

//Defines
#define MPI_TAG_SEND_COMMON                0
#define MPI_TAG_REQ_CHECK         0x00FFFFFD
#define MPI_TAG_ACC_REQ           0x00FFFFFE
#define MPI_TAG_KILL_RECEIVER     0x00FFFFFF



//Includes
#ifndef WIN32
#include <sched.h>
#endif

#include <chrono>
#include <mpi.h>

#include "include/PAStar.h"
#include "include/backtrace.h"
#include "include/Coord.h"
#include "include/Node.h"
#include "include/TimeCounter.h"
#include "include/lz4sup.h"
#include "include/Sequences.h"


//Inclusion of some of MPI-PAStar functions and threads
#include "pastar_functions/PAStarSender.cpp"
#include "pastar_functions/PAStarReceiver.cpp"
#include "pastar_functions/PAStarMessageProcesser.cpp"
#include "pastar_functions/PAStarSyncData.cpp"
#include "pastar_functions/PAStarDistributedBacktrace.cpp"





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
	set_affinity(tid);
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
    pastar_instance.distributed_backtrace_n_print();

    return 0;
}

#define PASTAR_DECLARE_TEMPLATE( X ) \
template class PAStar< X >; \

MAX_NUM_SEQ_HELPER(PASTAR_DECLARE_TEMPLATE);
