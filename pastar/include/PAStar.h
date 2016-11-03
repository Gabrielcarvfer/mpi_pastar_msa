#ifndef _PSTAR_H
#define _PSTAR_H
/*!
 * \class PAStar
 * \author Daniel Sundfeld
 * \copyright MIT License
 *
 * \brief Do a multiple sequence alignment reducing the search space
 * with parallel a-star algorithm
 */

#ifndef WIN32
#define linux //mac for thread port
#endif

#include <atomic>
#include <condition_variable>
#include <iostream>
#include <map>
#include <string>
#include <thread>
#include <vector>
#include <tuple>
#include <queue>
#include "AStar.h"
#include "Coord.h"
#include "CoordHash.h"
#include "Node.h"
#include "PriorityList.h"

#include <sstream>
// include input and output archivers

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

// include this header to serialize vectors, maps and other types
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/collection_size_type.hpp>
#include <boost/serialization/list.hpp>



#ifndef THREADS_NUM
    #define THREADS_NUM std::thread::hardware_concurrency()
#endif

#ifndef WIN32

#ifdef MAC
    //http://yyshen.github.io/2015/01/18/binding_threads_to_cores_osx.html
    //https://gist.github.com/Coneko/4234842
    #define SYSCTL_CORE_COUNT   "machdep.cpu.core_count"
    #import <mach/thread_act.h>
    #include <mach/mach_types.h>
    #include <pthread.h>


    typedef struct cpu_set {
        uint32_t    count;
    } cpu_set_t;

    static inline void
    CPU_ZERO(cpu_set_t *cs) { cs->count = 0; }

    static inline void
    CPU_SET(int num, cpu_set_t *cs) { cs->count |= (1 << num); }

    static inline int
    CPU_ISSET(int num, cpu_set_t *cs) { return (cs->count & (1 << num)); }

    kern_return_t	thread_policy_set(thread_t                      thread,
                                      thread_policy_flavor_t		flavor,
                                      thread_policy_t			policy_info,
                                      mach_msg_type_number_t		count);


    int sched_setaffinity(pthread_t thread, size_t cpu_size, cpu_set_t *cpu_set);
#endif

#endif

/*!
 * \brief Arguments for PAStar class
 */
struct PAStarOpt {
    AStarOpt common_options;
    hashType hash_type;
    int hash_shift = 0;
    int threads_num = 0;

    int mpiRank = 0;
    int mpiCommSize = 0;
    int mpiMin = 0;
    int mpiMax = 0;
    int totalThreads =0;

    PAStarOpt()
    {
        hash_type = HashFZorder;
        hash_shift = HASH_SHIFT;
        threads_num = THREADS_NUM;
    }
    PAStarOpt(AStarOpt &common, hashType type, int shift, int th)
    {
        common_options = common;
        hash_type = type;
        hash_shift = shift;
        threads_num = th;
    }
};

template < int N >
class PAStar {
    public:
        static int pa_star(const Node<N> &node_zero, const Coord<N> &coord_final, const PAStarOpt &options);

        std::vector<std::thread> proc_threads;

    private:
        // Members
        const PAStarOpt m_options;
        PriorityList<N> *OpenList;
        std::map< Coord<N>, Node<N> > *ClosedList;

        //Accountability data
        long long int *nodes_count;
        long long int *nodes_reopen;

        long long int *nodes_countFinal;
        long long int *nodes_reopenFinal;

        long long int * nodes_openListSizeFinal;
        long long int * nodes_closedListSizeFinal;

        //Control access variables
        std::mutex *queue_mutex;
        std::condition_variable *queue_condition;
        std::vector< Node<N> > *queue_nodes;
        std::atomic<bool> end_cond;
		std::atomic<bool> end_condLocal;
		std::mutex * check_stop_mutex;

		//Structure that hold
		std::vector< Node<N> > **send_queue; //In that way, we can have a vector of vectors of nodes

		//Structures for syncing final node and other data
        std::mutex final_node_mutex;
        Node<N> final_node;
        
        std::atomic<int> final_node_count;
	
		bool remoteLowerVal;

		//Synchronization variables
        std::mutex sync_mutex;
        std::atomic<int> sync_count;
        std::condition_variable sync_condition;

		//Sender mutex and condition variables
		std::condition_variable sender_condition;
		std::mutex sender_mutex;
		std::mutex squeues_mutex;
		std::mutex *squeue_mutex;
		bool sender_empty;
		int * threadLookupTable;

        //Receiver mutex and condition variables
        std::condition_variable receiver_condition;
        std::mutex receiver_mutex;
        std::condition_variable *processing_condition;
        std::vector<char*> * processing_queue;
        std::mutex * processing_mutexes;
        bool recv_goodbye = false;
        long long int recv_cnt = 0;

		//Node pairwise costs variables
#ifdef WIN32
		int ** pairwise_costs;
#endif


        // Constructor
        PAStar(const Node<N> &node_zero, const PAStarOpt &opt);
        ~PAStar();

        // Misc functions
#ifndef WIN32
        int set_affinity(int tid);
#endif
        void sync_threads(bool flushing);
		void sync_threads_local();
        void print_nodes_count();

        // Queue functions
        void enqueue(int tid, std::vector< Node<N> > &nodes);
        void consume_queue(int tid);
        void wait_queue(int tid);
        void wake_all_queue();

        // End functions
        void process_final_node(int tid, const Node<N> &n);
        bool check_stop(int tid);
		bool check_stop_global(int tid);

        // Worker Functions
        void worker_inner(int tid, const Coord<N> &coord_final);
        int worker(int tid, const Coord<N> &coord_final);

        // Receiver Function
        int receiver(PAStar<N> * pastar_inst);
		int sender(void);
		int process_message(int tid);
		std::mutex processing_mutex;
		std::mutex sync_mutex_global;
		std::condition_variable sync_condition_global;
		void flush_sender();
		void flush_receiver();
		bool end_of_transmission;

        // Backtrack
        void sync_pastar_data(void);
        void print_answer(void);


};
#endif
