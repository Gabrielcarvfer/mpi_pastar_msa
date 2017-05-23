/*!
 * \authors {Daniel Sundfeld, Gabriel Ferreira}
 * \copyright MIT License
 */

//Included in pastar.cpp


//! Execute a message processing thread
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
