/*!
 * \authors {Daniel Sundfeld, Gabriel Ferreira}
 * \copyright MIT License
 */

//Included in pastar.cpp


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

                send_queue[i] = ptr; //prealocate first to reduce time blocked

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

//! Flush sender thread message qeueu
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
