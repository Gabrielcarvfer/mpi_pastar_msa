/*!
 * \authors {Daniel Sundfeld, Gabriel Ferreira}
 * \copyright MIT License
 */

//Included in pastar.cpp


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