/*!
 * \authors {Daniel Sundfeld, Gabriel Ferreira}
 * \copyright MIT License
 */

//Included in pastar.cpp


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