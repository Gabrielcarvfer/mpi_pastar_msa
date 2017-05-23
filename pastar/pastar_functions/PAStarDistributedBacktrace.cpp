/*!
 * \authors {Daniel Sundfeld, Gabriel Ferreira}
 * \copyright MIT License
 */

//Included in pastar.cpp


#define BACKTRACE_WORKING                  0
#define BACKTRACE_LISTENING                1
#define BACKTRACE_NEXTNODE                30
#define BACKTRACE_PARTIAL_ALIGNMENT       31
#define BACKTRACE_END                     32

/*!
*  Execute distributed backtrace and print results
*/
template < int N >
void PAStar<N>::distributed_backtrace_n_print()
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