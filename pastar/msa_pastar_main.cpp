/*!
 * \authors {Daniel Sundfeld, Gabriel Ferreira}
 * \copyright MIT License
 *
 * \brief The main function for msa_pastar project
 */

#include <iostream>
#include <mpi.h>
#include "include/HeuristicHPair.h"
#include "include/max_seq_helper.h"
#include "include/msa_options.h"
#include "include/PAStar.h"
#include "include/Sequences.h"
#include "include/read_fasta.h"
#include "include/lz4sup.h"



int pa_star_run_core(const PAStarOpt &opt)
{
    HeuristicHPair::getInstance()->init();

    // This macro is expanded to every supported number of sequences
    #define RUN_PASTAR(X) \
        case X : \
            return PAStar< X >::pa_star(Sequences::get_initial_node< X >(), Sequences::get_final_coord< X >(), opt);

    std::cout << "Performing search with Parallel A-Star.\n";
    switch (Sequences::get_seq_num())
    {
        MAX_NUM_SEQ_HELPER( RUN_PASTAR );
        default:
            std::cerr << "Fatal error: Invalid number of sequences: " << Sequences::get_seq_num() << std::endl;
    }
    return -1;
}

int pa_star_run(const PAStarOpt &opt)
{
    try
    {
        return pa_star_run_core(opt);
    }
    catch (std::exception &e)
    {
        std::cerr << "Running fatal error: " << e.what() << std::endl;
    }
    catch (...)
    {
        std::cerr << "Unknown fatal error while running!\n";
    }
    return -1;
}

int main(int argc, char *argv[])
{
	//Inicia MPI e comunicador global
    //MPI_Init(&argc, &argv);
	int provided = 0;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

	if (provided != MPI_THREAD_MULTIPLE)
	{
		std::cout << "asked " << MPI_THREAD_MULTIPLE << " and " << provided << " was provided" << std::endl;
		//MPI_Abort(MPI_COMM_WORLD, 1);
		//MPI_Finalize();
		//exit(1);
	}

    MPI_Barrier(MPI_COMM_WORLD);
    PAStarOpt opt;
    std::string filename;

    MPI_Comm_rank(MPI_COMM_WORLD, &opt.mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &opt.mpiCommSize);

    //std::cout << opt.mpiRank << opt.mpiCommSize << std::endl;

    //std::cout << opt.mpiRank << ": fase 1" << std::endl;
    //Todos os processos carregam argumentos de configuração
    if (msa_pastar_options(argc, argv, filename, opt) != 0)
    {
        //no caso de erro nas opcoes
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    opt.mpiMin = 0 + opt.mpiRank * opt.threads_num;
    opt.mpiMax = opt.mpiMin + opt.threads_num;
    opt.totalThreads = opt.mpiCommSize * opt.threads_num;

    int numSeq, seqLen, i;

    Sequences *sequences = Sequences::getInstance();

    //std::cout << opt.mpiRank << ": fase 2" << std::endl;
    //Se noh de rank 0
    if (opt.mpiRank == 0)
    {
        //read sequences from file
        if (read_fasta_file(filename) != 0)
        {
            //erro no arquivo, aborta execucao
            MPI_Abort(MPI_COMM_WORLD, 1);
            MPI_Finalize();
            exit(1);
        }

        //retrieve number of sequences
        numSeq = Sequences::get_seq_num();

        //Prepare structures for serialization
        std::ostringstream ss (std::ios_base::binary);
        boost::archive::binary_oarchive oa{ ss };

        oa & numSeq;

        //for each read sequence
        for (i = 0; i < numSeq; i++)
        {
	       //load sequence
	       std::string seq = sequences->get_seq(i);
           oa & seq;
        }


        char ** lz4Data = NULL;
        lz4Data = (char**)calloc(1, sizeof(char*));

        int lz4Size = pastar_lz4_en( ss.str().c_str(), lz4Data, ss.str().length() );
        //Sending message
        
        //Sending message
        for (i = 1; i < opt.mpiCommSize; i++)
            MPI_Send(*lz4Data, lz4Size, MPI_CHAR, i, 0, MPI_COMM_WORLD);

        //Clearing buffer
        delete[] *lz4Data;
        free(lz4Data);
    }
    else
    {
            int i, sender, n_bytes = 0;
            MPI_Status status;

            MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_CHAR, &n_bytes);
            sender = status.MPI_SOURCE;

            char ** buffer = NULL;
            buffer = (char**)calloc(1, sizeof(char*));
            *buffer = new char[n_bytes]();

            MPI_Recv(*buffer, n_bytes, MPI_CHAR, sender, 0, MPI_COMM_WORLD, &status);

            int len = pastar_lz4_dec(buffer);

            //Unserialize data into place
            std::istringstream ss(std::string(*buffer, *buffer + len), std::ios_base::binary);
            
            //Free buffer
            delete[] *buffer;
            boost::archive::binary_iarchive ia{ ss };

            std::string seq;
            int numSeq = 0;

            ia & numSeq;
            //Recover stuff in same order as saved
            for (int j = 0; j < numSeq; j++)
            {
                ia & seq;
                //save sequence to local process
                sequences->set_seq(seq);
            }   

            free(buffer);

        
    }

    //std::cout << opt.mpiRank << ": fase 3" << std::endl;

    //Inicia rotina principal
    int ret = pa_star_run(opt);

    //Espera todos os nos terminarem o processo e threads para finalizar
    MPI_Barrier(MPI_COMM_WORLD);

    //Termina graciosamente
    MPI_Finalize();

    return ret;
}
