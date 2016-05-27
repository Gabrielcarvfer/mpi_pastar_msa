/*!
 * \author Daniel Sundfeld
 * \copyright MIT License
 *
 * \brief The main function for msa_pastar project
 */
					  
#include <iostream>
#include "HeuristicHPair.h"
#include "max_seq_helper.h"
#include "mpi_dependencies.h"
#include "msa_options.h"
#include "PAStar.h"
#include "Sequences.h"
#include "read_fasta.h"



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
		MPI_Abort(MPI_COMM_WORLD, 1);
		MPI_Finalize();
		exit(1);
	}

    PAStarOpt opt;
    std::string filename;

    MPI_Comm_rank(MPI_COMM_WORLD, &opt.mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &opt.mpiCommSize);

    opt.mpiMin = 0 + opt.mpiRank * opt.threads_num;
    opt.mpiMax = opt.mpiMin + opt.threads_num;
    opt.totalThreads = opt.mpiCommSize * opt.threads_num;

    //std::cout << mpiRank << mpiCommSize << std::endl;
    //Todos os processos carregam argumentos de configuração
	if (msa_pastar_options(argc, argv, filename, opt) != 0)
	{
        //no caso de erro nas opcoes
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int numSeq, seqLen, i;

    Sequences *sequences = Sequences::getInstance();

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
	//std::cout << numSeq << std::endl;

        //send sequences to other processes
	MPI_Bcast(&numSeq, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//for each read sequence
	for (i = 0; i < numSeq; i++)
	{
	    //load sequence
	    std::string seq = sequences->get_seq(i);

	    //broadcast size and content
	    seqLen = seq.size()+1;
	    MPI_Bcast(&seqLen, 1, MPI_INT, 0, MPI_COMM_WORLD);
	    MPI_Bcast((void*)seq.c_str(), seqLen, MPI_CHAR, 0, MPI_COMM_WORLD);
	    }
    }
    else
    {
        //TODO: receive fasta from node 0
        MPI_Bcast(&numSeq, 1, MPI_INT, 0, MPI_COMM_WORLD);
        //std::cout << numSeq << std::endl;
	
	char * seqBuff = NULL;
        //for each announced sequence
        for (i = 0; i < numSeq; i++)
        {
            //receive size of next sequence
            MPI_Bcast(&seqLen, 1, MPI_INT, 0, MPI_COMM_WORLD);
            //std::cout << seqLen << std::endl;

            //prepare buffer
            seqBuff = new char[seqLen];

            //receive sequence
            MPI_Bcast(seqBuff, seqLen, MPI_CHAR, 0, MPI_COMM_WORLD);

            std::string seq (seqBuff);

            //save sequence to local process
            sequences->set_seq(seq);

	    delete[] seqBuff;
        }

        //std::string seq = sequences->get_seq(0);
	}

    //Espera todos os nos finalizarem sincronizacao para iniciar
    //MPI_Barrier(MPI_COMM_WORLD);

    int ret = pa_star_run(opt);

    //Espera todos os nos terminarem o processo e threads para finalizar
    MPI_Barrier(MPI_COMM_WORLD);
	//std::cout << "pastar returned " << ret << std::endl;
    //Termina graciosamente
    MPI_Finalize();

    return ret;
}
