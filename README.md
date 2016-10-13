Parallel A-Star MSA with MPI:

This project is an extension of Sundfeld's Parallel A-Star MSA, enabling the usage of two or more machines by using MPI to exchange data between different nodes. You can check out his amazing work in the following website https://bitbucket.org/danielsundfeld/astar_msa

#Pre-requisites
You will need Boost C++, LZ4 and Mpich to run that program. If you don't have them installed, it's a good time to do that.

#Building
To build, enter the cloned folder and execute make. A folder ./bin will be created with pastar executable inside.

#Running
You can run the program with mpiexec command, just like any other program with MPI. You don't really need to have NFS running to run that program, only make sure that the executable path is available in every MPI node of the cluster. The rank 0 will read the input file and send to all other nodes.

#More info
Check Sundfeld's page for more info. 

#Citing
Check Sundfeld's page for more info, until I publish my own paper.