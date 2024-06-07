#include <mpi.h>
#include <iostream>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size); // total number of CPU cores (threads) on the system

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); // this is a CPU core identifier, could be threads if --use-hwthread-cpus flag is used

    std::cout << "Hello from processor " << world_rank;
    std::cout << " out of " << world_size << " processors." << std::endl;

    if(world_rank == 0) {
        std::cout << "All processors reached the barrier. Proceeding..." << std::endl;
    }

    MPI_Finalize();
    return 0;
}