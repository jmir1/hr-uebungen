#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

int *init(int rank, int alloc_size, int chunk_size) {

    int *buf = (int *)malloc(sizeof(int) * alloc_size);

    srand(time(NULL) + rank);

    for (int i = 0; i < chunk_size; i++) {
        // Do not modify "% 25"
        buf[i] = rand() % 25;
    }

    return buf;
}

int circle(int *buf) {
    // TODO
    return 0;
}

void target_communication(int rank, int nprocs, int *buf) {
    if (rank == 0) {
        MPI_Ssend(buf, 1, MPI_INT, nprocs - 1, 0, MPI_COMM_WORLD);
    } else if (rank == nprocs - 1) {
        MPI_Status status;
        MPI_Recv(buf, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        int target = *buf;
        printf("\nTARGET: %d\n", target);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    sleep(0);
}

void print_array(int N, int rank, int nprocs, int chunk_size, int alloc_size, int *buf,
                 int quotient, int remainder) {
    int proc_chunk_size;

    if (rank == 0) {
        for (int i = 0; i < chunk_size; i++) {
            printf("rank %d: %d\n", rank, buf[i]);
        }

        MPI_Status status;
        for (int p = 1; p < nprocs; p++) {
            if (p < remainder) {
                proc_chunk_size = quotient + 1;
            } else {
                proc_chunk_size = quotient;
            }
            int proc_buf[proc_chunk_size];

            MPI_Recv(proc_buf, proc_chunk_size, MPI_INT, p, 0, MPI_COMM_WORLD, &status);
            for (int i = 0; i < proc_chunk_size; i++) {
                printf("rank %d: %d\n", p, proc_buf[i]);
            }
        }

    } else {
        MPI_Ssend(buf, chunk_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}

int main(int argc, char **argv) {
    int N;
    int rank;
    int nprocs;
    int *buf;
    int ret;
    int alloc_size;
    int chunk_size;

    if (argc < 2) {
        printf("Arguments error!\nPlease specify a buffer size.\n");
        return EXIT_FAILURE;
    }

    // Array length
    N = atoi(argv[1]);

    if (nprocs > N) {
        // TODO
        printf("Array lässt sich nicht auf alle Prozesse aufteilen\n");
        return EXIT_FAILURE;
    }

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int quotient = N / nprocs;
    int remainder = N % nprocs;

    if (remainder != 0) {
        alloc_size = quotient + 1;
        if (rank < remainder) {
            chunk_size = alloc_size;
        } else {
            chunk_size = quotient;
        }
    } else {
        alloc_size = quotient;
        chunk_size = quotient;
    }

    buf = init(rank, alloc_size, chunk_size);

    target_communication(rank, nprocs, buf);

    if (rank == 0) {
        printf("\nBEFORE\n");
    }
    print_array(N, rank, nprocs, chunk_size, alloc_size, buf, quotient, remainder);

    ret = circle(buf);

    if (rank == 0) {
        printf("\nAFTER\n");
    }
    print_array(N, rank, nprocs, chunk_size, alloc_size, buf, quotient, remainder);

    return EXIT_SUCCESS;
}
