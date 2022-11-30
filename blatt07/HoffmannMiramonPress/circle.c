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

int get_chunk_size(int rank, int remainder, int quotient) {
    if (rank < remainder) {
        return quotient + 1;
    } else {
        return quotient;
    }
}

int circle(int *buf, int rank, int nprocs, int chunk_size, int remainder, int quotient,
           int target) {
    int iteration = 0;
    int start = 0;
    int end = remainder - 1;
    int abort;

    while (1) {
        int next = rank + 1;
        int previous = rank - 1;
        start += 1;
        end += 1;

        if (rank == nprocs - 1) {
            next = 0;
        } else if (rank == 0) {
            previous = nprocs - 1;
        }

        MPI_Request request;
        MPI_Issend(buf, chunk_size, MPI_INT, next, 0, MPI_COMM_WORLD, &request);
        fflush(stdout);
        MPI_Status status;

        if (remainder != 0) {
            if (end >= nprocs) {
                if (rank >= start || rank <= end % nprocs) {
                    chunk_size = quotient + 1;
                } else {
                    chunk_size = quotient;
                }
            } else {
                if (rank >= start && rank <= end) {
                    chunk_size = quotient + 1;
                } else {
                    chunk_size = quotient;
                }
            }
        } else {
            chunk_size = quotient;
        }
        // fprintf(stdout, "\n%d: rank: %d, chunk_size: %d, (%d-%d)\n", iteration, rank,
        //         chunk_size, start, end);
        MPI_Recv(buf, chunk_size, MPI_INT, previous, 0, MPI_COMM_WORLD, &status);
        MPI_Wait(&request, &status);
        // fprintf(stdout, "\n%d: rank: %d, buf[0]: %d, chunk_size: %d, (%d-%d)\n",
        //         iteration, rank, buf[0], proc_chunk_size, start, end);
        if (rank == nprocs - 1) {
            if (target == buf[0]) {
                abort = 1;
            } else {
                abort = 0;
            }
        }

        MPI_Bcast(&abort, 1, MPI_INT, nprocs - 1, MPI_COMM_WORLD);
        if (abort == 1) {
            break;
        }

        iteration += 1;
        sleep(1);
    }
    printf("\nTARGET: %d\n", target);
    printf("\nITERATION: %d\n", iteration);

    return 0;
}

int target_communication(int rank, int nprocs, int *buf) {
    if (rank == 0) {
        MPI_Ssend(buf, 1, MPI_INT, nprocs - 1, 0, MPI_COMM_WORLD);
    } else if (rank == nprocs - 1) {
        MPI_Status status;
        MPI_Recv(buf, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        int target = *buf;
        printf("\nTARGET: %d\n", target);

        MPI_Barrier(MPI_COMM_WORLD);
        sleep(0);

        return target;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    sleep(0);
}

void print_array(int N, int rank, int nprocs, int chunk_size, int alloc_size, int *buf,
                 int quotient, int remainder) {
    if (rank == 0) {
        for (int i = 0; i < chunk_size; i++) {
            printf("rank %d: %d\n", rank, buf[i]);
        }

        MPI_Status status;
        for (int p = 1; p < nprocs; p++) {
            int proc_chunk_size = get_chunk_size(p, remainder, quotient);
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
    int target;

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
        chunk_size = get_chunk_size(rank, remainder, quotient);
    } else {
        alloc_size = quotient;
        chunk_size = quotient;
    }

    buf = init(rank, alloc_size, chunk_size);

    target = target_communication(rank, nprocs, buf);

    if (rank == 0) {
        printf("\nBEFORE\n");
    }
    print_array(N, rank, nprocs, chunk_size, alloc_size, buf, quotient, remainder);

    ret = circle(buf, rank, nprocs, chunk_size, remainder, quotient, target);

    if (rank == 0) {
        printf("\nAFTER\n");
    }
    print_array(N, rank, nprocs, chunk_size, alloc_size, buf, quotient, remainder);

    return EXIT_SUCCESS;
}
