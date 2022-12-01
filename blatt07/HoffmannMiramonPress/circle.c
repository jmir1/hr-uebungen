#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// This function initializes the buffer values.
int *init(int rank, int alloc_size, int chunk_size) {

    int *buf = (int *)malloc(sizeof(int) * alloc_size);

    // Add the rank to generate different seeds
    srand(time(NULL) + rank);

    for (int i = 0; i < chunk_size; i++) {
        // Do not modify "% 25"
        buf[i] = rand() % 25;
    }

    // Mark the unused space in the buffer via -1
    if (chunk_size < alloc_size) {
        buf[alloc_size - 1] = -1;
    }

    return buf;
}

// In this function, a process passes the buffer to its neighboring process until the
// target value is reached in the process with last rank.
int circle(int *buf, int rank, int nprocs, int alloc_size, int target, int *iteration) {
    int abort;

    while (1) {
        int next = rank + 1;
        int previous = rank - 1;

        if (rank == nprocs - 1) {
            next = 0;
        } else if (rank == 0) {
            previous = nprocs - 1;
        }

        MPI_Request request;
        MPI_Issend(buf, alloc_size, MPI_INT, next, 0, MPI_COMM_WORLD, &request);
        fflush(stdout);
        MPI_Status status;

        MPI_Recv(buf, alloc_size, MPI_INT, previous, 0, MPI_COMM_WORLD, &status);
        MPI_Wait(&request, &status);

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

        *iteration += 1;
    }
    return 0;
}

// This function sets the first value of the process with rank 0 as the target value
// and sends it to the process with last rank.
void target_communication(int rank, int nprocs, int *buf, int *target) {
    if (rank == 0) {
        MPI_Ssend(buf, 1, MPI_INT, nprocs - 1, 0, MPI_COMM_WORLD);
        *target = buf[0];
    } else if (rank == nprocs - 1) {
        MPI_Status status;
        MPI_Recv(buf, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        *target = *buf;
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

// This function prints the buffer of every process in order of its ranks.
void print_array(int rank, int nprocs, int alloc_size, int *buf) {
    if (rank == 0) {
        for (int i = 0; i < alloc_size; i++) {
            // Don't print surplus values in the buffer.
            if (buf[i] == -1) {
                continue;
            }
            printf("rank %d: %d\n", rank, buf[i]);
        }

        MPI_Status status;
        for (int p = 1; p < nprocs; p++) {
            int proc_buf[alloc_size];

            MPI_Recv(proc_buf, alloc_size, MPI_INT, p, 0, MPI_COMM_WORLD, &status);
            for (int i = 0; i < alloc_size; i++) {
                // Don't print surplus values in the buffer.
                if (proc_buf[i] == -1) {
                    continue;
                }
                printf("rank %d: %d\n", p, proc_buf[i]);
            }
        }

    } else {
        MPI_Ssend(buf, alloc_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}

int main(int argc, char **argv) {
    int N;
    int rank;
    int nprocs;
    int *buf;
    int alloc_size; // size of allocated memory
    int chunk_size; // size of relevant memory
    int target;
    int iteration = 0;

    if (argc < 2) {
        printf("Arguments error!\nPlease specify a buffer size.\n");
        return EXIT_FAILURE;
    }

    // Array length
    N = atoi(argv[1]);

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (nprocs > N) {
        // TODO
        printf("Array mit Größe %d lässt sich nicht auf alle Prozesse (%d) aufteilen\n",
               N, nprocs);
        return EXIT_FAILURE;
    }

    int quotient = N / nprocs;
    int remainder = N % nprocs;

    // Set alloc and chunk size depending on whether the array can be split evenly.
    if (remainder != 0) {
        alloc_size = quotient + 1;
        if (rank < remainder) {
            chunk_size = quotient + 1;
        } else {
            chunk_size = quotient;
        }
    } else {
        alloc_size = quotient;
        chunk_size = quotient;
    }

    buf = init(rank, alloc_size, chunk_size);

    target_communication(rank, nprocs, buf, &target);

    if (rank == 0) {
        printf("\nBEFORE\n");
    }
    print_array(rank, nprocs, alloc_size, buf);

    circle(buf, rank, nprocs, alloc_size, target, &iteration);

    if (rank == 0) {
        printf("\nAFTER\n");
    }
    print_array(rank, nprocs, alloc_size, buf);

    if (rank == 0) {
        printf("\nTARGET: %d\n", target);
        printf("ITERATION: %d\n", iteration);
    }

    free(buf);

    return EXIT_SUCCESS;
}

// TODO: 1 process
// TODO: nprocs > N
// TODO: nprocs == N

/*
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

int circle(int *buf, int rank, int nprocs, int *chunk_size, int remainder, int quotient,
           int target, int *iteration) {
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
        MPI_Issend(buf, *chunk_size, MPI_INT, next, 0, MPI_COMM_WORLD, &request);
        fflush(stdout);
        MPI_Status status;

        if (remainder != 0) {
            if (end >= nprocs) {
                if (rank >= start || rank <= end % nprocs) {
                    *chunk_size = quotient + 1;
                } else {
                    *chunk_size = quotient;
                }
            } else {
                if (rank >= start && rank <= end) {
                    *chunk_size = quotient + 1;
                } else {
                    *chunk_size = quotient;
                }
            }
        } else {
            *chunk_size = quotient;
        }
        // fprintf(stdout, "\n%d: rank: %d, chunk_size: %d, (%d-%d)\n", *iteration,
        // rank,
        //         chunk_size, start, end);
        MPI_Recv(buf, *chunk_size, MPI_INT, previous, 0, MPI_COMM_WORLD, &status);
        MPI_Wait(&request, &status);
        // fprintf(stdout, "\n%d: rank: %d, buf[0]: %d, chunk_size: %d, (%d-%d)\n",
        //         *iteration, rank, buf[0], proc_chunk_size, start, end);
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

        *iteration += 1;
    }
    return 0;
}

int target_communication(int rank, int nprocs, int *buf) {
    int target = -1;
    if (rank == 0) {
        MPI_Ssend(buf, 1, MPI_INT, nprocs - 1, 0, MPI_COMM_WORLD);
        target = buf[0];
    } else if (rank == nprocs - 1) {
        MPI_Status status;
        MPI_Recv(buf, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        target = *buf;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return target;
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

void print_array2(int N, int rank, int nprocs, int chunk_size, int alloc_size, int *buf,
                  int quotient, int remainder, int iteration) {
    if (rank == 0) {
        for (int i = 0; i < chunk_size; i++) {
            printf("rank %d: %d\n", rank, buf[i]);
        }

        MPI_Status status;
        for (int p = 1; p < nprocs; p++) {
            int proc_chunk_size;

            int start = 0 + iteration + 1;
            int end = remainder + iteration;
            if (remainder != 0) {
                if (end >= nprocs) {
                    if (p >= start || p <= end % nprocs) {
                        proc_chunk_size = quotient + 1;
                    } else {
                        proc_chunk_size = quotient;
                    }
                } else {
                    if (p >= start && p <= end) {
                        proc_chunk_size = quotient + 1;
                    } else {
                        proc_chunk_size = quotient;
                    }
                }
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
    int target;
    int iteration = 0;

    if (argc < 2) {
        printf("Arguments error!\nPlease specify a buffer size.\n");
        return EXIT_FAILURE;
    }

    // Array length
    N = atoi(argv[1]);

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (nprocs > N) {
        // TODO
        printf("Array mit Größe %d lässt sich nicht auf alle Prozesse (%d) aufteilen\n",
               N, nprocs);
        return EXIT_FAILURE;
    }

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

    ret =
        circle(buf, rank, nprocs, &chunk_size, remainder, quotient, target, &iteration);

    if (rank == 0) {
        printf("\nAFTER\n");
    }
    print_array2(N, rank, nprocs, chunk_size, alloc_size, buf, quotient, remainder,
                 iteration);

    if (rank == 0) {
        printf("\nTARGET: %d\n", target);
        printf("ITERATION: %d\n", iteration);
    }

    free(buf);

    return EXIT_SUCCESS;
}
*/