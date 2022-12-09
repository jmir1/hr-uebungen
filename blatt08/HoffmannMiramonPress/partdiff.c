/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                 TU München - Institut für Informatik                   **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff.c                                                  **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauß-Seidel and    **/
/**            Jacobi method.                                              **/
/**                                                                        **/
/****************************************************************************/
/****************************************************************************/

/* ************************************************************************ */
/* Include standard header file.                                            */
/* ************************************************************************ */
#define _POSIX_C_SOURCE 200809L

#include <inttypes.h>
#include <malloc.h>
#include <math.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "partdiff.h"

struct calculation_arguments {
    uint64_t N;            /* number of spaces between lines (lines=N+1)     */
    uint64_t num_matrices; /* number of matrices                             */
    double h;              /* length of a space between two lines            */
    double ***Matrix;      /* index matrix used for addressing M             */
    double *M;             /* two matrices with real values                  */
};

struct calculation_results {
    uint64_t m;
    uint64_t stat_iteration; /* number of current iteration                    */
    double stat_precision;   /* actual precision of all slaves in iteration    */
};

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time; /* time when program started                      */
struct timeval comp_time;  /* time when calculation completed                */

/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static void initVariables(struct calculation_arguments *arguments,
                          struct calculation_results *results,
                          struct options const *options) {
    arguments->N = (options->interlines * 8) + 9 - 1;
    arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
    arguments->h = 1.0 / arguments->N;

    results->m = 0;
    results->stat_iteration = 0;
    results->stat_precision = 0;
}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static void freeMatrices(struct calculation_arguments *arguments) {
    uint64_t i;

    for (i = 0; i < arguments->num_matrices; i++) {
        free(arguments->Matrix[i]);
    }

    free(arguments->Matrix);
    free(arguments->M);
}

/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
static void *allocateMemory(size_t size) {
    void *p;

    if ((p = malloc(size)) == NULL) {
        printf("Speicherprobleme! (%" PRIu64 " Bytes angefordert)\n", size);
        exit(1);
    }

    return p;
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static void allocateMatrices(struct calculation_arguments *arguments) {
    uint64_t i, j;

    uint64_t const N = arguments->N;

    arguments->M =
        allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double));
    arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double **));

    for (i = 0; i < arguments->num_matrices; i++) {
        arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double *));

        for (j = 0; j <= N; j++) {
            arguments->Matrix[i][j] =
                arguments->M + (i * (N + 1) * (N + 1)) + (j * (N + 1));
        }
    }
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static void allocateMatrices_MPI(struct calculation_arguments *arguments, int lines) {
    uint64_t i, j;

    uint64_t const N = arguments->N;
    uint64_t p_N = lines + 2;

    arguments->M =
        allocateMemory(arguments->num_matrices * p_N * (N + 1) * sizeof(double));
    arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double **));

    for (i = 0; i < arguments->num_matrices; i++) {
        arguments->Matrix[i] = allocateMemory((p_N + 1) * sizeof(double *));

        for (j = 0; j < p_N; j++) {
            arguments->Matrix[i][j] =
                arguments->M + (i * p_N * (N + 1)) + (j * (N + 1));
        }
    }
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static void initMatrices(struct calculation_arguments *arguments,
                         struct options const *options) {
    uint64_t g, i, j; /* local variables for loops */

    uint64_t const N = arguments->N;
    double const h = arguments->h;
    double ***Matrix = arguments->Matrix;

    /* initialize matrix/matrices with zeros */
    for (g = 0; g < arguments->num_matrices; g++) {
        for (i = 0; i <= N; i++) {
            for (j = 0; j <= N; j++) {
                Matrix[g][i][j] = 0.0;
            }
        }
    }

    /* initialize borders, depending on function (function 2: nothing to do) */
    if (options->inf_func == FUNC_F0) {
        for (i = 0; i < N; i++) {
            for (j = 0; j < arguments->num_matrices; j++) {
                Matrix[j][i][0] = 1 + (1 - (h * i)); // Linke Kante
                Matrix[j][N][i] = 1 - (h * i);       // Untere Kante
                Matrix[j][N - i][N] = h * i;         // Rechte Kante
                Matrix[j][0][N - i] = 1 + h * i;     // Obere Kante
            }
        }
    }
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static void initMatrices_MPI(struct calculation_arguments *arguments,
                             struct options const *options, int size, int rank,
                             uint64_t start, uint64_t end) {
    uint64_t g, i, j; /* local variables for loops */

    uint64_t const N = arguments->N;
    double const h = arguments->h;
    double ***Matrix = arguments->Matrix;

    uint64_t p_N = end - start + 1;

    /* initialize matrix/matrices with zeros */
    for (g = 0; g < arguments->num_matrices; g++) {
        for (i = 0; i <= p_N; i++) {
            for (j = 0; j <= N; j++) {
                Matrix[g][i][j] = 0.0;
            }
        }
    }

    /* initialize borders, depending on function (function 2: nothing to do) */
    if (options->inf_func == FUNC_F0) {
        for (g = 0; g < arguments->num_matrices; g++) {
            for (i = start - 1; i < end; i++) {
                Matrix[g][i - (start - 1)][0] = 1 + (1 - (h * i));     // Linke Kante
                Matrix[g][(i - (start - 1) + 1)][N] = h * (N - i - 1); // Rechte Kante
            }

            if (rank == 0) {
                for (j = 0; j < N; j++) {
                    Matrix[g][0][N - j] = 1 + h * j; // Obere Kante
                }
            }

            if (rank == size - 1 || end == N) {
                for (j = 0; j < N; j++) {
                    Matrix[g][p_N][j] = 1 - (h * j); // Untere Kante
                }
            }
        }
    }
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static void calculate(struct calculation_arguments const *arguments,
                      struct calculation_results *results,
                      struct options const *options) {
    int i, j;           /* local variables for loops */
    int m1, m2;         /* used as indices for old and new matrices */
    double star;        /* four times center value minus 4 neigh.b values */
    double residuum;    /* residuum of current iteration */
    double maxResiduum; /* maximum residuum value of a slave in iteration */

    int const N = arguments->N;
    double const h = arguments->h;

    double pih = 0.0;
    double fpisin = 0.0;

    int term_iteration = options->term_iteration;

    /* initialize m1 and m2 depending on algorithm */
    if (options->method == METH_JACOBI) {
        m1 = 0;
        m2 = 1;
    } else {
        m1 = 0;
        m2 = 0;
    }

    if (options->inf_func == FUNC_FPISIN) {
        pih = PI * h;
        fpisin = 0.25 * TWO_PI_SQUARE * h * h;
    }

    while (term_iteration > 0) {
        double **Matrix_Out = arguments->Matrix[m1];
        double **Matrix_In = arguments->Matrix[m2];

        maxResiduum = 0;

        /* over all rows */
        for (i = 1; i < N; i++) {
            double fpisin_i = 0.0;

            if (options->inf_func == FUNC_FPISIN) {
                fpisin_i = fpisin * sin(pih * (double)i);
            }

            /* over all columns */
            for (j = 1; j < N; j++) {
                star = 0.25 * (Matrix_In[i - 1][j] + Matrix_In[i][j - 1] +
                               Matrix_In[i][j + 1] + Matrix_In[i + 1][j]);

                if (options->inf_func == FUNC_FPISIN) {
                    star += fpisin_i * sin(pih * (double)j);
                }

                if (options->termination == TERM_PREC || term_iteration == 1) {
                    residuum = Matrix_In[i][j] - star;
                    residuum = (residuum < 0) ? -residuum : residuum;
                    maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
                }

                Matrix_Out[i][j] = star;
            }
        }

        results->stat_iteration++;
        results->stat_precision = maxResiduum;

        /* exchange m1 and m2 */
        i = m1;
        m1 = m2;
        m2 = i;

        /* check for stopping calculation depending on termination method */
        if (options->termination == TERM_PREC) {
            if (maxResiduum < options->term_precision) {
                term_iteration = 0;
            }
        } else if (options->termination == TERM_ITER) {
            term_iteration--;
        }
    }

    results->m = m2;
}

/* ************************************************************************ */
/* calculate: solves the equation with MPI parallelization */
/* ************************************************************************ */
static void calculate_MPI(struct calculation_arguments const *arguments,
                          struct calculation_results *results,
                          struct options const *options, int rank, int size, int start,
                          int end, MPI_Comm calc_comm) {
    int i, j;           /* local variables for loops */
    int m1, m2;         /* used as indices for old and new matrices */
    double star;        /* four times center value minus 4 neigh.b values */
    double residuum;    /* residuum of current iteration */
    double maxResiduum; /* maximum residuum value of a slave in iteration */

    int const N = arguments->N;
    double const h = arguments->h;
    int p_N = end - start + 1;

    double pih = 0.0;
    double fpisin = 0.0;

    int term_iteration = options->term_iteration;

    /* initialize m1 and m2 */
    m1 = 0;
    m2 = 1;

    if (options->inf_func == FUNC_FPISIN) {
        pih = PI * h;
        fpisin = 0.25 * TWO_PI_SQUARE * h * h;
    }

    while (term_iteration > 0) {
        double **Matrix_Out = arguments->Matrix[m1];
        double **Matrix_In = arguments->Matrix[m2];

        maxResiduum = 0;

        /* over all rows */
        for (i = 1; i < p_N; i++) {
            double fpisin_i = 0.0;

            if (options->inf_func == FUNC_FPISIN) {
                fpisin_i = fpisin * sin(pih * (double)i);
            }

            /* over all columns */
            for (j = 1; j < N; j++) {
                star = 0.25 * (Matrix_In[i - 1][j] + Matrix_In[i][j - 1] +
                               Matrix_In[i][j + 1] + Matrix_In[i + 1][j]);

                if (options->inf_func == FUNC_FPISIN) {
                    star += fpisin_i * sin(pih * (double)j);
                }

                if (options->termination == TERM_PREC || term_iteration == 1) {
                    residuum = Matrix_In[i][j] - star;
                    residuum = (residuum < 0) ? -residuum : residuum;
                    maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
                }

                Matrix_Out[i][j] = star;
            }
        }

        int next = rank + 1;
        int previous = rank - 1;

        MPI_Status status;
        if (rank == 0) {
            MPI_Request request_last_row;
            MPI_Issend(Matrix_Out[p_N - 1], (N + 1), MPI_DOUBLE, next, 0,
                       MPI_COMM_WORLD, &request_last_row);

            MPI_Recv(Matrix_Out[p_N], (N + 1), MPI_DOUBLE, next, 0, MPI_COMM_WORLD,
                     &status);

            MPI_Wait(&request_last_row, &status);

        } else if (rank == size - 1 || end == N) {
            MPI_Request request_first_row;
            MPI_Issend(Matrix_Out[1], (N + 1), MPI_DOUBLE, previous, 0, MPI_COMM_WORLD,
                       &request_first_row);

            MPI_Recv(Matrix_Out[0], (N + 1), MPI_DOUBLE, previous, 0, MPI_COMM_WORLD,
                     &status);

            MPI_Wait(&request_first_row, &status);

        } else {
            // Send first row
            MPI_Request request_first_row, request_last_row;
            MPI_Issend(Matrix_Out[1], (N + 1), MPI_DOUBLE, previous, 0, MPI_COMM_WORLD,
                       &request_first_row);
            // Send last row
            MPI_Issend(Matrix_Out[p_N - 1], (N + 1), MPI_DOUBLE, next, 0,
                       MPI_COMM_WORLD, &request_last_row);

            MPI_Recv(Matrix_Out[0], (N + 1), MPI_DOUBLE, previous, 0, MPI_COMM_WORLD,
                     &status);
            MPI_Recv(Matrix_Out[p_N], (N + 1), MPI_DOUBLE, next, 0, MPI_COMM_WORLD,
                     &status);

            MPI_Wait(&request_first_row, &status);
            MPI_Wait(&request_last_row, &status);
        }

        /* exchange m1 and m2 */
        i = m1;
        m1 = m2;
        m2 = i;

        double all_maxResiduum;
        /* check for stopping calculation depending on termination method */
        if (options->termination == TERM_PREC) {
            MPI_Allreduce(&maxResiduum, &all_maxResiduum, 1, MPI_DOUBLE, MPI_MAX,
                          MPI_COMM_WORLD);
            maxResiduum = all_maxResiduum;
            if (maxResiduum < options->term_precision) {
                term_iteration = 0;
            }
        } else if (options->termination == TERM_ITER) {
            term_iteration--;
            if (term_iteration == 0) {
                MPI_Allreduce(&maxResiduum, &all_maxResiduum, 1, MPI_DOUBLE, MPI_MAX,
                              calc_comm);
                maxResiduum = all_maxResiduum;
            }
        }

        results->stat_iteration++;
        results->stat_precision = maxResiduum;
    }

    results->m = m2;
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static void displayStatistics(struct calculation_arguments const *arguments,
                              struct calculation_results const *results,
                              struct options const *options) {
    int N = arguments->N;
    double time = (comp_time.tv_sec - start_time.tv_sec) +
                  (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

    printf("Berechnungszeit:    %f s \n", time);
    printf("Speicherbedarf:     %f MiB\n", (N + 1) * (N + 1) * sizeof(double) *
                                               arguments->num_matrices / 1024.0 /
                                               1024.0);
    printf("Berechnungsmethode: ");

    if (options->method == METH_GAUSS_SEIDEL) {
        printf("Gauß-Seidel");
    } else if (options->method == METH_JACOBI) {
        printf("Jacobi");
    }

    printf("\n");
    printf("Interlines:         %" PRIu64 "\n", options->interlines);
    printf("Stoerfunktion:      ");

    if (options->inf_func == FUNC_F0) {
        printf("f(x,y) = 0");
    } else if (options->inf_func == FUNC_FPISIN) {
        printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
    }

    printf("\n");
    printf("Terminierung:       ");

    if (options->termination == TERM_PREC) {
        printf("Hinreichende Genaugkeit");
    } else if (options->termination == TERM_ITER) {
        printf("Anzahl der Iterationen");
    }

    printf("\n");
    printf("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration);
    printf("Norm des Fehlers:   %e\n", results->stat_precision);
    printf("\n");
}

/****************************************************************************/
/** Beschreibung der Funktion displayMatrix:                               **/
/**                                                                        **/
/** Die Funktion displayMatrix gibt eine Matrix                            **/
/** in einer "ubersichtlichen Art und Weise auf die Standardausgabe aus.   **/
/**                                                                        **/
/** Die "Ubersichtlichkeit wird erreicht, indem nur ein Teil der Matrix    **/
/** ausgegeben wird. Aus der Matrix werden die Randzeilen/-spalten sowie   **/
/** sieben Zwischenzeilen ausgegeben.                                      **/
/****************************************************************************/
static void displayMatrix(struct calculation_arguments *arguments,
                          struct calculation_results *results,
                          struct options *options) {
    int x, y;

    double **Matrix = arguments->Matrix[results->m];

    int const interlines = options->interlines;

    printf("Matrix:\n");

    for (y = 0; y < 9; y++) {
        for (x = 0; x < 9; x++) {
            printf("%7.4f", Matrix[y * (interlines + 1)][x * (interlines + 1)]);
        }

        printf("\n");
    }

    fflush(stdout);
}

/**
 * rank and size are the MPI rank and size, respectively.
 * from and to denote the global(!) range of lines that this process is responsible for.
 *
 * Example with 9 matrix lines and 4 processes:
 * - rank 0 is responsible for 1-2, rank 1 for 3-4, rank 2 for 5-6 and rank 3 for 7.
 *   Lines 0 and 8 are not included because they are not calculated.
 * - Each process stores two halo lines in its matrix (except for ranks 0 and 3 that
 * only store one).
 * - For instance: Rank 2 has four lines 0-3 but only calculates 1-2 because 0 and 3 are
 * halo lines for other processes. It is responsible for (global) lines 5-6.
 */
static void displayMatrix_MPI(struct calculation_arguments *arguments,
                              struct calculation_results *results,
                              struct options *options, int rank, int size, int start,
                              int end) {
    int const elements = 8 * options->interlines + 9;
    int N = arguments->N;

    int x, y;
    double **Matrix = arguments->Matrix[results->m];
    MPI_Status status;

    /* first line belongs to rank 0 */
    if (rank == 0) {
        start--;
        printf("Matrix:\n");
    }

    /* last line belongs to rank size - 1 */
    if (rank == size - 1 || end == N) {
        end++;
    }

    for (y = 0; y < 9; y++) {
        int line = y * (options->interlines + 1);

        if (rank == 0) {
            /* check whether this line belongs to rank 0 */
            if (line >= end) {
                /* use the tag to receive the lines in the correct order
                 * the line is stored in Matrix[0], because we do not need it anymore */
                MPI_Recv(Matrix[0], elements, MPI_DOUBLE, MPI_ANY_SOURCE, 42 + y,
                         MPI_COMM_WORLD, &status);
            }
        } else {
            if (line >= start && line < end) {
                /* if the line belongs to this process, send it to rank 0
                 * (line - from + 1) is used to calculate the correct local address */
                MPI_Send(Matrix[line - start + 1], elements, MPI_DOUBLE, 0, 42 + y,
                         MPI_COMM_WORLD);
            }
        }

        if (rank == 0) {
            for (x = 0; x < 9; x++) {
                int col = x * (options->interlines + 1);

                if (line < end) {
                    /* this line belongs to rank 0 */
                    printf("%7.4f", Matrix[line][col]);
                } else {
                    /* this line belongs to another rank and was received above */
                    printf("%7.4f", Matrix[0][col]);
                }
            }

            printf("\n");
        }
    }

    fflush(stdout);
}

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int main(int argc, char **argv) {
    struct options options;
    struct calculation_arguments arguments;
    struct calculation_results results;
    int rank, size, start, end;
    MPI_Comm calc_comm;

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (options.method == METH_GAUSS_SEIDEL && rank != 0) {
        return 0;
    }

    if (rank == 0) {
        askParams(&options, argc, argv);
    }

    MPI_Bcast(&options, sizeof(options), MPI_CHAR, 0, MPI_COMM_WORLD);

    initVariables(&arguments, &results, &options);

    if (options.method == METH_JACOBI && size > 1) {
        int N = arguments.N;

        int lines;
        int quotient = (N - 1) / size;
        int remainder = (N - 1) % size;

        if (rank < remainder) {
            lines = quotient + 1;
            start = rank * lines + 1;
        } else {
            lines = quotient;
            start = rank * lines + 1 + remainder;
        }
        end = start + lines;

        if (start == end) {
            MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, rank, &calc_comm);
            MPI_Barrier(MPI_COMM_WORLD);
            return 0;
        }

        MPI_Comm_split(MPI_COMM_WORLD, 0, rank, &calc_comm);

        allocateMatrices_MPI(&arguments, lines);
        initMatrices_MPI(&arguments, &options, size, rank, start, end);

        gettimeofday(&start_time, NULL);
        calculate_MPI(&arguments, &results, &options, rank, size, start, end,
                      calc_comm);
        gettimeofday(&comp_time, NULL);

        MPI_Comm_free(&calc_comm);

        if (rank == 0) {
            displayStatistics(&arguments, &results, &options);
        }
        displayMatrix_MPI(&arguments, &results, &options, rank, size, start, end);

    } else {
        allocateMatrices(&arguments);
        initMatrices(&arguments, &options);

        gettimeofday(&start_time, NULL);
        calculate(&arguments, &results, &options);
        gettimeofday(&comp_time, NULL);

        displayStatistics(&arguments, &results, &options);
        displayMatrix(&arguments, &results, &options);
    }

    freeMatrices(&arguments);
    MPI_Barrier(MPI_COMM_WORLD);

    return 0;
}
