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
    uint64_t p_N = lines + 2; // The number of rows of the process

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

    uint64_t p_N = end - start + 1; // The max index of rows of the process

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
static void calculate_MPI_Jacobi(struct calculation_arguments const *arguments,
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
        for (i = start; i < end; i++) {
            double fpisin_i = 0.0;
            int x = i - (start - 1); // Matrix index for i

            if (options->inf_func == FUNC_FPISIN) {
                fpisin_i = fpisin * sin(pih * (double)i);
            }

            /* over all columns */
            for (j = 1; j < N; j++) {
                star = 0.25 * (Matrix_In[x - 1][j] + Matrix_In[x][j - 1] +
                               Matrix_In[x][j + 1] + Matrix_In[x + 1][j]);

                if (options->inf_func == FUNC_FPISIN) {
                    star += fpisin_i * sin(pih * (double)j);
                }

                if (options->termination == TERM_PREC || term_iteration == 1) {
                    residuum = Matrix_In[x][j] - star;
                    residuum = (residuum < 0) ? -residuum : residuum;
                    maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
                }

                Matrix_Out[x][j] = star;
            }
        }

        int next = rank + 1;
        int previous = rank - 1;

        MPI_Status status;
        if (rank == 0) {
            // The first process only sends its last row and and recieves the lower halo
            // row.
            MPI_Request request_last_row;
            MPI_Issend(Matrix_Out[p_N - 1], (N + 1), MPI_DOUBLE, next, 0,
                       MPI_COMM_WORLD, &request_last_row);

            MPI_Recv(Matrix_Out[p_N], (N + 1), MPI_DOUBLE, next, 0, MPI_COMM_WORLD,
                     &status);

            MPI_Wait(&request_last_row, &status);

        } else if (rank == size - 1 || end == N) {
            // The last process only sends its first row and and recieves the upper halo
            // row.
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

            // Recieve upper halo row
            MPI_Recv(Matrix_Out[0], (N + 1), MPI_DOUBLE, previous, 0, MPI_COMM_WORLD,
                     &status);
            // Recieve lower halo row
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
            // Get the maxResiduum from all processes in calc_comm and send the largest
            // maxResiduum back.
            MPI_Allreduce(&maxResiduum, &all_maxResiduum, 1, MPI_DOUBLE, MPI_MAX,
                          calc_comm);
            maxResiduum = all_maxResiduum;
            if (maxResiduum < options->term_precision) {
                term_iteration = 0;
            }
        } else if (options->termination == TERM_ITER) {
            term_iteration--;
            if (term_iteration == 0) {
                // Get the maxResiduum from all processes in calc_comm and send the
                // largest maxResiduum back.
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
/* calculate: solves the equation with MPI parallelization */
/* ************************************************************************ */
static void calculate_MPI(struct calculation_arguments const *arguments,
                          struct calculation_results *results,
                          struct options const *options, int rank, int size, int start,
                          int end, MPI_Comm calc_comm) {
    int i, j;           /* local variables for loops */
    double star;        /* four times center value minus 4 neigh.b values */
    double residuum;    /* residuum of current iteration */
    double maxResiduum; /* maximum residuum value of a slave in iteration */
	int terminate = options->termination == TERM_ITER;  // termination signal. default: do not terminate

    int const N = arguments->N;
    double const h = arguments->h;
    int p_N = end - start + 1;

    double pih = 0.0;
    double fpisin = 0.0;

    int term_iteration = options->term_iteration;

    if (options->inf_func == FUNC_FPISIN) {
        pih = PI * h;
        fpisin = 0.25 * TWO_PI_SQUARE * h * h;
    }

    while (term_iteration > 0) {
        double **Matrix = arguments->Matrix[0];

        int next = rank + 1;
        int previous = rank - 1;


        MPI_Status status;
		if (rank > 0 && rank < size - 1 && end != N) {
            // Recieve upper halo row
			if (options->termination == TERM_ITER || term_iteration > 1) {
				MPI_Recv(Matrix[0], (N + 1), MPI_DOUBLE, previous, 0, calc_comm,
						&status);
			}
			
		}

        maxResiduum = 0;

        /* over all rows */
        for (i = start; i < end; i++) {
            double fpisin_i = 0.0;
            int x = i - (start - 1); // Matrix index for i

            if (options->inf_func == FUNC_FPISIN) {
                fpisin_i = fpisin * sin(pih * (double)i);
            }

            /* over all columns */
            for (j = 1; j < N; j++) {
                star = 0.25 * (Matrix[x - 1][j] + Matrix[x][j - 1] +
                               Matrix[x][j + 1] + Matrix[x + 1][j]);

                if (options->inf_func == FUNC_FPISIN) {
                    star += fpisin_i * sin(pih * (double)j);
                }

                if (options->termination == TERM_PREC || term_iteration == 1) {
                    residuum = Matrix[x][j] - star;
                    residuum = (residuum < 0) ? -residuum : residuum;
                    maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
                }

                Matrix[x][j] = star;
            }
        }

        if (rank == 0) {
            // The first process only sends its last row and and recieves the lower halo
            // row.
            MPI_Request request_last_row, request_max_residuum;
			MPI_Issend(Matrix[p_N - 1], (N + 1), MPI_DOUBLE, next, 0,
					calc_comm, &request_last_row);
			if (!terminate) {
				MPI_Issend(&maxResiduum, 1, MPI_DOUBLE, next, 0,
							calc_comm, &request_max_residuum);
			}

            MPI_Recv(Matrix[p_N], (N + 1), MPI_DOUBLE, next, 0, calc_comm,
                     &status);

			MPI_Wait(&request_last_row, &status);
			if (!terminate) {
				MPI_Wait(&request_max_residuum, &status);
				MPI_Recv(&terminate, 1, MPI_INT, next, 0, calc_comm, &status);
				if (terminate) {
					int receivedIteration;
					MPI_Allreduce(&results->stat_iteration, &receivedIteration, 1, MPI_INT, MPI_MAX, calc_comm);
					term_iteration = receivedIteration - results->stat_iteration + 1;
				}
			}
        } else if (rank == size - 1 || end == N) {
            // The last process only sends its first row and and recieves the upper halo
            // row.
            MPI_Request request_first_row;
			double receivedMaxResiduum = 0;
			if (terminate || term_iteration > 1) {
				MPI_Issend(Matrix[1], (N + 1), MPI_DOUBLE, previous, 0, calc_comm,
						&request_first_row);
			}
			MPI_Recv(Matrix[0], (N + 1), MPI_DOUBLE, previous, 0, calc_comm,
					&status);
			if (!terminate) {
				MPI_Recv(&receivedMaxResiduum, 1, MPI_DOUBLE, previous, 0, calc_comm,
				&status);

				maxResiduum = (receivedMaxResiduum < maxResiduum) ? maxResiduum : receivedMaxResiduum;
			}

			if (terminate || term_iteration > 1) {
            	MPI_Wait(&request_first_row, &status);
			}
			if (!terminate) {
				if (maxResiduum < options->term_precision) {
					terminate = 1; // terminate
					term_iteration = rank + 1;
            	}
				MPI_Request request_terminate;
				MPI_Issend(&terminate, 1, MPI_INT, previous, 0, calc_comm, &request_terminate);
				if (terminate) {
					int receivedIteration;
					MPI_Allreduce(&results->stat_iteration, &receivedIteration, 1, MPI_INT, MPI_MAX, calc_comm);
					term_iteration = receivedIteration - results->stat_iteration + 1;
				}
				MPI_Wait(&request_terminate, &status);
			}
        } else {
            // Send first row
            MPI_Request request_first_row, request_last_row, request_max_residuum, request_lower_halo;
			if (terminate || term_iteration > 1) {
            	MPI_Issend(Matrix[1], (N + 1), MPI_DOUBLE, previous, 0, calc_comm,
                       		&request_first_row);
			}
			// Send last row
			MPI_Issend(Matrix[p_N - 1], (N + 1), MPI_DOUBLE, next, 0,
					calc_comm, &request_last_row);
            // Recieve lower halo row
            MPI_Irecv(Matrix[p_N], (N + 1), MPI_DOUBLE, next, 0, calc_comm,
                     &request_lower_halo);
			if (!terminate) {
				double receivedMaxResiduum = 0;
				MPI_Recv(&receivedMaxResiduum, 1, MPI_DOUBLE, previous, 0, calc_comm,
						&status);
				maxResiduum = (receivedMaxResiduum < maxResiduum) ? maxResiduum : receivedMaxResiduum;
				MPI_Issend(&maxResiduum, 1, MPI_DOUBLE, next, 0,
							calc_comm, &request_max_residuum);
			}
			if (terminate || term_iteration > 1) {
				MPI_Wait(&request_first_row, &status);
			}
			MPI_Wait(&request_last_row, &status);
			MPI_Wait(&request_lower_halo, &status);
			if (!terminate) {
				MPI_Wait(&request_max_residuum, &status);
				MPI_Recv(&terminate, 1, MPI_INT, next, 0, calc_comm, &status);
				MPI_Request request_terminate;
				MPI_Issend(&terminate, 1, MPI_INT, previous, 0, calc_comm, &request_terminate);
				if (terminate) {
					int receivedIteration;
					MPI_Allreduce(&results->stat_iteration, &receivedIteration, 1, MPI_INT, MPI_MAX, calc_comm);
					term_iteration = receivedIteration - results->stat_iteration + 1;
				}
				MPI_Wait(&request_terminate, &status);
			}
        }

		if (terminate) {
			term_iteration--;
		}
		if (terminate && term_iteration == 0) {
			double all_maxResiduum = 0;
			MPI_Allreduce(&maxResiduum, &all_maxResiduum, 1, MPI_DOUBLE, MPI_MAX,
						calc_comm);
			maxResiduum = all_maxResiduum;
		}
		results->stat_iteration++;
        results->stat_precision = maxResiduum;
    }

    results->m = 0;
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

#include "displaymatrix-mpi.c"

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

    // The first process gets the input and sends it to the other processes.
    if (rank == 0) {
        askParams(&options, argc, argv);
    }
    MPI_Bcast(&options, sizeof(options), MPI_CHAR, 0, MPI_COMM_WORLD);

    initVariables(&arguments, &results, &options);

    if (options.method == METH_GAUSS_SEIDEL && size > 1) {
        // MPI parallelized computation if the method is Gauss-Seidel
        // and there is more than 1 process.

        int lines;
        int N = arguments.N;
        int quotient = (N - 1) / size;
        int remainder = (N - 1) % size;

        // Set the start and end indices for a process
        if (rank < remainder) {
            lines = quotient + 1;
            start = rank * lines + 1;
        } else {
            lines = quotient;
            start = rank * lines + 1 + remainder;
        }
        end = start + lines;

        // Returns unused processes if there are more processes than rows in the matrix.
        if (start == end) {
            MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, rank, &calc_comm);
            MPI_Finalize();
            return 0;
        }

        // Create a communicator with all processes that are involved in the
        // calculation.
        MPI_Comm_split(MPI_COMM_WORLD, 0, rank, &calc_comm);

        // Allocate and initialize the chunk of the matrix assigned to the process.
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
        DisplayMatrix(&arguments, &results, &options, rank, size, start, end);

    } else if (options.method == METH_JACOBI && size > 1) {
		
        // MPI parallelized computation if the method is Jacobi
        // and there is more than 1 process.

        int lines;
        int N = arguments.N;
        int quotient = (N - 1) / size;
        int remainder = (N - 1) % size;

        // Set the start and end indices for a process
        if (rank < remainder) {
            lines = quotient + 1;
            start = rank * lines + 1;
        } else {
            lines = quotient;
            start = rank * lines + 1 + remainder;
        }
        end = start + lines;

        // Returns unused processes if there are more processes than rows in the matrix.
        if (start == end) {
            MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, rank, &calc_comm);
            MPI_Finalize();
            return 0;
        }

        // Create a communicator with all processes that are involved in the
        // calculation.
        MPI_Comm_split(MPI_COMM_WORLD, 0, rank, &calc_comm);

        // Allocate and initialize the chunk of the matrix assigned to the process.
        allocateMatrices_MPI(&arguments, lines);
        initMatrices_MPI(&arguments, &options, size, rank, start, end);

        gettimeofday(&start_time, NULL);
        calculate_MPI_Jacobi(&arguments, &results, &options, rank, size, start, end,
                      calc_comm);
        gettimeofday(&comp_time, NULL);

        MPI_Comm_free(&calc_comm);

        if (rank == 0) {
            displayStatistics(&arguments, &results, &options);
        }
        DisplayMatrix(&arguments, &results, &options, rank, size, start, end);

	} else {
        // Sequential computation if the method is Jacobi and/or
        // there are less than 2 processes.

        allocateMatrices(&arguments);
        initMatrices(&arguments, &options);

        gettimeofday(&start_time, NULL);
        calculate(&arguments, &results, &options);
        gettimeofday(&comp_time, NULL);

        displayStatistics(&arguments, &results, &options);
        displayMatrix(&arguments, &results, &options);
    }

    freeMatrices(&arguments);
    MPI_Finalize();

    return 0;
}
