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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#ifdef OPENMP
#include <omp.h>
#endif

#include "partdiff.h"

struct calculation_arguments
{
    uint64_t N;            /* number of spaces between lines (lines=N+1)     */
    uint64_t num_matrices; /* number of matrices                             */
    double h;              /* length of a space between two lines            */
    double ***Matrix;      /* index matrix used for addressing M             */
    double *M;             /* two matrices with real values                  */
};

struct calculation_results
{
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
static void
initVariables(struct calculation_arguments *arguments, struct calculation_results *results, struct options const *options)
{
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
static void
freeMatrices(struct calculation_arguments *arguments)
{
    uint64_t i;

    for (i = 0; i < arguments->num_matrices; i++)
    {
        free(arguments->Matrix[i]);
    }

    free(arguments->Matrix);
    free(arguments->M);
}

/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
static void *
allocateMemory(size_t size)
{
    void *p;

    if ((p = malloc(size)) == NULL)
    {
        printf("Speicherprobleme! (%" PRIu64 " Bytes angefordert)\n", size);
        exit(1);
    }

    return p;
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static void
allocateMatrices(struct calculation_arguments *arguments)
{
    uint64_t i, j;

    uint64_t const N = arguments->N;

    arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double));
    arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double **));

    for (i = 0; i < arguments->num_matrices; i++)
    {
        arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double *));

        for (j = 0; j <= N; j++)
        {
            arguments->Matrix[i][j] = arguments->M + (i * (N + 1) * (N + 1)) + (j * (N + 1));
        }
    }
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static void
initMatrices(struct calculation_arguments *arguments, struct options const *options)
{
    uint64_t g, i, j; /* local variables for loops */

    uint64_t const N = arguments->N;
    double const h = arguments->h;
    double ***Matrix = arguments->Matrix;

    /* initialize matrix/matrices with zeros */
    for (g = 0; g < arguments->num_matrices; g++)
    {
        for (i = 0; i <= N; i++)
        {
            for (j = 0; j <= N; j++)
            {
                Matrix[g][i][j] = 0.0;
            }
        }
    }

    /* initialize borders, depending on function (function 2: nothing to do) */
    if (options->inf_func == FUNC_F0)
    {
        for (g = 0; g < arguments->num_matrices; g++)
        {
            for (i = 0; i <= N; i++)
            {
                Matrix[g][i][0] = 1 + (1 - (h * i)); // Linke Kante
                Matrix[g][N][i] = 1 - (h * i);       // Untere Kante
                Matrix[g][N - i][N] = h * i;         // Rechte Kante
                Matrix[g][0][N - i] = 1 + h * i;     // Obere Kante
            }
        }
    }
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static void
calculate_seq(struct calculation_arguments const *arguments,
			  struct calculation_results *results,
			  struct options const *options,
			  int first_row,
			  int last_row,
			  int first_column,
			  int last_column,
			  int elements)
{
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
    if (options->method == METH_JACOBI)
    {
        m1 = 0;
        m2 = 1;
    }
    else
    {
        m1 = 0;
        m2 = 0;
    }

    if (options->inf_func == FUNC_FPISIN)
    {
        pih = PI * h;
        fpisin = 0.25 * TWO_PI_SQUARE * h * h;
    }

    while (term_iteration > 0)
    {
        double **Matrix_Out = arguments->Matrix[m1];
        double **Matrix_In = arguments->Matrix[m2];

        maxResiduum = 0;

        /* over all rows */
        for (i = 1; i < N; i++)
        {
			if (!elements && i < first_row)
			{
				continue;
			}
			else if (!elements && i > last_row)
			{
				break;
			}
            double fpisin_i = 0.0;

            if (options->inf_func == FUNC_FPISIN)
            {
                fpisin_i = fpisin * sin(pih * (double)i);
            }

            /* over all columns */
            for (j = 1; j < N; j++)
            {
				if ((!elements && j < first_column) || (elements && i == first_row && j < first_column))
				{
					continue;
				}
				else if ((!elements && j > last_column) || (elements && i == last_row && j > last_column))
				{
					break;
				}
                star = 0.25 * (Matrix_In[i - 1][j] + Matrix_In[i][j - 1] + Matrix_In[i][j + 1] + Matrix_In[i + 1][j]);

                if (options->inf_func == FUNC_FPISIN)
                {
                    star += fpisin_i * sin(pih * (double)j);
                }

                if (options->termination == TERM_PREC || term_iteration == 1)
                {
                    residuum = Matrix_In[i][j] - star;
                    residuum = (residuum < 0) ? -residuum : residuum;
                    maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
                }
                Matrix_Out[i][j] = star;
            }
        }
        #ifdef OPENMP
		#pragma omp critical
		#endif
        {
            results->stat_iteration++;
            results->stat_precision = (results->stat_precision < maxResiduum) ? maxResiduum : results->stat_precision;
        }

        /* exchange m1 and m2 */
        i = m1;
        m1 = m2;
        m2 = i;

        /* check for stopping calculation depending on termination method */
        if (options->termination == TERM_PREC)
        {
            if (maxResiduum < options->term_precision)
            {
                term_iteration = 0;
            }
        }
        else if (options->termination == TERM_ITER)
        {
            term_iteration--;
        }
		#ifdef OPENMP
        #pragma omp barrier
		#endif
    }

    results->m = m2;
}

#ifdef OPENMP
static void calculate_omp(struct calculation_arguments const *arguments,
                          struct calculation_results *results,
                          struct options const *options)
{
    int const N = arguments->N;
    int const n_threads = options->number;
    omp_set_num_threads(n_threads);
    #pragma omp parallel default(none) shared(arguments, results, N, options, n_threads)
    {
        int k = omp_get_thread_num();
		#ifdef SPALTE
    	int cells_per_thread = (N - 1) / n_threads;
		int first_column = k * cells_per_thread + 1;
		int last_column = (k == n_threads - 1) ? N - 1
											: (k + 1) * cells_per_thread;
		calculate_seq(arguments, results, options, 1, N - 1, first_column, last_column, 0); /* go over every row of specific columns */
		#elif ELEMENT
		int inner_elements = (N - 2) * (N - 2);
		int interval = inner_elements / n_threads;
		int chunksize;
		int first_index;
		int rest = inner_elements % n_threads;
		if (k < rest)
		{
			first_index = k * (interval + 1);
			chunksize = interval + 1;
		}
		else
		{
			first_index = k * interval + rest;
			chunksize = interval;
		}
		int last_index = first_index + chunksize - 1;
		int first_row = first_index % (N - 1);
		int last_row = last_index % (N - 1);
		int first_column = first_index / (N - 1);
		int last_column = last_index / (N - 1);
		calculate_seq(arguments, results, options, first_row, last_row, first_column, last_column, 1); /* go over a range of elements */
		#else /* ZEILE */
    	int cells_per_thread = (N - 1) / n_threads;
		int first_row = k * cells_per_thread + 1;
		int last_row = (k == n_threads - 1) ? N - 1
											: (k + 1) * cells_per_thread;
		calculate_seq(arguments, results, options, first_row, last_row, 1, N - 1, 0); /* go over every column of specific rows */
		#endif
    }
}
#endif

static void calculate(struct calculation_arguments const *arguments,
                      struct calculation_results *results,
                      struct options const *options)
{
	#ifdef OPENMP
    /* when using Gauss-Seidel or 1 thread, we still use the sequential version without openmp. */
    if (options->number == 1 || options->method == METH_GAUSS_SEIDEL)
    {
        calculate_seq(arguments, results, options, 1, arguments->N - 1, 1, arguments-> N - 1, 0);
    }
    /* otherwise, use the parallel version */
    else
    {
        calculate_omp(arguments, results, options);
    }
	#else
	calculate_seq(arguments, results, options, 1, arguments->N - 1, 1, arguments-> N - 1, 0);
	#endif
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static void
displayStatistics(struct calculation_arguments const *arguments, struct calculation_results const *results, struct options const *options)
{
    int N = arguments->N;
    double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

    printf("Berechnungszeit:    %f s \n", time);
    printf("Speicherbedarf:     %f MiB\n", (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
    printf("Berechnungsmethode: ");

    if (options->method == METH_GAUSS_SEIDEL)
    {
        printf("Gauß-Seidel");
    }
    else if (options->method == METH_JACOBI)
    {
        printf("Jacobi");
    }

    printf("\n");
    printf("Interlines:         %" PRIu64 "\n", options->interlines);
    printf("Stoerfunktion:      ");

    if (options->inf_func == FUNC_F0)
    {
        printf("f(x,y) = 0");
    }
    else if (options->inf_func == FUNC_FPISIN)
    {
        printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
    }

    printf("\n");
    printf("Terminierung:       ");

    if (options->termination == TERM_PREC)
    {
        printf("Hinreichende Genaugkeit");
    }
    else if (options->termination == TERM_ITER)
    {
        printf("Anzahl der Iterationen");
    }

    printf("\n");
	long unsigned int iterations =  /* iterations need to be divided because each thread increments this value */
	#ifdef OPENMP
	results->stat_iteration / options->number;
	#else
	results->stat_iteration;
	#endif
    printf("Anzahl Iterationen: %" PRIu64 "\n", iterations);
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
static void
displayMatrix(struct calculation_arguments *arguments, struct calculation_results *results, struct options *options)
{
    int x, y;

    double **Matrix = arguments->Matrix[results->m];

    int const interlines = options->interlines;

    printf("Matrix:\n");

    for (y = 0; y < 9; y++)
    {
        for (x = 0; x < 9; x++)
        {
            printf("%7.4f", Matrix[y * (interlines + 1)][x * (interlines + 1)]);
        }

        printf("\n");
    }

    fflush(stdout);
}

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int main(int argc, char **argv)
{
    struct options options;
    struct calculation_arguments arguments;
    struct calculation_results results;

    askParams(&options, argc, argv);

    initVariables(&arguments, &results, &options);

    allocateMatrices(&arguments);
    initMatrices(&arguments, &options);

    gettimeofday(&start_time, NULL);
    calculate(&arguments, &results, &options);
    gettimeofday(&comp_time, NULL);

    displayStatistics(&arguments, &results, &options);
    displayMatrix(&arguments, &results, &options);

    freeMatrices(&arguments);

    return 0;
}
