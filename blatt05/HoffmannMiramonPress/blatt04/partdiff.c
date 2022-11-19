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
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

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
static void initVariables(struct calculation_arguments *arguments,
                          struct calculation_results *results,
                          struct options const *options)
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
static void freeMatrices(struct calculation_arguments *arguments)
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
static void *allocateMemory(size_t size)
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
static void allocateMatrices(struct calculation_arguments *arguments)
{
    uint64_t i, j;

    uint64_t const N = arguments->N;

    arguments->M =
        allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double));
    arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double **));

    for (i = 0; i < arguments->num_matrices; i++)
    {
        arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double *));

        for (j = 0; j <= N; j++)
        {
            arguments->Matrix[i][j] =
                arguments->M + (i * (N + 1) * (N + 1)) + (j * (N + 1));
        }
    }
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static void initMatrices(struct calculation_arguments *arguments,
                         struct options const *options)
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
                Matrix[g][i][0] = 1.0 - (h * i);
                Matrix[g][i][N] = h * i;
                Matrix[g][0][i] = 1.0 - (h * i);
                Matrix[g][N][i] = h * i;
            }

            Matrix[g][N][0] = 0.0;
            Matrix[g][0][N] = 0.0;
        }
    }
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */ 
/* basically unchanged from the original version, only some logic for the   */
/* barrier and mutexes is added                                             */
/* ************************************************************************ */
static void calculate_seq(struct calculation_arguments const *arguments,
                          struct calculation_results *results,
                          struct options const options,
                          unsigned int i_start,
                          unsigned int i_end,
                          pthread_barrier_t *barrier,
                          pthread_mutex_t *mutex)
{
    unsigned int i, j;  /* local variables for loops */
    int m1, m2;         /* used as indices for old and new matrices */
    double star;        /* four times center value minus 4 neigh.b values */
    double residuum;    /* residuum of current iteration */
    double maxResiduum; /* maximum residuum value of a slave in iteration */

    unsigned int const N = arguments->N;
    double const h = arguments->h;

    double pih = 0.0;
    double fpisin = 0.0;

    int term_iteration = options.term_iteration;

    /* initialize m1 and m2 depending on algorithm */
    if (options.method == METH_JACOBI)
    {
        m1 = 0;
        m2 = 1;
    }
    else
    {
        m1 = 0;
        m2 = 0;
    }

    if (options.inf_func == FUNC_FPISIN)
    {
        pih = PI * h;
        fpisin = 0.25 * TWO_PI_SQUARE * h * h;
    }

    while (term_iteration > 0)
    {
        double **Matrix_Out = arguments->Matrix[m1];
        double **Matrix_In = arguments->Matrix[m2];

        maxResiduum = 0;

        /* over specific rows */
        for (i = i_start; i <= i_end; i++)
        {
            double fpisin_i = 0.0;

            if (options.inf_func == FUNC_FPISIN)
            {
                fpisin_i = fpisin * sin(pih * (double)i);
            }

            /* over all columns */
            for (j = 1; j < N; j++)
            {
                star = 0.25 * (Matrix_In[i - 1][j] + Matrix_In[i][j - 1] +
                               Matrix_In[i][j + 1] + Matrix_In[i + 1][j]);

                if (options.inf_func == FUNC_FPISIN)
                {
                    star += fpisin_i * sin(pih * (double)j);
                }

                if (options.termination == TERM_PREC || term_iteration == 1)
                {
                    residuum = Matrix_In[i][j] - star;
                    residuum = (residuum < 0) ? -residuum : residuum;
                    maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
                }

                Matrix_Out[i][j] = star;
            }
        }

        /* here we change some global variables so we need to use a mutex */
        if (mutex)
        {
            pthread_mutex_lock(mutex);
        }
        results->stat_iteration++;
        results->stat_precision = (results->stat_precision < maxResiduum) ? maxResiduum : results->stat_precision;
        if (mutex)
        {
            pthread_mutex_unlock(mutex);
        }

        /* exchange m1 and m2 */
        i = m1;
        m1 = m2;
        m2 = i;

        /* check for stopping calculation depending on termination method */
        if (options.termination == TERM_PREC)
        {
            if (maxResiduum < options.term_precision)
            {
                term_iteration = 0;
            }
        }
        else if (options.termination == TERM_ITER)
        {
            term_iteration--;
        }
        if (barrier)
        {
            pthread_barrier_wait(barrier);  /* wait for all other threads */
        }
    }
    results->m = m2;
}

/* this function is executed by every thread */
void *posix_fun(void *_args)
{
    struct calc_posix_args *args = (struct calc_posix_args *)_args;
    const struct calculation_arguments *arguments = args->arguments;
    const struct options *options = args->opts;
    struct calculation_results *results = args->results;
    pthread_barrier_t *barrier = args->barrier;
    pthread_mutex_t *mutex = args->mutex;
    calculate_seq(arguments, results, *options, args->first_i, args->last_i, barrier, mutex);

    return NULL;
}

static void calculate_posix(struct calculation_arguments const *arguments,
                            struct calculation_results *results,
                            struct options const *options)
{
    unsigned int const N = arguments->N;
    pthread_t threads[options->number];
    pthread_barrier_t barrier;
    pthread_mutex_t mutex;

    unsigned int cells_per_thread = (N - 1) / options->number;
    struct calc_posix_args args[options->number];
    pthread_barrier_init(&barrier, NULL, options->number);                          /* barrier for waiting at the end of every while loop iteration */
    pthread_mutex_init(&mutex, NULL);                                               /* mutex for waiting before changing "global" variables         */

    /* here we create the necessary arguments and the threads */
    for (unsigned int k = 0; k < options->number; k++)
    {
        args[k].first_i = k * cells_per_thread + 1;                                 /* first line index to calculate */
        args[k].last_i = (k == options->number - 1) ? N - 1
                                                    : (k + 1) * cells_per_thread;   /* last line index to calculate  */
        args[k].arguments = arguments;
        args[k].opts = options;
        args[k].barrier = &barrier;
        args[k].results = results;
        args[k].mutex = &mutex;
        pthread_create(&threads[k], NULL, &posix_fun, &args[k]);
    }

    for (unsigned int k = 0; k < options->number; k++)
    {
        pthread_join(threads[k], NULL);                 /* we wait for the threads to stop */
    }
}

static void calculate(struct calculation_arguments const *arguments,
                      struct calculation_results *results,
                      struct options const *options)
{
    /* when using Gauss-Seidel or 1 thread, we use the sequential version without mutex or barrier. */
    if (options->number == 1 || options->method == METH_GAUSS_SEIDEL)
    {
        calculate_seq(arguments, results, *options, 1, arguments->N - 1, NULL, NULL);
    }
    /* otherwise, use the parallel version */
    else
    {
        calculate_posix(arguments, results, options);
    }
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static void displayStatistics(struct calculation_arguments const *arguments,
                              struct calculation_results const *results,
                              struct options const *options)
{
    int N = arguments->N;
    double time = (comp_time.tv_sec - start_time.tv_sec) +
                  (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

    printf("Berechnungszeit:    %f s \n", time);
    printf("Speicherbedarf:     %f MiB\n", (N + 1) * (N + 1) * sizeof(double) *
                                               arguments->num_matrices / 1024.0 /
                                               1024.0);
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
    printf("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration / options->number); /* iterations need to be divided because each thread increments this value */
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
                          struct options *options)
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
