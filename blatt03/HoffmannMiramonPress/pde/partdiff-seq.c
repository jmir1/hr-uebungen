/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                TU München - Institut fuer Informatik                   **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff-seq.c                                              **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauss-Seidel and   **/
/**            Jacobi methods.                                             **/
/**                                                                        **/
/****************************************************************************/
/****************************************************************************/

/* ************************************************************************ */
/* Include standard header file.                                            */
/* ************************************************************************ */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include "partdiff-seq.h"

struct calculation_arguments
{
	int     N;              /* number of spaces between lines (lines=N+1)     */
	int     num_matrices;   /* number of matrices                             */
	double  ***Matrix;      /* index matrix used for addressing M             */
	double  *M;             /* two matrices with real values                  */
	double  h;              /* length of a space between two lines            */
};

struct calculation_results
{
	int     m;
	int     stat_iteration; /* number of current iteration                    */
	double  stat_precision; /* actual precision of all slaves in iteration    */
};

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time;       /* time when program started                      */
struct timeval comp_time;        /* time when calculation completed                */


/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static
void
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
{
	arguments->N = options->interlines * 8 + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = (float)( ( (float)(1) ) / (arguments->N));

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;
}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static
void
freeMatrices (struct calculation_arguments* arguments)
{
	int i;

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
static
void*
allocateMemory (size_t size)
{
	void *p;

	if ((p = malloc(size)) == NULL)
	{
		printf("\n\nSpeicherprobleme!\n");
		/* exit program */
		exit(1);
	}

	return p;
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static
void
allocateMatrices (struct calculation_arguments* arguments)
{
	int i, j;

	int N = arguments->N;

	arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double*));

		for (j = 0; j <= N; j++)
		{
			arguments->Matrix[i][j] = (double*)(arguments->M + (i * (N + 1) * (N + 1)) + (j * (N + 1)));
		}
	}
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices (struct calculation_arguments* arguments, struct options* options)
{
	int g, i, j;                                /*  local variables for loops   */

	int N = arguments->N;
	double h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i <= N; i++)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[g][i][j] = 0;
			}
		}
	}

	/* initialize borders, depending on function (function 2: nothing to do) */
	if (options->inf_func == FUNC_F0)
	{
		for(i = 0; i < N; i++)
		{
			for (j = 0; j < arguments->num_matrices; j++)
			{
				Matrix[j][i][0] = 1 + (1 - (h * i)); // Linke Kante
				Matrix[j][N][i] = 1 - (h * i); // Untere Kante
				Matrix[j][N - i][N] = h * i; // Rechte Kante
				Matrix[j][0][N - i] = 1 + h * i; // Obere Kante
			}
		}
	}
}

/* ************************************************************************ */
/* getResiduum: calculates residuum                                         */
/* Input: x,y - actual column and row                                       */
/* ************************************************************************ */
double
getResiduum (struct calculation_arguments* arguments, struct options* options, int x, int y, double star)
{
	if (options->inf_func == FUNC_F0)
	{
		return ((-star) / 4.0);
	}
	else
	{
		return ((TWO_PI_SQUARE * sin((double)(y) * PI * arguments->h) * sin((double)(x) * PI * arguments->h) * arguments->h * arguments->h - star) / 4.0);
	}
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate (struct calculation_arguments* arguments, struct calculation_results *results, struct options* options)
{
	int i, j;                                   /* local variables for loops  */
	int m1, m2;                                 /* used as indices for old and new matrices       */
	double star;                                /* four times center value minus 4 neigh.b values */
	double korrektur;
	double residuum;                            /* residuum of current iteration                  */

	int N = arguments->N;
	double*** Matrix = arguments->Matrix;

	/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_GAUSS_SEIDEL)
	{
		m1=0; m2=0;
	}
	else
	{
		m1=0; m2=1;
	}

	while (options->term_iteration > 0)
	{
		results->stat_precision = 0;

		/* over all rows */
		for (j = 1; j < N; j++)
		{
			/* over all columns */
			for (i = 1; i < N; i++)
			{
				star = -Matrix[m2][i-1][j] - Matrix[m2][i][j-1] - Matrix[m2][i][j+1] - Matrix[m2][i+1][j] + 4.0 * Matrix[m2][i][j];

				residuum = getResiduum(arguments, options, i, j, star);
				korrektur = residuum;
				residuum = (residuum < 0) ? -residuum : residuum;
				results->stat_precision = (residuum < results->stat_precision) ? results->stat_precision : residuum;

				Matrix[m1][i][j] = Matrix[m2][i][j] + korrektur;
			}
		}

		results->stat_iteration++;

		/* exchange m1 and m2 */
		i=m1; m1=m2; m2=i;

		/* check for stopping calculation, depending on termination method */
		if (options->termination == TERM_PREC)
		{
			if (results->stat_precision < options->term_precision)
			{
				options->term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			options->term_iteration--;
		}
	}

	results->m = m2;
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static
void
displayStatistics (struct calculation_arguments* arguments, struct calculation_results *results, struct options* options)
{
	(void)arguments;

	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;
	printf("Berechnungszeit:    %f s \n", time);

	printf("Berechnungsmethode: ");

	if (options->method == METH_GAUSS_SEIDEL)
	{
		printf("Gauss-Seidel");
	}
	else if (options->method == METH_JACOBI)
	{
		printf("Jacobi");
	}

	printf("\n");
	printf("Interlines:         %d\n",options->interlines);
	printf("Stoerfunktion:      ");

	if (options->inf_func == FUNC_F0)
	{
		printf("f(x,y)=0");
	}
	else if (options->inf_func == FUNC_FPISIN)
	{
		printf("f(x,y)=2pi^2*sin(pi*x)sin(pi*y)");
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
	printf("Anzahl Iterationen: %d\n", results->stat_iteration);
	printf("Norm des Fehlers:   %e\n", results->stat_precision);
}

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{
	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;

	/* get parameters */
	AskParams(&options, argc, argv);              /* ************************* */

	initVariables(&arguments, &results, &options);           /* ******************************************* */

	allocateMatrices(&arguments);        /*  get and initialize variables and matrices  */
	initMatrices(&arguments, &options);            /* ******************************************* */

	gettimeofday(&start_time, NULL);                   /*  start timer         */
	calculate(&arguments, &results, &options);                                      /*  solve the equation  */
	gettimeofday(&comp_time, NULL);                   /*  stop timer          */

	displayStatistics(&arguments, &results, &options);                                  /* **************** */
	DisplayMatrix("Matrix:",                              /*  display some    */
			arguments.Matrix[results.m][0], options.interlines);            /*  statistics and  */

	freeMatrices(&arguments);                                       /*  free memory     */

	return 0;
}
