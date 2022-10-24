/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                TU Muenchen - Institut fuer Informatik                  **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff-seq.c                                              **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauss-Seidel and   **/
/**            Jacobi method.                                              **/
/**                                                                        **/
/** $Id: partdiff-seq.c,v 2.4 2008-04-17 14:34:21 kunkel Exp $             **/
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
#include <time.h>
#include "partdiff-seq.h"


/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */
/* input variables */
int method = 0;			/* Gauss Seidel or Jacobi method of iteration     */
int interlines = 0;		/* matrix size = interlines*8+9                   */
int inf_func = 0;		/* inference function                             */
int termination = 0;		/* termination condition                          */
int term_iteration = 0;		/* terminate if iteration number reached          */
double term_precision = 0;	/* terminate if precision reached                 */

/* numerical calculation variables */
int N = 0;			/* number of spaces between lines (lines=N+1)     */
int m1 = 0, m2 = 0;		/* used as indices for old and new matrices       */
int stat_iteration = 0;		/* number of current iteration                    */
double ***Matrix;		/* index matrix used for addressing M             */
double *M;			/* two matrices with real values                  */
double star = 0;		/* four times center value minus 4 neigh.b values */
double h = 0;			/* length of a space between two lines            */
double residuum = 0;		/* residuum of current iteration                  */
double korrektur = 0;		/* necessary correction value for star            */
double maxresiduum = 0;		/* maximum residuum value of a slave in iteration */
double stat_precision = 0;	/* actual precision of all slaves in iteration    */

/* time measurement variables */
time_t start_time;		/* time when program started                      */
time_t comp_time;		/* time when calculation completed                */


/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
void
initVariables (void)
{
  N = interlines * 8 + 9 - 1;
  h = (float) (((float) (1)) / (N));
}


/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
void
allocateMatrices (void)
{
  int i, j;			/* local variables */

  Matrix = (double ***) calloc (2, sizeof (double **));	/* allocate memory */
  if (Matrix == 0)
    {
      errorQuit ();
    }				/* quit if error   */

  Matrix[0] = (double **) calloc ((N + 1), sizeof (double *));	/* allocate memory */
  if (Matrix[0] == 0)
    {
      errorQuit ();
    }				/* quit if error   */

  Matrix[1] = (double **) calloc ((N + 1), sizeof (double *));	/* allocate memory */
  if (Matrix[1] == 0)
    {
      errorQuit ();
    }				/* quit if error   */

  M = malloc (sizeof (double) * (N + 1) * (N - 1) * 2);	/* allocate memory */
  if (M == 0)
    {
      errorQuit ();
    }				/* quit if error   */

  for (i = 0; i <= 1; i++)
    for (j = 0; j <= N; j++)
      Matrix[i][j] = (double *) (M + (i * (N + 1) * (N + 1)) + (j * (N + 1)));
}


/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
void
initMatrices (void)
{
  int g, i, j;			/*  local variables for loops   */

  if (method == GAUSS_SEIDEL)	/* **************************** */
    {
      m1 = 0;
      m2 = 0;
    }				/*  initialize m1 and m2        */
  else				/*  depending on algorithm       */
    {
      m1 = 0;
      m2 = 1;
    }				/* **************************** */

  for (g = 0; g <= 1; g++)	/* **************************** */
    {
      for (i = 0; i <= N; i++)	/*  initialize matrix/matrices  */
	{
	  for (j = 0; j <= N; j++)	/*  with zeros                  */
	    {
	      Matrix[g][i][j] = 0;
	    }
	}
    }				/* **************************** */

  if (inf_func == FUNC_1)	/* ******************************* */
    {
      for (i = 0; i <= N; i++)	/*  initialize borders, depending  */
	{
	  Matrix[0][i][0] = 1 - (h * i);	/*  on function                    */
	  Matrix[0][i][N] = h * i;	/*  (function 2: nothing to do)    */
	  Matrix[0][0][i] = 1 - (h * i);	/* ******************************* */
	  Matrix[0][N][i] = h * i;
	  Matrix[1][i][0] = 1 - (h * i);
	  Matrix[1][i][N] = h * i;
	  Matrix[1][0][i] = 1 - (h * i);
	  Matrix[1][N][i] = h * i;
	}
      Matrix[0][N][0] = 0;
      Matrix[0][0][N] = 0;
      Matrix[1][N][0] = 0;
      Matrix[1][0][N] = 0;
    }
}


/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
void
freeMatrices (void)
{
  free (Matrix);
  if (Matrix[1] != 0)
    free (Matrix[1]);
  if (Matrix[0] != 0)
    free (Matrix[0]);
}


/* ************************************************************************ */
/* getResiduum: calculates residuum                                         */
/* Input: x,y - actual column and row                                       */
/* ************************************************************************ */
double
getResiduum (int x, int y)
{
  if (inf_func == FUNC_1)	/* ************************ */
    return ((-1 * star) / 4.0);	/*  calculate the residuum  */
  else				/*  depending on function   */
    return ((TWO_PI_SQUARE *	/* ************************ */
	     sin ((double) (y) * PI * h) *
	     sin ((double) (x) * PI * h) * h * h - star) / 4.0);
}


/* ************************************************************************ */
/* checkQuit: exchanges matrices and decides, whether calculation is         */
/*            ready or not.                                                 */
/* ************************************************************************ */
void
checkQuit (void)
{
  int i;			/*  local variable          */

  i = m1;
  m1 = m2;
  m2 = i;			/*  exchange m1 and m2      */

  if (termination == TERM_1)	/* ************************ */
    if (maxresiduum < term_precision)	/*  check for stopping       */
      {
	term_iteration = 0;
      }				/*  calculation, depending  */
  if (termination == TERM_2)	/*  on termination method   */
    {
      term_iteration = term_iteration - 1;
    }				/* ************************ */
}


/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
void
calculate (void)
{
  int i, j;			/* local variables for loops  */

  while (term_iteration > 0)
    {
      maxresiduum = 0;
      for (j = 1; j < N; j++)	/* over all columns     */
	{			/*                   */
	  for (i = 1; i < N; i++)	/* over all rows  */
	    {
	      star = -Matrix[m2][i - 1][j]
		- Matrix[j - 1][m2][i] + 4 * Matrix[m2][i][j] -
		Matrix[m2][i][j + 1] - Matrix[m2][i + 1][j];

	      residuum = getResiduum (i, j);
	      korrektur = residuum;
	      if (residuum < 0)
		residuum = residuum * (-1);
	      if (residuum >= maxresiduum)
		maxresiduum = residuum;
	      Matrix[m1][i][j] = Matrix[m2][i][j] + korrektur;
	    }
	}
      stat_iteration = stat_iteration + 1;
      stat_precision = maxresiduum;
      checkQuit ();
    }
}


/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
void
displayStatistics (void)
{
  printf ("\n\nErgebnis:\n\n");
  printf ("Berechnungszeit:    ");
  printf ("%d", (int) (comp_time - start_time));
  printf (" Sekunden\n");
  printf ("Berechnungsmethode: ");
  if (method == GAUSS_SEIDEL)
    printf ("Gauss-Seidel");
  if (method == JACOBI)
    printf ("Jacobi");
  printf ("\n");
  printf ("Interlines:         %d\n", interlines);
  printf ("Stoerfunktion:      ");
  if (inf_func == FUNC_1)
    printf ("f(x,y)=0");
  if (inf_func == FUNC_2)
    printf ("f(x,y)=2pi^2*sin(pi*x)sin(pi*y)");
  printf ("\n");
  printf ("Terminierung:       ");
  if (termination == TERM_1)
    printf ("Hinreichende Genaugkeit");
  if (termination == TERM_2)
    printf ("Anzahl der Iterationen");
  printf ("\n");
  printf ("Anzahl Iterationen: %d\n", stat_iteration);
  printf ("Norm des Fehlers:   %e\n", stat_precision);
}


/* ************************************************************************ */
/* errorQuit ()                                                             */
/* frees memory for matrices and quits if there was a memory allocation-    */
/* problem                                                                  */
/* ************************************************************************ */
int
errorQuit (void)
{
  freeMatrices ();
  printf ("\n\nSpeicherprobleme!\n");
  exit (1);			/* exit program */
}



/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argC, char **argV)
{
  AskParams (&method, &interlines, &inf_func,	/* ************************* */
	     &termination, &term_precision,	/*  get parameters           */
	     &term_iteration, argC, argV);	/* ************************* */

  initVariables ();		/* ******************************************* */
  allocateMatrices ();		/*  get and initialize variables and matrices  */
  initMatrices ();		/* ******************************************* */

  start_time = time (NULL);	/*  start timer         */
  calculate ();			/*  solve the equation  */
  comp_time = time (NULL);	/*  stop timer          */

  displayStatistics ();		/* **************** */
  DisplayMatrix ("Matrix: ",	/*  display some    */
		 Matrix[m2][0], interlines);	/*  statistics and  */
  freeMatrices ();		/*  free memory     */
  exit (0);			/* **************** */
}
