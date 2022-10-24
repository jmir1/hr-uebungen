/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**  	      	   TU Muenchen - Institut fuer Informatik                  **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Thomas A. Zochler, Andreas C. Schmidt                       **/
/**                                                                        **/
/** File:      partdiff-seq.h                                              **/
/**                                                                        **/
/****************************************************************************/
/****************************************************************************/

/* *********************************** */
/* Include some standard header files. */
/* *********************************** */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <malloc.h>

/* ************* */
/* Some defines. */
/* ************* */
#ifndef PI
#define PI 			3.141592653589793
#endif
#ifndef FALSE
#define FALSE			0
#define TRUE			1
#endif
#define TWO_PI_SQUARE 		(2*PI*PI)
#define MAX_ITERATION  		200000
#define GAUSS_SEIDEL 		1
#define JACOBI 			2
#define FUNC_1			1
#define FUNC_2			2
#define TERM_1			1
#define TERM_2			2

/* *************************** */
/* Some function declarations. */
/* *************************** */
/* Documentation in files      */
/* - askparams.c               */
/* - displaymatrix.c           */
/* *************************** */
void AskParams (int *method,
		int *interlines,
		int *func,
		int *termination,
		double *term_precision,
		int *term_iteration, int argC, char **argV);

void DisplayMatrix (char *s, double *v, int interlines);

void DisplayMatrixAddr (char *s, double ***v, int interlines, int matrixnum);

void initVariables (void);

void allocateMatrices (void);

void initMatrices (void);

void freeMatrices (void);

double getResiduum (int, int);

void checkQuit (void);

void calculate (void);

void displayStatistics (void);

int errorQuit (void);
