/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**  	      	   TU Muenchen - Institut fuer Informatik                  **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Thomas A. Zochler, Andreas C. Schmidt                       **/
/**                                                                        **/
/** File:      partdiff.h                                                  **/
/**                                                                        **/
/****************************************************************************/
/****************************************************************************/

/* *********************************** */
/* Include some standard header files. */
/* *********************************** */
#include <math.h>
#include <stdint.h>
#include <pthread.h>

/* ************* */
/* Some defines. */
/* ************* */
#ifndef PI
#define PI 3.141592653589793
#endif
#define TWO_PI_SQUARE (2 * PI * PI)
#define MAX_INTERLINES 10240
#define MAX_ITERATION 200000
#define MAX_THREADS 1024
#define METH_GAUSS_SEIDEL 1
#define METH_JACOBI 2
#define FUNC_F0 1
#define FUNC_FPISIN 2
#define TERM_PREC 1
#define TERM_ITER 2

struct options
{
	uint64_t number;		 /* Number of threads                              */
	uint64_t method;		 /* Gauss Seidel or Jacobi method of iteration     */
	uint64_t interlines;	 /* matrix size = interlines*8+9                   */
	uint64_t inf_func;		 /* inference function                             */
	uint64_t termination;	 /* termination condition                          */
	uint64_t term_iteration; /* terminate if iteration number reached          */
	double term_precision;	 /* terminate if precision reached                 */
};

struct calc_posix_args
{
	unsigned int first_i;						   /* first line index to calculate              */
	unsigned int last_i;						   /* last line index to calculate               */
	const struct calculation_arguments *arguments; /* calculation arguments                      */
	const struct options *opts;					   /* options, see above                         */
	struct calculation_results *results;		   /* calculation results                        */
	pthread_barrier_t *barrier;					   /* barrier for waiting on other threads       */
	pthread_mutex_t *mutex;						   /* mutex for safely changing global variables */
};

/* *************************** */
/* Some function declarations. */
/* *************************** */
/* Documentation in files      */
/* - askparams.c               */
/* - displaymatrix.c           */
/* *************************** */
void askParams(struct options *, int, char **);
