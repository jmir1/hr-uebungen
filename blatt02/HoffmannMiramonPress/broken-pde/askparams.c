/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**  	      	   TU Muenchen - Institut fuer Informatik                  **/
/**                                                                        **/
/** Copyright: Dr. Thomas Ludwig                                           **/
/**            Thomas A. Zochler                                           **/
/**                                                                        **/
/** File:      askparams.c                                                 **/
/** Version 23                                                             **/
/****************************************************************************/
/****************************************************************************/

/****************************************************************************/
/** Beschreibung der Funktion AskParams():                                 **/
/**                                                                        **/
/** Die Funktion AskParams liest sechs Parameter (Erkl"arung siehe unten)  **/
/** entweder von der Standardeingabe oder von Kommandozeilenoptionen ein.  **/
/**                                                                        **/
/** Ziel dieser Funktion ist es, die Eingabe der Parameter sowohl inter-   **/
/** aktiv als auch als Kommandozeilenparameter zu erm"oglichen.            **/
/**                                                                        **/
/** F"ur die Parameter argC und argV k"onnen direkt die vom System         **/
/** gelieferten Variablen der Funktion main verwendet werden.              **/
/**                                                                        **/
/** Beispiel:                                                              **/
/**                                                                        **/
/** int main ( int argC, char **argV )                                     **/
/** {                                                                      **/
/**   ...                                                                  **/
/**   AskParams( ..., argC, argV );                                        **/
/**   ...                                                                  **/
/** }                                                                      **/
/**                                                                        **/
/** Dabei wird argV[0] ignoriert und weiter eingegebene Parameter der      **/
/** Reihe nach verwendet.                                                  **/
/**                                                                        **/
/** Falls bei Aufruf von AskParams() argC < 2 "ubergeben wird, werden      **/
/** die Parameter statt dessen von der Standardeingabe gelesen.            **/
/****************************************************************************/
/** int *method;                                                           **/
/**         Bezeichnet das bei der L"osung der Poissongleichung zu         **/
/**         verwendende Verfahren ( Gauss-Seidel oder Jacobi ).            **/
/** Werte:  GAUSS_SEIDEL  oder JACOBI (definierte Konstanten)              **/
/****************************************************************************/
/** int *interlines:                                                       **/
/**         Gibt die Zwischenzeilen zwischen den auszugebenden             **/
/**         neun Zeilen an. Die Gesamtanzahl der Zeilen ergibt sich als    **/
/**         lines = 8 * (*interlines) + 9. Diese Art der Berechnung der    **/
/**         Problemgr"o"se (auf dem Aufgabenblatt mit N bezeichnet)        **/
/**         wird benutzt, um mittels der Ausgaberoutine DisplayMatrix()    **/
/**         immer eine "ubersichtliche Ausgabe zu erhalten.                **/
/** Werte:  0 < *interlines                                                **/
/****************************************************************************/
/** int *func:                                                             **/
/**         Bezeichnet die St"orfunktion ( I oder II ) und damit auch      **/
/**         die Randbedingungen.                                           **/
/** Werte:  FUNC_1: f(x,y)=0, 0<x<1, 0<y<1                                 **/
/**         FUNC_2: f(x,y)=2pi^2*sin(pi*x)sin(pi*y), 0<x<1, 0<y<1          **/
/****************************************************************************/
/** int *termination:                                                      **/
/**         Gibt die Art der Abbruchbedingung an.                          **/
/** Werte:  TERM_1: Abbruchbedingung ist die Genauigkeit der bereits       **/
/**                 berechneten N"aherung. Diese soll unter die            **/
/**                 Grenze term_precision kommen.                          **/
/**         TERM_2: Abbruchbedingung ist die Anzahl der Iterationen. Diese **/
/**                 soll gr"o"ser als term_iteration sein.                 **/
/****************************************************************************/
/** double *term_precision:                                                **/
/** int t*erm_iteration:                                                   **/
/**         Es wird jeweils nur einer der beiden Parameter f"ur die        **/
/**         Abbruchbedingung eingelesen.                                   **/
/****************************************************************************/

#include "partdiff-seq.h"
#include <string.h>

void
AskParams (int *method,
	   int *interlines,
	   int *func,
	   int *termination,
	   double *term_precision, int *term_iteration, int argC, char **argV)
{
  printf ("\n");
  printf ("============================================================\n");
  printf ("Program for calculation of partial differential equations.  \n");
  printf ("============================================================\n");
  printf ("(c) Dr. Thomas Ludwig, TU M\"unchen.\n");
  printf ("    Thomas A. Zochler, TU M\"unchen.\n");
  printf ("    Andreas C. Schmidt, TU M\"unchen.\n");
  printf ("============================================================\n");

  if (argC < 2)
    {
      /* ----------------------------------------------- */
      /* Get input: method, interlines, func, precision. */
      /* ----------------------------------------------- */
      do
	{
	  printf ("\n");
	  printf ("Select calculation method:\n");
	  printf ("  %1d: Gauss-Seidel.\n", GAUSS_SEIDEL);
	  printf ("  %1d: Jacobi.\n", JACOBI);
	  printf ("method> ");
	  fflush (stdout);
	  scanf ("%d", method);
	}
      while ((*method < GAUSS_SEIDEL) || (*method > JACOBI));
      do
	{
	  printf ("\n");
	  printf ("Matrixsize = Interlines*8+9\n");
	  printf ("Interlines> ");
	  fflush (stdout);
	  scanf ("%d", interlines);
	}
      while ((*interlines < 0) || (*interlines > 1000));
      do
	{
	  printf ("\n");
	  printf ("Select interference function:\n");
	  printf (" %1d: f(x,y)=0.\n", FUNC_1);
	  printf (" %1d: f(x,y)=2pi^2*sin(pi*x)sin(pi*y).\n", FUNC_2);
	  printf ("interference function> ");
	  fflush (stdout);
	  scanf ("%d", func);
	}
      while ((*func < FUNC_1) || (*func > FUNC_2));
      do
	{
	  printf ("\n");
	  printf ("Select termination:\n");
	  printf (" %1d: sufficient precision.\n", TERM_1);
	  printf (" %1d: number of iterations.\n", TERM_2);
	  printf ("termination> ");
	  fflush (stdout);
	  scanf ("%d", termination);
	}
      while ((*termination < TERM_1) || (*termination > TERM_2));
      if (*termination == TERM_1)
	{
	  do
	    {
	      printf ("\n");
	      printf ("Select precision:\n");
	      printf ("  Range: 1e-4 .. 1e-20.\n");
	      printf ("precision> ");
	      fflush (stdout);
	      scanf ("%lf", term_precision);
	    }
	  while ((*term_precision < 1e-20) || (*term_precision > 1e-4));
	  *term_iteration = MAX_ITERATION;
	}
      if (*termination == TERM_2)
	{
	  do
	    {
	      printf ("\n");
	      printf ("Select number of iterations:\n");
	      printf ("  Range: 1 .. %d.\n", MAX_ITERATION);
	      printf ("Iterations> ");
	      fflush (stdout);
	      scanf ("%d", term_iteration);
	    }
	  while ((*term_iteration < 1) || (*term_iteration > MAX_ITERATION));
	  *term_precision = 0;
	}
    }
  else
    {
      if (strcmp ("help", (char *) argV[1]) == 0
	  || strcmp ("-h", (char *) argV[1]) == 0
	  || strcmp ("-?", (char *) argV[1]) == 0
	  || strcmp ("?", (char *) argV[1]) == 0 || argC < 6)
	{
	  printf ("\nUsage:\n");
	  printf
	    ("partdiff [num] [method] [lines] [func] [term] [prec/iter]\n");
	  printf ("  - num:    number of process/threads to use\n");
	  printf ("  - method: %1d: Gauss-Seidel.\n", GAUSS_SEIDEL);
	  printf ("            %1d: Jacobi.\n", JACOBI);
	  printf
	    ("  - lines:  (lines=interlines) matrixsize = interlines*8+9\n");
	  printf ("  - func:   %1d: f(x,y)=0.\n", FUNC_1);
	  printf ("            %1d: f(x,y)=2pi^2*sin(pi*x)sin(pi*y).\n",
		  FUNC_2);
	  printf ("  - term:   %1d: sufficient precision.\n", TERM_1);
	  printf ("            %1d: number of iterations.\n", TERM_2);
	  printf ("  - prec/iter: depending on term:\n");
	  printf ("            precision:  Range: 1e-4 .. 1e-20.\n");
	  printf ("            iterations: Range: 1 .. %d.\n", MAX_ITERATION);
	  printf ("\n");
	  exit (0);
	}
      sscanf (argV[2], "%d", method);
      sscanf (argV[3], "%d", interlines);
      sscanf (argV[4], "%d", func);
      sscanf (argV[5], "%d", termination);
      if (*termination == 1)
	{
	  sscanf (argV[6], "%lf", term_precision);
	  *term_iteration = MAX_ITERATION;
	}
      else
	{
	  sscanf (argV[6], "%d", term_iteration);
	  *term_precision = 0;
	}
    }
}
