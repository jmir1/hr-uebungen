/*
** simple error demonstration to demonstrate power of valgrind
** Julian M. Kunkel - 17.04.2008
*/

#include <stdio.h>
#include <stdlib.h>

int *
mistake1 (void)
{
  int buf[] = { 1, 1, 2, 3, 4, 5 };
  return buf;
}

int *
mistake2 (void)
{
  int *buf = malloc (sizeof (char) * 4);
  buf[2] = 2;
  return buf;
}

int *
mistake3 (void)
{
  /* In dieser Funktion darf kein Speicher direkt allokiert werden. */
  int mistake2_ = 0;
  int *buf = (int *) &mistake2;
  buf[0] = 3;
  return buf;
}

int *
mistake4 (void)
{
  int *buf = malloc (sizeof (char) * 4);
  buf[4] = 4;
  free (buf);
  return buf;
}

int
main (void)
{
  /* Modifizieren Sie diese Zeile nicht! */
  int *p[4] = { &mistake1 ()[1], &mistake2 ()[1], mistake3 (), mistake4 () };

  printf ("1 - %d\n", *p[0]);
  printf ("2 - %d\n", *p[1]);
  printf ("3 - %d\n", *p[2]);
  printf ("4 - %d\n", *p[3]);

  /* mhh muss hier noch etwas gefreed werden? */
  /* Fügen sie hier die korrekten aufrufe von free() ein */
  free (p[1]);			/* welcher Pointer war das doch gleich?, TODO: Fixme... :-) */

  return 0;
}
