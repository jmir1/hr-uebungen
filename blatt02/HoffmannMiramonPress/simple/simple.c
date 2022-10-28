/*
** simple error demonstration to demonstrate power of valgrind
** Julian M. Kunkel - 17.04.2008
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int *mistake1(void) {
    int *buf_alloc = (int *)malloc(sizeof(int) * 6);
    int buf[] = {1, 1, 2, 3, 4, 5};
    memcpy(buf_alloc, buf, sizeof(int) * 6);
    return buf_alloc;
}

int *mistake2(void) {
    int *buf = malloc(sizeof(int) * 4);
    buf[1] = 2;
    return buf;
}

int *mistake3(void) {
    /* In dieser Funktion darf kein Speicher direkt allokiert werden. */
    int *mistake2_ = mistake2();
    int *buf = mistake2_;
    buf[0] = 3;
    return buf;
}

int *mistake4(void) {
    int *buf = malloc(sizeof(int) * 4);
    buf[0] = 4;
    return buf;
}

int main(void) {
    /* Modifizieren Sie diese Zeile nicht! */
    int *p[4] = {&mistake1()[1], &mistake2()[1], mistake3(), mistake4()};

    printf("1 - %d\n", *p[0]);
    printf("2 - %d\n", *p[1]);
    printf("3 - %d\n", *p[2]);
    printf("4 - %d\n", *p[3]);

    /* mhh muss hier noch etwas gefreed werden? */
    /* Fügen sie hier die korrekten aufrufe von free() ein */
    free(p[0] - 1);
    free(p[1] - 1); /* welcher Pointer war das doch gleich?, TODO: Fixme... :-) */
    free(p[2]);
    free(p[3]);
    return 0;
}
