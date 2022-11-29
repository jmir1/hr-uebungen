#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int *init(int N) {
    // TODO
    int *buf = (int *)malloc(sizeof(int) * N);

    srand(time(NULL));

    for (int i = 0; i < N; i++) {
        // Do not modify "% 25"
        buf[i] = rand() % 25;
    }

    return buf;
}

int circle(int *buf) {
    // TODO
    return 0;
}

int main(int argc, char **argv) {
    int N;
    int rank;
    int *buf;
    int ret;

    if (argc < 2) {
        printf("Arguments error!\nPlease specify a buffer size.\n");
        return EXIT_FAILURE;
    }

    // Array length
    N = atoi(argv[1]);
    buf = init(N);

    // TODO
    rank = 0;

    printf("\nBEFORE\n");

    for (int i = 0; i < N; i++) {
        printf("rank %d: %d\n", rank, buf[i]);
    }

    ret = circle(buf);

    printf("\nAFTER\n");

    for (int j = 0; j < N; j++) {
        printf("rank %d: %d\n", rank, buf[j]);
    }

    return EXIT_SUCCESS;
}
