#define _DEFAULT_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

int main(void) {

    struct timeval tv;
    time_t time;
    int micro_sec;
    char time_string[30];
    char output[80];
    char hostname[30];

    gettimeofday(&tv, NULL);
    gethostname(hostname, 30);

    time = tv.tv_sec;
    micro_sec = tv.tv_usec;

    strftime(time_string, 30, "%Y-%m-%d %T", localtime(&time));
    snprintf(output, 80, "%s : %s.%d", hostname, time_string, (int)micro_sec);

    printf("%s\n", output);
    printf("%d\n", (int)micro_sec);

    return 0;
}