#define _DEFAULT_SOURCE

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

void get_output(char *output, int p_rank) {

    struct timeval tv;
    time_t time;
    int micro_sec;
    char time_string[30];
    char hostname[30];

    gettimeofday(&tv, NULL);
    gethostname(hostname, 30);

    time = tv.tv_sec;
    micro_sec = tv.tv_usec;

    strftime(time_string, 30, "%Y-%m-%d %T", localtime(&time));
    snprintf(output, 80, "[%d] %s : %s.%d", p_rank, hostname, time_string,
             (int)micro_sec);
}

int main(int argc, char **argv) {

    int p_rank;
    int p_num;
    char output[80];

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &p_num);

    if (p_rank == p_num - 1) {
        if (p_num == 1) {
            // Print own output if the last process is the only existing process
            get_output(output, p_rank);
            printf("%s\n", output);
        }

        MPI_Status status;
        for (int p = 0; p < p_num - 1; p++) {
            // Recieve and print the output of the processes with ranks 0 to n-2
            MPI_Recv(output, 80, MPI_CHAR, p, 0, MPI_COMM_WORLD, &status);
            // We use fprintf instead of printf because of following issue:
            // https://stackoverflow.com/questions/42743100/mpi-barrier-doesnt-seem-to-work-reordering-printf-stdout-messages
            fprintf(stdout, "%s\n", output);
        }
    } else {
        // Generate the output for the processes with ranks 0 to n-2
        get_output(output, p_rank);

        // Send the output to the process with the last rank (n-1)
        MPI_Send(output, 80, MPI_CHAR, p_num - 1, 0, MPI_COMM_WORLD);
    }

    // Wait until all output is printed
    MPI_Barrier(MPI_COMM_WORLD);
    // Without sleep(0) the order of the outputs would be non-deterministic.
    // There is a problem with MPI_Barrier and the buffering in libc/glibc printing.
    // The best case would be to let the last process print all "[p] beendet jetzt!"
    // but the exercise requires that each process should print the text.
    sleep(0);

    printf("[%d] beendet jetzt!\n", p_rank);

    // Finalize the MPI environment
    MPI_Finalize();

    return 0;
}
