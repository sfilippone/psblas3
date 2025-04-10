#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void compare_files(const char *file1, const char *file2) {
    FILE *file_s = fopen(file1, "r");
    FILE *file_p = fopen(file2, "r");

    if (!file_s || !file_p) {
        perror("Error opening files");
        exit(EXIT_FAILURE);
    }

    int n_s, n_p;
    fscanf(file_s, "%d", &n_s);
    fscanf(file_p, "%d", &n_p); // Assuming both files have the same number of lines

    if (n_s != n_p) {
        fprintf(stderr, "Error, differnet file sizes $d, $d", n_s, n_p);
        exit(EXIT_FAILURE);
    }

    double value_s, value_p;
    

    for (int i = 0; i < n_s; i++) {
        fscanf(file_s, "%lf", &value_s);
        fscanf(file_p, "%lf", &value_p);
        
        if (value_s != value_p) {
            printf("Index %d: %.2lf != %.2lf\n", i, value_s, value_p);
        }
    }

    fclose(file_s);
    fclose(file_p);
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <partial_filename>\n", argv[0]);
        return EXIT_FAILURE;
    }

    char file_s[256], file_p[256];
    snprintf(file_s, sizeof(file_s), "%s_serial.txt", argv[1]);
    snprintf(file_p, sizeof(file_p), "%s_parallel.txt", argv[1]);

    compare_files(file_s, file_p);

    return EXIT_SUCCESS;
}
