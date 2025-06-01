#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>

#define MAX_COLUMNS 100
#define MAX_LINE_LENGTH 1024

// Structure to hold a column for sorting
typedef struct {
    char content[256]; // Column content
} Column;

// Natural comparison function for qsort
int natural_compare(const void *a, const void *b) {
    const char *str_a = *(const char **)a;
    const char *str_b = *(const char **)b;

    while (*str_a && *str_b) {
        if (isdigit(*str_a) && isdigit(*str_b)) {
            // Compare numeric parts
            char *end_a, *end_b;
            long num_a = strtol(str_a, &end_a, 10);
            long num_b = strtol(str_b, &end_b, 10);

            if (num_a != num_b) return num_a - num_b;

            // Move pointers past the numeric part
            str_a = end_a;
            str_b = end_b;
        } else {
            // Lexicographical comparison for non-numeric parts
            if (*str_a != *str_b) return *str_a - *str_b;
            str_a++;
            str_b++;
        }
    }
    return *str_a - *str_b;
}

// Function to process the file
void process_file(const char *input_file, const char *output_file, const char *cols, const char *sort_cols) {
    FILE *infile = fopen(input_file, "r");
    FILE *outfile = fopen(output_file, "w");
    if (!infile || !outfile) {
        perror("Error opening files");
        exit(EXIT_FAILURE);
    }

    int selected_columns[MAX_COLUMNS], sort_columns[MAX_COLUMNS];
    int num_selected = 0, num_sort = 0;

    // Parse the `cols` parameter
    char cols_copy[256];
    strncpy(cols_copy, cols, sizeof(cols_copy));
    char *token = strtok(cols_copy, ",");
    while (token) {
        selected_columns[num_selected++] = atoi(token) - 1; // Convert to 0-based index
        token = strtok(NULL, ",");
    }

    // Parse the `sort_cols` parameter, if it is provided
    if (strlen(sort_cols) > 0) {
        char sort_cols_copy[256];
        strncpy(sort_cols_copy, sort_cols, sizeof(sort_cols_copy));
        token = strtok(sort_cols_copy, ",");
        while (token) {
            sort_columns[num_sort++] = atoi(token) - 1; // Convert to 0-based index
            token = strtok(NULL, ",");
        }
    }

    char line[MAX_LINE_LENGTH];
    while (fgets(line, sizeof(line), infile)) {
        char *columns[MAX_COLUMNS];
        int col_count = 0;

        // Split the line into columns (handling both space and tab delimiters)
        token = strtok(line, "\t ");
        while (token) {
            columns[col_count++] = token;
            token = strtok(NULL, "\t ");
        }

        // Prepare output fields
        char *output_fields[MAX_COLUMNS];
        int output_count = 0;

        char *to_sort[MAX_COLUMNS];
        int sort_count = 0;

        // Collect selected columns and separate sortable ones
        for (int i = 0; i < num_selected; i++) {
            int col_idx = selected_columns[i];
            if (col_idx < col_count) {
                int is_sort_col = 0;

                for (int j = 0; j < num_sort; j++) {
                    if (sort_columns[j] == col_idx) {
                        to_sort[sort_count++] = columns[col_idx];
                        is_sort_col = 1;
                        break;
                    }
                }

                if (!is_sort_col) {
                    output_fields[output_count++] = columns[col_idx];
                }
            }
        }

        // If sort_columns is provided, sort the sortable columns
        if (num_sort > 0) {
            qsort(to_sort, sort_count, sizeof(char *), natural_compare);
        }

        // Write the selected and sorted fields to output
        for (int i = 0; i < output_count; i++) {
            fprintf(outfile, "%s", output_fields[i]);
            if (i < output_count - 1 || sort_count > 0) {
                fprintf(outfile, ":");
            }
        }

        // Write sorted fields to output
        for (int i = 0; i < sort_count; i++) {
            fprintf(outfile, "%s", to_sort[i]);
            if (i < sort_count - 1) {
                fprintf(outfile, ":");
            }
        }

        // End the line
        fprintf(outfile, "\n");
    }

    fclose(infile);
    fclose(outfile);
}

// Function to check if a file is a gzip file
int is_gzip_file(const char *filename) {
    FILE *file = fopen(filename, "rb");
    if (!file) return 0;

    unsigned char buffer[2];
    fread(buffer, 1, 2, file);
    fclose(file);

    return buffer[0] == 0x1f && buffer[1] == 0x8b;
}

// Function to decompress a gzip file
void decompress_gzip(const char *input_file, const char *temp_file) {
    gzFile gz = gzopen(input_file, "rb");
    FILE *out = fopen(temp_file, "wb");
    if (!gz || !out) {
        perror("Failed to open files for decompression");
        exit(EXIT_FAILURE);
    }

    char buffer[1024];
    int bytes_read;
    while ((bytes_read = gzread(gz, buffer, sizeof(buffer))) > 0) {
        fwrite(buffer, 1, bytes_read, out);
    }

    gzclose(gz);
    fclose(out);
}

// Main function to handle arguments and process the files
int main(int argc, char *argv[]) {
    if (argc < 4 || argc > 5) {
        fprintf(stderr, "Usage: %s <input_file> <suffix> <cols> [sort_cols]\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char *input_file = argv[1];
    const char *suffix = argv[2];
    const char *cols = argv[3];
    const char *sort_cols = argc == 5 ? argv[4] : "";  // Use empty string if sort_cols is not provided

    // Extract base name and extension
    char base_name[MAX_LINE_LENGTH], extension[32], output_file[MAX_LINE_LENGTH], temp_file[MAX_LINE_LENGTH];
    strncpy(base_name, input_file, sizeof(base_name));
    char *dot = strrchr(base_name, '.');
    if (dot) {
        strncpy(extension, dot + 1, sizeof(extension));
        *dot = '\0';
    } else {
        strncpy(extension, "txt", sizeof(extension));
    }

    // Construct output file name by adding the suffix
    if (snprintf(output_file, sizeof(output_file), "%s_%s.%s", base_name, suffix, extension) >= sizeof(output_file)) {
        fprintf(stderr, "Error: Output file name truncated.\n");
        return EXIT_FAILURE;
    }

    // Construct temporary file name
    if (snprintf(temp_file, sizeof(temp_file), "%s_temp.txt", base_name) >= sizeof(temp_file)) {
        fprintf(stderr, "Error: Temporary file name truncated.\n");
        return EXIT_FAILURE;
    }

    // Check if the file is gzipped, if so, decompress
    if (is_gzip_file(input_file)) {
        decompress_gzip(input_file, temp_file);
        process_file(temp_file, output_file, cols, sort_cols);
        remove(temp_file);
    } else {
        process_file(input_file, output_file, cols, sort_cols);
    }

    return EXIT_SUCCESS;
}
