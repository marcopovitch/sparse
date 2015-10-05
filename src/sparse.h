#include <stdio.h>
#include <stdlib.h>
/* #include <malloc.h> */
#include <assert.h>
#include <math.h>

#include "matrice.h"

#ifndef __SPARSE_H__
#define __SPARSE_H__

#define EPS_SPARSE 1.0e-6
/* #define NOVALUE -999 */
#define NOVALUE 0

enum {
    SPARSE_COL_LINK = 1
};

struct sparse_item_t {
    long int col_index;
    long int line_index;
    double val;
    struct sparse_item_t *next_in_line;
    struct sparse_item_t *next_in_col;
};

/*
 * col_link_status=SPARSE_COL_LINK, speeds up importation of sparse matrix,
 * ONLY IF the data are ordered  in the file such as :   for (l=0;
 * l<line<l++) { for (c=0; c<col; c++) read_from_file (element[l][c]); }
 */
struct sparse_matrix_t {
    long int nb_line;
    long int nb_col;
    long int nb_item;
    struct sparse_item_t **line;
    struct sparse_item_t **col;
    int col_link_status;
    struct sparse_item_t **last_col;
};

char *libsparseversion();
struct sparse_matrix_t *new_sparse_matrix(long int nb_line,
                                          long int nb_col,
                                          int col_link_status);
void free_sparse_matrix(struct sparse_matrix_t *m);

double sparse_get_value(struct sparse_matrix_t *m, long int i, long int j);
struct sparse_item_t *sparse_set_value(struct sparse_matrix_t *m,
                                       long int i, long int j, double val,
                                       struct sparse_item_t *last_item);
void sparse_update_col_link(struct sparse_matrix_t *m,
                            struct sparse_item_t *item);

struct sparse_matrix_t *read_sparse_matrix(char *filename,
                                           int col_link_status);
struct sparse_matrix_t *read_ijk_sparse_matrix(char *filename,
                                               int col_link_status);

struct sparse_matrix_t *import_sparse_matrix(struct sparse_matrix_t *a,
                                             char *filename);
void write_sparse_matrix(struct sparse_matrix_t *A, char *filename);
void write_sparse_matrix_with_line_offset(struct sparse_matrix_t *A,
                                          long int offset, char *filename);

void dump_sparse_matrix(struct sparse_matrix_t *m);
void dump_sparse_matrix_to_scilab(struct sparse_matrix_t *m);

struct vector_t *sparse_extract_col(struct sparse_matrix_t *A, long int c);
struct vector_t *sparse_extract_line(struct sparse_matrix_t *A,
                                     long int l);

struct sparse_matrix_t *sparsify(struct matrix_t *M, int col_link_status);

struct sparse_matrix_t *sparse_matrix_resize(struct sparse_matrix_t *m,
                                             long int nbline,
                                             long int nbcol);

int check_sparse_matrix(struct sparse_matrix_t *m);
void sparse_compute_length(struct sparse_matrix_t *m, char *filename);
struct sparse_matrix_t *AtransA(struct sparse_matrix_t *A);
double mean_diag_AtA(struct sparse_matrix_t *A);
void show_sparse_stats(struct sparse_matrix_t *A);
#endif
