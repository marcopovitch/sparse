#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "sparse.h"

/** \brief Library information **/
char *libsparseversion()
{
    char *s;

    s = malloc(strlen(PACKAGE) + strlen(VERSION) + 3);
    sprintf(s, "%s-%s", PACKAGE, VERSION);
    return (s);
}

/** \brief Create a sparse matrix **/
struct sparse_matrix_t *new_sparse_matrix(long int nb_line,
                                          long int nb_col,
                                          int col_link_status)
{
    struct sparse_matrix_t *matrix;

    matrix = (struct sparse_matrix_t *)
        malloc(sizeof(struct sparse_matrix_t));
    assert(matrix);

    matrix->nb_line = nb_line;
    matrix->nb_col = nb_col;
    matrix->line = (struct sparse_item_t **)
        calloc(nb_line, sizeof(struct sparse_item_t *));
    assert(matrix->line);
    matrix->col = (struct sparse_item_t **)
        calloc(nb_col, sizeof(struct sparse_item_t *));
    assert(matrix->col);

    matrix->col_link_status = col_link_status;
    if (col_link_status == SPARSE_COL_LINK) {
        matrix->last_col = (struct sparse_item_t **)
            calloc(nb_col, sizeof(struct sparse_item_t *));
        assert(matrix->last_col);
    } else {
        matrix->last_col = NULL;
    }

    matrix->nb_item = 0;

    return (matrix);
}

/** \brief return value at (i,j) position in the sparse matrix **/
double sparse_get_value(struct sparse_matrix_t *m, long int i, long int j)
{
    struct sparse_item_t *cur_item;

    assert(!(i > m->nb_line));
    assert(!(j > m->nb_col));

    cur_item = m->line[i];

    if (!cur_item) {
        return (0);
    }
    while (cur_item && cur_item->col_index < j) {
        cur_item = cur_item->next_in_line;
    }

    if (cur_item && cur_item->col_index == j) {
        return (cur_item->val);
    }
    return (0);
}

void
sparse_update_col_link(struct sparse_matrix_t *m,
                       struct sparse_item_t *item)
{
    struct sparse_item_t *cur_item;

    struct sparse_item_t *last_item = NULL;

    long int i, j;

    i = item->line_index;
    j = item->col_index;

    if (m->col_link_status == SPARSE_COL_LINK) {
        cur_item = m->last_col[j];
    } else {
        cur_item = m->col[j];
    }
    assert(!(cur_item && cur_item->line_index >= item->line_index));

    if (cur_item && cur_item->line_index >= item->line_index) {
        fprintf(stderr, "?");
    }
    /* first item */
    if (!cur_item) {
        m->col[j] = item;
        if (m->col_link_status == SPARSE_COL_LINK) {
            m->last_col[j] = item;
        }
        return;
    }
    if (m->col_link_status == SPARSE_COL_LINK) {
        assert(j == cur_item->col_index);
    }
    /* find item place */
    while (cur_item && cur_item->line_index < i) {
        last_item = cur_item;
        cur_item = cur_item->next_in_col;
    }

    /* the last one */
    if (!cur_item) {
        last_item->next_in_col = item;
        if (m->col_link_status == SPARSE_COL_LINK) {
            m->last_col[j] = item;
        }
        return;
    }
    if (m->col_link_status == SPARSE_COL_LINK) {
        /*
         * we should not be there, the item to insert is always at
         * the end of the list
         */
        assert(0);
    }
    /* item found */
    if (cur_item->line_index == i) {
        /* already done */
        return;
    }
    /* new item to be inserted here */
    if (!last_item) {
        /* become the first item */

        item->next_in_col = cur_item;
        m->col[j] = item;
    } else {
        last_item->next_in_col = item;
        item->next_in_col = cur_item;
    }
}

/** \brief set value at (i,j) position in the sparse matrix m **/
struct sparse_item_t *sparse_set_value(struct sparse_matrix_t *m,
                                       long int i, long int j, double val,
                                       struct sparse_item_t *previous)
{
    struct sparse_item_t *cur_item;

    struct sparse_item_t *new_item, *last_item = NULL;

    assert(!(i >= m->nb_line));
    assert(!(j >= m->nb_col));

    if (previous) {
        assert(!(previous->line_index != i && previous->col_index > j));
        cur_item = previous;
    } else {
        cur_item = m->line[i];

        /* first item */
        if (!cur_item) {
            /*
             * fprintf(stderr, "set: first item [%ld,%ld]\n",
             * i,j);
             */
            new_item = (struct sparse_item_t *)
                calloc(1, sizeof(struct sparse_item_t));
            assert(new_item);

            new_item->line_index = i;
            new_item->col_index = j;
            new_item->val = val;

            m->line[i] = new_item;
            m->nb_item++;

            if (m->col_link_status == SPARSE_COL_LINK) {
                sparse_update_col_link(m, new_item);
            }
            return (new_item);
        }
    }

    /* find item place */
    while (cur_item && cur_item->col_index < j) {
        last_item = cur_item;
        cur_item = cur_item->next_in_line;
    }

    /* the last one */
    if (!cur_item) {
        /* fprintf(stderr, "set: last item [%ld,%ld]\n", i,j); */
        new_item = (struct sparse_item_t *)
            calloc(1, sizeof(struct sparse_item_t));
        assert(new_item);

        last_item->next_in_line = new_item;

        new_item->line_index = i;
        new_item->col_index = j;
        new_item->val = val;
        new_item->next_in_line = NULL;

        m->nb_item++;
        if (m->col_link_status == SPARSE_COL_LINK) {
            sparse_update_col_link(m, new_item);
        }
        return (new_item);
    }
    /* item found */
    if (cur_item->col_index == j) {
        /* replace the current item */
        /* fprintf(stderr, "set: found item [%ld,%ld]\n", i,j); */

        fprintf(stderr,
                "sparse_set_value: duplicate (%ld,%ld) old=%f/new=%f\n", i,
                j, cur_item->val, val);

        cur_item->val += val;
        /* m->nb_item ++; */
        /* sparse_update_col_link (m, cur_item); */
        return (cur_item);
    }
    /* new item to be inserted here */
    /* fprintf(stderr, "set: inserted item [%ld,%ld]\n", i,j); */
    new_item = (struct sparse_item_t *)
        calloc(1, sizeof(struct sparse_item_t));
    assert(new_item);

    if (!last_item) {
        /* become the first item */
        m->line[i] = new_item;
    } else {
        last_item->next_in_line = new_item;
    }

    new_item->next_in_line = cur_item;

    new_item->line_index = i;
    new_item->col_index = j;
    new_item->val = val;

    m->nb_item++;
    if (m->col_link_status == SPARSE_COL_LINK) {
        sparse_update_col_link(m, new_item);
    }
    return (new_item);
}

/** \brief Read sparse matrix file formated such as :
	ni nj
	i j value
	...

   col_link_status=SPARSE_COL_LINK, speeds up importation of sparse matrix,
   ONLY IF the data are ordered  in the file such as :
	for (l=0; l<line<l++) {
 		for (c=0; c<col; c++)
 			get element[l][c];
	
**/
struct sparse_matrix_t *read_ijk_sparse_matrix(char *filename,
                                               int col_link_status)
{
    struct sparse_matrix_t *a;

    long int m, n, j, i;

    double val;

    FILE *fd;

    int nb_read;

    long int nb_item;

    long int cpt = 0;

    fprintf(stdout, "reading sparse matrix from '%s' ... ", filename);
    fflush(stdout);

    if (!(fd = fopen(filename, "r"))) {
        perror(filename);
        exit(1);
    }
    nb_read = fscanf(fd, "%ld %ld", &m, &n);
    if (nb_read != 2) {
        fprintf(stdout, "\n");
        fprintf(stderr,
                "read_ijk_sparse_matrix: error reading (m,n) in '%s'\n",
                filename);
        exit(1);
    }
    fprintf(stdout, "(%ldx%ld) ", m, n);
    a = new_sparse_matrix(m, n, col_link_status);

    while (1) {

        nb_read = fscanf(fd, "%ld %ld %lf", &i, &j, &val);

        if (feof(fd)) {
            break;
        }
        if (nb_read != 3) {
            fprintf(stdout, "\n");
            fprintf(stderr,
                    "read_ijk_sparse_matrix: file '%s' corrupted nread=%d\n",
                    filename, nb_read);
            exit(1);
        }
        //fprintf(stderr, "%ld %ld %lf\n", i, j, val);
        sparse_set_value(a, i, j, val, NULL);
        cpt++;
    }

    fprintf(stdout, "%ld lines\n", cpt);
    fflush(stdout);

    fclose(fd);
    return (a);
}

/** \brief Readme sparse matrix from file

 col_link_status=SPARSE_COL_LINK, speeds up importation of sparse matrix,
 ONLY IF the data are ordered  in the file such as :
 	for (l=0; l<line<l++) {
  		for (c=0; c<col; c++)
 			get element[l][c];
**/
struct sparse_matrix_t *read_sparse_matrix(char *filename,
                                           int col_link_status)
{
    struct sparse_matrix_t *a;

    long int m, n, j, index;

    double val;

    FILE *fd;

    int nb_read;

    long int nb_item;

    long int rayid;

    long int cpt = 0;

    struct sparse_item_t *last_item = NULL;     /* to speed up things */

    fprintf(stdout, "reading sparse matrix from '%s' ... ", filename);
    fflush(stdout);

    if (!(fd = fopen(filename, "r"))) {
        perror(filename);
        exit(1);
    }
    nb_read = fscanf(fd, "%ld %ld", &m, &n);
    if (nb_read != 2) {
        fprintf(stdout, "\n");
        fprintf(stderr,
                "read_sparse_matrix: error reading (m,n) in '%s'\n",
                filename);
        exit(1);
    }
    fprintf(stdout, "(%ldx%ld) ", m, n);
    a = new_sparse_matrix(m, n, col_link_status);

    while (1) {

        nb_read = fscanf(fd, "%ld %ld", &rayid, &nb_item);

        if (feof(fd)) {
            break;
        }
        if (nb_read != 2) {
            fprintf(stdout, "\n");
            fprintf(stderr,
                    "read_sparse_matrix: file '%s' corrupted nread=%d\n",
                    filename, nb_read);
            exit(1);
        }
        if (a->col_link_status == SPARSE_COL_LINK) {
            last_item = NULL;
        }
        for (j = 0; j < nb_item; j++) {

            nb_read = fscanf(fd, "%ld %lf", &index, &val);

            if (nb_read != 2 && !feof(fd)) {
                fprintf(stdout, "\n");
                fprintf(stderr,
                        "read_sparse_matrix: error reading item (%ld,%ld) in '%s' nread=%d\n",
                        rayid, j, filename, nb_read);
                exit(1);
            }
            if (a->col_link_status == SPARSE_COL_LINK) {
                last_item =
                    sparse_set_value(a, rayid, index, val, last_item);
            } else {
                sparse_set_value(a, rayid, index, val, NULL);
            }
        }
        cpt++;
    }

    fprintf(stdout, "%ld lines\n", cpt);
    fflush(stdout);

    fclose(fd);
    return (a);
}

struct sparse_matrix_t *import_sparse_matrix(struct sparse_matrix_t *a,
                                             char *filename)
{
    struct sparse_matrix_t *sparse;

    struct sparse_item_t *last_item;

    long int m, n, j, index;

    double val;

    FILE *fd;

    int nb_read;

    long int nb_item;

    long int rayid;

    long int cpt = 0;

    if (!a) {
        sparse = read_sparse_matrix(filename, SPARSE_COL_LINK);
        return (sparse);
    }
    fprintf(stdout, "importing sparse matrix from '%s' into (%p) ... ",
            filename, a);
    fflush(stdout);
    if (!(fd = fopen(filename, "r"))) {
        perror(filename);
        exit(1);
    }
    nb_read = fscanf(fd, "%ld %ld\n", &m, &n);
    if (nb_read != 2) {
        fprintf(stdout, "\n");
        fprintf(stderr, "Error reading (m,n) in '%s'\n", filename);
        exit(1);
    }
    if (a->nb_line != m && a->nb_col != n) {
        fprintf(stdout, "\n");
        fprintf(stderr,
                "Error importing '%s' into (%p), nb_line=%ld/%ld nb_col=%ld/%ld\n",
                filename, a, a->nb_line, m, a->nb_col, n);
        exit(1);
    }
    while (1) {
        nb_read = fscanf(fd, "%ld %ld\n", &rayid, &nb_item);

        if (feof(fd)) {
            break;
        }
        if (nb_read != 2) {
            fprintf(stdout, "\n");
            fprintf(stderr, "import_sparse_matrix: file '%s' corrupted\n",
                    filename);
            exit(1);
        }
        last_item = NULL;
        for (j = 0; j < nb_item; j++) {
            nb_read = fscanf(fd, "%ld %lf", &index, &val);
            if (nb_read != 2 && !feof(fd)) {
                fprintf(stdout, "\n");
                fprintf(stderr,
                        "import_sparse_matrix: error reading item (%ld,%ld) in '%s'\n",
                        rayid, j, filename);
                exit(1);
            }
            if (a->col_link_status == SPARSE_COL_LINK) {
                last_item =
                    sparse_set_value(a, rayid, index, val, last_item);
            } else {
                sparse_set_value(a, rayid, index, val, NULL);
            }
        }
        cpt++;
    }
    fclose(fd);
    fprintf(stdout, "%ld lines\n", cpt);
    fflush(stdout);
    return (a);
}

/** \brief Write sparse matrix A to file **/
void write_sparse_matrix(struct sparse_matrix_t *A, char *filename)
{
    long int i, nb_item;

    FILE *fd;

    struct sparse_item_t *cur_item;

    if (!(fd = fopen(filename, "w"))) {
        perror(filename);
        exit(1);
    }
    fprintf(stdout, "writing sparse matrix (%p) to '%s' %ld items\n",
            A, filename, A->nb_item);

    fprintf(fd, "%ld %ld\n", A->nb_line, A->nb_col);

    for (i = 0; i < A->nb_line; i++) {

        /* how many item in the current line */
        cur_item = A->line[i];
        nb_item = 0;
        while (cur_item) {
            nb_item++;
            cur_item = cur_item->next_in_line;
        }

        if (!nb_item) {
            continue;
        }
        /* write items value */
        cur_item = A->line[i];
        fprintf(fd, "%ld %ld\n", cur_item->line_index, nb_item);
        while (cur_item) {
            fprintf(fd, "%ld %lf ", cur_item->col_index, cur_item->val);
            cur_item = cur_item->next_in_line;
        }
        fprintf(fd, "\n");
    }
    fclose(fd);
}

void
write_sparse_matrix_with_line_offset(struct sparse_matrix_t *A,
                                     long int offset, char *filename)
{
    long int i, nb_item;

    FILE *fd;

    struct sparse_item_t *cur_item;

    if (!(fd = fopen(filename, "w"))) {
        perror(filename);
        exit(1);
    }
    fprintf(stdout,
            "writing sparse matrix (%p) offset=%ld to '%s' %ld items\n", A,
            offset, filename, A->nb_item);

    fprintf(fd, "%ld %ld\n", A->nb_line + offset, A->nb_col);

    for (i = 0; i < A->nb_line; i++) {

        /* how many item in the current line */
        cur_item = A->line[i];
        nb_item = 0;
        while (cur_item) {
            nb_item++;
            cur_item = cur_item->next_in_line;
        }

        if (!nb_item) {
            continue;
        }
        /* write items value */
        cur_item = A->line[i];
        fprintf(fd, "%ld %ld\n", cur_item->line_index + offset, nb_item);
        while (cur_item) {
            fprintf(fd, "%ld %lf ", cur_item->col_index, cur_item->val);
            cur_item = cur_item->next_in_line;
        }
        fprintf(fd, "\n");
    }
    fclose(fd);
}

struct vector_t *sparse_extract_line(struct sparse_matrix_t *A, long int l)
{
    struct vector_t *tmp;

    struct sparse_item_t *cur_item;

    tmp = new_vector(A->nb_col);
    cur_item = A->line[l];

    while (cur_item) {
        tmp->mat[cur_item->col_index] = cur_item->val;
        cur_item = cur_item->next_in_line;
    }

    return (tmp);
}

struct vector_t *sparse_extract_col(struct sparse_matrix_t *A, long int c)
{
    struct vector_t *tmp;

    struct sparse_item_t *cur_item;

    tmp = new_vector(A->nb_line);
    cur_item = A->col[c];

    while (cur_item) {
        tmp->mat[cur_item->line_index] = cur_item->val;
        cur_item = cur_item->next_in_col;
    }

    return (tmp);
}

void dump_sparse_matrix(struct sparse_matrix_t *m)
{
    long int i, j;

    double val;

    struct vector_t *quick_access;

    for (i = 0; i < m->nb_line; i++) {
        quick_access = sparse_extract_line(m, i);
        for (j = 0; j < m->nb_col; j++) {
            val = quick_access->mat[j];
            if (fabs(val) > EPS_SPARSE)
                fprintf(stderr, "%f\n", val);
        }
        free_vector(quick_access);
        /* fprintf(stderr,"\n"); */
    }

}

void dump_sparse_matrix_to_scilab(struct sparse_matrix_t *m)
{
    long int i, j;

    double val;

    struct vector_t *quick_access;

    for (i = 0; i < m->nb_line; i++) {
        if (i == 0) {
            fprintf(stderr, "[");
        }
        quick_access = sparse_extract_line(m, i);
        for (j = 0; j < m->nb_col; j++) {
            val = quick_access->mat[j];
            fprintf(stderr, "%f", val);
            if (j != m->nb_col - 1) {
                fprintf(stderr, ",");
            }
        }
        free_vector(quick_access);
        if (i != m->nb_line - 1) {
            fprintf(stderr, ";");
        }
    }
    fprintf(stderr, "]\n");

}

void free_sparse_matrix(struct sparse_matrix_t *m)
{
    long int i;

    struct sparse_item_t *cur_item, *last_item;

    long int n = 0;

    for (i = 0; i < m->nb_line; i++) {
        cur_item = m->line[i];
        while (cur_item) {
            last_item = cur_item;
            cur_item = cur_item->next_in_line;
            free(last_item);
            n++;
        }
    }

    if (m->col_link_status == SPARSE_COL_LINK) {
        free(m->last_col);
    }
    fprintf(stdout, "free sparse matrix (%p): %ld items\n", m, n);
    fflush(stdout);

    free(m->line);
    free(m->col);
    free(m);
}

struct sparse_matrix_t *sparsify(struct matrix_t *M, int col_link_status)
{
    struct sparse_matrix_t *spm;

    long int i, j;

    double val;

    spm = new_sparse_matrix(M->nb_line, M->nb_col, col_link_status);

    for (i = 0; i < M->nb_line; i++) {
        for (j = 0; j < M->nb_col; j++) {
            val = M->mat[i][j];
            if (fabs(val) > EPS_SPARSE) {
                sparse_set_value(spm, i, j, val, NULL);
            }
        }
    }

    return (spm);
}

struct sparse_matrix_t *sparse_matrix_resize(struct sparse_matrix_t *m,
                                             long int nbline,
                                             long int nbcol)
{
    long int i;

    assert(m);

    fprintf(stdout,
            "sparse_matrix_resize (%p) from (%ldx%ld) to (%ldx%ld)\n", m,
            m->nb_line, m->nb_col, nbline, nbcol);

    if (nbline != m->nb_line) {
        if (nbline > m->nb_line) {
            m->line = (struct sparse_item_t **)
                realloc(m->line, nbline * sizeof(struct sparse_item_t *));
            assert(m->line);
            for (i = m->nb_line; i < nbline; i++) {
                m->line[i] = NULL;
            }
            m->nb_line = nbline;
        } else {
            for (i = nbline; i < m->nb_line; i++) {
                if (m->line[i] != NULL) {
                    fprintf(stderr,
                            "sparse_matrix_resize: [line] some items will be lost !\n");
                    /* fixme free those items */
                    break;
                }
                m->line =
                    (struct sparse_item_t **) realloc(m->line,
                                                      nbline *
                                                      sizeof(struct
                                                             sparse_item_t
                                                             *));
                assert(m->line);
                m->nb_line = nbline;
            }
        }
    }
    if (nbcol != m->nb_col) {
        if (nbcol > m->nb_col) {
            m->col = (struct sparse_item_t **)
                realloc(m->col, nbcol * sizeof(struct sparse_item_t *));
            assert(m->col);
            if (m->col_link_status == SPARSE_COL_LINK) {
                m->last_col = (struct sparse_item_t **)
                    realloc(m->col,
                            nbcol * sizeof(struct sparse_item_t *));
                assert(m->last_col);
            }
            for (i = m->nb_col; i < nbcol; i++) {
                m->col[i] = NULL;
                if (m->col_link_status == SPARSE_COL_LINK) {
                    m->last_col[i] = NULL;
                }
            }
            m->nb_col = nbcol;
        } else {
            for (i = nbcol; i < m->nb_col; i++) {
                if (m->col[i] != NULL) {
                    fprintf(stderr,
                            "sparse_matrix_resize: [col] some items will be lost !\n");
                    /* fixme free those items */
                    break;
                }
                m->col =
                    (struct sparse_item_t **) realloc(m->col,
                                                      nbcol *
                                                      sizeof(struct
                                                             sparse_item_t
                                                             *));
                assert(m->col);
                m->nb_col = nbcol;

                if (m->col_link_status == SPARSE_COL_LINK) {
                    m->last_col = (struct sparse_item_t **)
                        realloc(m->col,
                                nbcol * sizeof(struct sparse_item_t *));
                    assert(m->last_col);
                }
            }
        }
    }
    return (m);
}

/************************* debug stuff ***********************************/
int check_sparse_matrix(struct sparse_matrix_t *m)
{
    long int i, j, cpt_col, cpt_line;

    struct sparse_item_t *cur_item, *last_item;

    fprintf(stderr, "Checking sparse matrix (%p) ... ", m);

    cpt_col = 0;
    for (j = 0; j < m->nb_col; j++) {
        if (!m->col[j])
            continue;

        cur_item = m->col[j];
        while (cur_item) {
            cpt_col++;
            last_item = cur_item;
            cur_item = cur_item->next_in_col;
            assert(!(cur_item && cur_item->col_index != j));
            assert(!
                   (cur_item
                    && last_item->line_index >= cur_item->line_index));
        }
    }

    cpt_line = 0;
    for (i = 0; i < m->nb_line; i++) {
        if (!m->line[i])
            continue;

        cur_item = m->line[i];
        while (cur_item) {
            cpt_line++;
            last_item = cur_item;
            cur_item = cur_item->next_in_line;

            assert(!(cur_item && cur_item->line_index != i));
            /*
             * if (cur_item && last_item->col_index >=
             * cur_item->col_index) { fprintf(stderr,
             * "last-item[%ld,%ld] = %f, item[%ld,%ld] = %f\n",
             * last_item->line_index, last_item->col_index,
             * last_item->val, cur_item->line_index,
             * cur_item->col_index, cur_item->val); }
             */
            assert(!
                   (cur_item
                    && last_item->col_index >= cur_item->col_index));

        }
    }

    if (cpt_col != cpt_line) {
        fprintf(stderr, "item found using col = %ld, using line = %ld\n",
                cpt_col, cpt_line);
        return (0);
    }
    fprintf(stderr, "(%ld items) ok\n", cpt_col);
    return (0);
}

void sparse_compute_length(struct sparse_matrix_t *m, char *filename)
{
    long int i;

    double length;

    struct sparse_item_t *cur_item;

    FILE *fd;

    if (!(fd = fopen(filename, "w"))) {
        perror(filename);
        exit(1);
    }
    fprintf(fd, "%ld %ld\n", m->nb_line, m->nb_col);

    for (i = 0; i < m->nb_line; i++) {
        if (!m->line[i])
            continue;

        cur_item = m->line[i];
        length = 0;
        while (cur_item) {
            length += cur_item->val;
            cur_item = cur_item->next_in_line;
        }
        fprintf(fd, "%ld %f\n", i, length);

    }

    fclose(fd);
}

struct sparse_matrix_t *AtransA(struct sparse_matrix_t *A)
{
    struct sparse_matrix_t *AtA;

    long int i, j;

    double sum;

    int sum_updated;

    struct sparse_item_t *last_item = NULL;

    struct sparse_item_t *item1, *item2;

    AtA = new_sparse_matrix(A->nb_col, A->nb_col, SPARSE_COL_LINK);

    for (i = 0; i < A->nb_col; i++) {
        last_item = NULL;
        for (j = 0; j < A->nb_col; j++) {
            item1 = A->col[i];
            item2 = A->col[j];

            sum = 0.;
            sum_updated = 0;
            while (item1 && item2) {
                if (item1->line_index == item2->line_index) {
                    sum_updated = 1;
                    sum += item1->val * item2->val;
                    item2 = item2->next_in_col;
                    item1 = item1->next_in_col;
                    continue;
                }
                if (item1->line_index > item2->line_index) {
                    item2 = item2->next_in_col;
                } else {
                    item1 = item1->next_in_col;
                }
            }
            if (sum_updated) {
                last_item = sparse_set_value(AtA, i, j, sum, last_item);
            }
        }
    }

    return (AtA);
}

double mean_diag_AtA(struct sparse_matrix_t *A)
{
    long int i, j;

    double sum, diag_sum = 0;

    int sum_updated;

    struct sparse_item_t *item1, *item2;

    for (i = 0; i < A->nb_col; i++) {
        j = i;

        item1 = A->col[i];
        item2 = A->col[j];

        sum = 0.;
        sum_updated = 0;
        while (item1 && item2) {
            if (item1->line_index == item2->line_index) {
                sum_updated = 1;
                sum += item1->val * item2->val;
                item2 = item2->next_in_col;
                item1 = item1->next_in_col;
                continue;
            }
            if (item1->line_index > item2->line_index) {
                item2 = item2->next_in_col;
            } else {
                item1 = item1->next_in_col;
            }
        }
        if (sum_updated) {
            diag_sum += sum;
        }
    }
    return (diag_sum / (A->nb_col));
}

void show_sparse_stats(struct sparse_matrix_t *A)
{
    float density;

    density = 100.0 * (double) A->nb_item /
        ((double) A->nb_line * (double) A->nb_col);

    fprintf(stdout, "sparse stats (%p):\n", A);
    fprintf(stdout, "\tsize: %ldx%ld, nb items: %ld\n",
            A->nb_line, A->nb_col, A->nb_item);
    fprintf(stdout, "\tdensity: %f\n", density);
}
