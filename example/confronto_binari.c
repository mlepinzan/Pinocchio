/*
 * compare_binary.c
 *
 * Confronta due file binari interpretandoli come array di double (float64).
 * Calcola statistiche di differenza: max assoluto, RMSE, MAE, e relative.
 *
 * Compilazione:
 *   gcc -O2 -o compare_binary compare_binary.c -lm
 *
 * Uso:
 *   ./compare_binary file1.bin file2.bin
 *
 * Opzionale: puoi specificare un offset in bytes con cui iniziare la lettura:
 *   ./compare_binary file1.bin file2.bin [offset_bytes]
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <errno.h>

/* ------------------------------------------------------------------ */
/* Struttura per le statistiche di confronto                           */
/* ------------------------------------------------------------------ */
typedef struct {
    size_t  n_elements;         /* numero di double confrontati         */
    size_t  n_exact;            /* coppie identiche bit-a-bit           */
    size_t  n_nan_mismatch;     /* uno NaN, l'altro no                  */
    size_t  n_inf_mismatch;     /* uno Inf, l'altro no                  */

    double  max_abs_diff;       /* max |a - b|                          */
    double  max_rel_diff;       /* max |a - b| / max(|a|,|b|, eps)      */
    double  sum_abs_diff;       /* somma |a - b|                        */
    double  sum_sq_diff;        /* somma (a - b)^2                      */
    double  sum_sq_ref;         /* somma b^2  (per errore relativo RMS) */

    size_t  idx_max_abs;        /* indice del massimo scarto assoluto   */
    double  val_a_max_abs;      /* valore in file1 a quell'indice       */
    double  val_b_max_abs;      /* valore in file2 a quell'indice       */
} Stats;

/* ------------------------------------------------------------------ */
/* Lettura di un intero file binario in memoria                        */
/* ------------------------------------------------------------------ */
static double *read_doubles(const char *path, size_t *out_n, size_t offset_bytes)
{
    FILE *fp = fopen(path, "rb");
    if (!fp) {
        fprintf(stderr, "Errore apertura '%s': %s\n", path, strerror(errno));
        return NULL;
    }

    /* dimensione del file */
    if (fseek(fp, 0, SEEK_END) != 0) {
        perror("fseek"); fclose(fp); return NULL;
    }
    long file_size = ftell(fp);
    if (file_size < 0) { perror("ftell"); fclose(fp); return NULL; }

    if ((size_t)file_size < offset_bytes) {
        fprintf(stderr, "Offset (%zu) maggiore della dimensione del file (%ld) per '%s'\n",
                offset_bytes, file_size, path);
        fclose(fp); return NULL;
    }

    size_t data_bytes = (size_t)file_size - offset_bytes;
    size_t n = data_bytes / sizeof(double);

    if (n == 0) {
        fprintf(stderr, "Nessun double leggibile in '%s' (dati: %zu bytes)\n", path, data_bytes);
        fclose(fp); return NULL;
    }

    if (fseek(fp, (long)offset_bytes, SEEK_SET) != 0) {
        perror("fseek"); fclose(fp); return NULL;
    }

    double *buf = (double *)malloc(n * sizeof(double));
    if (!buf) {
        fprintf(stderr, "Allocazione fallita per %zu doubles\n", n);
        fclose(fp); return NULL;
    }

    size_t read_n = fread(buf, sizeof(double), n, fp);
    if (read_n != n) {
        fprintf(stderr, "Lettura parziale da '%s': attesi %zu, letti %zu\n",
                path, n, read_n);
        free(buf); fclose(fp); return NULL;
    }

    fclose(fp);
    *out_n = n;
    return buf;
}

/* ------------------------------------------------------------------ */
/* Confronto elemento per elemento                                     */
/* ------------------------------------------------------------------ */
static void compute_stats(const double *a, const double *b, size_t n, Stats *s)
{
    memset(s, 0, sizeof(*s));
    s->n_elements    = n;
    s->max_abs_diff  = 0.0;
    s->max_rel_diff  = 0.0;
    s->idx_max_abs   = 0;

    const double eps = 2.2204460492503131e-16; /* DBL_EPSILON */

    for (size_t i = 0; i < n; i++) {
        double ai = a[i], bi = b[i];

        /* identici bit-a-bit */
        if (ai == bi) { s->n_exact++; }

        /* NaN mismatch */
        if (isnan(ai) != isnan(bi)) { s->n_nan_mismatch++; continue; }
        if (isnan(ai) && isnan(bi)) { s->n_exact++;         continue; } /* entrambi NaN */

        /* Inf mismatch */
        if (isinf(ai) != isinf(bi)) { s->n_inf_mismatch++; continue; }

        double diff    = ai - bi;
        double abs_diff = fabs(diff);
        double denom   = fmax(fabs(ai), fabs(bi));
        double rel_diff = abs_diff / fmax(denom, eps);

        s->sum_abs_diff += abs_diff;
        s->sum_sq_diff  += diff * diff;
        s->sum_sq_ref   += bi  * bi;

        if (abs_diff > s->max_abs_diff) {
            s->max_abs_diff  = abs_diff;
            s->idx_max_abs   = i;
            s->val_a_max_abs = ai;
            s->val_b_max_abs = bi;
        }
        if (rel_diff > s->max_rel_diff) {
            s->max_rel_diff = rel_diff;
        }
    }
}

/* ------------------------------------------------------------------ */
/* Stampa risultati                                                    */
/* ------------------------------------------------------------------ */
static void print_stats(const Stats *s, const char *f1, const char *f2)
{
    size_t n = s->n_elements;
    double mae  = (n > 0) ? s->sum_abs_diff / (double)n : 0.0;
    double rmse = (n > 0) ? sqrt(s->sum_sq_diff / (double)n) : 0.0;
    double nrmse = (s->sum_sq_ref > 0.0)
                 ? sqrt(s->sum_sq_diff / s->sum_sq_ref)
                 : 0.0;

    printf("\n");
    printf("========================================================\n");
    printf("  Confronto binario a double-precision\n");
    printf("  FILE 1 : %s\n", f1);
    printf("  FILE 2 : %s\n", f2);
    printf("========================================================\n");
    printf("  Elementi confrontati      : %zu\n",   n);
    printf("  Identici (bit-a-bit)      : %zu  (%.4f%%)\n",
           s->n_exact, 100.0 * s->n_exact / (double)n);
    printf("  NaN mismatch              : %zu\n",   s->n_nan_mismatch);
    printf("  Inf mismatch              : %zu\n",   s->n_inf_mismatch);
    printf("--------------------------------------------------------\n");
    printf("  Max scarto assoluto       : %24.17e\n", s->max_abs_diff);
    printf("    @ indice                : %zu\n",   s->idx_max_abs);
    printf("    file1[i]                : %24.17e\n", s->val_a_max_abs);
    printf("    file2[i]                : %24.17e\n", s->val_b_max_abs);
    printf("  Max scarto relativo       : %24.17e\n", s->max_rel_diff);
    printf("--------------------------------------------------------\n");
    printf("  MAE  (Mean Abs Error)     : %24.17e\n", mae);
    printf("  RMSE (Root Mean Sq Error) : %24.17e\n", rmse);
    printf("  NRMSE (norm. su ||b||_2)  : %24.17e\n", nrmse);
    printf("========================================================\n\n");

    /* Verdetto rapido */
    if (s->n_exact == n) {
        printf("  >> I due file sono IDENTICI bit-a-bit.\n\n");
    } else if (s->max_abs_diff < 1e-15) {
        printf("  >> Differenze entro la precisione macchina (< 1e-15).\n\n");
    } else if (s->max_abs_diff < 1e-10) {
        printf("  >> Differenze molto piccole (< 1e-10), probabilmente numeriche.\n\n");
    } else {
        printf("  >> Differenze significative rilevate.\n\n");
    }
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */
int main(int argc, char *argv[])
{
    if (argc < 3 || argc > 4) {
        fprintf(stderr, "Uso: %s <file1.bin> <file2.bin> [offset_bytes]\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char *path1 = argv[1];
    const char *path2 = argv[2];
    size_t offset = 0;

    if (argc == 4) {
        char *end;
        offset = (size_t)strtoull(argv[3], &end, 10);
        if (*end != '\0') {
            fprintf(stderr, "Offset non valido: '%s'\n", argv[3]);
            return EXIT_FAILURE;
        }
    }

    size_t n1 = 0, n2 = 0;
    double *buf1 = read_doubles(path1, &n1, offset);
    double *buf2 = read_doubles(path2, &n2, offset);

    if (!buf1 || !buf2) {
        free(buf1); free(buf2);
        return EXIT_FAILURE;
    }

    size_t n = n1;
    if (n1 != n2) {
        fprintf(stderr,
            "Attenzione: file di dimensione diversa (%zu vs %zu doubles).\n"
            "Il confronto si limita ai primi %zu elementi.\n",
            n1, n2, (n1 < n2 ? n1 : n2));
        n = (n1 < n2) ? n1 : n2;
    }

    Stats s;
    compute_stats(buf1, buf2, n, &s);
    print_stats(&s, path1, path2);

    free(buf1);
    free(buf2);
    return EXIT_SUCCESS;
}
