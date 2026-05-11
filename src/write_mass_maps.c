/*
 * MASS_MAPS: per-segment PLC entries -> HEALPix maps (products-only)
 * ------------------------------------------------------------------
 * Detect entries of mass elements into a past light-cone (PLC) between adjacent
 * fragmentation outputs and accumulate counts on one HEALPix map per segment.
 *
 * Pipeline (high level):
 *   - Cheap, conservative block-level culls (radial shell; conservative angular bound) before touching particles.
 *   - For passing blocks/particles, compute previous/current Eulerian positions (q + LPT displacements),
 *     detect PLC boundary crossings via H = chi(z) - r sign change, interpolate the entry position,
 *     perform pixel-based angular selection (see below), and accumulate into the map.
 *
 * Angular selection (pixel-based):
 *   - Convert the entry position to a HEALPix pixel at the configured NSIDE using the PLC-oriented basis.
 *   - For partial sky (aperture < 180 deg), accept iff the pixel index lies in the cap prefix [0 .. Ncap-1]
 *     in RING ordering. This matches downstream "bins populated" semantics (mm!=UNSEEN).
 *   - For full sky, all pixels are accepted.
 *   - The previous geometric angle test is retained only for diagnostics (see mass_maps_entry_inside_aperture).
 *
 * Configuration overview:
 *   - PLC                       Enable past-light-cone logic (required).
 *   - MASS_MAPS                 Enable mass maps (segment accumulation and FITS output).
 *   - PLCAperture [deg]         Angular radius of the PLC cone; <180 selects a spherical cap map.
 *   - PLC axis/center           Provided by PLC setup; used for pixelization basis and origin.
 *   - NSIDE                     HEALPix resolution for maps (RING ordering).
 *   - MASS_MAPS_* toggles:
 *       MASS_MAPS_FULLSKY_OUTPUT   If nonzero, also write full-sky IMAGE HDU even when aperture<180 (default 0).
 *       MASS_MAPS_CAP_DIAG         Emit cap-boundary diagnostics (prefix misses and angle extrema).
 *       MASS_MAPS_THREAD_ACCUM     Use per-thread HEALPix accumulators and merge at the end to reduce atomics.
 *
 * Output geometry:
 *   - Full-sky (aperture=180): a full-sphere RING map with NPIX = 12*NSIDE^2.
 *   - Partial-sky (aperture<180): a RING-ordered spherical cap with Ncap pixels (prefix 0..Ncap-1).
 *     If the computed Ncap is 0 (too small aperture for the chosen NSIDE), maps are disabled with a warning.
 *   - FITS headers include basic HEALPix metadata, PLC axis (unit vector), FILTER state, and map geometry.
 */
#include "pinocchio.h"
#ifdef MASS_MAPS
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifndef PLC
#error "MASS_MAPS requires PLC"
#endif

#ifndef RECOMPUTE_DISPLACEMENTS
#error "MASS_MAPS requires RECOMPUTE_DISPLACEMENTS"
#endif

/* Ensure variable pad macros are available before any use (histogram header, etc.) */
#ifndef MASS_MAPS_VARPAD_FACTOR
#define MASS_MAPS_VARPAD_FACTOR 20.0
#endif
#ifndef MASS_MAPS_VARPAD_WARN_FRAC
#define MASS_MAPS_VARPAD_WARN_FRAC (0.05)
#endif

#ifndef HAVE_CFITSIO
#error "MASS_MAPS requires HAVE_CFITSIO"
#endif
#include <fitsio.h>
#ifdef MASS_MAPS
#include <chealpix.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

/* Forward declaration: needed because mass_maps_init_sheets() calls it on failure */
void mass_maps_free_sheets(void);

/* ---------------------------------------------------------- */
/* Optional displacement histogram (compile-time, default OFF) */
/* ---------------------------------------------------------- */
#ifdef DISP_HISTOGRAM
#ifndef DISP_HIST_NBINS
#define DISP_HIST_NBINS 100
#endif
#ifndef DISP_HIST_DMIN
#define DISP_HIST_DMIN 0.0
#endif
#ifndef DISP_HIST_DMAX
#define DISP_HIST_DMAX 100.0
#endif
#ifndef DISP_HIST_STRIDE
#define DISP_HIST_STRIDE 1
#endif
#ifndef DISP_HIST_AT_PREV
#define DISP_HIST_AT_PREV 1
#endif
#ifndef DISP_HIST_AT_CURR
#define DISP_HIST_AT_CURR 1
#endif
/** Write an ASCII histogram file on rank 0
 *
 * mass_maps_disp_hist_write_file
 * ------------------------------
 * Purpose
 *   Write the global displacement histogram |x(q,z)-q| (measured in InterPartDist
 *   units) to an ASCII file. Only MPI rank 0 performs I/O.
 *
 * Inputs
 *   - z_target         Redshift at which the endpoint histogram was sampled.
 *   - hist_global      Array of length DISP_HIST_NBINS with global bin counts
 *                      (already MPI-reduced across tasks).
 *   - underflow_global Global count of samples with d < DMIN.
 *   - overflow_global  Global count of samples with d >= DMAX.
 *   - sampled_global   Global number of samples considered (including under/overflow).
 *
 * Behavior and file format
 *   - Filename: pinocchio.&lt;z&gt;.&lt;RunFlag&gt;.disp_hist.out where &lt;z&gt; is formatted with
 *     4 decimals for stable sorting/comparisons.
 *   - Header: prints grid/analysis provenance and the variable AABB pad used by
 *     mass-maps at this redshift, reported in Mpc, Mpc/h, and PAD_CELLS
 *     (PAD_PHYS / InterPartDist) for diagnostics.
 *   - Data rows (one per bin):
 *       bin_lo  bin_hi  count  fraction  cum_fraction
 *
 * Units
 *   Histogram d values are in InterPartDist units.
 *
 * Thread-safety / MPI
 *   Inputs must be globally reduced; only rank 0 writes.
 */
static void mass_maps_disp_hist_write_file(double z_target,
                                           const unsigned long long *hist_global,
                                           unsigned long long underflow_global,
                                           unsigned long long overflow_global,
                                           unsigned long long sampled_global)
{
  /* Only root performs I/O; other ranks return immediately */
  if (ThisTask != 0)
    return;
  /* Build filename: pinocchio.<z>.<RunFlag>.disp_hist.out with fixed-width z */
  char fname[4 * LBLENGTH];
  snprintf(fname, sizeof(fname), "pinocchio.%.4f.%s.disp_hist.out", z_target, params.RunFlag);
  FILE *fd = fopen(fname, "w");
  if (!fd)
  {
    fprintf(stderr, "Task 0 ERROR: cannot open %s for writing.\n", fname);
    return;
  }
  double dmin = (double)DISP_HIST_DMIN;
  double dmax = (double)DISP_HIST_DMAX;
  int nb = (int)DISP_HIST_NBINS;
  int stride = (int)DISP_HIST_STRIDE;
  /* Compute variable buffer pad at this redshift for header context */
  double sigma8 = params.Sigma8;
  double factor = MASS_MAPS_VARPAD_FACTOR;
  double res = params.InterPartDist; /* true Mpc */
  double h = params.Hubble100;
  double pad_phys_h = factor * pow(sigma8 / 0.8, 1.5) / (1.0 + z_target); /* Mpc/h */
  double pad_phys = (h > 0.0) ? (pad_phys_h / h) : pad_phys_h;            /* true Mpc */
  double pad_cells = (res > 0.0) ? (pad_phys / res) : 0.0;                /* dimensionless */
  /* Header with provenance and units */
  fprintf(fd, "# Displacement histogram |x(q,z)-q| in InterPartDist units\n");
  fprintf(fd, "# z=%.8f NBINS=%d DMIN=%.6g DMAX=%.6g SPACING=linear STRIDE=%d\n", z_target, nb, dmin, dmax, stride);
  fprintf(fd, "# GRID=(%d,%d,%d) InterPartDist(Mpc)=%.9g LPT_ORDER=%d MPI_TASKS=%d RUN=%s\n",
          (int)MyGrids[0].GSglobal[_x_], (int)MyGrids[0].GSglobal[_y_], (int)MyGrids[0].GSglobal[_z_],
          params.InterPartDist,
#ifdef THREE_LPT
          3,
#else
#ifdef TWO_LPT
          2,
#else
          1,
#endif
#endif
          NTasks, params.RunFlag);
  fprintf(fd, "# VARPAD: FACTOR=%.6g SIGMA8=%.6g PAD_PHYS(Mpc)=%.9g PAD_PHYS(Mpc/h)=%.9g PAD_CELLS=%.9g\n",
          factor, sigma8, pad_phys, pad_phys_h, pad_cells);
  fprintf(fd, "# SAMPLED_GLOBAL=%llu UNDERFLOW=%llu OVERFLOW=%llu\n",
          (unsigned long long)sampled_global,
          (unsigned long long)underflow_global,
          (unsigned long long)overflow_global);
  fprintf(fd, "# Columns: bin_lo bin_hi count fraction cum_fraction\n");
  /* Emit rows */
  double width = (nb > 0) ? ((dmax - dmin) / (double)nb) : 0.0;
  unsigned long long cum = 0ULL;
  for (int b = 0; b < nb; ++b)
  {
    double blo = dmin + width * (double)b;
    double bhi = blo + width;
    unsigned long long c = hist_global[b];
    cum += c;
    double frac = (sampled_global > 0ULL) ? ((double)c / (double)sampled_global) : 0.0;
    double cfrac = (sampled_global > 0ULL) ? ((double)cum / (double)sampled_global) : 0.0;
    fprintf(fd, "%.9g %.9g %llu %.9g %.9g\n", blo, bhi, (unsigned long long)c, frac, cfrac);
  }
  fclose(fd);
  if (internal.verbose_level >= VDIAG)
  {
    printf("[%s] DISP_HIST: wrote %s (z=%.6f nbins=%d sampled=%llu uf=%llu of=%llu)\n",
           fdate(), fname, z_target, nb,
           (unsigned long long)sampled_global,
           (unsigned long long)underflow_global,
           (unsigned long long)overflow_global);
  }
}

/**
 * mass_maps_disp_hist_accumulate_endpoint
 * ---------------------------------------
 * Purpose
 *   Accumulate a local displacement histogram for one endpoint selector:
 *   use_curr=0 for the previous endpoint ("prev"), use_curr=1 for the current
 *   endpoint ("curr"). Displacements are |x(q,z)-q| in InterPartDist units
 *   (i.e., multiples of the InterPartDist). Underflow/overflow and
 *   total sampled counters are tracked alongside per-bin counts.
 *
 * Parameters
 *   - use_curr         0: use products.*_prev fields; 1: use current products.* fields.
 *   - hist_local       Output array of length DISP_HIST_NBINS (must be zero-initialized by caller).
 *   - underflow_local  Output: incremented by samples with d < DMIN.
 *   - overflow_local   Output: incremented by samples with d >= DMAX.
 *   - sampled_local    Output: incremented by total samples considered (after stride filter).
 *
 * Behavior
 *   - Iterates over this task's local FFT tile; optional sub-sampling via DISP_HIST_STRIDE
 *     by applying a deterministic (i+j+k) % stride == 0 filter on global indices.
 *   - Builds endpoint displacement magnitude per lattice site including compiled LPT orders.
 *   - Uses per-thread private histograms and reduces into hist_local to minimize contention.
 *
 * Threading
 *   OpenMP is supported; per-thread histograms are accumulated with a final reduction.
 *   No MPI calls happen here; caller performs cross-task reductions/writes.
 **/
static void mass_maps_disp_hist_accumulate_endpoint(int use_curr,
                                                    unsigned long long *hist_local,
                                                    unsigned long long *underflow_local,
                                                    unsigned long long *overflow_local,
                                                    unsigned long long *sampled_local)
{
  const int nb = (int)DISP_HIST_NBINS;
  const double dmin = (double)DISP_HIST_DMIN;
  const double dmax = (double)DISP_HIST_DMAX;
  const int stride = (int)DISP_HIST_STRIDE;
  if (nb <= 0 || dmax <= dmin)
    return;
  int nx = (int)MyGrids[0].GSlocal[_x_];
  int ny = (int)MyGrids[0].GSlocal[_y_];
  int nz = (int)MyGrids[0].GSlocal[_z_];
  int gx0 = (int)MyGrids[0].GSstart[_x_];
  int gy0 = (int)MyGrids[0].GSstart[_y_];
  int gz0 = (int)MyGrids[0].GSstart[_z_];
  double invw = (double)nb / (dmax - dmin);

  unsigned long long uf = 0ULL, of = 0ULL, tot = 0ULL;

#ifdef _OPENMP
  int nthreads = omp_get_max_threads();
#else
  int nthreads = 1;
#endif
  unsigned long long *thread_hist = (unsigned long long *)calloc((size_t)nthreads * (size_t)nb, sizeof(unsigned long long));
  if (!thread_hist)
    return;

  /* Parallel over local tile */
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
#ifdef _OPENMP
    int tid = omp_get_thread_num();
#else
    int tid = 0;
#endif
    unsigned long long *h = thread_hist + (size_t)tid * (size_t)nb;

#ifdef _OPENMP
#pragma omp for collapse(3) schedule(static) reduction(+ : uf, of, tot)
#endif
    for (int li = 0; li < nx; ++li)
      for (int lj = 0; lj < ny; ++lj)
        for (int lk = 0; lk < nz; ++lk)
        {
          if (stride > 1)
          {
            int gi = gx0 + li, gj = gy0 + lj, gk = gz0 + lk;
            /* Simple deterministic sub-sampling */
            if (((gi + gj + gk) % stride) != 0)
              continue;
          }
          unsigned int idx = COORD_TO_INDEX(li, lj, lk, MyGrids[0].GSlocal);
          /* Build displacement vector in grid units (InterPartDist units) */
          double vx = 0.0, vy = 0.0, vz = 0.0;
          if (use_curr)
          {
            vx = products[idx].Vel[0];
            vy = products[idx].Vel[1];
            vz = products[idx].Vel[2];
#ifdef TWO_LPT
            vx += products[idx].Vel_2LPT[0];
            vy += products[idx].Vel_2LPT[1];
            vz += products[idx].Vel_2LPT[2];
#ifdef THREE_LPT
            vx += products[idx].Vel_3LPT_1[0] + products[idx].Vel_3LPT_2[0];
            vy += products[idx].Vel_3LPT_1[1] + products[idx].Vel_3LPT_2[1];
            vz += products[idx].Vel_3LPT_1[2] + products[idx].Vel_3LPT_2[2];
#endif
#endif
          }
          else
          {
            vx = products[idx].Vel_prev[0];
            vy = products[idx].Vel_prev[1];
            vz = products[idx].Vel_prev[2];
#ifdef TWO_LPT
            vx += products[idx].Vel_2LPT_prev[0];
            vy += products[idx].Vel_2LPT_prev[1];
            vz += products[idx].Vel_2LPT_prev[2];
#ifdef THREE_LPT
            vx += products[idx].Vel_3LPT_1_prev[0] + products[idx].Vel_3LPT_2_prev[0];
            vy += products[idx].Vel_3LPT_1_prev[1] + products[idx].Vel_3LPT_2_prev[1];
            vz += products[idx].Vel_3LPT_1_prev[2] + products[idx].Vel_3LPT_2_prev[2];
#endif
#endif
          }
          double d = sqrt(vx * vx + vy * vy + vz * vz); /* in InterPartDist units */
          ++tot;
          if (d < dmin)
          {
            ++uf;
            continue;
          }
          if (!(d < dmax))
          {
            ++of;
            continue;
          }
          int bin = (int)floor((d - dmin) * invw);
          if (bin < 0)
            bin = 0;
          if (bin >= nb)
            bin = nb - 1;
          h[bin]++;
        }
  }

  /* Reduce thread hist to task hist */
  for (int t = 0; t < nthreads; ++t)
  {
    unsigned long long *h = thread_hist + (size_t)t * (size_t)nb;
    for (int b = 0; b < nb; ++b)
      hist_local[b] += h[b];
  }
  free(thread_hist);
  *underflow_local += uf;
  *overflow_local += of;
  *sampled_local += tot;
}

/**
 * mass_maps_disp_hist_on_segment
 * ------------------------------
 * Purpose
 *   Orchestrate per-segment displacement histograms at the two segment endpoints
 *   (previous and current). Optionally computes for each endpoint depending on
 *   DISP_HIST_AT_PREV / DISP_HIST_AT_CURR, performs MPI reductions to rank 0,
 *   and writes ASCII histogram files via mass_maps_disp_hist_write_file.
 *
 * Parameters
 *   - segment_index  Index of the current fragmentation segment (not used for I/O;
 *                    retained for context and future use).
 *   - z_prev         Redshift at previous endpoint of the segment.
 *   - z_curr         Redshift at current endpoint of the segment.
 *
 * Behavior
 *   - Allocates per-endpoint task-local histograms (length DISP_HIST_NBINS) when enabled.
 *   - Invokes mass_maps_disp_hist_accumulate_endpoint for prev/curr to fill local bins and
 *     counters (underflow/overflow/sampled) in InterPartDist units.
 *   - Reduces bins and counters across MPI tasks to rank 0 and writes one file per
 *     enabled endpoint using z_prev/z_curr in the filename and header.
 *   - Frees temporary buffers.
 *
 * Notes
 *   - Units: The accumulator operates in InterPartDist (true Mpc) units; the writer prints
 *     both InterPartDist and variable pad context in headers. No unit conversion happens here.
 *   - Concurrency: This function uses MPI_Reduce but no OpenMP loops directly.
 **/
/* Internal helper implementing the documented per-segment histogram workflow */
static void mass_maps_disp_hist_on_segment(int segment_index, double z_prev, double z_curr)
{
  (void)segment_index;
  const int nb = (int)DISP_HIST_NBINS;
  if (nb <= 0 || DISP_HIST_DMAX <= DISP_HIST_DMIN)
    return;

  /* Allocate task-local buffers (conditional on compile-time toggles) */
  unsigned long long *hist_prev_local = NULL, *hist_curr_local = NULL;
  unsigned long long uf_prev_local = 0ULL, of_prev_local = 0ULL, tot_prev_local = 0ULL;
  unsigned long long uf_curr_local = 0ULL, of_curr_local = 0ULL, tot_curr_local = 0ULL;

  if (DISP_HIST_AT_PREV)
    hist_prev_local = (unsigned long long *)calloc((size_t)nb, sizeof(unsigned long long));
  if (DISP_HIST_AT_CURR)
    hist_curr_local = (unsigned long long *)calloc((size_t)nb, sizeof(unsigned long long));

  if (hist_prev_local)
    mass_maps_disp_hist_accumulate_endpoint(0, hist_prev_local, &uf_prev_local, &of_prev_local, &tot_prev_local);
  if (hist_curr_local)
    mass_maps_disp_hist_accumulate_endpoint(1, hist_curr_local, &uf_curr_local, &of_curr_local, &tot_curr_local);

  /* Global reductions (MPI) and writes on rank 0 */
  if (hist_prev_local)
  {
    unsigned long long *hist_prev_global = NULL;
    if (!ThisTask)
      hist_prev_global = (unsigned long long *)calloc((size_t)nb, sizeof(unsigned long long));
    MPI_Reduce(hist_prev_local, (!ThisTask ? hist_prev_global : NULL), nb, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    unsigned long long uf_prev_global = 0ULL, of_prev_global = 0ULL, tot_prev_global = 0ULL;
    MPI_Reduce(&uf_prev_local, &uf_prev_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&of_prev_local, &of_prev_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tot_prev_local, &tot_prev_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    mass_maps_disp_hist_write_file(z_prev, hist_prev_global, uf_prev_global, of_prev_global, tot_prev_global);
    if (!ThisTask && hist_prev_global)
      free(hist_prev_global);
  }
  if (hist_curr_local)
  {
    unsigned long long *hist_curr_global = NULL;
    if (!ThisTask)
      hist_curr_global = (unsigned long long *)calloc((size_t)nb, sizeof(unsigned long long));
    MPI_Reduce(hist_curr_local, (!ThisTask ? hist_curr_global : NULL), nb, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    unsigned long long uf_curr_global = 0ULL, of_curr_global = 0ULL, tot_curr_global = 0ULL;
    MPI_Reduce(&uf_curr_local, &uf_curr_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&of_curr_local, &of_curr_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tot_curr_local, &tot_curr_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    mass_maps_disp_hist_write_file(z_curr, hist_curr_global, uf_curr_global, of_curr_global, tot_curr_global);
    if (!ThisTask && hist_curr_global)
      free(hist_curr_global);
  }
  if (hist_prev_local)
    free(hist_prev_local);
  if (hist_curr_local)
    free(hist_curr_local);
}
#endif /* DISP_HISTOGRAM */

/* ---------------------------------------------------------- */
/* Optional: carry ZPLC>ZACC filter diagnostics into FITS     */
/* ---------------------------------------------------------- */
typedef struct
{
  /* Totals */
  long long considered;
  long long excluded;
  long long included;
  double frac_excluded;
  /* Ranges and consistency */
  double zplc_min, zplc_max;
  double zacc_min, zacc_max;
  long long zplc_out_of_segment;
  /* Segment id for context */
  int segment_index;
  /* Flag */
  int present;
} MassMapsFilterDiag;

static MassMapsFilterDiag g_massmaps_filter_diag = {0};

/**
 * mass_maps_write_filter_keywords
 * -------------------------------
 * Purpose
 *   If MASS_MAPS_FILTER_UNCOLLAPSED is enabled and diagnostics are present for
 *   the given segment, write ZPLC>ZACC filter metadata into the segment FITS
 *   header. This makes filtering decisions visible and reproducible in outputs.
 *
 * Parameters
 *   - fptr           Open CFITSIO file handle for the segment map being written.
 *   - segment_index  Segment index that the metadata refers to; used to guard
 *                    against stale/mismatched diagnostics.
 *
 * Behavior
 *   - Validates that diagnostics are present and match segment_index.
 *   - Writes the following header keys on rank 0 (caller writes on root only):
 *       FILTER   = 'ZPLC>ZACC'              (filter condition)
 *       ZF_CONS  = considered entries       (TLONGLONG)
 *       ZF_EXCL  = excluded entries         (TLONGLONG)
 *       ZF_INCL  = included entries         (TLONGLONG)
 *       ZF_FEXCL = excluded fraction [0..1] (TDOUBLE)
 *       ZPLC_MIN / ZPLC_MAX                 (TDOUBLE)
 *       ZACC_MIN / ZACC_MAX                 (TDOUBLE)
 *       ZPLC_OOS  = z_plc outside seg (tol) (TLONGLONG)
 *     Also appends a concise HISTORY line summarizing counts.
 *
 * Notes
 *   - No MPI calls occur here; inputs are expected to be globally reduced
 *     beforehand during segment processing.
 *   - If MASS_MAPS_FILTER_UNCOLLAPSED is not defined, this is a no-op.
 **/
static void mass_maps_write_filter_keywords(fitsfile *fptr, int segment_index)
{
#ifdef MASS_MAPS_FILTER_UNCOLLAPSED
  /* Basic guards: null file, no diagnostics, or mismatched segment */
  if (!fptr)
    return;
  if (!g_massmaps_filter_diag.present)
    return;
  if (g_massmaps_filter_diag.segment_index != segment_index)
  {
    /* Stale or mismatched diagnostics; skip to avoid wrong metadata */
    return;
  }
  int status = 0;
  {
    char filter_desc[] = "ZPLC>ZACC";
    fits_update_key(fptr, TSTRING, "FILTER", filter_desc, "Mass-maps inclusion filter", &status);
  }
  {
    long long considered = g_massmaps_filter_diag.considered;
    long long excluded = g_massmaps_filter_diag.excluded;
    long long included = g_massmaps_filter_diag.included;
    double frac_ex = g_massmaps_filter_diag.frac_excluded;
    fits_update_key(fptr, TLONGLONG, "ZF_CONS", &considered, "Entries considered by filter", &status);
    fits_update_key(fptr, TLONGLONG, "ZF_EXCL", &excluded, "Entries excluded by filter", &status);
    fits_update_key(fptr, TLONGLONG, "ZF_INCL", &included, "Entries included after filter", &status);
    fits_update_key(fptr, TDOUBLE, "ZF_FEXCL", &frac_ex, "Excluded fraction (0..1)", &status);
  }
  {
    double zplc_min = g_massmaps_filter_diag.zplc_min;
    double zplc_max = g_massmaps_filter_diag.zplc_max;
    double zacc_min = g_massmaps_filter_diag.zacc_min;
    double zacc_max = g_massmaps_filter_diag.zacc_max;
    long long outseg = g_massmaps_filter_diag.zplc_out_of_segment;
    fits_update_key(fptr, TDOUBLE, "ZPLC_MIN", &zplc_min, "Min PLC redshift encountered", &status);
    fits_update_key(fptr, TDOUBLE, "ZPLC_MAX", &zplc_max, "Max PLC redshift encountered", &status);
    fits_update_key(fptr, TDOUBLE, "ZACC_MIN", &zacc_min, "Min collapse redshift (ZACC)", &status);
    fits_update_key(fptr, TDOUBLE, "ZACC_MAX", &zacc_max, "Max collapse redshift (ZACC)", &status);
    fits_update_key(fptr, TLONGLONG, "ZPLC_OOS", &outseg, "Count of z_plc outside segment bounds (tolerance)", &status);
  }
  {
    /* Add concise HISTORY lines */
    char hist[128];
    snprintf(hist, sizeof(hist), "FILTER: seg=%d ZPLC>ZACC considered=%lld excluded=%lld included=%lld",
             segment_index,
             (long long)g_massmaps_filter_diag.considered,
             (long long)g_massmaps_filter_diag.excluded,
             (long long)g_massmaps_filter_diag.included);
    fits_write_history(fptr, hist, &status);
  }
  if (status)
  {
    /* Non-fatal: report but continue */
    fits_report_error(stderr, status);
  }
#else
  (void)fptr;
  (void)segment_index;
#endif
}

/**
 * mass_maps_write_cosmology_keywords
 * -----------------------------------
 * Purpose
 *   Write the main cosmological parameters from `params` into the current HDU
 *   of the given FITS file so that downstream consumers can read the assumed
 *   cosmology directly from the map header.
 *
 * Parameters
 *   - fptr  Open CFITSIO file handle positioned at the HDU to annotate.
 *
 * Notes
 *   - Non-fatal on error: reports via fits_report_error but does not abort.
 **/
static void mass_maps_write_cosmology_keywords(fitsfile *fptr)
{
  if (!fptr)
    return;
  int status = 0;
  double val;

  val = params.Omega0;
  fits_update_key(fptr, TDOUBLE, "COS_OM0", &val, "Omega_m at z=0", &status);
  val = params.OmegaLambda;
  fits_update_key(fptr, TDOUBLE, "COS_OL0", &val, "Omega_Lambda at z=0", &status);
  val = params.OmegaBaryon;
  fits_update_key(fptr, TDOUBLE, "COS_OB0", &val, "Omega_baryon at z=0", &status);
  val = params.Hubble100;
  fits_update_key(fptr, TDOUBLE, "COS_H100", &val, "H0 / (100 km/s/Mpc)", &status);
  val = params.Sigma8;
  fits_update_key(fptr, TDOUBLE, "COS_S8", &val, "sigma8 normalization", &status);
  val = params.PrimordialIndex;
  fits_update_key(fptr, TDOUBLE, "COS_NS", &val, "Primordial spectral index", &status);
  val = params.DEw0;
  fits_update_key(fptr, TDOUBLE, "COS_W0", &val, "Dark-energy w0", &status);
  val = params.DEwa;
  fits_update_key(fptr, TDOUBLE, "COS_WA", &val, "Dark-energy wa", &status);

  {
    char hist[128];
    snprintf(hist, sizeof(hist),
             "COSMOLOGY: Om0=%.4f OL0=%.4f Ob0=%.4f h=%.4f s8=%.4f ns=%.4f w0=%.4f wa=%.4f",
             params.Omega0, params.OmegaLambda, params.OmegaBaryon,
             params.Hubble100, params.Sigma8, params.PrimordialIndex,
             params.DEw0, params.DEwa);
    fits_write_history(fptr, hist, &status);
  }

  if (status)
  {
    /* Non-fatal: report but continue */
    fits_report_error(stderr, status);
  }
}

#ifdef SNAPSHOT
/* Forward declaration: back-distribute ZACC (and group_ID) from fragment sub-boxes to products FFT tiles */
int distribute_back(void);
int my_distribute_back(void);
int distribute_back_alltoall(void);
#endif

/* Local constants */
#ifndef MASS_MAPS_Z_EPS
#define MASS_MAPS_Z_EPS 1e-8
#endif
/* Crossing epsilon for alpha denominator */
#ifndef MASS_MAPS_F_EPS
#define MASS_MAPS_F_EPS 1e-7
#endif
/* Culling defaults  */
#ifndef MASS_MAPS_BLOCKS
#define MASS_MAPS_BLOCKS 8
#endif
/* Default per-axis block counts fall back to MASS_MAPS_BLOCKS if not provided */
#ifndef MASS_MAPS_BLOCKS_X
#define MASS_MAPS_BLOCKS_X MASS_MAPS_BLOCKS
#endif
#ifndef MASS_MAPS_BLOCKS_Y
#define MASS_MAPS_BLOCKS_Y MASS_MAPS_BLOCKS
#endif
#ifndef MASS_MAPS_BLOCKS_Z
#define MASS_MAPS_BLOCKS_Z MASS_MAPS_BLOCKS
#endif

/* Forward decls for helpers provided elsewhere */
void set_weight(pos_data *);

/* Global sheet array */
MassSheet *MassSheets = NULL;
int NMassSheets = 0;
double *MassMapBoundaryZ = NULL;
double *MassMapBoundaryChi = NULL;
double *MassMapBoundaryDA = NULL;
/* HEALPix maps: one map per mass sheet (segment). */
static int MassMapNSIDE_current = 0;
static long MassMapNPIX = 0;
static double *MassMapSegmentMaps = NULL; /* length = NMassSheets * MassMapNPIX */
/* Optional diagnostics for cap prefix boundary investigation */
#ifndef MASS_MAPS_CAP_DIAG
#define MASS_MAPS_CAP_DIAG 0
#endif
#if MASS_MAPS_CAP_DIAG
static unsigned long long MassMapCapPrefixMiss = 0ULL; /* entries inside aperture but ipix >= Ncap */
static double MassMapCapMissAngleMax = 0.0;            /* max angle (deg) among missed pixels */
static double MassMapCapAcceptAngleMax = 0.0;          /* max angle (deg) among accepted pixels */
static int MassMapCapDiagActive = 1;                   /* enable instrumentation (set 0 to disable) */
#else
/* Stubs to satisfy references when diagnostics are disabled; optimized out at runtime */
#define MassMapCapDiagActive 0
static unsigned long long MassMapCapPrefixMiss = 0ULL;
static double MassMapCapMissAngleMax = 0.0;
static double MassMapCapAcceptAngleMax = 0.0;
#endif

static double *MassMapReduceBuffer = NULL; /* root-only temporary buffer for reductions */
/* Per-segment replication candidate list (radial overlap with segment chi span) */
static int *MassMapSegmentReplications = NULL;
static int MassMapSegmentReplicationCount = 0;
static int MassMapSegmentReplicationsCapacity = 0;

/* Cached PLC center (phys) and sky basis (e1,e2,e3) for pixelization */
static double MassMapPLC_CenterPhys[3] = {0.0, 0.0, 0.0};
static double MassMapPLC_e1[3] = {1.0, 0.0, 0.0};
static double MassMapPLC_e2[3] = {0.0, 1.0, 0.0};
static double MassMapPLC_e3[3] = {0.0, 0.0, 1.0};
static int MassMapPLC_BasisReady = 0;
/* Cache of cap size for current NSIDE (when using partial-sky) */
static long MassMapNCAP = 0;

/* Optional: keep full-sky output (IMAGE HDU) even if aperture < 180 */
#ifndef MASS_MAPS_FULLSKY_OUTPUT
/* Undefine or set to 0 to use compact partial-sky by default */
#define MASS_MAPS_FULLSKY_OUTPUT 0
#endif

/**
 * mass_maps_cache_plc_basis_and_center
 * ------------------------------------
 * Purpose
 *   Lazily cache the PLC center in physical units (true Mpc) and the
 *   orthonormal sky basis vectors (e1,e2,e3) used for pixelization and
 *   aperture tests. Avoids recomputing these per call in hot paths.
 *
 * Behavior
 *   - If the cache is already initialized (MassMapPLC_BasisReady!=0), returns.
 *   - Otherwise, converts plc.center (grid units) to physical units via
 *     InterPartDist, and copies the precomputed unit basis vectors
 *     (plc.xvers, plc.yvers, plc.zvers) into local statics.
 *
 * Notes
 *   - Units: positions scaled by params.InterPartDist (true Mpc).
 *   - Threading: benign double-init is possible if called concurrently, but
 *     values are idempotent and identical, so no functional impact.
 */
static inline void mass_maps_cache_plc_basis_and_center(void)
{
  if (MassMapPLC_BasisReady)
    return;
  MassMapPLC_CenterPhys[0] = plc.center[0] * params.InterPartDist;
  MassMapPLC_CenterPhys[1] = plc.center[1] * params.InterPartDist;
  MassMapPLC_CenterPhys[2] = plc.center[2] * params.InterPartDist;
  /* Use PLC basis computed at initialization (orthonormal and unit) */
  MassMapPLC_e1[0] = plc.xvers[0];
  MassMapPLC_e1[1] = plc.xvers[1];
  MassMapPLC_e1[2] = plc.xvers[2];
  MassMapPLC_e2[0] = plc.yvers[0];
  MassMapPLC_e2[1] = plc.yvers[1];
  MassMapPLC_e2[2] = plc.yvers[2];
  MassMapPLC_e3[0] = plc.zvers[0];
  MassMapPLC_e3[1] = plc.zvers[1];
  MassMapPLC_e3[2] = plc.zvers[2];
  MassMapPLC_BasisReady = 1;
}

/**
 * mass_maps_entry_inside_aperture
 * -------------------------------
 * Purpose
 *   Geometric angle check for whether an absolute position lies inside the PLC
 *   spherical-cap aperture (diagnostics only).
 *
 * Current selection semantics
 *   - Mass-map accumulation now uses pixel-based selection: entries are kept
 *     iff their HEALPix pixel is within the cap prefix [0..Ncap-1]. This function
 *     is retained for diagnostics and optional instrumentation.
 *
 * Parameters
 *   - pos  Absolute position in true Mpc (same units as InterPartDist-scaled
 *           positions elsewhere in mass-maps).
 *
 * Returns
 *   1 if the geometric angle test passes (angle <= PLCAperture); otherwise 0.
 *   Returns 0 for the degenerate case where |pos − center| == 0.
 *
 * Notes
 *   - Uses the cached center/basis via mass_maps_cache_plc_basis_and_center().
 *   - Clamps numerical round-off on cosines to [-1,1] before acos.
 */
static inline int mass_maps_entry_inside_aperture(const double pos[3])
{
  mass_maps_cache_plc_basis_and_center();
  double vx = pos[0] - MassMapPLC_CenterPhys[0];
  double vy = pos[1] - MassMapPLC_CenterPhys[1];
  double vz = pos[2] - MassMapPLC_CenterPhys[2];
  double vnorm = sqrt(vx * vx + vy * vy + vz * vz);
  if (vnorm == 0.0)
    return 0;
  double invv = 1.0 / vnorm;
  double vhat_z = (vx * MassMapPLC_e3[0] + vy * MassMapPLC_e3[1] + vz * MassMapPLC_e3[2]) * invv;
  if (vhat_z > 1.0)
    vhat_z = 1.0;
  else if (vhat_z < -1.0)
    vhat_z = -1.0;
  double angle = acos(vhat_z) * (180.0 / acos(-1.0));
  double A = params.PLCAperture;
  return angle <= A;
}

/**
 * mass_maps_angle_deg_from_pos
 * ----------------------------
 * Purpose
 *   Compute the angular separation θ(pos) in degrees between (pos − PLC_center)
 *   and the PLC axis (plc.zvers), using the cached PLC basis and center.
 *   Intended for diagnostics (e.g., cap prefix miss tracking).
 *
 * Parameters
 *   - pos  Absolute position in true Mpc (same units as InterPartDist-scaled positions).
 *
 * Returns
 *   Angle in degrees in [0,180]. Returns 0 for |pos − center| == 0.
 *
 * Notes
 *   - Clamps cosine to [-1,1] before acos to stabilize numerics.
 *   - Thread-safe: relies on idempotent cache init.
 */
static inline double mass_maps_angle_deg_from_pos(const double pos[3])
{
  mass_maps_cache_plc_basis_and_center();
  double vx = pos[0] - MassMapPLC_CenterPhys[0];
  double vy = pos[1] - MassMapPLC_CenterPhys[1];
  double vz = pos[2] - MassMapPLC_CenterPhys[2];
  double vnorm = sqrt(vx * vx + vy * vy + vz * vz);
  if (vnorm == 0.0)
    return 0.0;
  double invv = 1.0 / vnorm;
  double vhat_z = (vx * MassMapPLC_e3[0] + vy * MassMapPLC_e3[1] + vz * MassMapPLC_e3[2]) * invv;
  if (vhat_z > 1.0)
    vhat_z = 1.0;
  else if (vhat_z < -1.0)
    vhat_z = -1.0;
  return acos(vhat_z) * (180.0 / acos(-1.0));
}

/**
 * mass_maps_init_sheets
 * ---------------------
 * Build the array of MassSheet structures from the global output redshift list.
 * Validates:
 *   - At least 2 output redshifts.
 *   - Strictly descending ordering (z[i] > z[i+1]).
 *   - First / last redshifts match PLC StartingzForPLC / LastzForPLC.
 *   - No adjacent duplicates within MASS_MAPS_Z_EPS.
 * Allocates:
 *   - MassSheets (NMassSheets = outputs.n - 1)
 *   - Boundary arrays (Z, Chi, DA) of length outputs.n.
 * Precomputes per-sheet: chi / da bounds, deltas, inverse delta, chi^3 difference.
 * Returns 0 on success, 1 on validation / allocation failure (after emitting message on task 0).
 */
int mass_maps_init_sheets(void)
{
  /* Preconditions: outputs.z[] filled, descending order, outputs.n >= 2 */
  if (outputs.n < 2)
  {
    if (!ThisTask)
      fprintf(stderr, "MASS_MAPS ERROR: Need at least 2 output redshifts (found %d) to define mass sheets.\n", outputs.n);
    return 1;
  }

  /* Ensure coverage matches PLC redshift range: first output == StartingzForPLC, last output == LastzForPLC */
  {
    double z_first = outputs.z[0];
    double z_last = outputs.z[outputs.n - 1];
    if (fabs(z_first - params.StartingzForPLC) > MASS_MAPS_Z_EPS)
    {
      if (!ThisTask)
        fprintf(stderr, "MASS_MAPS ERROR: First output redshift %g does not match StartingzForPLC=%g (|Δ|=%g > %g).\n",
                z_first, params.StartingzForPLC, fabs(z_first - params.StartingzForPLC), MASS_MAPS_Z_EPS);
      return 1;
    }
    if (fabs(z_last - params.LastzForPLC) > MASS_MAPS_Z_EPS)
    {
      if (!ThisTask)
        fprintf(stderr, "MASS_MAPS ERROR: Last output redshift %g does not match LastzForPLC=%g (|Δ|=%g > %g).\n",
                z_last, params.LastzForPLC, fabs(z_last - params.LastzForPLC), MASS_MAPS_Z_EPS);
      return 1;
    }
  }

  /* Verify strictly descending and no duplicates within eps */
  for (int i = 0; i < outputs.n - 1; ++i)
  {
    double z_hi = outputs.z[i];
    double z_lo = outputs.z[i + 1];
    if (!(z_hi > z_lo))
    {
      if (!ThisTask)
        fprintf(stderr, "MASS_MAPS ERROR: Output redshifts must be in strictly descending order (z[%d]=%g, z[%d]=%g).\n", i, z_hi, i + 1, z_lo);
      return 1;
    }
    if (fabs(z_hi - z_lo) < MASS_MAPS_Z_EPS)
    {
      if (!ThisTask)
        fprintf(stderr, "MASS_MAPS ERROR: Adjacent output redshifts too close (|%g-%g| < %g).\n", z_hi, z_lo, MASS_MAPS_Z_EPS);
      return 1;
    }
  }
  /*
   * Number of mass sheets equals the number of adjacent redshift pairs
   * in the output list: sheets s = 0..(n-2) are bounded by z[s] and z[s+1].
   * We therefore allocate:
   *  - MassSheets: one entry per sheet (outputs.n - 1)
   *  - Boundary arrays (Z/Chi/DA): one value per boundary, length outputs.n
   *    (there are outputs.n boundaries for outputs.n-1 sheets).
   */
  NMassSheets = outputs.n - 1;
  MassSheets = NULL;
  MassMapBoundaryZ = NULL;
  MassMapBoundaryChi = NULL;
  MassMapBoundaryDA = NULL;

  MassSheets = (MassSheet *)malloc(sizeof(MassSheet) * NMassSheets);
  MassMapBoundaryZ = (double *)malloc(sizeof(double) * outputs.n);
  MassMapBoundaryChi = (double *)malloc(sizeof(double) * outputs.n);
  MassMapBoundaryDA = (double *)malloc(sizeof(double) * outputs.n);
  if (!MassSheets || !MassMapBoundaryZ || !MassMapBoundaryChi || !MassMapBoundaryDA)
  {
    if (!ThisTask)
      fprintf(stderr, "MASS_MAPS ERROR: Cannot allocate MassSheets (%d entries) or boundary arrays (n=%d).\n", NMassSheets, outputs.n);
    goto fail;
  }

  /* Fill boundary arrays */
  for (int k = 0; k < outputs.n; ++k)
  {
    double z = outputs.z[k];
    MassMapBoundaryZ[k] = z;
    MassMapBoundaryChi[k] = ComovingDistance(z);
    MassMapBoundaryDA[k] = DiameterDistance(z);
  }

  /*
   * For each mass sheet s (bounded by output boundaries s and s+1),
   * cache the upper/lower redshift (z), comoving distance (chi), and
   * angular-diameter distance (da). Compute the spans Δz and Δχ, check
   * that Δχ>0 (monotonic chi with descending z), store its inverse, and
   * precompute χ_hi^3 − χ_lo^3 (useful for volume-weighted quantities).
   */
  for (int s = 0; s < NMassSheets; ++s)
  {
    MassSheet *ms = &MassSheets[s];
    ms->z_hi = MassMapBoundaryZ[s];
    ms->z_lo = MassMapBoundaryZ[s + 1];
    ms->delta_z = ms->z_hi - ms->z_lo;
    ms->chi_hi = MassMapBoundaryChi[s];
    ms->chi_lo = MassMapBoundaryChi[s + 1];
    ms->delta_chi = ms->chi_hi - ms->chi_lo;
    if (!(ms->delta_chi > 0.0))
    {
      if (!ThisTask)
        fprintf(stderr, "MASS_MAPS ERROR: Non-positive comoving distance span for sheet %d (chi_hi=%g chi_lo=%g).\n", s, ms->chi_hi, ms->chi_lo);
      goto fail;
    }
    ms->inv_dchi = 1.0 / ms->delta_chi;
    ms->da_hi = MassMapBoundaryDA[s];
    ms->da_lo = MassMapBoundaryDA[s + 1];
    double chi_hi3 = ms->chi_hi * ms->chi_hi * ms->chi_hi;
    double chi_lo3 = ms->chi_lo * ms->chi_lo * ms->chi_lo;
    ms->chi3_diff = chi_hi3 - chi_lo3;
  }

  if (!ThisTask && internal.verbose_level >= VDIAG)
  {
    printf("[%s] MASS_MAPS: initialized %d mass sheets\n", fdate(), NMassSheets);
    for (int s = 0; s < NMassSheets; ++s)
    {
      MassSheet *ms = &MassSheets[s];

      /* After traversal, print cap prefix diagnostics once per segment (root only) */
#if MASS_MAPS_CAP_DIAG
      if (!ThisTask && MassMapCapDiagActive && internal.verbose_level >= VDIAG && params.PLCAperture < 180.0 && MassMapNSIDE_current > 0)
      {
        double A = params.PLCAperture;
        printf("[%s] MASS_MAPS DIAG: aperture=%.2f deg accept_max_angle=%.3f deg miss_max_angle=%.3f deg prefix_miss=%llu\n",
               fdate(), A, MassMapCapAcceptAngleMax, MassMapCapMissAngleMax, MassMapCapPrefixMiss);
      }
#endif
      printf("  sheet %3d: z_hi=%8.4f z_lo=%8.4f Δz=%8.4f chi_hi=%10.4f chi_lo=%10.4f Δchi=%10.4f\n", s, ms->z_hi, ms->z_lo, ms->delta_z, ms->chi_hi, ms->chi_lo, ms->delta_chi);
    }
  }

  return 0;

fail:
  mass_maps_free_sheets();
  return 1;
}

/**
 * mass_maps_free_sheets
 * ---------------------
 * Free MassSheets and the boundary arrays; safe to call multiple times.
 */
void mass_maps_free_sheets(void)
{
  if (MassSheets)
    free(MassSheets);
  MassSheets = NULL;
  NMassSheets = 0;
  if (MassMapBoundaryZ)
    free(MassMapBoundaryZ);
  if (MassMapBoundaryChi)
    free(MassMapBoundaryChi);
  if (MassMapBoundaryDA)
    free(MassMapBoundaryDA);
  MassMapBoundaryZ = MassMapBoundaryChi = MassMapBoundaryDA = NULL;
}

/* ---------------------------------------------------------- */
/* HEALPix map allocation (one map per mass sheet)            */
/* ---------------------------------------------------------- */

/**
 * mass_maps_npix_from_nside
 * -------------------------
 * Purpose
 *   Return the total number of HEALPix pixels for a given NSIDE without
 *   requiring any external HEALPix library calls.
 *
 * Parameters
 *   - nside  HEALPix NSIDE (power-of-two, >0 expected).
 *
 * Returns
 *   12 * nside^2 as a long integer. For non-positive nside, the arithmetic
 *   still evaluates consistently to 0 or negative; callers are expected to
 *   validate nside beforehand.
 */
static inline long mass_maps_npix_from_nside(int nside)
{
  return 12L * (long)nside * (long)nside;
}

/**
 * mass_maps_ncap_from_aperture
 * ----------------------------
 * Purpose
 *   Compute the number of HEALPix RING pixels contained in a spherical cap of
 *   angular radius 'aperture_deg' centered at theta=0 (i.e., around the north pole
 *   aligned with the PLC axis). The result corresponds to the contiguous prefix
 *   [0 .. Ncap-1] in RING ordering for a cap centered at the pole.
 *
 * Parameters
 *   - nside         HEALPix NSIDE (power-of-two, >0).
 *   - aperture_deg  Spherical cap radius in degrees; clamped to [0, 180].
 *
 * Returns
 *   The number of pixels Ncap in the cap (0 <= Ncap <= 12*Nside^2). For aperture=0,
 *   returns 0. For aperture>=180, returns full-sky count (12*Nside^2).
 *
 * Notes
 *   - Uses analytic HEALPix ring geometry without calling chealpix:
 *       • North polar rings i=1..N-1: z_i = 1 - i^2/(3N^2), npix_i = 4i
 *       • Equatorial rings i=N..3N:   z_i = (4N - 2i)/(3N), npix_i = 4N
 *       • South polar rings i=3N+1..4N-1 with m=4N-i: z_i = -1 + m^2/(3N^2), npix_i = 4m
 *     Rings with z >= z0, where z0=cos(theta_cap), are fully included.
 *   - O(1) cost; avoids per-pixel work and chealpix dependency.
 *   - The cap is symmetric about the pole, so the included pixels form a RING prefix.
 */
static long mass_maps_ncap_from_aperture(int nside, double aperture_deg)
{
  if (nside <= 0)
    return 0;
  long npix_full = mass_maps_npix_from_nside(nside);
  /* Clamp aperture to [0,180] */
  if (aperture_deg <= 0.0)
    return 0;
  if (aperture_deg >= 180.0)
    return npix_full;
  double pi = acos(-1.0);
  double theta_cap = aperture_deg * (pi / 180.0);
  double z0 = cos(theta_cap); /* include rings with z >= z0 */

  long ncap = 0;
  int N = nside;

  /* North polar rings: i = 1..N-1, z_i = 1 - i^2/(3 N^2), npix_i = 4 i */
  double tmp = 3.0 * (double)N * (double)N * (1.0 - z0);
  long ipolar_max = (long)floor(sqrt(fmax(0.0, tmp)) + 1e-12);
  if (ipolar_max > (long)N - 1)
    ipolar_max = (long)N - 1;
  if (ipolar_max > 0)
  {
    /* sum_{i=1}^{ipolar_max} 4 i = 2 * ipolar_max * (ipolar_max + 1) */
    ncap += 2L * ipolar_max * (ipolar_max + 1L);
  }

  /* Equatorial rings: i = N..3N, z_i = (4N - 2i)/(3N), npix_i = 4N */
  double i2max_d = ((4.0 * N) - 3.0 * N * z0) / 2.0;
  long i2_max = (long)floor(i2max_d + 1e-12);
  if (i2_max > 3L * N)
    i2_max = 3L * N;
  if (i2_max >= (long)N)
  {
    long neq_rings = i2_max - (long)N + 1L;
    if (neq_rings > 0)
      ncap += (long)(4L * (long)N) * neq_rings;
  }

  /* South polar rings: i = 3N+1..4N-1, z_i = -1 + m^2/(3N^2), where m = 4N - i, npix_i = 4 m */
  if (i2_max >= 3L * N)
  {
    /* Only if cap extends beyond equatorial band */
    double tmp_s = 3.0 * (double)N * (double)N * (1.0 + z0);
    long m_min = (long)ceil(sqrt(fmax(0.0, tmp_s)) - 1e-12);
    if (m_min < 1)
      m_min = 1;
    if (m_min <= (long)N - 1)
    {
      long m_max = (long)N - 1;
      long count = m_max - m_min + 1L;
      /* sum_{m=m_min}^{m_max} 4 m = 4 * ( (m_min + m_max) * count / 2 ) */
      long sum_m = (m_min + m_max) * count;
      ncap += 2L * sum_m;
    }
  }
  if (ncap > npix_full)
    ncap = npix_full;
  /* For very large apertures (> 90 deg), the subset defined by z>=z0 spans
     all north polar + equatorial + part of south polar rings, still forming
     a RING prefix. However, callers using this as a cap representation should
     treat the map as a large cap (not its complement). Diagnostics will help
     verify boundary coverage. */
  return ncap;
}

/**
 * mass_maps_segment_ptr
 * ---------------------
 * Purpose
 *   Accessor for the per-segment HEALPix map storage.
 *
 * Parameters
 *   - s  Segment index (0 .. NMassSheets-1).
 *
 * Returns
 *   Pointer to the beginning of the map buffer for segment s (length
 *   MassMapNPIX), or NULL if maps are not allocated or s is out of range.
 *
 * Notes
 *   - The memory layout is a flat array of size NMassSheets * MassMapNPIX,
 *     with segment s at offset s * MassMapNPIX.
 *   - Caller must not write past MassMapNPIX elements.
 */
static inline double *mass_maps_segment_ptr(int s)
{
  if (!MassMapSegmentMaps || s < 0 || s >= NMassSheets)
    return NULL;
  return MassMapSegmentMaps + ((long)s) * MassMapNPIX;
}

/* Forward decl: returns pole pixel (theta=0) in RING ordering */
static inline long mass_maps_axis_pixel_ring(int nside);

/**
 * mass_maps_zero_all_maps
 * -----------------------
 * Purpose
 *   Zero-initialize all per-segment HEALPix map buffers.
 *
 * Behavior
 *   - No-op if maps are not allocated (MassMapSegmentMaps == NULL).
 *   - Iterates over the flat buffer of size NMassSheets * MassMapNPIX and
 *     sets all entries to 0.0.
 *
 * Notes
 *   - Single-threaded loop; safe to call from any context before concurrent
 *     accumulation begins. If parallelization is desired, this loop can be
 *     guarded with OpenMP pragmas as needed.
 */
static void mass_maps_zero_all_maps(void)
{
  if (!MassMapSegmentMaps)
    return;
  long total = (long)NMassSheets * MassMapNPIX;
  for (long i = 0; i < total; ++i)
    MassMapSegmentMaps[i] = 0.0;
}

/**
 * mass_maps_free_maps
 * -------------------
 * Purpose
 *   Release HEALPix map memory owned by this module and reset associated
 *   state so maps are considered uninitialized.
 *
 * Behavior
 *   - Frees MassMapSegmentMaps (per-segment map buffer) if allocated and
 *     sets the pointer to NULL.
 *   - Frees MassMapReduceBuffer (temporary MPI reduction buffer) if allocated
 *     and sets the pointer to NULL.
 *   - Resets MassMapNSIDE_current and MassMapNPIX to 0 to mark maps as absent.
 *
 * Notes
 *   - Idempotent: safe to call multiple times.
 *   - Does not modify MassMapNCAP; callers that change NSIDE/aperture should
 *     recompute NCAP during re-initialization (via mass_maps_init_healpix_maps).
 */
void mass_maps_free_maps(void)
{
  if (MassMapSegmentMaps)
    free(MassMapSegmentMaps);
  MassMapSegmentMaps = NULL;
  if (MassMapReduceBuffer)
    free(MassMapReduceBuffer);
  MassMapReduceBuffer = NULL;
  MassMapNSIDE_current = 0;
  MassMapNPIX = 0;
}

/**
 * mass_maps_init_healpix_maps
 * ---------------------------
 * Purpose
 *   Ensure per-segment HEALPix maps are allocated and initialized for the
 *   requested NSIDE, honoring the configured PLC aperture and output mode.
 *
 * Parameters
 *   - nside  Desired HEALPix NSIDE (power-of-two, >0).
 *
 * Behavior
 *   - Frees any existing maps, then computes:
 *       npix_full = 12 * nside^2
 *       Ncap      = mass_maps_ncap_from_aperture(nside, params.PLCAperture)
 *     and selects MassMapNPIX based on MASS_MAPS_FULLSKY_OUTPUT:
 *       • FULLSKY: NPIX = npix_full (entire sphere)
 *       • CAP:     NPIX = Ncap      (RING prefix 0..Ncap-1 around north pole)
 *   - Allocates a flat buffer of doubles of size NMassSheets * NPIX and zeros it.
 *   - On rank 0 and at diagnostic verbosity, prints a summary (ORDERING=RING).
 *
 * Returns
 *   0 on success; 1 if inputs are invalid, memory allocation fails, or the
 *   computed NPIX is zero (maps disabled). In failure cases, internal state
 *   (MassMapNSIDE_current, MassMapNPIX, MassMapNCAP) is reset to 0.
 *
 * Notes
 *   - The map memory is owned by this module and must be released with
 *     mass_maps_free_maps().
 *   - Aperture is in degrees; cap computation is analytic in RING geometry.
 */
int mass_maps_init_healpix_maps(int nside)
{
  if (nside <= 0)
    return 1;
  if (NMassSheets <= 0)
    return 1;
  if (MassMapSegmentMaps && MassMapNSIDE_current == nside)
    return 0; /* already ready */

  /* (Re)allocate */
  mass_maps_free_maps();
  MassMapNSIDE_current = nside;
  long npix_full = mass_maps_npix_from_nside(nside);
  MassMapNCAP = mass_maps_ncap_from_aperture(nside, params.PLCAperture);
#if MASS_MAPS_FULLSKY_OUTPUT
  MassMapNPIX = npix_full;
#else
  MassMapNPIX = MassMapNCAP;
#endif
  if (MassMapNPIX <= 0)
  {
    /* Ensure we still have a valid array; keep zero-sized map disabled */
    if (!ThisTask)
      fprintf(stderr, "MASS_MAPS WARNING: computed Ncap=0 for NSIDE=%d, aperture=%.6g deg. Maps disabled.\n", nside, params.PLCAperture);
    MassMapNSIDE_current = 0;
    MassMapNPIX = 0;
    MassMapNCAP = 0;
    return 1;
  }
  long total = (long)NMassSheets * MassMapNPIX;
  MassMapSegmentMaps = (double *)malloc(sizeof(double) * (size_t)total);
  if (!MassMapSegmentMaps)
  {
    if (!ThisTask)
      fprintf(stderr, "MASS_MAPS ERROR: cannot allocate HEALPix maps (segments=%d npix=%ld total=%ld).\n", NMassSheets, MassMapNPIX, total);
    MassMapNSIDE_current = 0;
    MassMapNPIX = 0;
    MassMapNCAP = 0;
    return 1;
  }
  mass_maps_zero_all_maps();

  if (!ThisTask && internal.verbose_level >= VDIAG)
  {
    long axis_pix = mass_maps_axis_pixel_ring(nside);
    printf("[%s] MASS_MAPS: allocated HEALPix maps ORDERING=RING NSIDE=%d fullNPIX=%ld aperture=%.6gdeg Ncap=%ld using_%s npix=%ld segments=%d axis_pix=%ld\n",
           fdate(), nside, npix_full, params.PLCAperture, MassMapNCAP,
#if MASS_MAPS_FULLSKY_OUTPUT
           "FULLSKY",
#else
           "CAP",
#endif
           MassMapNPIX, NMassSheets, axis_pix);
  }
  return 0;
}

/* ---------------------------------------------------------- */
/* PLC-oriented basis and pixelization helpers                */
/* ---------------------------------------------------------- */

/**
 * mass_maps_map_stats
 * -------------------
 * Purpose
 *   Compute simple statistics over a map buffer: sum, minimum, maximum, and
 *   the count of non-zero entries (nnz).
 *
 * Parameters
 *   - map      Pointer to the map array (length npix). Must be non-NULL if
 *              npix > 0. Not dereferenced when npix <= 0.
 *   - npix     Number of pixels/elements in the map (may be 0).
 *   - sum_out  If non-NULL, receives the sum of all elements (0.0 if npix==0).
 *   - min_out  If non-NULL, receives the minimum value (0.0 if npix==0).
 *   - max_out  If non-NULL, receives the maximum value (0.0 if npix==0).
 *   - nnz_out  If non-NULL, receives the count of elements with value != 0.0.
 *
 * Notes
 *   - Non-zero detection uses exact floating-point comparison (v != 0.0), with
 *     no epsilon tolerance. Callers needing a threshold should apply it
 *     externally.
 *   - When npix <= 0, outputs (if requested) are set to their neutral defaults:
 *     sum=0, min=0, max=0, nnz=0.
 */
static inline void mass_maps_map_stats(const double *map, long npix,
                                       double *sum_out, double *min_out, double *max_out, long *nnz_out)
{
  double s = 0.0;
  double mn = 0.0;
  double mx = 0.0;
  long nnz = 0;
  if (npix > 0)
  {
    mn = mx = map[0];
    for (long i = 0; i < npix; ++i)
    {
      double v = map[i];
      s += v;
      if (v < mn)
        mn = v;
      if (v > mx)
        mx = v;
      if (v != 0.0)
        nnz++;
    }
  }
  if (sum_out)
    *sum_out = s;
  if (min_out)
    *min_out = mn;
  if (max_out)
    *max_out = mx;
  if (nnz_out)
    *nnz_out = nnz;
}

/**
 * mass_maps_axis_pixel_ring
 * -------------------------
 * Purpose
 *   Return the HEALPix pixel index (RING ordering) at polar angle theta=0 and
 *   phi=0 for a given NSIDE. In this module, theta=0 corresponds to the PLC
 *   axis direction used to orient maps.
 *
 * Parameters
 *   - nside  HEALPix NSIDE (power-of-two, >0 expected).
 *
 * Returns
 *   Pixel index in [0, 12*nside^2 - 1] on success; -1 if the mapping fails or
 *   nside is invalid.
 *
 * Notes
 *   - Uses RING ordering (ang2pix_ring). Not valid for NESTED ordering.
 *   - This is a convenience helper for diagnostics and header metadata.
 */
static inline long mass_maps_axis_pixel_ring(int nside)
{
  long pix = -1;
  ang2pix_ring((long)nside, 0.0, 0.0, &pix); /* theta=0 on PLC axis */
  return pix;
}

/**
 * mass_maps_crossing_from_shifted_positions
 * ----------------------------------------
 * Purpose
 *   Detect and locate a radial crossing of the PLC spherical boundary along
 *   a particle segment, using absolute positions that already include the
 *   replication shift, and precomputed chi values at the segment endpoints.
 *   The boundary is defined by r(alpha) = chi(alpha), with
 *     H = chi - r, and linear interpolation in alpha between endpoints.
 *
 * Parameters
 *   - Pprev      Absolute position at the previous endpoint (true Mpc),
 *                replication shift already applied.
 *   - Pcurr      Absolute position at the current endpoint (true Mpc),
 *                replication shift already applied.
 *   - c_phys     PLC center in true Mpc (same frame/units as Pprev/Pcurr).
 *   - chi_prev   Comoving distance chi at the previous endpoint (Mpc). In the
 *                present-epoch frame used here, this equals the corresponding
 *                physical distance, so it is directly comparable to r.
 *   - chi_curr   Comoving distance chi at the current endpoint (Mpc).
 *   - alpha_out  Optional; on success, receives the segment interpolation
 *                parameter alpha in [0,1] locating the crossing.
 *   - entry_pos  Optional; on success, receives the 3D crossing position
 *                Pprev + alpha * (Pcurr - Pprev).
 *   - chi_cross_out Optional; on success, receives the interpolated comoving
 *                distance at the crossing: chi_prev + alpha * (chi_curr - chi_prev).
 *                Preferred over |entry_pos - c| for deriving the crossing redshift.
 *
 * Returns
 *   1 if a crossing is detected and the intersection is computed; otherwise 0.
 *
 * Behavior
 *   - Let H_prev = chi_prev - |Pprev - c| and H_curr = chi_curr - |Pcurr - c|.
 *     A crossing is reported only for the inside→outside case
 *       (H_prev > 0 && H_curr <= 0).
 *   - Uses linear interpolation of H in alpha and solves H(alpha)=0 via
 *       alpha0 = H_prev / (H_prev - H_curr) as a secant predictor, then
 *       applies one Newton corrector step on F(alpha) = chi(alpha) - |P(alpha) - c|
 *       to refine the crossing fraction. The corrector is skipped when
 *       |P(alpha0)-c| == 0 or the derivative is too small.
 *       A small-denominator guard using MASS_MAPS_F_EPS protects the predictor.
 *       The resulting alpha is clamped to [0,1] after both steps.
 *
 * Notes
 *   - Purely local computation; no shared state is modified.
 *   - Caller must ensure units are consistent: P*, c_phys in true Mpc; chi_* in
 *     Mpc (comoving), which here are comparable to present-epoch physical Mpc.
 */
static inline int mass_maps_crossing_from_shifted_positions(const double Pprev[3],
                                                            const double Pcurr[3],
                                                            const double c_phys[3],
                                                            double chi_prev,
                                                            double chi_curr,
                                                            double *alpha_out,
                                                            double entry_pos[3],
                                                            double *chi_cross_out)
{
  /* Distances to PLC center */
  double dx0 = Pprev[0] - c_phys[0];
  double dy0 = Pprev[1] - c_phys[1];
  double dz0 = Pprev[2] - c_phys[2];
  double dx1 = Pcurr[0] - c_phys[0];
  double dy1 = Pcurr[1] - c_phys[1];
  double dz1 = Pcurr[2] - c_phys[2];
  double r_prev = sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0);
  double r_curr = sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
  double H_prev = chi_prev - r_prev;
  double H_curr = chi_curr - r_curr;
  if (!(H_prev > 0.0 && H_curr <= 0.0))
    return 0;
  double denom = H_prev - H_curr;
  if (fabs(denom) < MASS_MAPS_F_EPS)
    return 0;
  /* Secant predictor */
  double alpha = H_prev / denom;
  if (alpha < 0.0)
    alpha = 0.0;
  else if (alpha > 1.0)
    alpha = 1.0;

  /* One Newton corrector step on F(alpha) = chi(alpha) - |P(alpha) - c_phys| */
  {
    const double newton_eps = 1e-12;
    double dP[3] = {Pcurr[0] - Pprev[0], Pcurr[1] - Pprev[1], Pcurr[2] - Pprev[2]};
    double dchi = chi_curr - chi_prev;
    /* Evaluate at alpha0 */
    double P0[3] = {Pprev[0] + alpha * dP[0] - c_phys[0],
                    Pprev[1] + alpha * dP[1] - c_phys[1],
                    Pprev[2] + alpha * dP[2] - c_phys[2]};
    double r0 = sqrt(P0[0] * P0[0] + P0[1] * P0[1] + P0[2] * P0[2]);
    if (r0 > 0.0)
    {
      double F0 = chi_prev + alpha * dchi - r0;
      double Fp0 = dchi - (P0[0] * dP[0] + P0[1] * dP[1] + P0[2] * dP[2]) / r0;
      if (fabs(Fp0) > newton_eps)
      {
        double alpha1 = alpha - F0 / Fp0;
        if (alpha1 < 0.0)
          alpha1 = 0.0;
        else if (alpha1 > 1.0)
          alpha1 = 1.0;
        alpha = alpha1;
      }
    }
  }

  if (alpha_out)
    *alpha_out = alpha;
  if (entry_pos)
  {
    entry_pos[0] = Pprev[0] + alpha * (Pcurr[0] - Pprev[0]);
    entry_pos[1] = Pprev[1] + alpha * (Pcurr[1] - Pprev[1]);
    entry_pos[2] = Pprev[2] + alpha * (Pcurr[2] - Pprev[2]);
  }
  if (chi_cross_out)
    *chi_cross_out = chi_prev + alpha * (chi_curr - chi_prev);
  return 1;
}

/* ---------------------------------------------------------- */
/* Sub-volume culling helpers (AABB vs PLC shell)             */
/* ---------------------------------------------------------- */

/**
 * aabb_distance_bounds_to_point
 * -----------------------------
 * Purpose
 *   Compute tight lower and upper bounds for the Euclidean distance from a
 *   point c to an axis-aligned bounding box (AABB) defined by per-axis
 *   minima and maxima.
 *
 * Parameters
 *   - minv      3-vector of AABB minimum coordinates (x_min, y_min, z_min).
 *   - maxv      3-vector of AABB maximum coordinates (x_max, y_max, z_max).
 *   - c         3-vector of point coordinates.
 *   - dmin_out  If non-NULL, receives the minimum possible distance from c to
 *               any point in the AABB (0 if c lies inside the box).
 *   - dmax_out  If non-NULL, receives the maximum possible distance from c to
 *               any point in the AABB (distance to the farthest corner).
 *
 * Behavior
 *   - d_min is computed as the length of the per-axis clamp residuals (zero
 *     when c is inside the interval [minv,maxv] on that axis).
 *   - d_max is computed as the distance to the farthest AABB corner from c,
 *     found by taking the larger absolute delta on each axis.
 *
 * Invariants
 *   - For any inputs, d_min <= d_max, and both are >= 0.
 *   - Units follow the inputs (e.g., true Mpc when min/max are given in Mpc).
 *
 * Notes
 *   - Used for quick culling against the PLC spherical shell: compare these
 *     bounds with the segment's chi interval to early-accept or early-reject
 *     replications before finer tests.
 */
static inline void aabb_distance_bounds_to_point(const double minv[3], const double maxv[3], const double c[3], double *dmin_out, double *dmax_out)
{
  /* d_min: point-to-box distance */
  double d2min = 0.0;
  for (int a = 0; a < 3; ++a)
  {
    double da = 0.0;
    if (c[a] < minv[a])
      da = minv[a] - c[a];
    else if (c[a] > maxv[a])
      da = c[a] - maxv[a];
    d2min += da * da;
  }
  /* d_max: farthest corner distance */
  double d2max = 0.0;
  for (int a = 0; a < 3; ++a)
  {
    double da0 = fabs(minv[a] - c[a]);
    double da1 = fabs(maxv[a] - c[a]);
    double da = (da0 > da1) ? da0 : da1;
    d2max += da * da;
  }
  if (dmin_out)
    *dmin_out = sqrt(d2min);
  if (dmax_out)
    *dmax_out = sqrt(d2max);
}

/**
 * mass_maps_compute_pixel_from_pos
 * --------------------------------
 * Purpose
 *   Compute the HEALPix RING pixel index corresponding to an absolute
 *   3D position expressed in physical units (true Mpc), using the PLC
 *   orientation (e1,e2,e3) and center previously cached.
 *
 * Parameters
 *   - pos      Absolute position (true Mpc) in the same frame as the cached
 *              PLC center.
 *   - nside    HEALPix NSIDE (power-of-two, >0).
 *   - pix_out  Output pointer; on success receives the pixel index (>=0).
 *
 * Returns
 *   1 on success; 0 on invalid inputs (pix_out==NULL or nside<=0), when the
 *   position coincides with the PLC center (undefined direction), or if the
 *   RING mapping fails.
 *
 * Behavior
 *   - Translates pos by the cached PLC center and normalizes to obtain vhat.
 *   - Computes spherical angles relative to the PLC basis:
 *       cos(theta) = vhat·e3,   phi = atan2(vhat·e2, vhat·e1) in [0,2π).
 *   - Calls ang2pix_ring(NSIDE, theta, phi) and stores the result in pix_out.
 *
 * Notes
 *   - Uses RING ordering only. For partial-sky outputs, callers must also
 *     check that the returned pixel lies within the active prefix (0..Ncap-1).
 *   - Invokes mass_maps_cache_plc_basis_and_center() to ensure the PLC center
 *     and basis are available.
 */
int mass_maps_compute_pixel_from_pos(const double pos[3], int nside, long *pix_out)
{
  if (!pix_out || nside <= 0)
    return 0;
  /* Cache center/basis once per process to avoid recomputing in hot paths */
  mass_maps_cache_plc_basis_and_center();
  double vx = pos[0] - MassMapPLC_CenterPhys[0];
  double vy = pos[1] - MassMapPLC_CenterPhys[1];
  double vz = pos[2] - MassMapPLC_CenterPhys[2];
  double vnorm = sqrt(vx * vx + vy * vy + vz * vz);
  if (vnorm == 0.0)
    return 0;
  double invv = 1.0 / vnorm;
  double vhat[3] = {vx * invv, vy * invv, vz * invv};
  const double *e1 = MassMapPLC_e1;
  const double *e2 = MassMapPLC_e2;
  const double *e3 = MassMapPLC_e3;
  double cos_theta = vhat[0] * e3[0] + vhat[1] * e3[1] + vhat[2] * e3[2];
  if (cos_theta > 1.0)
    cos_theta = 1.0;
  else if (cos_theta < -1.0)
    cos_theta = -1.0;
  double theta = acos(cos_theta);
  double x1 = vhat[0] * e1[0] + vhat[1] * e1[1] + vhat[2] * e1[2];
  double y1 = vhat[0] * e2[0] + vhat[1] * e2[1] + vhat[2] * e2[2];
  double phi = atan2(y1, x1);
  if (phi < 0.0)
  {
    double twopi = 2.0 * acos(-1.0);
    phi += twopi;
  }
  long ipix = -1;
  /* RING-only ordering */
  ang2pix_ring((long)nside, theta, phi, &ipix);
  if (ipix < 0)
    return 0;
  *pix_out = ipix;
  return 1;
}

/**
 * mass_maps_accumulate_entry_pos
 * ------------------------------
 * Purpose
 *   Accumulate a scalar weight into the HEALPix map for a given mass sheet
 *   segment at the pixel corresponding to an absolute entry position.
 *
 * Parameters
 *   - s          Segment index (0 .. NMassSheets-1).
 *   - entry_pos  Absolute entry position in true Mpc (same frame as PLC).
 *   - weight     Scalar to add to the corresponding pixel (e.g., 1.0).
 *
 * Behavior
 *   - Validates segment and NSIDE; computes the pixel via
 *     mass_maps_compute_pixel_from_pos(entry_pos, MassMapNSIDE_current, &ipix).
 *   - If the segment map exists and ipix lies within [0, MassMapNPIX), adds
 *     'weight' to map[ipix]. Out-of-range pixels are ignored.
 *
 * Threading
 *   - This helper does not perform any synchronization. Callers must ensure
 *     either single-threaded use, per-thread private maps (later reduced), or
 *     protect the update with an atomic/lock at the call site when sharing a
 *     global map among threads.
 */
static inline void mass_maps_accumulate_entry_pos(int s,
                                                  const double entry_pos[3],
                                                  double weight)
{
  if (s < 0 || s >= NMassSheets || MassMapNSIDE_current <= 0)
    return;
  long ipix;
  if (!mass_maps_compute_pixel_from_pos(entry_pos, MassMapNSIDE_current, &ipix))
    return;
  double *map = mass_maps_segment_ptr(s);
  if (!map)
    return;
  if ((unsigned long)ipix < (unsigned long)MassMapNPIX)
    map[ipix] += weight;
}

/**
 * mass_maps_write_segment_map
 * ---------------------------
 * Purpose
 *   Rank 0 only: write the HEALPix map for a single mass sheet segment to a
 *   FITS file named pinocchio.<RunFlag>.massmap.seg%03d.fits.
 *
 * Parameters
 *   - s  Segment index (0 .. NMassSheets-1).
 *
 * Behavior
 *   - Precondition checks: only executes on task 0, with maps allocated and
 *     MassMapNSIDE_current > 0 and MassMapNPIX > 0.
 *   - Overwrites any existing file (CFITSIO "!" prefix).
 *   - Two output modes controlled by MASS_MAPS_FULLSKY_OUTPUT:
 *       • FULLSKY: Dual-format output for maximum compatibility:
 *           - Primary IMAGE HDU (length NPIX = npix_full, RING ordering)
 *           - Secondary BINTABLE extension (EXTNAME='HEALPIX') with columns:
 *             PIXEL (0..npix_full-1, EXPLICIT indexing) and TEMPERATURE (map values)
 *         Healpy can read either the IMAGE or the BINTABLE path; providing both
 *         avoids version-dependent ingestion issues.
 *       • CAP:     BINTABLE with explicit indexing (columns: PIXEL, TEMPERATURE)
 *                  where PIXEL runs 0..Ncap-1 and TEMPERATURE holds map values.
 *   - Adds HEALPix metadata (PIXTYPE=HEALPIX, ORDERING=RING, NSIDE, FIRSTPIX,
 *     LASTPIX) and PLC context (APERTURE, AXISV1/2/3). If available, also
 *     writes filter diagnostics via mass_maps_write_filter_keywords().
 *
 * Returns
 *   0 on success; 1 on invalid input or CFITSIO error (error is reported to
 *   stderr via fits_report_error). On success, prints a diagnostic line when
 *   verbosity >= VDIAG.
 */
int mass_maps_write_segment_map(int s)
{
  if (ThisTask != 0)
    return 0;
  if (s < 0 || s >= NMassSheets)
    return 1;
  if (!MassMapSegmentMaps || MassMapNSIDE_current <= 0 || MassMapNPIX <= 0)
    return 1;

  char fname[3 * LBLENGTH];
  snprintf(fname, sizeof(fname), "pinocchio.%s.massmap.seg%03d.fits", params.RunFlag, s);
  char fitsname[3 * LBLENGTH + 2];
  snprintf(fitsname, sizeof(fitsname), "!%s", fname); /* overwrite if exists */

  fitsfile *fptr = NULL;
  int status = 0;
  long naxes[1];
  naxes[0] = MassMapNPIX;
  double *map = mass_maps_segment_ptr(s);
  if (!map)
    return 1;

  fits_create_file(&fptr, fitsname, &status);
  if (status)
  {
    fits_report_error(stderr, status);
    return 1;
  }
  /* Choose output format: full-sky IMAGE or partial-sky BINTABLE (IMPLICIT indexing) */
#if MASS_MAPS_FULLSKY_OUTPUT
  /* Write full-sky map directly as primary IMAGE HDU (simplest for many healpy versions) */
  fits_create_img(fptr, DOUBLE_IMG, 1, naxes, &status);
  if (status)
  {
    fits_report_error(stderr, status);
    fits_close_file(fptr, &status);
    return 1;
  }
  {
    char pixtype[] = "HEALPIX";
    char ordering[] = "RING";
    fits_update_key(fptr, TSTRING, "PIXTYPE", pixtype, "HEALPix pixelisation", &status);
    fits_update_key(fptr, TSTRING, "ORDERING", ordering, "Pixel ordering scheme (RING/NESTED)", &status);
    fits_update_key(fptr, TINT, "NSIDE", &MassMapNSIDE_current, "HEALPix NSIDE", &status);
    long firstpix = 0;
    long lastpix = MassMapNPIX - 1;
    fits_update_key(fptr, TLONG, "FIRSTPIX", &firstpix, "First pixel number (0-based)", &status);
    fits_update_key(fptr, TLONG, "LASTPIX", &lastpix, "Last pixel number (0-based)", &status);
    fits_update_key(fptr, TDOUBLE, "APERTURE", &params.PLCAperture, "PLC cone aperture (deg)", &status);
    {
      char seltype[32] = "PIXEL_CAP_PREFIX";
      fits_update_key(fptr, TSTRING, "SELTYPE", seltype, "Angular selection mode", &status);
    }
    fits_update_key(fptr, TDOUBLE, "AXISV1", &plc.zvers[0], "PLC axis x-component", &status);
    fits_update_key(fptr, TDOUBLE, "AXISV2", &plc.zvers[1], "PLC axis y-component", &status);
    fits_update_key(fptr, TDOUBLE, "AXISV3", &plc.zvers[2], "PLC axis z-component", &status);
    mass_maps_write_filter_keywords(fptr, s);
    mass_maps_write_cosmology_keywords(fptr);
  }
  {
    long fpixel = 1;
    long nelem = MassMapNPIX;
    fits_write_img(fptr, TDOUBLE, fpixel, nelem, map, &status);
  }
  /* Also provide a HEALPIX BINTABLE extension with explicit indexing and TEMPERATURE column */
  if (!status)
  {
    int tfields = 2;
    char *ttype[] = {"PIXEL", "TEMPERATURE"};
    char *tform[] = {"1J", "1D"};
    char *tunit[] = {"", ""};
    long nrows = MassMapNPIX;
    fits_create_tbl(fptr, BINARY_TBL, nrows, tfields, ttype, tform, tunit, "HEALPIX", &status);
    if (!status)
    {
      char pixtype[] = "HEALPIX";
      char ordering[] = "RING";
      char idxschm[] = "EXPLICIT";
      fits_update_key(fptr, TSTRING, "PIXTYPE", pixtype, "HEALPix pixelisation", &status);
      fits_update_key(fptr, TSTRING, "ORDERING", ordering, "Pixel ordering scheme (RING)", &status);
      fits_update_key(fptr, TINT, "NSIDE", &MassMapNSIDE_current, "HEALPix NSIDE", &status);
      fits_update_key(fptr, TSTRING, "INDXSCHM", idxschm, "Indexing: EXPLICIT (PIXEL column)", &status);
      long firstpix = 0;
      long lastpix = MassMapNPIX - 1;
      fits_update_key(fptr, TLONG, "FIRSTPIX", &firstpix, "First pixel number (0-based)", &status);
      fits_update_key(fptr, TLONG, "LASTPIX", &lastpix, "Last pixel number (0-based)", &status);
      fits_update_key(fptr, TDOUBLE, "APERTURE", &params.PLCAperture, "PLC cone aperture (deg)", &status);
      {
        char seltype[32] = "PIXEL_CAP_PREFIX";
        fits_update_key(fptr, TSTRING, "SELTYPE", seltype, "Angular selection mode", &status);
      }
      fits_update_key(fptr, TDOUBLE, "AXISV1", &plc.zvers[0], "PLC axis x-component", &status);
      fits_update_key(fptr, TDOUBLE, "AXISV2", &plc.zvers[1], "PLC axis y-component", &status);
      fits_update_key(fptr, TDOUBLE, "AXISV3", &plc.zvers[2], "PLC axis z-component", &status);
      mass_maps_write_filter_keywords(fptr, s);
      mass_maps_write_cosmology_keywords(fptr);
      /* Write columns */
      long firstrow = 1;  /* 1-based */
      long firstelem = 1; /* scalar */
      long nelements = MassMapNPIX;
      int *pix_index = (int *)malloc(sizeof(int) * (size_t)MassMapNPIX);
      if (!pix_index)
      {
        fits_report_error(stderr, status);
        fits_close_file(fptr, &status);
        return 1;
      }
      for (long i = 0; i < MassMapNPIX; ++i)
        pix_index[i] = (int)i;
      fits_write_col(fptr, TINT, 1, firstrow, firstelem, nelements, pix_index, &status);
      free(pix_index);
      fits_write_col(fptr, TDOUBLE, 2, firstrow, firstelem, nelements, map, &status);
    }
  }
#else
  /* Partial-sky in a binary table; choose EXPLICIT indexing for broad compatibility */
  int tfields = 2;
  char *ttype[] = {"PIXEL", "TEMPERATURE"};
  char *tform[] = {"1J", "1D"};
  char *tunit[] = {"", ""};
  long nrows = MassMapNPIX;
  fits_create_tbl(fptr, BINARY_TBL, nrows, tfields, ttype, tform, tunit, "HEALPIX", &status);
  if (status)
  {
    fits_report_error(stderr, status);
    fits_close_file(fptr, &status);
    return 1;
  }
  {
    /* HEALPix metadata */
    char pixtype[] = "HEALPIX";
    char ordering[] = "RING";
    char idxschm[] = "EXPLICIT";
    fits_update_key(fptr, TSTRING, "PIXTYPE", pixtype, "HEALPix pixelisation", &status);
    fits_update_key(fptr, TSTRING, "ORDERING", ordering, "Pixel ordering scheme (RING)", &status);
    fits_update_key(fptr, TINT, "NSIDE", &MassMapNSIDE_current, "HEALPix NSIDE", &status);
    fits_update_key(fptr, TSTRING, "INDXSCHM", idxschm, "Indexing: EXPLICIT (PIXEL column)", &status);
    /* Explicit indexing: we also include FIRSTPIX/LASTPIX for convenience */
    long firstpix = 0;
    long lastpix = MassMapNPIX - 1; /* 0..Ncap-1 */
    fits_update_key(fptr, TLONG, "FIRSTPIX", &firstpix, "First pixel number (0-based)", &status);
    fits_update_key(fptr, TLONG, "LASTPIX", &lastpix, "Last pixel number (0-based)", &status);
    /* Context: PLC geometry */
    fits_update_key(fptr, TDOUBLE, "APERTURE", &params.PLCAperture, "PLC cone aperture (deg)", &status);
    {
      char seltype[32] = "PIXEL_CAP_PREFIX";
      fits_update_key(fptr, TSTRING, "SELTYPE", seltype, "Angular selection mode", &status);
    }
    fits_update_key(fptr, TDOUBLE, "AXISV1", &plc.zvers[0], "PLC axis x-component", &status);
    fits_update_key(fptr, TDOUBLE, "AXISV2", &plc.zvers[1], "PLC axis y-component", &status);
    fits_update_key(fptr, TDOUBLE, "AXISV3", &plc.zvers[2], "PLC axis z-component", &status);
    /* Filter metadata (if available) */
    mass_maps_write_filter_keywords(fptr, s);
    mass_maps_write_cosmology_keywords(fptr);
  }
  {
    /* Write columns: PIXEL (0..Ncap-1) and TEMPERATURE */
    long firstrow = 1;  /* 1-based */
    long firstelem = 1; /* scalar */
    long nelements = MassMapNPIX;
    /* Build pixel index array */
    int *pix_index = (int *)malloc(sizeof(int) * (size_t)MassMapNPIX);
    if (!pix_index)
    {
      fits_report_error(stderr, status);
      fits_close_file(fptr, &status);
      return 1;
    }
    for (long i = 0; i < MassMapNPIX; ++i)
      pix_index[i] = (int)i;
    fits_write_col(fptr, TINT, 1, firstrow, firstelem, nelements, pix_index, &status);
    free(pix_index);
    fits_write_col(fptr, TDOUBLE, 2, firstrow, firstelem, nelements, map, &status);
  }
#endif

  fits_close_file(fptr, &status);
  if (status)
  {
    fits_report_error(stderr, status);
    return 1;
  }
  if (internal.verbose_level >= VDIAG)
    printf("[%s] MASS_MAPS: wrote %s (NSIDE=%d NPIX=%ld, format=%s)\n", fdate(), fname, MassMapNSIDE_current, MassMapNPIX,
#if MASS_MAPS_FULLSKY_OUTPUT
           "IMAGE"
#else
           "BINTABLE"
#endif
    );
  return 0;
}

/**
 * mass_maps_write_sheet_table
 * ---------------------------
 * Purpose
 *   Rank 0 only: write an ASCII table describing each mass sheet (segment)
 *   to pinocchio.<RunFlag>.sheets.out for reproducibility and external
 *   analysis.
 *
 * Format
 *   - Plain text, whitespace-separated values.
 *   - Comment lines start with '#', including a column legend.
 *   - File is overwritten if it already exists.
 *
 * Columns and units
 *   1) id             Integer segment id (0..NMassSheets-1)
 *   2) z_hi           Upper redshift bound (dimensionless)
 *   3) z_lo           Lower redshift bound (dimensionless), z_hi > z_lo
 *   4) delta_z        z_hi - z_lo
 *   5) chi_hi         Upper comoving distance (Mpc)
 *   6) chi_lo         Lower comoving distance (Mpc)
 *   7) delta_chi      chi_hi - chi_lo (Mpc)
 *   8) inv_delta_chi  1 / (chi_hi - chi_lo) (1/Mpc)
 *   9) da_hi          Angular diameter distance at z_hi (Mpc)
 *  10) da_lo          Angular diameter distance at z_lo (Mpc)
 *  11) chi3_diff      chi_hi^3 - chi_lo^3 (Mpc^3)
 *
 * Returns
 *   0 on success (including the no-op case when there is nothing to write),
 *   1 on I/O error opening or writing the file.
 *
 * Notes
 *   - Only task 0 writes this file; other ranks return 0 immediately.
 *   - No MPI synchronization is performed here. The table reflects the
 *     current in-memory MassSheets array.
 */
int mass_maps_write_sheet_table(void)
{
  if (ThisTask != 0)
    return 0;
  if (!MassSheets || NMassSheets <= 0)
    return 0;

  char fname[2 * LBLENGTH];
  snprintf(fname, sizeof(fname), "pinocchio.%s.sheets.out", params.RunFlag);
  FILE *fd = fopen(fname, "w");
  if (!fd)
  {
    fprintf(stderr, "Task 0 ERROR: cannot open %s for writing.\n", fname);
    return 1;
  }
  fprintf(fd, "# Mass sheet definitions\n");
  fprintf(fd, "# Columns:\n");
  fprintf(fd, "#  %-3s %-10s %-10s %-10s %-12s %-12s %-12s %-14s %-12s %-12s %-12s\n",
          "id", "z_hi", "z_lo", "delta_z", "chi_hi", "chi_lo", "delta_chi", "inv_delta_chi", "da_hi", "da_lo", "chi3_diff");
  for (int s = 0; s < NMassSheets; ++s)
  {
    MassSheet *ms = &MassSheets[s];
    fprintf(fd, "%3d %-10.8f %-10.8f %-10.8f %-12.8f %-12.8f %-12.8f %-14.12g %-12.8f %-12.8f %-12.8g\n",
            s, ms->z_hi, ms->z_lo, ms->delta_z, ms->chi_hi, ms->chi_lo, ms->delta_chi, ms->inv_dchi, ms->da_hi, ms->da_lo, ms->chi3_diff);
  }
  fclose(fd);
  if (internal.verbose_level >= VDIAG)
    printf("[%s] MASS_MAPS: wrote %s (N=%d)\n", fdate(), fname, NMassSheets);
  return 0;
}

/* ---------------------------------------------------------- */
/* Auxiliary geometry / crossing utilities (minimal version) */
/* ---------------------------------------------------------- */

/**
 * mass_maps_replication_chi_bounds
 * --------------------------------
 * Purpose
 *   Compute the radial comoving distance interval [rmin, rmax] (in Mpc)
 *   covered by a replication volume, using its stored scale factors F1, F2
 *   (where F = 1 + z) from plc.repls[].
 *
 * Parameters
 *   - rep_id  Replication index (0 .. plc.Nreplications-1). Must be valid.
 *   - rmin    Output pointer; receives the smaller of the two chi values.
 *   - rmax    Output pointer; receives the larger of the two chi values.
 *
 * Behavior
 *   - Converts F1,F2 to redshifts z1=F1-1 and z2=F2-1, then clamps negatives
 *     to zero (z<0 -> z=0).
 *   - Maps each z to comoving distance via ComovingDistance(z) and orders the
 *     results so that rmin <= rmax regardless of input ordering.
 *
 * Returns
 *   void. rmin and rmax must be non-NULL.
 *
 * Notes
 *   - Units: returns distances in Mpc (comoving).
 *   - If z1 == z2, then rmin == rmax.
 *   - This helper performs no bounds checks on rep_id and assumes plc.repls
 *     has been initialized.
 */
static inline void mass_maps_replication_chi_bounds(int rep_id, double *rmin, double *rmax)
{
  double F1 = plc.repls[rep_id].F1;
  double F2 = plc.repls[rep_id].F2;
  double z1 = F1 - 1.0;
  if (z1 < 0)
    z1 = 0;
  double z2 = F2 - 1.0;
  if (z2 < 0)
    z2 = 0;
  double chi1 = ComovingDistance(z1);
  double chi2 = ComovingDistance(z2);
  if (chi1 < chi2)
  {
    *rmin = chi1;
    *rmax = chi2;
  }
  else
  {
    *rmin = chi2;
    *rmax = chi1;
  }
}

/**
 * mass_maps_select_segment_replications
 * -------------------------------------
 * Purpose
 *   Build the list of replication indices whose radial comoving-distance
 *   interval [rmin, rmax] overlaps the current segment's chi span. This is a
 *   coarse radial prefilter; angular aperture culling is performed later per
 *   block using tighter bounds.
 *
 * Parameters
 *   - segment_index  Index into the ScaleDep outputs defining segment bounds.
 *                    Segment spans between ScaleDep.z[segment_index-1] (prev)
 *                    and ScaleDep.z[segment_index] (curr). For segment 0 or
 *                    invalid indices, the function returns early.
 *
 * Behavior
 *   - Computes chi_prev=ComovingDistance(z_prev) and chi_curr likewise, then
 *     sets seg_chi_min/max as the ordered pair.
 *   - Ensures the replication list buffer has capacity plc.Nreplications
 *     (realloc as needed). On allocation failure, logs a warning (task 0) and
 *     returns without populating the list.
 *   - For each replication r, obtains [rmin,rmax] via
 *     mass_maps_replication_chi_bounds(r, &rmin, &rmax) and appends r if the
 *     intervals overlap (rmax >= seg_chi_min && rmin <= seg_chi_max).
 *   - Updates the globals MassMapSegmentReplications (array of indices) and
 *     MassMapSegmentReplicationCount.
 *
 * Returns
 *   void. The replication candidate list reflects radial overlap only.
 *
 * Notes
 *   - Units: chi values are in Mpc (comoving). ComovingDistance increases with
 *     z, so when ScaleDep.z is in descending order, chi_prev > chi_curr.
 *   - Angular aperture filtering and per-block AABB tests are applied later
 *     in the block-processing loop to further reduce candidates.
 */
static void mass_maps_select_segment_replications(int segment_index)
{
  MassMapSegmentReplicationCount = 0;
  if (segment_index <= 0)
    return; /* first segment has no crossings */
  /* Segment spans between previous and current fragmentation redshifts */
  double z_prev = ScaleDep.z[segment_index - 1];
  double z_curr = ScaleDep.z[segment_index];
  /* Convert to chi; ComovingDistance grows with z so chi_prev > chi_curr given descending z array. */
  double chi_prev = ComovingDistance(z_prev);
  double chi_curr = ComovingDistance(z_curr);
  /* Segment radial span is [chi_curr, chi_prev] (chi_prev >= chi_curr). */
  double seg_chi_min = (chi_curr < chi_prev) ? chi_curr : chi_prev;
  double seg_chi_max = (chi_curr > chi_prev) ? chi_curr : chi_prev;
  /* Ensure capacity */
  if (MassMapSegmentReplicationsCapacity < plc.Nreplications)
  {
    int newcap = plc.Nreplications;
    int *tmp = (int *)realloc(MassMapSegmentReplications, sizeof(int) * newcap);
    if (!tmp)
    {
      if (!ThisTask)
        fprintf(stderr, "MASS_MAPS WARNING: cannot allocate segment replication list (N=%d).\n", plc.Nreplications);
      return;
    }
    MassMapSegmentReplications = tmp;
    MassMapSegmentReplicationsCapacity = newcap;
  }
  for (int r = 0; r < plc.Nreplications; ++r)
  {
    double rmin, rmax;
    mass_maps_replication_chi_bounds(r, &rmin, &rmax);
    if (rmax < seg_chi_min || rmin > seg_chi_max)
      continue; /* no overlap */
    MassMapSegmentReplications[MassMapSegmentReplicationCount++] = r;
  }
  if (!ThisTask && internal.verbose_level >= VDBG)
  {
    printf("[%s] MASS_MAPS: segment %d chi-span=[%g,%g] replication candidates=%d/%d\n",
           fdate(), segment_index, seg_chi_min, seg_chi_max,
           MassMapSegmentReplicationCount, plc.Nreplications);
  }
}

/* ---------------------------------------------------------- */
/* Displacement / position helpers                             */
/* ---------------------------------------------------------- */
/**
 * mass_maps_compute_prev_curr_positions
 * -------------------------------------
 * Purpose
 *   Build Eulerian positions for a particle at the previous and current
 *   segment endpoints from its stored LPT displacements, producing two
 *   absolute positions in physical units (true Mpc).
 *
 * Parameters
 *   - p          Pointer to pos_data with LPT displacement components.
 *   - pos_prev   Output array [3]; receives position at previous endpoint.
 *   - pos_curr   Output array [3]; receives position at current endpoint.
 *   - order      LPT order to include: 1, 2 (requires TWO_LPT), or 3
 *                (requires TWO_LPT and THREE_LPT). Higher terms are summed
 *                into both prev and curr as available.
 *
 * Behavior
 *   - Forms endpoint displacements by selecting per-endpoint components:
 *       prev: q + (v_prev [+ v2_prev [+ v31_prev + v32_prev]])
 *       curr: q + (v      [+ v2      [+ v31      + v32     ]])
 *   - Scales both positions from lattice units (cells) to physical box units
 *     using params.InterPartDist (true Mpc).
 *
 * Returns
 *   void. Outputs are written to pos_prev and pos_curr.
 *
 * Notes
 *   - Requires RECOMPUTE_DISPLACEMENTS (enforced for MASS_MAPS at compile time).
 *   - Replication shifts are NOT applied here; callers add them per replication.
 *   - Periodic wrapping, if desired, is also the caller's responsibility.
 */
static inline void mass_maps_compute_prev_curr_positions(const pos_data *p,
                                                         double pos_prev[3],
                                                         double pos_curr[3],
                                                         int order)
{
  /* Sum displacements per order (endpoints: weights collapse to selecting
     either *_prev or current components, so no interpolation weights used). */
  for (int i = 0; i < 3; ++i)
  {
    double dp_prev = p->v_prev[i];
    double dp_curr = p->v[i];
#ifdef TWO_LPT
    if (order > 1)
    {
      dp_prev += p->v2_prev[i];
      dp_curr += p->v2[i];
    }
#ifdef THREE_LPT
    if (order > 2)
    {
      dp_prev += p->v31_prev[i] + p->v32_prev[i];
      dp_curr += p->v31[i] + p->v32[i];
    }
#endif /* THREE_LPT */
#endif /* TWO_LPT */
    pos_prev[i] = p->q[i] + dp_prev;
    pos_curr[i] = p->q[i] + dp_curr;
  }
  /* Rescale from lattice units (grid cells) to physical box units using InterPartDist. */
  double scale = params.InterPartDist; /* typical: BoxSize / GridSize */
  pos_prev[0] *= scale;
  pos_prev[1] *= scale;
  pos_prev[2] *= scale;
  pos_curr[0] *= scale;
  pos_curr[1] *= scale;
  pos_curr[2] *= scale;
}

/* ---------------------------------------------------------- */
/* Build pos_data from products (FFT local tile, global coords) */
/* ---------------------------------------------------------- */
/**
 * mass_maps_set_point_from_products_global
 * ----------------------------------------
 * Purpose
 *   Populate a pos_data object for a particle referenced by its GLOBAL lattice
 *   coordinates (ig, jg, kg), loading displacement components from the
 *   products[] array on this task's FFT tile.
 *
 * Parameters
 *   - ig, jg, kg  Global lattice coordinates of the particle.
 *   - F           Scale factor (F = 1 + z) for the current segment endpoint.
 *   - myobj       Output pos_data to fill.
 *
 * Behavior
 *   - Translates global (ig,jg,kg) to local tile indices using MyGrids[0].
 *     If out of range, marks myobj->M=0 and returns (no data).
 *   - Sets myobj fields:
 *       z = F - 1;
 *       myk = Smoothing.k_GM_displ[last] (SCALE_DEPENDENT) or params.k_for_GM;
 *       weight via set_weight(myobj);
 *       M = 1;
 *       q[] = (ig+SHIFT, jg+SHIFT, kg+SHIFT) in GLOBAL lattice units;
 *       v[], and when enabled v2[], v31[], v32[] from products[idx];
 *       if RECOMPUTE_DISPLACEMENTS: v_prev[] and optionally v2_prev[],
 *       v31_prev[], v32_prev[].
 *
 * Returns
 *   void. On out-of-range indices, myobj is marked empty (M=0).
 *
 * Notes
 *   - q[] remains in lattice units; conversion to physical units is performed
 *     later when building Eulerian positions. Replication shifts are not
 *     applied here.
 *   - Caller is responsible for supplying a valid (ig,jg,kg) within this
 *     task's FFT tile when data is expected.
 */
static inline void mass_maps_set_point_from_products_global(int ig, int jg, int kg,
                                                            PRODFLOAT F,
                                                            pos_data *myobj)
{
  /* Translate global to local FFT indices */
  int li = ig - (int)MyGrids[0].GSstart[_x_];
  int lj = jg - (int)MyGrids[0].GSstart[_y_];
  int lk = kg - (int)MyGrids[0].GSstart[_z_];

  /* Safety: bound check (should not trigger if caller iterates tile properly) */
  if (li < 0 || lj < 0 || lk < 0 ||
      li >= (int)MyGrids[0].GSlocal[_x_] ||
      lj >= (int)MyGrids[0].GSlocal[_y_] ||
      lk >= (int)MyGrids[0].GSlocal[_z_])
  {
    /* Mark as empty */
    myobj->M = 0;
    return;
  }

  unsigned int idx = COORD_TO_INDEX(li, lj, lk, MyGrids[0].GSlocal);

  myobj->z = F - 1.0;
#ifdef SCALE_DEPENDENT
  int S = Smoothing.Nsmooth - 1;
  myobj->myk = Smoothing.k_GM_displ[S];
#else
  myobj->myk = params.k_for_GM;
#endif
  set_weight(myobj);

  myobj->M = 1;
  myobj->q[0] = ig + SHIFT;
  myobj->q[1] = jg + SHIFT;
  myobj->q[2] = kg + SHIFT;
  myobj->v[0] = products[idx].Vel[0];
  myobj->v[1] = products[idx].Vel[1];
  myobj->v[2] = products[idx].Vel[2];
#ifdef TWO_LPT
  myobj->v2[0] = products[idx].Vel_2LPT[0];
  myobj->v2[1] = products[idx].Vel_2LPT[1];
  myobj->v2[2] = products[idx].Vel_2LPT[2];
#ifdef THREE_LPT
  myobj->v31[0] = products[idx].Vel_3LPT_1[0];
  myobj->v31[1] = products[idx].Vel_3LPT_1[1];
  myobj->v31[2] = products[idx].Vel_3LPT_1[2];
  myobj->v32[0] = products[idx].Vel_3LPT_2[0];
  myobj->v32[1] = products[idx].Vel_3LPT_2[1];
  myobj->v32[2] = products[idx].Vel_3LPT_2[2];
#endif
#endif

#ifdef RECOMPUTE_DISPLACEMENTS
  myobj->v_prev[0] = products[idx].Vel_prev[0];
  myobj->v_prev[1] = products[idx].Vel_prev[1];
  myobj->v_prev[2] = products[idx].Vel_prev[2];
#ifdef TWO_LPT
  myobj->v2_prev[0] = products[idx].Vel_2LPT_prev[0];
  myobj->v2_prev[1] = products[idx].Vel_2LPT_prev[1];
  myobj->v2_prev[2] = products[idx].Vel_2LPT_prev[2];
#ifdef THREE_LPT
  myobj->v31_prev[0] = products[idx].Vel_3LPT_1_prev[0];
  myobj->v31_prev[1] = products[idx].Vel_3LPT_1_prev[1];
  myobj->v31_prev[2] = products[idx].Vel_3LPT_1_prev[2];
  myobj->v32_prev[0] = products[idx].Vel_3LPT_2_prev[0];
  myobj->v32_prev[1] = products[idx].Vel_3LPT_2_prev[1];
  myobj->v32_prev[2] = products[idx].Vel_3LPT_2_prev[2];
#endif
#endif
#endif
}

/**
 * mass_maps_get_zacc_global
 * -------------------------
 * Purpose
 *   Access the per-particle collapse redshift (ZACC) from the products array
 *   using GLOBAL lattice coordinates (ig, jg, kg) on this task's FFT tile.
 *
 * Parameters
 *   - ig, jg, kg  Global lattice indices of the particle.
 *
 * Returns
 *   The collapse redshift zacc as a double if the coordinates map inside this
 *   task's local tile; otherwise returns a large negative sentinel (-1e30)
 *   indicating "never collapsed / out of range".
 *
 * Notes
 *   - No MPI: caller must ensure the queried site resides on this rank if a
 *     valid value is expected.
 *   - Available only when SNAPSHOT is defined.
 */
#if defined(SNAPSHOT)
static inline double mass_maps_get_zacc_global(int ig, int jg, int kg)
{
  int li = ig - (int)MyGrids[0].GSstart[_x_];
  int lj = jg - (int)MyGrids[0].GSstart[_y_];
  int lk = kg - (int)MyGrids[0].GSstart[_z_];
  if (li < 0 || lj < 0 || lk < 0 || li >= (int)MyGrids[0].GSlocal[_x_] || lj >= (int)MyGrids[0].GSlocal[_y_] || lk >= (int)MyGrids[0].GSlocal[_z_])
    return -1e30; /* out of range; treat as never collapsed */
  unsigned int idx = COORD_TO_INDEX(li, lj, lk, MyGrids[0].GSlocal);
  return (double)products[idx].zacc;
}

/* Accessor for per-particle group_ID using global lattice coords (0 if none) */
/**
 * mass_maps_get_group_id_global
 * -----------------------------
 * Purpose
 *   Return the group_ID assigned to the particle at GLOBAL lattice
 *   coordinates (ig, jg, kg) on this task's FFT tile. If the coordinates are
 *   outside this tile, returns 0.
 *
 * Parameters
 *   - ig, jg, kg  Global lattice indices of the particle.
 *
 * Returns
 *   Integer group_ID for the particle (0 if none or if out of range).
 *
 * Behavior
 *   - Translates global to local indices using MyGrids[0]. If any index is
 *     outside the local tile bounds, returns 0.
 *   - Otherwise, loads products[idx].group_ID and returns it.
 *
 * Notes
 *   - Available only when SNAPSHOT is defined.
 *   - Caller must ensure that the queried (ig,jg,kg) lies on this task's tile
 *     when a non-zero ID is expected; no MPI communication occurs here.
 */
static inline int mass_maps_get_group_id_global(int ig, int jg, int kg)
{
  int li = ig - (int)MyGrids[0].GSstart[_x_];
  int lj = jg - (int)MyGrids[0].GSstart[_y_];
  int lk = kg - (int)MyGrids[0].GSstart[_z_];
  if (li < 0 || lj < 0 || lk < 0 || li >= (int)MyGrids[0].GSlocal[_x_] || lj >= (int)MyGrids[0].GSlocal[_y_] || lk >= (int)MyGrids[0].GSlocal[_z_])
    return 0;
  unsigned int idx = COORD_TO_INDEX(li, lj, lk, MyGrids[0].GSlocal);
  return products[idx].group_ID;
}
#endif

/**
 * mass_maps_particle_sign_change
 * ------------------------------
 * Test whether a particle leaves the PLC cone during the latest segment for a
 * specific replication, using H(z, r) = chi(z) - r with chi the comoving
 * distance to redshift and r the distance to the PLC center. A crossing is
 * detected when H_prev > 0 and H_curr <= 0 (inside -> outside across the
 * segment).
 *
 * Inputs:
 *   rep_id      - replication index
 *   q[3]        - base (Lagrangian) particle position
 *   disp_prev   - previous total Eulerian displacement vector
 *   disp_curr   - current total Eulerian displacement vector
 *   z_prev      - previous redshift (used to compute chi_prev)
 *   z_curr      - current redshift  (used to compute chi_curr)
 *   alpha_out   - (optional) on success, fraction in [0,1] from previous to current state
 *   entry_pos   - (optional) interpolated Eulerian position at crossing
 *   chi_cross_out - (optional) on success, interpolated comoving distance at crossing:
 *                  chi_prev + alpha * (chi_curr - chi_prev). Preferred over
 *                  |entry_pos - c| for deriving the crossing redshift.
 *
 * Notes:
 *   - If q == (0,0,0), disp_* are interpreted as absolute Eulerian positions.
 *     Otherwise positions are q + disp_*.
 *   - Replication shift equals (i,j,k) times the global box size in each axis.
 *
 * Logic:
 *   1. Build previous and current replicated positions.
 *   2. Compute r_prev, r_curr and chi_prev=ComovingDistance(z_prev), chi_curr.
 *   3. H_prev = chi_prev - r_prev, H_curr = chi_curr - r_curr.
 *   4. Crossing when H_prev > 0 and H_curr <= 0.
 *   5. Secant predictor: alpha0 = H_prev / (H_prev - H_curr); clamp to [0,1].
 *      Skip ambiguous cases where |H_prev - H_curr| < MASS_MAPS_F_EPS.
 *   6. One Newton corrector step on F(alpha) = chi(alpha) - |P(alpha) - c|
 *      to refine alpha. Skipped if |P(alpha0)-c| == 0 or derivative too small.
 *
 * Returns 1 if a crossing is detected per the above condition, else 0.
 */
int mass_maps_particle_sign_change(int rep_id,
                                   const double q[3],
                                   const double disp_prev[3],
                                   const double disp_curr[3],
                                   double z_prev,
                                   double z_curr,
                                   double *alpha_out,
                                   double entry_pos[3],
                                   double *chi_cross_out)
{
  if (rep_id < 0 || rep_id >= plc.Nreplications)
    return 0;
  int ii = plc.repls[rep_id].i;
  int jj = plc.repls[rep_id].j;
  int kk = plc.repls[rep_id].k;
  /* Replication shift in physical units: integer global box counts only (no subbox offsets) */
  double Bx = MyGrids[0].GSglobal[_x_] * params.InterPartDist;
  double By = MyGrids[0].GSglobal[_y_] * params.InterPartDist;
  double Bz = MyGrids[0].GSglobal[_z_] * params.InterPartDist;
  double shift[3] = {ii * Bx, jj * By, kk * Bz};

  /* If q==0 vector, treat disp_* as absolute positions (already q+dp). Otherwise treat as displacements relative to q. */
  int q_is_zero = (q[0] == 0.0 && q[1] == 0.0 && q[2] == 0.0);
  double pos_prev[3], pos_curr[3];
  if (q_is_zero)
  {
    pos_prev[0] = disp_prev[0] + shift[0];
    pos_prev[1] = disp_prev[1] + shift[1];
    pos_prev[2] = disp_prev[2] + shift[2];
    pos_curr[0] = disp_curr[0] + shift[0];
    pos_curr[1] = disp_curr[1] + shift[1];
    pos_curr[2] = disp_curr[2] + shift[2];
  }
  else
  {
    pos_prev[0] = q[0] + disp_prev[0] + shift[0];
    pos_prev[1] = q[1] + disp_prev[1] + shift[1];
    pos_prev[2] = q[2] + disp_prev[2] + shift[2];
    pos_curr[0] = q[0] + disp_curr[0] + shift[0];
    pos_curr[1] = q[1] + disp_curr[1] + shift[1];
    pos_curr[2] = q[2] + disp_curr[2] + shift[2];
  }

  /* Crossing logic:
     Define H_prev = chi_prev - r_prev and H_curr = chi_curr - r_curr.
     Crossing occurs if H_prev > 0 and H_curr <= 0. Estimate alpha by
     linear interpolation of H between the two endpoints. */
  double c0 = plc.center[0] * params.InterPartDist;
  double c1 = plc.center[1] * params.InterPartDist;
  double c2 = plc.center[2] * params.InterPartDist;
  double dxp0 = pos_prev[0] - c0, dyp0 = pos_prev[1] - c1, dzp0 = pos_prev[2] - c2;
  double dxp1 = pos_curr[0] - c0, dyp1 = pos_curr[1] - c1, dzp1 = pos_curr[2] - c2;
  double r_prev = sqrt(dxp0 * dxp0 + dyp0 * dyp0 + dzp0 * dzp0);
  double r_curr = sqrt(dxp1 * dxp1 + dyp1 * dyp1 + dzp1 * dzp1);
  double chi_prev = ComovingDistance(z_prev);
  double chi_curr = ComovingDistance(z_curr);
  double H_prev = chi_prev - r_prev;
  double H_curr = chi_curr - r_curr;

  /* Early exit if no crossing according to crossing condition */
  if (!(H_prev > 0.0 && H_curr <= 0.0))
  {
    static int sanity_no_cross_prints = 0;
    if (!ThisTask && internal.verbose_level >= VDIAG && sanity_no_cross_prints < 10)
    {
      printf("MASS_MAPS SANITY(no-entry): rep=%d z_prev=%.3f z_curr=%.3f H_prev=% .3e H_curr=% .3e r_prev=%.2f r_curr=%.2f chi_prev=%.2f chi_curr=%.2f\n",
             rep_id, z_prev, z_curr, H_prev, H_curr, r_prev, r_curr, chi_prev, chi_curr);
      sanity_no_cross_prints++;
    }
    return 0;
  }

  double denom = H_prev - H_curr;
  if (fabs(denom) < MASS_MAPS_F_EPS)
    return 0;                    /* ambiguous, skip */
  double alpha = H_prev / denom; /* secant predictor */
  /* Clamp to [0,1] just in case */
  if (alpha < 0.0)
    alpha = 0.0;
  else if (alpha > 1.0)
    alpha = 1.0;

  /* One Newton corrector step on F(alpha) = chi(alpha) - |P(alpha) - c| */
  {
    const double newton_eps = 1e-12;
    double dP[3] = {pos_curr[0] - pos_prev[0], pos_curr[1] - pos_prev[1], pos_curr[2] - pos_prev[2]};
    double dchi = chi_curr - chi_prev;
    double P0[3] = {pos_prev[0] + alpha * dP[0] - c0,
                    pos_prev[1] + alpha * dP[1] - c1,
                    pos_prev[2] + alpha * dP[2] - c2};
    double r0 = sqrt(P0[0] * P0[0] + P0[1] * P0[1] + P0[2] * P0[2]);
    if (r0 > 0.0)
    {
      double F0 = chi_prev + alpha * dchi - r0;
      double Fp0 = dchi - (P0[0] * dP[0] + P0[1] * dP[1] + P0[2] * dP[2]) / r0;
      if (fabs(Fp0) > newton_eps)
      {
        double alpha1 = alpha - F0 / Fp0;
        if (alpha1 < 0.0)
          alpha1 = 0.0;
        else if (alpha1 > 1.0)
          alpha1 = 1.0;
        alpha = alpha1;
      }
    }
  }

  /* Optional sanity print for a few crossings */
  {
    static int sanity_cross_prints = 0;
    if (!ThisTask && internal.verbose_level >= VDIAG && sanity_cross_prints < 10)
    {
      printf("MASS_MAPS SANITY(entry): rep=%d z_prev=%.3f z_curr=%.3f H_prev=% .3e H_curr=% .3e alpha=%.3f r_prev=%.2f r_curr=%.2f\n",
             rep_id, z_prev, z_curr, H_prev, H_curr, alpha, r_prev, r_curr);
      sanity_cross_prints++;
    }
  }
  if (alpha_out)
    *alpha_out = alpha;
  if (entry_pos)
  {
    entry_pos[0] = pos_prev[0] + alpha * (pos_curr[0] - pos_prev[0]);
    entry_pos[1] = pos_prev[1] + alpha * (pos_curr[1] - pos_prev[1]);
    entry_pos[2] = pos_prev[2] + alpha * (pos_curr[2] - pos_prev[2]);
  }
  if (chi_cross_out)
    *chi_cross_out = chi_prev + alpha * (chi_curr - chi_prev);
  return 1;
}

/* ---------------------------------------------------------- */
/* Orchestrator                                               */
/* ---------------------------------------------------------- */
/**
 * mass_maps_process_segment
 * -------------------------
 * Orchestrate mass-map accumulation for one fragmentation segment using the
 * "products" path (positions reconstructed from products and LPT displacements).
 *
 * Parameters:
 *   segment_index     - current segment index in [1..NMassSheets]; map index is s_map=segment_index-1
 *   z_segment         - redshift of the current segment (typically ScaleDep.z[segment_index])
 *   is_first_segment  - non-zero for the first segment; no crossings are computed and the function returns early
 *
 * Behavior:
 *   - Lazy allocation: on first invocation, allocate HEALPix maps for all sheets
 *     if params.MassMapNSIDE>0 (RING ordering).
 *   - Optional filter refresh: if MASS_MAPS_FILTER_UNCOLLAPSED and SNAPSHOT are defined,
 *     call distribute_back() so products[].zacc and products[].group_ID reflect current
 *     memberships before filtering/diagnostics; prints coverage diagnostics at VDBG.
 *   - Select replication candidates whose shifted Lagrangian AABBs overlap the segment
 *     in comoving distance (radial prune only).
 *   - Compute segment chi bounds and precompute per-candidate replication shifts
 *     in physical units (true Mpc).
 *   - Optional Lagrangian culling: if MASS_MAPS_BLOCKS_* > 0, split the local tile into
 *     blocks; build a padded AABB (variable pad in true Mpc via MASS_MAPS_VARPAD_FACTOR)
 *     and cull by radial [dmin,dmax] vs [chi_min,chi_max] plus an exact angular test:
 *       • Axis-piercing check in the PLC transverse plane (x/y spans contain 0) to early-accept.
 *       • Otherwise, enumerate the 8 AABB corners relative to the PLC center and compute the
 *         maximum corner cosine with the PLC axis; reject the block only if max_cos <= cos(A).
 *     Reuse cached prev/curr positions per block; optionally use per-thread LocalMap
 *     (MASS_MAPS_THREAD_ACCUM) to reduce atomic contention.
 *   - For each particle and passing replication:
 *       • Build previous and current Eulerian positions (mass_maps_compute_prev_curr_positions),
 *         apply replication shift.
 *       • Detect PLC boundary crossing inside the segment using
 *         mass_maps_crossing_from_shifted_positions with H=chi−r predicate.
 *       • Check angular aperture; optionally filter for uncollapsed-only (ZPLC>ZACC).
 *       • Accumulate one count into the HEALPix pixel of the entry position for sheet s_map.
 *   - Diagnostics: track per-rep and total entries; optional endpoint displacement vs buffer
 *     stats on the last segment; filter diagnostics and optional matching diagnostics.
 *   - MPI reductions: reduce per-rep counts and the segment map to rank 0, which writes the
 *     FITS file via mass_maps_write_segment_map(s_map) and prints map stats at VDBG.
 *
 * Threading and MPI:
 *   - OpenMP over blocks (when enabled); per-pixel updates use atomics or thread-local maps
 *     merged once per block. All ranks compute local contributions; rank 0 performs I/O.
 *
 * Units and conventions:
 *   - Positions and shifts are in true Mpc (InterPartDist units). ComovingDistance maps z→chi.
 *   - HEALPix uses RING ordering; aperture is in degrees. Map index s_map = segment_index-1.
 *
 * Early exits:
 *   - Feature disabled or not initialized (NMassSheets<=0 or params.MassMapNSIDE<=0).
 *   - First segment (no previous state to form crossings): returns immediately.
 *
 * Side effects:
 *   - Allocates global map storage on first use; persists across segments until freed by
 *     mass_maps_free_maps(). Updates global diagnostics used in FITS headers when filtering.
 **/
void mass_maps_process_segment(int segment_index, double z_segment, int is_first_segment)
{
  /* Reset cap diagnostics per segment */
#if MASS_MAPS_CAP_DIAG
  MassMapCapPrefixMiss = 0ULL;
  MassMapCapMissAngleMax = 0.0;
  MassMapCapAcceptAngleMax = 0.0;
#endif

  /* Ensure ZACC values are available in products when filtering uncollapsed-only */
#if defined(MASS_MAPS_FILTER_UNCOLLAPSED)
#if defined(SNAPSHOT)
  {
    /* Do back-distribution at EACH segment, as fragmentation updates across segments
       can assign new halo memberships and zacc values. This ensures the filter and
       matching diagnostics see up-to-date products[].zacc and products[].group_ID. */
    const char *method = getenv("PINOCCHIO_BACKDIST");
    double t_start, t_stop;

    if (method && strcmp(method, "bench") == 0)
    {
      if (!ThisTask && internal.verbose_level >= VDBG)
        printf("[%s] MASS_MAPS(filter): benchmarking back-distribution variants for segment %d\n", fdate(), segment_index);

      int rc_orig = 0, rc_round = 0, rc_all = 0;
      double dt_orig = 0.0, dt_round = 0.0, dt_all = 0.0;

      MPI_Barrier(MPI_COMM_WORLD);
      t_start = MPI_Wtime();
      rc_orig = distribute_back();
      MPI_Barrier(MPI_COMM_WORLD);
      t_stop = MPI_Wtime();
      dt_orig = t_stop - t_start;

      MPI_Barrier(MPI_COMM_WORLD);
      t_start = MPI_Wtime();
      rc_round = my_distribute_back();
      MPI_Barrier(MPI_COMM_WORLD);
      t_stop = MPI_Wtime();
      dt_round = t_stop - t_start;

      MPI_Barrier(MPI_COMM_WORLD);
      t_start = MPI_Wtime();
      rc_all = distribute_back_alltoall();
      MPI_Barrier(MPI_COMM_WORLD);
      t_stop = MPI_Wtime();
      dt_all = t_stop - t_start;

      if (!ThisTask)
        printf("[%s] MASS_MAPS(filter): back-distribution timings (segment %d): original=%.6fs roundrobin=%.6fs alltoall=%.6fs\n",
               fdate(), segment_index, dt_orig, dt_round, dt_all);

      if (!ThisTask && (rc_orig || rc_round || rc_all))
        fprintf(stderr, "[%s] MASS_MAPS(filter) WARNING: one or more back-distribution variants failed (orig=%d round=%d all=%d); ZACC may be incomplete.\n",
                fdate(), rc_orig, rc_round, rc_all);
    }
    else
    {
      int rc = 0;
      const char *label = NULL;

      if (method && strcmp(method, "original") == 0)
      {
        label = "original";
        MPI_Barrier(MPI_COMM_WORLD);
        t_start = MPI_Wtime();
        rc = distribute_back();
        MPI_Barrier(MPI_COMM_WORLD);
        t_stop = MPI_Wtime();
      }
      else if (method && strcmp(method, "round") == 0)
      {
        label = "roundrobin";
        MPI_Barrier(MPI_COMM_WORLD);
        t_start = MPI_Wtime();
        rc = my_distribute_back();
        MPI_Barrier(MPI_COMM_WORLD);
        t_stop = MPI_Wtime();
      }
      else
      {
        /* default to alltoall */
        label = "alltoall";
        MPI_Barrier(MPI_COMM_WORLD);
        t_start = MPI_Wtime();
        rc = distribute_back_alltoall();
        MPI_Barrier(MPI_COMM_WORLD);
        t_stop = MPI_Wtime();
      }

      if (!ThisTask)
        printf("[%s] MASS_MAPS(filter): back-distribution timing (segment %d, %s)=%.6fs\n",
               fdate(), segment_index, label, (t_stop - t_start));
      if (rc && !ThisTask)
        fprintf(stderr, "[%s] MASS_MAPS(filter) WARNING: back-distribution (%s) failed; ZACC may be unset (=-1) and the filter will include all.\n",
                fdate(), label);
    }
#ifdef MASS_MAPS_MATCH_DIAG
    /* After back-distribution, estimate coverage of group_ID and zacc on this task */
    unsigned long long gid_pos_local = 0ULL, zacc_set_local = 0ULL;
    unsigned long long gid_pos_global = 0ULL, zacc_set_global = 0ULL;
    /* Mismatch diagnostics: cells with zacc set but no halo ID, or vice versa */
    unsigned long long zacc_no_group_local = 0ULL, group_no_zacc_local = 0ULL;
    unsigned long long zacc_no_group_global = 0ULL, group_no_zacc_global = 0ULL;
    int nx = (int)MyGrids[0].GSlocal[_x_];
    int ny = (int)MyGrids[0].GSlocal[_y_];
    int nz = (int)MyGrids[0].GSlocal[_z_];
    for (int li = 0; li < nx; ++li)
      for (int lj = 0; lj < ny; ++lj)
        for (int lk = 0; lk < nz; ++lk)
        {
          unsigned int idx = COORD_TO_INDEX(li, lj, lk, MyGrids[0].GSlocal);
          int in_group = (products[idx].group_ID > FILAMENT);
          int has_zacc = (products[idx].zacc > -0.5);
          if (in_group)
            gid_pos_local++;
          if (has_zacc)
            zacc_set_local++;
          if (has_zacc && !in_group)
            zacc_no_group_local++;
          if (in_group && !has_zacc)
            group_no_zacc_local++;
        }
    MPI_Reduce(&gid_pos_local, &gid_pos_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&zacc_set_local, &zacc_set_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&zacc_no_group_local, &zacc_no_group_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&group_no_zacc_local, &group_no_zacc_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if (!ThisTask && internal.verbose_level >= VDBG)
    {
      unsigned long long total_local_size = (unsigned long long)MyGrids[0].GSglobal[_x_] * (unsigned long long)MyGrids[0].GSglobal[_y_] * (unsigned long long)MyGrids[0].GSglobal[_z_];
      double f_gid = (total_local_size > 0ULL) ? ((double)gid_pos_global / (double)total_local_size) : 0.0;
      double f_zacc = (total_local_size > 0ULL) ? ((double)zacc_set_global / (double)total_local_size) : 0.0;
      double f_zacc_no_gid = (total_local_size > 0ULL) ? ((double)zacc_no_group_global / (double)total_local_size) : 0.0;
      double f_gid_no_zacc = (total_local_size > 0ULL) ? ((double)group_no_zacc_global / (double)total_local_size) : 0.0;
      printf("[%s] MASS_MAPS(match): products coverage after back-distribution: frac(group_ID>FILAMENT)=%.3f frac(zacc_set)=%.3f frac(zacc_with_no_groupID)=%.3f frac(groupID_with_no_zacc)=%.3f\n",
             fdate(), f_gid, f_zacc, f_zacc_no_gid, f_gid_no_zacc);
    }
#endif
  }
#else
#error "MASS_MAPS_FILTER_UNCOLLAPSED requires SNAPSHOT"
#endif
#endif
  /* One-time HEALPix map allocator */
  {
    static int hpix_init_done = 0;
    if (!hpix_init_done)
    {
      if (params.MassMapNSIDE > 0)
        mass_maps_init_healpix_maps(params.MassMapNSIDE);
      hpix_init_done = 1;
    }
  }
  /* One-time units diagnostic to stdout (debug-level) */
  if (!ThisTask && internal.verbose_level >= VDBG)
  {
    static int printed_units = 0;
    if (!printed_units)
    {
      double c0 = plc.center[0] * params.InterPartDist;
      double c1 = plc.center[1] * params.InterPartDist;
      double c2 = plc.center[2] * params.InterPartDist;
      double Bx = MyGrids[0].GSglobal[_x_] * params.InterPartDist;
      double By = MyGrids[0].GSglobal[_y_] * params.InterPartDist;
      double Bz = MyGrids[0].GSglobal[_z_] * params.InterPartDist;
      printf("[%s] MASS_MAPS units: PLC center (phys)=(%.9g, %.9g, %.9g)  InterPartDist=%.9g  BoxSize=(%.9g, %.9g, %.9g)\n",
             fdate(), c0, c1, c2, params.InterPartDist, Bx, By, Bz);
      printed_units = 1;
    }
  }

  if (NMassSheets <= 0 || params.MassMapNSIDE <= 0)
    return; /* feature disabled or not initialized */
  if (is_first_segment)
    return; /* need previous segment for crossings */

  /* STEP 1: select replication candidates overlapping this segment radial span */
  mass_maps_select_segment_replications(segment_index);

  if (!ThisTask && internal.verbose_level >= VDBG)
  {
    int show = MassMapSegmentReplicationCount < 10 ? MassMapSegmentReplicationCount : 10;
    printf("[%s] MASS_MAPS: segment %d z=%g replication candidates %d/%d :",
           fdate(), segment_index, z_segment,
           MassMapSegmentReplicationCount, plc.Nreplications);
    for (int i = 0; i < show; ++i)
      printf(" %d", MassMapSegmentReplications[i]);
    if (MassMapSegmentReplicationCount > show)
      printf(" ...");
    printf("\n");
  }

  /* STEP 2: loop over replications and accumulate entries using products over the FFT tile */
  if (MassMapSegmentReplicationCount <= 0)
    return;

  /* Segment endpoint redshifts */
  double z_prev = ScaleDep.z[segment_index - 1];
  double z_curr = z_segment; /* expected to be ScaleDep.z[segment_index] */
  /* Optional: displacement histogram at endpoints (compile-time, default OFF) */
#ifdef DISP_HISTOGRAM
  mass_maps_disp_hist_on_segment(segment_index, z_prev, z_curr);
#endif
  /* Segment chi bounds */
  double chi_prev_seg = ComovingDistance(z_prev);
  double chi_curr_seg = ComovingDistance(z_curr);
  double seg_chi_min = chi_prev_seg < chi_curr_seg ? chi_prev_seg : chi_curr_seg;
  double seg_chi_max = chi_prev_seg > chi_curr_seg ? chi_prev_seg : chi_curr_seg;
  /* PLC center (phys) */
  double c_phys[3] = {plc.center[0] * params.InterPartDist,
                      plc.center[1] * params.InterPartDist,
                      plc.center[2] * params.InterPartDist};

  /* Determine LPT order compiled in */
  int lpt_order = 1;
#ifdef TWO_LPT
  lpt_order = 2;
#ifdef THREE_LPT
  lpt_order = 3;
#endif
#endif

  /* We'll count local crossings per replication and reduce globally */
  unsigned long long total_crossings_all_reps = 0ULL;
  unsigned long long total_crossings_all_reps_global = 0ULL;
#ifdef MASS_MAPS_FILTER_UNCOLLAPSED
  /* Diagnostics for ZPLC>ZACC filter: count how many entries were considered and excluded */
  unsigned long long zfilter_considered_local = 0ULL;
  unsigned long long zfilter_excluded_local = 0ULL;
  unsigned long long zfilter_considered_global = 0ULL;
  unsigned long long zfilter_excluded_global = 0ULL;
  /* Extra diagnostics: track ranges of z_plc and z_acc and check z_plc segment consistency */
  double zplc_min_local = 1e99, zplc_max_local = -1e99;
  double zacc_min_local = 1e99, zacc_max_local = -1e99;
  double zplc_min_global = 0.0, zplc_max_global = 0.0;
  double zacc_min_global = 0.0, zacc_max_global = 0.0;
  unsigned long long zplc_out_of_segment_local = 0ULL;
  unsigned long long zplc_out_of_segment_global = 0ULL;
#endif
  /* Store local per-rep counts to reduce in one shot */
  unsigned long long *rep_crossings_local = NULL;
  unsigned long long *rep_crossings_global = NULL;
  if (MassMapSegmentReplicationCount > 0)
  {
    rep_crossings_local = (unsigned long long *)calloc(MassMapSegmentReplicationCount, sizeof(unsigned long long));
    if (!ThisTask)
      rep_crossings_global = (unsigned long long *)calloc(MassMapSegmentReplicationCount, sizeof(unsigned long long));
  }

  /* Precompute scale factor F for set_point (1+z at current segment) */
  PRODFLOAT Fseg = (PRODFLOAT)(z_segment + 1.0);

  /* Optional: Lagrangian sub-volume culling setup (compile-time controlled) */
  int blocks_x = MASS_MAPS_BLOCKS_X;
  int blocks_y = MASS_MAPS_BLOCKS_Y;
  int blocks_z = MASS_MAPS_BLOCKS_Z;
  /* Unified AABB pad (phys units, same units as InterPartDist = true Mpc):
     Variable buffer based on (sigma8, z). We use the empirical form
       aabb_pad_phys = FACTOR * (sigma8/0.8)^{1.5} / (1+z)
     which is independent of resolution and in true Mpc to match positions formed by
     multiplying lattice coordinates by InterPartDist (true Mpc).
   */
  double aabb_pad_phys = 0.0;
  {
    double sigma8 = params.Sigma8;
    double z_now = z_curr; /* current segment redshift */
    double factor = MASS_MAPS_VARPAD_FACTOR;
    double h = params.Hubble100;
    /* Model defined in Mpc/h: pad_h = factor * (sigma8/0.8)^{1.5} / (1+z).
       Convert to true Mpc for internal positions by dividing by h. */
    aabb_pad_phys = (factor * pow(sigma8 / 0.8, 1.5) / (1.0 + z_now)) / (h > 0.0 ? h : 1.0); /* true Mpc */
  }
  /* Per-thread HEALPix accumulators: compile-time toggle */
#ifdef MASS_MAPS_THREAD_ACCUM
  int thread_accum = 1;
#else
  int thread_accum = 0;
#endif
  /* Precompute aperture and its cosine for angular culling */
  double A = params.PLCAperture;
  double cosA = cos(A * (acos(-1.0) / 180.0));
  int aperture_is_fullsky = (A >= 180.0 - 1e-9);
  int nx_local = (int)MyGrids[0].GSlocal[_x_];
  int ny_local = (int)MyGrids[0].GSlocal[_y_];
  int nz_local = (int)MyGrids[0].GSlocal[_z_];
  int gx_start = (int)MyGrids[0].GSstart[_x_];
  int gy_start = (int)MyGrids[0].GSstart[_y_];
  int gz_start = (int)MyGrids[0].GSstart[_z_];
  /* Cap blocks to local sizes */
  if (blocks_x > nx_local)
    blocks_x = nx_local;
  if (blocks_y > ny_local)
    blocks_y = ny_local;
  if (blocks_z > nz_local)
    blocks_z = nz_local;
  int culling_enabled = (blocks_x > 0 && blocks_y > 0 && blocks_z > 0);
  if (!ThisTask && internal.verbose_level >= VDIAG)
  {
    double pad_cells = (params.InterPartDist > 0.0) ? (aabb_pad_phys / params.InterPartDist) : 0.0;
    printf("[%s] MASS_MAPS culling: blocks=(%d,%d,%d) aabb_pad=%.3g Mpc (%.3g cells) enabled=%d (Lagrangian-only)\n",
           fdate(), blocks_x, blocks_y, blocks_z, aabb_pad_phys, pad_cells, culling_enabled);
  }

  /* Last-segment displacement vs buffer diagnostics threshold */
  double warn_frac = MASS_MAPS_VARPAD_WARN_FRAC;
  int last_segment = (segment_index == NMassSheets);
  int collect_disp_stats = (last_segment && (aabb_pad_phys > 0.0));
  unsigned long long disp_total_local = 0ULL;  /* number of particles evaluated */
  unsigned long long disp_exceed_local = 0ULL; /* number with |Pcurr-Pprev| > AABB pad */
  double disp_max_local = 0.0;                 /* maximum |Pcurr-Pprev| observed */
  if (!ThisTask && internal.verbose_level >= VDIAG)
  {
    if (aabb_pad_phys > 0.0)
      printf("[%s] MASS_MAPS culling: using unified AABB pad phys=%.6g\n", fdate(), aabb_pad_phys);
    if (collect_disp_stats)
      printf("[%s] MASS_MAPS buffer-diag: last segment will collect |x(q,z)-q| stats vs buffer pad=%.6g Mpc (%.6g cells) (warn frac=%.3g)\n",
             fdate(), aabb_pad_phys, (aabb_pad_phys / params.InterPartDist), warn_frac);
  }
  /* Build block edge indices in local coordinates (inclusive ranges per block) */
  int *edge_x = NULL, *edge_y = NULL, *edge_z = NULL;
  int nbx = 0, nby = 0, nbz = 0;
  if (culling_enabled)
  {
    nbx = blocks_x;
    nby = blocks_y;
    nbz = blocks_z;
    edge_x = (int *)malloc(sizeof(int) * (size_t)(nbx + 1));
    edge_y = (int *)malloc(sizeof(int) * (size_t)(nby + 1));
    edge_z = (int *)malloc(sizeof(int) * (size_t)(nbz + 1));
    if (!edge_x || !edge_y || !edge_z)
    {
      culling_enabled = 0;
    }
    else
    {
      for (int b = 0; b <= nbx; ++b)
        edge_x[b] = (int)((long)b * (long)nx_local / (long)nbx);
      for (int b = 0; b <= nby; ++b)
        edge_y[b] = (int)((long)b * (long)ny_local / (long)nby);
      for (int b = 0; b <= nbz; ++b)
        edge_z[b] = (int)((long)b * (long)nz_local / (long)nbz);
    }
  }

  /* Precompute box shifts for all replication candidates once
   * Rationale: avoids recomputing (ii,jj,kk)*Box for every particle.
   */
  double Bx = MyGrids[0].GSglobal[_x_] * params.InterPartDist;
  double By = MyGrids[0].GSglobal[_y_] * params.InterPartDist;
  double Bz = MyGrids[0].GSglobal[_z_] * params.InterPartDist;
  double *RepShift = NULL;
  if (MassMapSegmentReplicationCount > 0)
  {
    RepShift = (double *)malloc(sizeof(double) * 3 * (size_t)MassMapSegmentReplicationCount);
    if (RepShift)
    {
      for (int ir = 0; ir < MassMapSegmentReplicationCount; ++ir)
      {
        int rep_id = MassMapSegmentReplications[ir];
        int ii = plc.repls[rep_id].i;
        int jj = plc.repls[rep_id].j;
        int kk = plc.repls[rep_id].k;
        RepShift[3 * ir + 0] = ii * Bx;
        RepShift[3 * ir + 1] = jj * By;
        RepShift[3 * ir + 2] = kk * Bz;
      }
    }
  }

  /* Culling-enabled path: iterate blocks first to reuse particle positions across replications
   * Hotspot: the inner loops over particles and replications are dominant.
   * Current approach uses a per-block cache of positions to reduce recomputation
   * and OpenMP over blocks. Map updates default to atomics; an optional per-thread
   * accumulator (MASS_MAPS_THREAD_ACCUM) further reduces contention with a final merge.
   */
  if (culling_enabled && RepShift)
  {
    int s_map = segment_index - 1;
    /* Diagnostics: track block visitation effectiveness */
    long long blocks_total = (long long)nbx * (long long)nby * (long long)nbz;
    long long blocks_init_count = 0;    /* blocks with any particles (AABB initialized) */
    long long blocks_skipped_shell = 0; /* initialized blocks rejected by shell test for all reps */
    long long blocks_visited = 0;       /* initialized blocks with at least one passing replication */
                                        /* Parallelize over blocks; each iteration handles one block fully */
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(dynamic)
#endif
    for (int bz = 0; bz < nbz; ++bz)
      for (int by = 0; by < nby; ++by)
        for (int bx_i = 0; bx_i < nbx; ++bx_i)
        {
          int kg0_local = edge_z[bz];
          int kg1_local = edge_z[bz + 1] - 1;
          if (kg1_local < kg0_local)
            continue;
          int kg0 = gz_start + kg0_local;
          int kg1 = gz_start + kg1_local;
          int jg0_local = edge_y[by];
          int jg1_local = edge_y[by + 1] - 1;
          if (jg1_local < jg0_local)
            continue;
          int jg0 = gy_start + jg0_local;
          int jg1 = gy_start + jg1_local;
          int ig0_local = edge_x[bx_i];
          int ig1_local = edge_x[bx_i + 1] - 1;
          if (ig1_local < ig0_local)
            continue;
          int ig0 = gx_start + ig0_local;
          int ig1 = gx_start + ig1_local;
          /* Determine whether this block has content to consider */
          int block_initialized = 1; /* Lagrangian-only: any non-empty index range is considered initialized */
          if (!block_initialized)
            continue;
#ifdef _OPENMP
#pragma omp atomic update
#endif
          blocks_init_count++;

          /* Per-block buffer threshold for diagnostics */
          double scale_block = params.InterPartDist;
          double thresh_block = aabb_pad_phys;

          /* Build replication pass list for this block */
          int *rep_pass = (int *)alloca(sizeof(int) * (size_t)MassMapSegmentReplicationCount);
          int pass_n = 0;
          for (int ir = 0; ir < MassMapSegmentReplicationCount; ++ir)
          {
            const double *shift = &RepShift[3 * ir];
            /* Lagrangian AABB expanded by unified pad (phys units) */
            double min_s[3] = {(ig0 + SHIFT) * scale_block, (jg0 + SHIFT) * scale_block, (kg0 + SHIFT) * scale_block};
            double max_s[3] = {(ig1 + SHIFT) * scale_block, (jg1 + SHIFT) * scale_block, (kg1 + SHIFT) * scale_block};
            /* Effective AABB padding (variable): use per-segment aabb_pad_phys only */
            double thresh_block = aabb_pad_phys;
            double pad_total = aabb_pad_phys;
            for (int a = 0; a < 3; ++a)
            {
              min_s[a] -= pad_total;
              max_s[a] += pad_total;
              min_s[a] += shift[a];
              max_s[a] += shift[a];
            }
            double dmin, dmax;
            aabb_distance_bounds_to_point(min_s, max_s, c_phys, &dmin, &dmax);
            if (dmin > seg_chi_max || dmax < seg_chi_min)
              continue; /* radial shell rejects */
            /* Conservative angular cull using cosine upper bound:
               If max_{box} (u·r/|r|) ≤ cosA, the whole AABB is outside aperture.
               Upper bound <= (max_u_dot) / (min_norm). */
            double u[3] = {plc.zvers[0], plc.zvers[1], plc.zvers[2]}; /* unit axis */
            /* Build relative bounds r = x - c */
            double rmin[3] = {min_s[0] - c_phys[0], min_s[1] - c_phys[1], min_s[2] - c_phys[2]};
            double rmax[3] = {max_s[0] - c_phys[0], max_s[1] - c_phys[1], max_s[2] - c_phys[2]};
            /* Robust angular culling: exact corner-based upper bound.
               Enumerate all 8 corners of the padded AABB relative to PLC center.
               Compute cos(theta) = (u·r)/|r| at each corner; if the maximum
               corner cosine <= cosA, the entire box lies outside the cap.
               If PLC axis pierces the transverse rectangle (x,y spans contain 0),
               an interior point with cos(theta)=1 exists, so skip rejection. */
            int angular_pass = 1;
            if (!aperture_is_fullsky)
            {
              /* Axis piercing test via projections onto PLC x/y basis */
              double ex[3] = {plc.xvers[0], plc.xvers[1], plc.xvers[2]};
              double ey[3] = {plc.yvers[0], plc.yvers[1], plc.yvers[2]};
              double x_min = 0.0, x_max = 0.0, y_min = 0.0, y_max = 0.0;
              /* Extremal dot with ex */
              for (int a = 0; a < 3; ++a)
              {
                double va = ex[a];
                double lo = rmin[a];
                double hi = rmax[a];
                double add_min = (va >= 0.0) ? va * lo : va * hi;
                double add_max = (va >= 0.0) ? va * hi : va * lo;
                x_min += add_min;
                x_max += add_max;
                va = ey[a];
                add_min = (va >= 0.0) ? va * lo : va * hi;
                add_max = (va >= 0.0) ? va * hi : va * lo;
                y_min += add_min;
                y_max += add_max;
              }
              int axis_pierces = (x_min <= 0.0 && x_max >= 0.0 && y_min <= 0.0 && y_max >= 0.0);
              if (!axis_pierces)
              {
                double max_corner_cos = -1.0; /* track maximum corner cosine */
                int early_accept = 0;
                for (int ax = 0; ax < 2 && !early_accept; ++ax)
                  for (int ay = 0; ay < 2 && !early_accept; ++ay)
                    for (int az = 0; az < 2; ++az)
                    {
                      double rc[3] = {
                          (ax ? rmax[0] : rmin[0]),
                          (ay ? rmax[1] : rmin[1]),
                          (az ? rmax[2] : rmin[2])};
                      double dot = u[0] * rc[0] + u[1] * rc[1] + u[2] * rc[2];
                      double norm = sqrt(rc[0] * rc[0] + rc[1] * rc[1] + rc[2] * rc[2]);
                      if (norm == 0.0)
                        continue;
                      double cosv = dot / norm;
                      if (cosv > max_corner_cos)
                        max_corner_cos = cosv;
                      if (max_corner_cos > cosA)
                      {
                        early_accept = 1; /* no rejection possible */
                        break;
                      }
                    }
                if (!early_accept && max_corner_cos <= cosA)
                  angular_pass = 0; /* all corners outside spherical cap */
              }
            }
            if (!angular_pass)
              continue;
            rep_pass[pass_n++] = ir;
          }
          if (pass_n == 0)
          {
#ifdef _OPENMP
#pragma omp atomic update
#endif
            blocks_skipped_shell++;
            continue; /* no replication needs this block */
          }
#ifdef _OPENMP
#pragma omp atomic update
#endif
          blocks_visited++;

          /* Prepare buffers for this block's particles (thread-local) */
          long nx = (long)(ig1 - ig0 + 1);
          long ny = (long)(jg1 - jg0 + 1);
          long nz = (long)(kg1 - kg0 + 1);
          long nblock = nx * ny * nz;
          if (nblock <= 0)
            continue;
          double *PosPrevBuf = (double *)malloc(sizeof(double) * 3 * (size_t)nblock);
          double *PosCurrBuf = (double *)malloc(sizeof(double) * 3 * (size_t)nblock);
          /* Optional per-thread accumulator for this block to reduce atomics */
          double *LocalMap = NULL;
          if (thread_accum && MassMapNSIDE_current > 0)
            LocalMap = (double *)calloc((size_t)MassMapNPIX, sizeof(double));
          unsigned long long ipix_oob_local = 0ULL; /* guard counter */
          int use_cache = (PosPrevBuf && PosCurrBuf);
          if (use_cache)
          {
            /* Fill cached positions once */
            long idx = 0;
            /* Per-block diagnostics */
            unsigned long long blk_total = 0ULL, blk_exceed = 0ULL;
            double blk_max = 0.0;
            for (int ig = ig0; ig <= ig1; ++ig)
            {
              for (int jg = jg0; jg <= jg1; ++jg)
              {
                for (int kg = kg0; kg <= kg1; ++kg, ++idx)
                {
                  pos_data obj;
                  mass_maps_set_point_from_products_global(ig, jg, kg, Fseg, &obj);
                  if (obj.M == 0)
                  {
                    /* Mark as zero-length by duplicating prev/curr */
                    PosPrevBuf[3 * idx + 0] = PosPrevBuf[3 * idx + 1] = PosPrevBuf[3 * idx + 2] = 0.0;
                    PosCurrBuf[3 * idx + 0] = PosCurrBuf[3 * idx + 1] = PosCurrBuf[3 * idx + 2] = 0.0;
                    continue;
                  }
                  double pprev[3], pcurr[3];
                  mass_maps_compute_prev_curr_positions(&obj, pprev, pcurr, lpt_order);
                  PosPrevBuf[3 * idx + 0] = pprev[0];
                  PosPrevBuf[3 * idx + 1] = pprev[1];
                  PosPrevBuf[3 * idx + 2] = pprev[2];
                  PosCurrBuf[3 * idx + 0] = pcurr[0];
                  PosCurrBuf[3 * idx + 1] = pcurr[1];
                  PosCurrBuf[3 * idx + 2] = pcurr[2];
                  if (collect_disp_stats)
                  {
                    /* Displacement relative to Lagrangian position q */
                    double qx = (ig + SHIFT) * scale_block;
                    double qy = (jg + SHIFT) * scale_block;
                    double qz = (kg + SHIFT) * scale_block;
                    double dx = pcurr[0] - qx;
                    double dy = pcurr[1] - qy;
                    double dz = pcurr[2] - qz;
                    double d = sqrt(dx * dx + dy * dy + dz * dz);
                    blk_total++;
                    if (d > thresh_block)
                      blk_exceed++;
                    if (d > blk_max)
                      blk_max = d;
                  }
                }
              }
            }
            if (collect_disp_stats)
            {
#ifdef _OPENMP
#pragma omp atomic update
#endif
              disp_total_local += blk_total;
#ifdef _OPENMP
#pragma omp atomic update
#endif
              disp_exceed_local += blk_exceed;
#ifdef _OPENMP
#pragma omp critical
#endif
              {
                if (blk_max > disp_max_local)
                  disp_max_local = blk_max;
              }
            }

            /* Use cached positions for each passing replication */
            for (int ip = 0; ip < pass_n; ++ip)
            {
              int ir = rep_pass[ip];
              const double *shift = &RepShift[3 * ir];
              unsigned long long local_cross = 0ULL;
              long idx2 = 0;
              for (int ig = ig0; ig <= ig1; ++ig)
              {
                for (int jg = jg0; jg <= jg1; ++jg)
                {
                  for (int kg = kg0; kg <= kg1; ++kg, ++idx2)
                  {
                    double Pprev[3] = {PosPrevBuf[3 * idx2 + 0] + shift[0],
                                       PosPrevBuf[3 * idx2 + 1] + shift[1],
                                       PosPrevBuf[3 * idx2 + 2] + shift[2]};
                    double Pcurr[3] = {PosCurrBuf[3 * idx2 + 0] + shift[0],
                                       PosCurrBuf[3 * idx2 + 1] + shift[1],
                                       PosCurrBuf[3 * idx2 + 2] + shift[2]};
                    /* Skip if empty marked */
                    if (Pprev[0] == 0.0 && Pprev[1] == 0.0 && Pprev[2] == 0.0 &&
                        Pcurr[0] == 0.0 && Pcurr[1] == 0.0 && Pcurr[2] == 0.0)
                      continue;
                    double entry_pos[3];
                    double chi_cross;
                    if (!mass_maps_crossing_from_shifted_positions(Pprev, Pcurr, c_phys, chi_prev_seg, chi_curr_seg, NULL, entry_pos, &chi_cross))
                      continue;
                    /* Pixel-based angular selection: compute target pixel and keep only if inside cap */
                    long ipix_sel;
                    if (!mass_maps_compute_pixel_from_pos(entry_pos, MassMapNSIDE_current, &ipix_sel))
                      continue;
                    {
                      int prefix_accept = ((unsigned long)ipix_sel < (unsigned long)MassMapNCAP);
#if MASS_MAPS_CAP_DIAG
                      if (MassMapCapDiagActive && params.PLCAperture < 180.0)
                      {
                        int geom_accept = mass_maps_entry_inside_aperture(entry_pos);
                        double ang = mass_maps_angle_deg_from_pos(entry_pos);
                        if (geom_accept && !prefix_accept)
                        {
#ifdef _OPENMP
#pragma omp atomic update
#endif
                          MassMapCapPrefixMiss++;
#ifdef _OPENMP
#pragma omp critical
#endif
                          {
                            if (ang > MassMapCapMissAngleMax)
                              MassMapCapMissAngleMax = ang;
                          }
                        }
                        if (prefix_accept)
                        {
#ifdef _OPENMP
#pragma omp critical
#endif
                          {
                            if (ang > MassMapCapAcceptAngleMax)
                              MassMapCapAcceptAngleMax = ang;
                          }
                        }
                      }
#endif
                      if (params.PLCAperture < 180.0 && !prefix_accept)
                        continue;
                    }
#ifdef MASS_MAPS_FILTER_UNCOLLAPSED
                    {
                      /* Filter: include only if ZPLC > ZACC */
                      double z_plc = InverseComovingDistance(chi_cross);
#ifdef SNAPSHOT
                      double z_acc = mass_maps_get_zacc_global(ig, jg, kg);
#ifdef _OPENMP
#pragma omp atomic update
#endif
                      zfilter_considered_local++;
                  /* Track min/max and consistency for debug */
#ifdef _OPENMP
#pragma omp critical
#endif
                      {
                        if (z_plc < zplc_min_local)
                          zplc_min_local = z_plc;
                        if (z_plc > zplc_max_local)
                          zplc_max_local = z_plc;
                        if (z_acc < zacc_min_local)
                          zacc_min_local = z_acc;
                        if (z_acc > zacc_max_local)
                          zacc_max_local = z_acc;
                      }
                      {
                        /* Check that z_plc lies within this segment bounds (tolerance at edges) */
                        double zmin_seg = (z_prev < z_curr) ? z_prev : z_curr;
                        double zmax_seg = (z_prev > z_curr) ? z_prev : z_curr;
                        double tol = 1e-3;
                        if (z_plc < zmin_seg - tol || z_plc > zmax_seg + tol)
                        {
#ifdef _OPENMP
#pragma omp atomic update
#endif
                          zplc_out_of_segment_local++;
                        }
                      }
                      if (!(z_plc > z_acc))
                      {
#ifdef _OPENMP
#pragma omp atomic update
#endif
                        zfilter_excluded_local++;
                        continue;
                      }
#else
#error "MASS_MAPS_FILTER_UNCOLLAPSED requires SNAPSHOT"
#endif
                    }
#endif
                    local_cross++;
                    if (s_map >= 0 && s_map < NMassSheets && MassMapNSIDE_current > 0)
                    {
                      if ((unsigned long)ipix_sel < (unsigned long)MassMapNPIX)
                      {
                        if (LocalMap)
                        {
                          LocalMap[ipix_sel] += 1.0;
                        }
                        else
                        {
                          double *map = mass_maps_segment_ptr(s_map);
#ifdef _OPENMP
#pragma omp atomic update
#endif
                          map[ipix_sel] += 1.0;
                        }
                      }
                    }
                  }
                }
              }
          /* Thread-safe global counters */
#ifdef _OPENMP
#pragma omp atomic update
#endif
              total_crossings_all_reps += local_cross;
              if (rep_crossings_local)
              {
#ifdef _OPENMP
#pragma omp atomic update
#endif
                rep_crossings_local[ir] += local_cross;
              }
            }
          }
          else
          {
            /* No cache available: compute per particle per passing replication */
            /* Per-block diagnostics */
            unsigned long long blk_total = 0ULL, blk_exceed = 0ULL;
            double blk_max = 0.0;
            for (int ip = 0; ip < pass_n; ++ip)
            {
              int ir = rep_pass[ip];
              const double *shift = &RepShift[3 * ir];
              unsigned long long local_cross = 0ULL;
              for (int ig = ig0; ig <= ig1; ++ig)
              {
                for (int jg = jg0; jg <= jg1; ++jg)
                {
                  for (int kg = kg0; kg <= kg1; ++kg)
                  {
                    pos_data obj;
                    mass_maps_set_point_from_products_global(ig, jg, kg, Fseg, &obj);
                    if (obj.M == 0)
                      continue;
                    double pos_prev[3], pos_curr[3];
                    mass_maps_compute_prev_curr_positions(&obj, pos_prev, pos_curr, lpt_order);
                    if (collect_disp_stats && ip == 0)
                    {
                      double qx = (ig + SHIFT) * scale_block;
                      double qy = (jg + SHIFT) * scale_block;
                      double qz = (kg + SHIFT) * scale_block;
                      double dx = pos_curr[0] - qx;
                      double dy = pos_curr[1] - qy;
                      double dz = pos_curr[2] - qz;
                      double d = sqrt(dx * dx + dy * dy + dz * dz);
                      blk_total++;
                      if (d > thresh_block)
                        blk_exceed++;
                      if (d > blk_max)
                        blk_max = d;
                    }
                    double Pprev[3] = {pos_prev[0] + shift[0], pos_prev[1] + shift[1], pos_prev[2] + shift[2]};
                    double Pcurr[3] = {pos_curr[0] + shift[0], pos_curr[1] + shift[1], pos_curr[2] + shift[2]};
                    double entry_pos[3];
                    double chi_cross;
                    if (!mass_maps_crossing_from_shifted_positions(Pprev, Pcurr, c_phys, chi_prev_seg, chi_curr_seg, NULL, entry_pos, &chi_cross))
                      continue;
                    /* Pixel-based angular selection: compute target pixel and keep only if inside cap */
                    long ipix_sel;
                    if (!mass_maps_compute_pixel_from_pos(entry_pos, MassMapNSIDE_current, &ipix_sel))
                      continue;
                    {
                      int prefix_accept = ((unsigned long)ipix_sel < (unsigned long)MassMapNCAP);
#if MASS_MAPS_CAP_DIAG
                      if (MassMapCapDiagActive && params.PLCAperture < 180.0)
                      {
                        int geom_accept = mass_maps_entry_inside_aperture(entry_pos);
                        double ang = mass_maps_angle_deg_from_pos(entry_pos);
                        if (geom_accept && !prefix_accept)
                        {
#ifdef _OPENMP
#pragma omp atomic update
#endif
                          MassMapCapPrefixMiss++;
#ifdef _OPENMP
#pragma omp critical
#endif
                          {
                            if (ang > MassMapCapMissAngleMax)
                              MassMapCapMissAngleMax = ang;
                          }
                        }
                        if (prefix_accept)
                        {
#ifdef _OPENMP
#pragma omp critical
#endif
                          {
                            if (ang > MassMapCapAcceptAngleMax)
                              MassMapCapAcceptAngleMax = ang;
                          }
                        }
                      }
#endif
                      if (params.PLCAperture < 180.0 && !prefix_accept)
                        continue;
                    }
#ifdef MASS_MAPS_FILTER_UNCOLLAPSED
                    {
                      double z_plc = InverseComovingDistance(chi_cross);
#ifdef SNAPSHOT
                      double z_acc = mass_maps_get_zacc_global(ig, jg, kg);
#ifdef _OPENMP
#pragma omp atomic update
#endif
                      zfilter_considered_local++;
                  /* Track min/max and consistency for debug */
#ifdef _OPENMP
#pragma omp critical
#endif
                      {
                        if (z_plc < zplc_min_local)
                          zplc_min_local = z_plc;
                        if (z_plc > zplc_max_local)
                          zplc_max_local = z_plc;
                        if (z_acc < zacc_min_local)
                          zacc_min_local = z_acc;
                        if (z_acc > zacc_max_local)
                          zacc_max_local = z_acc;
                      }
                      {
                        /* Check that z_plc lies within this segment bounds (tolerance at edges) */
                        double zmin_seg = (z_prev < z_curr) ? z_prev : z_curr;
                        double zmax_seg = (z_prev > z_curr) ? z_prev : z_curr;
                        double tol = 1e-3;
                        if (z_plc < zmin_seg - tol || z_plc > zmax_seg + tol)
                        {
#ifdef _OPENMP
#pragma omp atomic update
#endif
                          zplc_out_of_segment_local++;
                        }
                      }
                      if (!(z_plc > z_acc))
                      {
#ifdef _OPENMP
#pragma omp atomic update
#endif
                        zfilter_excluded_local++;
                        continue;
                      }
#else
#error "MASS_MAPS_FILTER_UNCOLLAPSED requires SNAPSHOT"
#endif
                    }
#endif
                    local_cross++;
                    if (s_map >= 0 && s_map < NMassSheets && MassMapNSIDE_current > 0)
                    {
                      if ((unsigned long)ipix_sel < (unsigned long)MassMapNPIX)
                      {
                        double *map = mass_maps_segment_ptr(s_map);
#ifdef _OPENMP
#pragma omp atomic update
#endif
                        map[ipix_sel] += 1.0;
                      }
                    }
                  }
                }
              }
          /* Thread-safe global counters */
#ifdef _OPENMP
#pragma omp atomic update
#endif
              total_crossings_all_reps += local_cross;
              if (rep_crossings_local)
              {
#ifdef _OPENMP
#pragma omp atomic update
#endif
                rep_crossings_local[ir] += local_cross;
              }
            }
            if (collect_disp_stats)
            {
#ifdef _OPENMP
#pragma omp atomic update
#endif
              disp_total_local += blk_total;
#ifdef _OPENMP
#pragma omp atomic update
#endif
              disp_exceed_local += blk_exceed;
#ifdef _OPENMP
#pragma omp critical
#endif
              {
                if (blk_max > disp_max_local)
                  disp_max_local = blk_max;
              }
            }
          }
          /* If we used a local map, merge it into the global segment map once per block */
          if (LocalMap)
          {
            double *map = mass_maps_segment_ptr(s_map);
            if (map)
            {
              for (long i = 0; i < MassMapNPIX; ++i)
              {
                double v = LocalMap[i];
                if (v == 0.0)
                  continue;
#ifdef _OPENMP
#pragma omp atomic update
#endif
                map[i] += v;
              }
            }
          }
          if (ipix_oob_local > 0ULL && internal.verbose_level >= VDBG)
          {
            /* Per-block debug: report any out-of-cap pixel indices (should be zero) */
            printf("[%s] MASS_MAPS WARN: block (%d,%d,%d) ipix>=Ncap occurrences: %llu\n", fdate(), bx_i, by, bz, ipix_oob_local);
          }
          if (PosPrevBuf)
            free(PosPrevBuf);
          if (PosCurrBuf)
            free(PosCurrBuf);
          if (LocalMap)
            free(LocalMap);
        }

    /* Print diagnostics on block visitation */
    if (!ThisTask && internal.verbose_level >= VDBG)
    {
      double init_pct = (blocks_total > 0) ? (100.0 * (double)blocks_init_count / (double)blocks_total) : 0.0;
      double visit_pct = (blocks_total > 0) ? (100.0 * (double)blocks_visited / (double)blocks_total) : 0.0;
      double skip_pct = (blocks_total > 0) ? (100.0 * (double)blocks_skipped_shell / (double)blocks_total) : 0.0;
      printf("[%s] MASS_MAPS culling: blocks total=%lld init=%lld (%.1f%%) visited=%lld (%.1f%%) skipped_by_shell=%lld (%.1f%%)\n",
             fdate(), blocks_total, blocks_init_count, init_pct, blocks_visited, visit_pct, blocks_skipped_shell, skip_pct);
    }
  }
  else
  {
    /* Previous per-rep traversal (fallback without culling or on allocation failure) */
    for (int ir = 0; ir < MassMapSegmentReplicationCount; ++ir)
    {
      int rep_id = MassMapSegmentReplications[ir];
      unsigned long long rep_crossings = 0ULL;
      double shift[3];
      if (RepShift)
      {
        shift[0] = RepShift[3 * ir + 0];
        shift[1] = RepShift[3 * ir + 1];
        shift[2] = RepShift[3 * ir + 2];
      }
      else
      {
        int ii = plc.repls[rep_id].i;
        int jj = plc.repls[rep_id].j;
        int kk = plc.repls[rep_id].k;
        shift[0] = ii * Bx;
        shift[1] = jj * By;
        shift[2] = kk * Bz;
      }
      unsigned long long local_total = 0ULL, local_exceed = 0ULL;
      double local_max = 0.0;
      /* Buffer for fallback path: treat the whole local tile as one block (variable pad only) */
      double thresh_fallback = aabb_pad_phys;
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(static) reduction(+ : rep_crossings, local_total, local_exceed) reduction(max : local_max)
#endif
      for (int ig = gx_start; ig < gx_start + nx_local; ++ig)
      {
        for (int jg = gy_start; jg < gy_start + ny_local; ++jg)
        {
          for (int kg = gz_start; kg < gz_start + nz_local; ++kg)
          {
            pos_data obj;
            mass_maps_set_point_from_products_global(ig, jg, kg, Fseg, &obj);
            if (obj.M == 0)
              continue; /* safety */
            double pos_prev[3], pos_curr[3];
            mass_maps_compute_prev_curr_positions(&obj, pos_prev, pos_curr, lpt_order);
            if (collect_disp_stats && ir == 0)
            {
              double qx = (ig + SHIFT) * params.InterPartDist;
              double qy = (jg + SHIFT) * params.InterPartDist;
              double qz = (kg + SHIFT) * params.InterPartDist;
              double dx = pos_curr[0] - qx;
              double dy = pos_curr[1] - qy;
              double dz = pos_curr[2] - qz;
              double d = sqrt(dx * dx + dy * dy + dz * dz);
              local_total++;
              if (d > thresh_fallback)
                local_exceed++;
              if (d > local_max)
                local_max = d;
            }
            double Pprev[3] = {pos_prev[0] + shift[0], pos_prev[1] + shift[1], pos_prev[2] + shift[2]};
            double Pcurr[3] = {pos_curr[0] + shift[0], pos_curr[1] + shift[1], pos_curr[2] + shift[2]};
            double entry_pos[3];
            double chi_cross;
            if (!mass_maps_crossing_from_shifted_positions(Pprev, Pcurr, c_phys, chi_prev_seg, chi_curr_seg, NULL, entry_pos, &chi_cross))
              continue;
            /* Pixel-based angular selection */
            long ipix_sel;
            if (!mass_maps_compute_pixel_from_pos(entry_pos, MassMapNSIDE_current, &ipix_sel))
              continue;
            {
              int prefix_accept = ((unsigned long)ipix_sel < (unsigned long)MassMapNCAP);
#if MASS_MAPS_CAP_DIAG
              if (MassMapCapDiagActive && params.PLCAperture < 180.0)
              {
                int geom_accept = mass_maps_entry_inside_aperture(entry_pos);
                double ang = mass_maps_angle_deg_from_pos(entry_pos);
                if (geom_accept && !prefix_accept)
                {
#ifdef _OPENMP
#pragma omp atomic update
#endif
                  MassMapCapPrefixMiss++;
#ifdef _OPENMP
#pragma omp critical
#endif
                  {
                    if (ang > MassMapCapMissAngleMax)
                      MassMapCapMissAngleMax = ang;
                  }
                }
                if (prefix_accept)
                {
#ifdef _OPENMP
#pragma omp critical
#endif
                  {
                    if (ang > MassMapCapAcceptAngleMax)
                      MassMapCapAcceptAngleMax = ang;
                  }
                }
              }
#endif
              if (params.PLCAperture < 180.0 && !prefix_accept)
                continue;
            }
#ifdef MASS_MAPS_FILTER_UNCOLLAPSED
            {
              double z_plc = InverseComovingDistance(chi_cross);
#ifdef SNAPSHOT
              double z_acc = mass_maps_get_zacc_global(ig, jg, kg);
#ifdef _OPENMP
#pragma omp atomic update
#endif
              zfilter_considered_local++;
              /* Track min/max and consistency for debug */
#ifdef _OPENMP
#pragma omp critical
#endif
              {
                if (z_plc < zplc_min_local)
                  zplc_min_local = z_plc;
                if (z_plc > zplc_max_local)
                  zplc_max_local = z_plc;
                if (z_acc < zacc_min_local)
                  zacc_min_local = z_acc;
                if (z_acc > zacc_max_local)
                  zacc_max_local = z_acc;
              }
              {
                /* Check that z_plc lies within this segment bounds (tolerance at edges) */
                double zmin_seg = (z_prev < z_curr) ? z_prev : z_curr;
                double zmax_seg = (z_prev > z_curr) ? z_prev : z_curr;
                double tol = 1e-3;
                if (z_plc < zmin_seg - tol || z_plc > zmax_seg + tol)
                {
#ifdef _OPENMP
#pragma omp atomic update
#endif
                  zplc_out_of_segment_local++;
                }
              }
              if (!(z_plc > z_acc))
              {
#ifdef _OPENMP
#pragma omp atomic update
#endif
                zfilter_excluded_local++;
                continue;
              }
#else
#error "MASS_MAPS_FILTER_UNCOLLAPSED requires SNAPSHOT"
#endif
            }
#endif
            ++rep_crossings;
            int s_map = segment_index - 1;
            if (s_map >= 0 && s_map < NMassSheets && MassMapNSIDE_current > 0)
            {
              if ((unsigned long)ipix_sel < (unsigned long)MassMapNPIX)
              {
                double *map = mass_maps_segment_ptr(s_map);
#ifdef _OPENMP
#pragma omp atomic update
#endif
                map[ipix_sel] += 1.0;
              }
            }
          }
        }
      }
      total_crossings_all_reps += rep_crossings;
      if (rep_crossings_local)
        rep_crossings_local[ir] = rep_crossings;
      if (collect_disp_stats && ir == 0)
      {
        disp_total_local += local_total;
        disp_exceed_local += local_exceed;
        if (local_max > disp_max_local)
          disp_max_local = local_max;
      }
      if (!ThisTask && internal.verbose_level >= VDBG)
      {
        printf("[%s] MASS_MAPS(products): segment %d rep %d local entries %llu\n", fdate(), segment_index, rep_id, rep_crossings);
      }
    }
  }

  if (RepShift)
    free(RepShift);
  /* Free culling structures */
  if (edge_x)
    free(edge_x);
  if (edge_y)
    free(edge_y);
  if (edge_z)
    free(edge_z);

  /* Reduce per-rep and total counts across tasks */
  if (MassMapSegmentReplicationCount > 0)
  {
    /* Per-rep reductions (only if allocations succeeded) */
    if (rep_crossings_local)
    {
      MPI_Reduce(rep_crossings_local, rep_crossings_global,
                 MassMapSegmentReplicationCount, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    /* Total across reps */
    unsigned long long local_total = total_crossings_all_reps;
    MPI_Reduce(&local_total, &total_crossings_all_reps_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  }

  if (!ThisTask && internal.verbose_level >= VDBG)
  {
    printf("[%s] MASS_MAPS(products): segment %d total local entries across %d reps: %llu\n",
           fdate(), segment_index, MassMapSegmentReplicationCount, total_crossings_all_reps);
    printf("[%s] MASS_MAPS(products): segment %d total global entries across %d reps: %llu\n",
           fdate(), segment_index, MassMapSegmentReplicationCount, total_crossings_all_reps_global);
  }

#ifdef MASS_MAPS_FILTER_UNCOLLAPSED
  /* Reduce and report ZPLC>ZACC filter diagnostics */
  MPI_Reduce(&zfilter_considered_local, &zfilter_considered_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&zfilter_excluded_local, &zfilter_excluded_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  /* Reduce extra diagnostics */
  MPI_Reduce(&zplc_min_local, &zplc_min_global, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&zplc_max_local, &zplc_max_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&zacc_min_local, &zacc_min_global, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&zacc_max_local, &zacc_max_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&zplc_out_of_segment_local, &zplc_out_of_segment_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  if (!ThisTask && internal.verbose_level >= VDIAG)
  {
    unsigned long long zf_included = (zfilter_considered_global >= zfilter_excluded_global)
                                         ? (zfilter_considered_global - zfilter_excluded_global)
                                         : 0ULL;
    double frac_excl = (zfilter_considered_global > 0ULL)
                           ? ((double)zfilter_excluded_global / (double)zfilter_considered_global)
                           : 0.0;
    /* Stash diagnostics for FITS header metadata */
    g_massmaps_filter_diag.considered = (long long)zfilter_considered_global;
    g_massmaps_filter_diag.excluded = (long long)zfilter_excluded_global;
    g_massmaps_filter_diag.included = (long long)zf_included;
    g_massmaps_filter_diag.frac_excluded = frac_excl;
    g_massmaps_filter_diag.zplc_min = zplc_min_global;
    g_massmaps_filter_diag.zplc_max = zplc_max_global;
    g_massmaps_filter_diag.zacc_min = zacc_min_global;
    g_massmaps_filter_diag.zacc_max = zacc_max_global;
    g_massmaps_filter_diag.zplc_out_of_segment = (long long)zplc_out_of_segment_global;
    g_massmaps_filter_diag.segment_index = segment_index - 1; /* map index s_map below */
    g_massmaps_filter_diag.present = 1;
    printf("[%s] MASS_MAPS(filter): segment %d considered=%llu excluded=%llu (%.2f%%) included=%llu (%.2f%%)\n",
           fdate(), segment_index,
           zfilter_considered_global,
           zfilter_excluded_global, 100.0 * frac_excl,
           zf_included, 100.0 * (1.0 - frac_excl));
    if (internal.verbose_level >= VDBG)
      printf("[%s] MASS_MAPS(filter-diag): seg %d z_plc[min,max]=[%.6g, %.6g]  z_acc[min,max]=[%.6g, %.6g]  z_plc_out_of_seg=%llu\n",
             fdate(), segment_index,
             zplc_min_global, zplc_max_global,
             zacc_min_global, zacc_max_global,
             zplc_out_of_segment_global);
  }
#endif

  /* Report last-segment displacement vs buffer stats */
  if (collect_disp_stats)
  {
    unsigned long long disp_total_global = 0ULL, disp_exceed_global = 0ULL;
    double disp_max_global = 0.0;
    MPI_Reduce(&disp_total_local, &disp_total_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&disp_exceed_local, &disp_exceed_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&disp_max_local, &disp_max_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (!ThisTask)
    {
      double frac = (disp_total_global > 0ULL) ? ((double)disp_exceed_global / (double)disp_total_global) : 0.0;
      if (internal.verbose_level >= VDBG)
        printf("[%s] MASS_MAPS buffer-diag: last segment |x(q,z)-q| vs buffer pad=%.6g Mpc (%.6g cells): total=%llu exceed=%llu (%.3f%%) max=%.6g\n",
               fdate(), aabb_pad_phys, (aabb_pad_phys / params.InterPartDist),
               disp_total_global, disp_exceed_global, 100.0 * frac, disp_max_global);
      if (disp_total_global > 0ULL && frac >= warn_frac)
      {
        fprintf(stderr,
                "[%s] MASS_MAPS WARNING: %.2f%% of particles exceeded the buffer pad=%.6g Mpc (%.6g cells). Consider increasing MASS_MAPS_VARPAD_FACTOR (warn threshold=%.2f%%).\n",
                fdate(), 100.0 * frac, aabb_pad_phys, (aabb_pad_phys / params.InterPartDist), 100.0 * warn_frac);
      }
    }
  }

  /* Reduce the segment map to root (task 0) */
  if (MassMapNSIDE_current > 0)
  {
    int s_map = segment_index - 1;
    if (s_map >= 0 && s_map < NMassSheets)
    {
      double *map_local = mass_maps_segment_ptr(s_map);
      if (!ThisTask)
      {
        if (!MassMapReduceBuffer)
          MassMapReduceBuffer = (double *)malloc(sizeof(double) * (size_t)MassMapNPIX);
      }
      MPI_Reduce(map_local,
                 (!ThisTask ? MassMapReduceBuffer : NULL),
                 (int)MassMapNPIX,
                 MPI_DOUBLE,
                 MPI_SUM,
                 0,
                 MPI_COMM_WORLD);
      if (!ThisTask && MassMapReduceBuffer)
      {
        /* Copy reduced result back into the segment map on root */
        memcpy(map_local, MassMapReduceBuffer, sizeof(double) * (size_t)MassMapNPIX);
        /* Write to FITS */
        mass_maps_write_segment_map(s_map);
        /* Clear saved filter diagnostics once written for this segment */
        g_massmaps_filter_diag.present = 0;
        if (internal.verbose_level >= VDBG)
        {
          double sum, vmin, vmax;
          long nnz;
          mass_maps_map_stats(map_local, MassMapNPIX, &sum, &vmin, &vmax, &nnz);
          long axis_pix = mass_maps_axis_pixel_ring(MassMapNSIDE_current);
          double axis_val = (axis_pix >= 0 && axis_pix < MassMapNPIX) ? map_local[axis_pix] : -1.0;
          printf("[%s] MASS_MAPS(diag): seg %d map stats: sum=%.0f nnz=%ld min=%.0f max=%.0f axis_pix=%ld axis_val=%.0f\n",
                 fdate(), s_map, sum, nnz, vmin, vmax, axis_pix, axis_val);
        }
      }
    }
  }

  if (rep_crossings_local)
    free(rep_crossings_local);
  if (rep_crossings_global)
    free(rep_crossings_global);
  return; /* products path handled this segment */
}

#endif /* MASS_MAPS */
