/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "pinocchio.h"
#include "def_splines.h"
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

/* Wrapper around gsl_spline_init that checks strict monotonicity of x
   and prints a diagnostic message identifying the offending spline. */
int checked_spline_init(gsl_spline *spline, const double xa[], const double ya[],
                        size_t size, const char *label)
{
  if (size < 2)
  {
    printf("ERROR on task %d: spline \"%s\" requires at least 2 points, got %zu\n",
           ThisTask, label, size);
    fflush(stdout);
    return 1;
  }
  for (size_t i = 1; i < size; i++)
    if (xa[i] <= xa[i - 1])
    {
      printf("ERROR on task %d: x values not strictly increasing in spline \"%s\" "
             "at index %zu: x[%zu]=%.15g  x[%zu]=%.15g\n",
             ThisTask, label, i, i - 1, xa[i - 1], i, xa[i]);
      fflush(stdout);
      return 1;
    }
  int gsl_ret = gsl_spline_init(spline, xa, ya, size);
  if (gsl_ret)
  {
    printf("ERROR on task %d: gsl_spline_init failed for spline \"%s\" (GSL error %d)\n",
           ThisTask, label, gsl_ret);
    fflush(stdout);
  }
  return gsl_ret;
}

/* Wrapper around gsl_spline2d_init that checks strict monotonicity of
   both x and y arrays and prints diagnostics identifying the offending spline. */
int checked_spline2d_init(gsl_spline2d *spline, const double xa[], const double ya[],
                          const double za[], size_t xsize, size_t ysize, const char *label)
{
  if (xsize < 2 || ysize < 2)
  {
    printf("ERROR on task %d: 2D spline \"%s\" requires at least 2 points per axis, "
           "got xsize=%zu ysize=%zu\n",
           ThisTask, label, xsize, ysize);
    fflush(stdout);
    return 1;
  }
  for (size_t i = 1; i < xsize; i++)
    if (xa[i] <= xa[i - 1])
    {
      printf("ERROR on task %d: x values not strictly increasing in 2D spline \"%s\" "
             "at index %zu: x[%zu]=%.15g  x[%zu]=%.15g\n",
             ThisTask, label, i, i - 1, xa[i - 1], i, xa[i]);
      fflush(stdout);
      return 1;
    }
  for (size_t i = 1; i < ysize; i++)
    if (ya[i] <= ya[i - 1])
    {
      printf("ERROR on task %d: y values not strictly increasing in 2D spline \"%s\" "
             "at index %zu: y[%zu]=%.15g  y[%zu]=%.15g\n",
             ThisTask, label, i, i - 1, ya[i - 1], i, ya[i]);
      fflush(stdout);
      return 1;
    }
  int gsl_ret = gsl_spline2d_init(spline, xa, ya, za, xsize, ysize);
  if (gsl_ret)
  {
    printf("ERROR on task %d: gsl_spline2d_init failed for 2D spline \"%s\" (GSL error %d)\n",
           ThisTask, label, gsl_ret);
    fflush(stdout);
  }
  return gsl_ret;
}

/* Lazy-built inverse comoving-distance spline: chi [Mpc] -> log10(a) */
static gsl_spline *SPLINE_INVCOMVDIST = 0x0;
static gsl_interp_accel *ACCEL_INVCOMVDIST = 0x0;
static int INVCOMVDIST_READY = 0;

#define NVAR (1 + NkBINS * 8)
#define NBB 10
#ifdef NORADIATION
#define OMEGARAD_H2 ((double)0.0)
#else
#define OMEGARAD_H2 ((double)4.2e-5)
#endif
#define UnitLength_in_cm ((double)3.085678e24)
#define HUBBLETIME_GYR ((double)3.085678e24 / (double)1.e7 / (double)3.1558150e16)
#define DELTA_C ((double)1.686)
#define SHAPE_EFST ((double)0.21)

// #define FOMEGA_GAMMA 0.554
// TUTTE LE CONDIZIONI SULLE DIRETTIVE DEVONO ESSERE MESSE INSIEME
#if defined(FOMEGA_GAMMA) && defined(SCALE_DEPENDENT)
#error Do not use FOMEGA_GAMMA with SCALE_DEPENDENT
#endif

static int Today;
static int WhichSpectrum, NPowerTable = 0, NtabEoS = 0;
static double PkNorm, MatterDensity, OmegaK, OmegaRad;

/* declaration of gsl quantities */

#ifdef SCALE_DEPENDENT
static double kmin, kmax;
#endif

int system_of_ODEs(double, const double[], double *, void *);
int system_of_ODEs_small(double, const double[], double *, void *);
int read_TabulatedEoS(void);
#ifdef READ_HUBBLE_TABLE
int read_TabulatedHubble(void);
#endif
int initialize_PowerSpectrum(void);
int normalize_PowerSpectrum(void);
int read_Pk_from_file(void);
double IntegrandForEoS(double, void *);
double DE_EquationOfState(double);
double IntegrandComovingDistance(double, void *);
double ComputeMassVariance(double);
double ComputeDisplVariance(double);
#ifdef READ_PK_TABLE
int read_Pk_table_from_CAMB(double *, double *, double *, double *, double *, double *, double *, double *, double *);
#endif

/**************************/
/* INITIALIZATION SECTION */
/**************************/

int initialize_cosmology()
{

  /*
    Computes the following functions:
    Scale factor, growth first, second and third-order LPT, cosmic time,
    comoving and diameter distance on a grid of values to be interpolated
  */

  double ode_param;
  double y[NVAR], x1, x2, hh, norm, result, error, SqrtOmegaK, R0, k;
  int status = GSL_SUCCESS, i, j;
#ifdef SCALE_DEPENDENT
  int ik;
#endif
  char filename[LBLENGTH];
  FILE *fd;
  double log_amin = -4., dloga = -log_amin / (double)(NBINS - NBB);

  double *scalef, *cosmtime, *grow1, *grow2, *IntEoS, *comvdist, *diamdist,
      *fomega1, *fomega2, *grow31, *grow32, *fomega31, *fomega32;

  gsl_function Function_cosmo;

  /* ODE objects: hoisted to function scope so the fail label can free them */
  gsl_odeiv2_step *ode_s = NULL;
  gsl_odeiv2_control *ode_c = NULL;
  gsl_odeiv2_evolve *ode_e = NULL;

#ifdef MOD_GRAV_FR
  H_over_c = 100. / SPEEDOFLIGHT;
#endif
  OmegaRad = OMEGARAD_H2 / params.Hubble100 / params.Hubble100;
  OmegaK = 1.0 - params.Omega0 - params.OmegaLambda - OmegaRad;
  SqrtOmegaK = sqrt(fabs(OmegaK));
  MatterDensity = 2.775499745e11 * params.Hubble100 * params.Hubble100 * params.Omega0;
  if (params.DEw0 == -1 && params.DEwa == 0 && !strcmp(params.TabulatedEoSfile, "no"))
    params.simpleLambda = 1;
  else
    params.simpleLambda = 0;

  /* allocation of (most) splines */
  SPLINE = (gsl_spline **)calloc(NSPLINES, sizeof(gsl_spline *));
  for (i = 0; i < NSPLINES - 3; i++)
    if (i != SP_COMVDIST && i != SP_DIAMDIST)
      SPLINE[i] = gsl_spline_alloc(gsl_interp_cspline, NBINS);
    else
      SPLINE[i] = gsl_spline_alloc(gsl_interp_cspline, NBINS - NBB);
  ACCEL = (gsl_interp_accel **)calloc(NSPLINES, sizeof(gsl_interp_accel *));
  for (i = 0; i < NSPLINES; i++)
    ACCEL[i] = gsl_interp_accel_alloc();

  /* if needed, read the tabulated Equation of State of the dark energy
     and initialize its spline, then compute the integrand for the DE EoS */
  if (!params.simpleLambda)
  {
    if (strcmp(params.TabulatedEoSfile, "no"))
    {
      if (read_TabulatedEoS()) /* this allocates and initializes the spline */
        return 1;
    }
    else
      NtabEoS = 0;

    SPLINE[SP_INTEOS] = gsl_spline_alloc(gsl_interp_cspline, NBINS);

    scalef = (double *)malloc(NBINS * sizeof(double));
    IntEoS = (double *)malloc(NBINS * sizeof(double));
    Function_cosmo.function = &IntegrandForEoS;
    for (i = 0; i < NBINS; i++)
    {
      x2 = pow(10., log_amin + (i + 1) * dloga);
      gsl_integration_qags(&Function_cosmo, x2, 1.0, 0.0, TOLERANCE, NWINT, workspace, &result, &error);
      scalef[i] = log_amin + (i + 1) * dloga;
      IntEoS[i] = result;
    }

    if (checked_spline_init(SPLINE[SP_INTEOS], scalef, IntEoS, NBINS, "SP_INTEOS"))
    {
      free(IntEoS);
      free(scalef);
      return 1;
    }

    free(IntEoS);
    free(scalef);
  }

#ifdef READ_HUBBLE_TABLE
  /* If requested, read tabulated E(z)=H(z)/H0 and create its spline over log10(a). */
  if (read_TabulatedHubble())
    return 1;
#endif

#ifdef SCALE_DEPENDENT
  kmin = pow(10., LOGKMIN);
  kmax = pow(10., LOGKMIN + (NkBINS - 1) * DELTALOGK);
#endif

  /* The power spectrum is initialized; in case, it is read from file(s) */
  if (initialize_PowerSpectrum())
    return 1;

  /* allocation of vectors for interpolation */
  scalef = (double *)malloc(NBINS * sizeof(double));
  cosmtime = (double *)malloc(NBINS * sizeof(double));
  comvdist = (double *)malloc((NBINS - NBB) * sizeof(double));
  diamdist = (double *)malloc((NBINS - NBB) * sizeof(double));
  grow1 = (double *)malloc(NBINS * NkBINS * sizeof(double));
  grow2 = (double *)malloc(NBINS * NkBINS * sizeof(double));
  grow31 = (double *)malloc(NBINS * NkBINS * sizeof(double));
  grow32 = (double *)malloc(NBINS * NkBINS * sizeof(double));
  fomega1 = (double *)malloc(NBINS * NkBINS * sizeof(double));
  fomega2 = (double *)malloc(NBINS * NkBINS * sizeof(double));
  fomega31 = (double *)malloc(NBINS * NkBINS * sizeof(double));
  fomega32 = (double *)malloc(NBINS * NkBINS * sizeof(double));

#ifdef READ_PK_TABLE
  if (WhichSpectrum != 5)
#endif
  {
    /* Runge-Kutta integration of cosmic time and growth rate */
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;
    ode_s = gsl_odeiv2_step_alloc(T, NVAR);
    ode_c = gsl_odeiv2_control_standard_new(1.0e-8, 1.0e-8, 1.0, 1.0);
    ode_e = gsl_odeiv2_evolve_alloc(NVAR);
    gsl_odeiv2_system ode_sys = {system_of_ODEs, jac, NVAR, (void *)&ode_param};

    /* ICs for the runge-kutta integration */
    x1 = pow(10., log_amin - 2.);  /*  initial value of scale factor a */
    y[0] = 2. / 3. * pow(x1, 1.5); /*  initial value of t(a)*Hubble0   */

    /* this is valid for scale-independent and scale-dependent functions */
    for (j = 0; j < NkBINS; j++)
    {
      y[1 + j * 8] = 1.0;                      /*  initial value of dD1/da   */
      y[2 + j * 8] = x1;                       /*  initial value of D1(a)    */
      y[3 + j * 8] = -6. / 7. * x1;            /*  initial value of dD2/da   */
      y[4 + j * 8] = -3. / 7. * x1 * x1;       /*  initial value of D2(a)    */
      y[5 + j * 8] = -x1 * x1;                 /*  initial value of dD^3a/da */
      y[6 + j * 8] = -1. / 3. * x1 * x1 * x1;  /*  initial value of D^3a     */
      y[7 + j * 8] = 10. / 7. * x1 * x1;       /*  initial value of dD^3b/da */
      y[8 + j * 8] = 10. / 21. * x1 * x1 * x1; /*  initial value of D^3b     */
    }

    hh = x1 / 10.; /*  initial guess of time-step */

    /* this function will be integrated within the loop */
    Function_cosmo.function = &IntegrandComovingDistance;

    /***********************************************/
    /* ODE integration of time-dependent functions */
    /***********************************************/
    for (i = 0, Today = 0; i < NBINS; i++)
    {
      x2 = pow(10., log_amin + i * dloga);
      if (fabs(log_amin + i * dloga) < dloga / 10.)
        x2 = 1.0;

      /* integration of ODE system */
      while (x1 < x2 && status == GSL_SUCCESS)
      {
        status = gsl_odeiv2_evolve_apply(ode_e, ode_c, ode_s, &ode_sys, &x1, x2, &hh, y);
        if (status != GSL_SUCCESS)
        {
          printf("ERROR on task %d: integration of cosmological quantities failed\n", ThisTask);
          fflush(stdout);
          goto fail;
        }
      }

      scalef[i] = x2;
      cosmtime[i] = log10(y[0] * HUBBLETIME_GYR / params.Hubble100);
      if (!Today && x2 >= 1.)
        Today = i;

      for (j = 0; j < NkBINS; j++) /* First-order growth rate */
        grow1[i + j * NBINS] = y[2 + j * 8];
      for (j = 0; j < NkBINS; j++) /* Second-order growth rate */
        grow2[i + j * NBINS] = -y[4 + j * 8];
      for (j = 0; j < NkBINS; j++) /* Third-order first growth rate */
        grow31[i + j * NBINS] = -y[6 + j * 8] / 3.;
      for (j = 0; j < NkBINS; j++) /* Third-order second growth rate */
        grow32[i + j * NBINS] = y[8 + j * 8] / 4.;

      for (j = 0; j < NkBINS; j++) /* First-order f(Omega) */
        fomega1[i + j * NBINS] = x2 * y[1 + j * 8] / y[2 + j * 8];
      for (j = 0; j < NkBINS; j++) /* Second-order f(Omega) */
        fomega2[i + j * NBINS] = x2 * y[3 + j * 8] / y[4 + j * 8];
      for (j = 0; j < NkBINS; j++) /* Third-order first f(Omega) */
        fomega31[i + j * NBINS] = x2 * y[5 + j * 8] / y[6 + j * 8];
      for (j = 0; j < NkBINS; j++) /* Third-order second f(Omega) */
        fomega32[i + j * NBINS] = x2 * y[7 + j * 8] / y[8 + j * 8];

      /* calculation of comoving distance (Mpc) for a generic cosmology (i.e. flat, open or closed)
         and generic equation of state for the DE component  */
      if (i < NBINS - NBB)
      {
        gsl_integration_qags(&Function_cosmo, 0.0, 1. / x2 - 1., 0.0, TOLERANCE, NWINT, workspace, &result, &error);
        comvdist[i] = SPEEDOFLIGHT * result;
        if (fabs(OmegaK) < 1.e-4)
          diamdist[i] = x2 * comvdist[i];
        else if (OmegaK < 0)
        {
          R0 = SPEEDOFLIGHT / params.Hubble100 / 100. / SqrtOmegaK;
          diamdist[i] = x2 * R0 * sin(comvdist[i] / R0);
        }
        else
        {
          R0 = SPEEDOFLIGHT / params.Hubble100 / 100. / SqrtOmegaK;
          diamdist[i] = x2 * R0 * sinh(comvdist[i] / R0);
        }
      }

      /* closing the loop on integrations */
      x1 = x2;
    }

    /* normalization of the first- and second- order growth rate */
    /* this is valid for LambdaCDM; for scale-dependent growth due to modified gravity,
       the power spectrum is given as the LambdaCDM P(k) extrapolated at z=0, but this
       is valid only at high redshift; this normalization is still correct
       when the k=0 growth rate (identical to LambdaCDM) is used at all scales */
    norm = grow1[Today];
    for (i = 0; i < NBINS * NkBINS; i++)
    {
      grow1[i] /= norm;
      grow2[i] /= norm * norm;
      grow31[i] /= norm * norm * norm;
      grow32[i] /= norm * norm * norm;
    }

    gsl_odeiv2_evolve_free(ode_e);
    gsl_odeiv2_control_free(ode_c);
    gsl_odeiv2_step_free(ode_s);
    ode_e = NULL;
    ode_c = NULL;
    ode_s = NULL;
  }
#ifdef READ_PK_TABLE
  else
  {
    if (!ThisTask)
      printf("Only the cosmic time is integrated, the growth rate is read from CAMB files\n");

    /* in this case the integration is limited only to the cosmic time */
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;
    ode_s = gsl_odeiv2_step_alloc(T, 1);
    ode_c = gsl_odeiv2_control_standard_new(1.0e-8, 1.0e-8, 1.0, 1.0);
    ode_e = gsl_odeiv2_evolve_alloc(1);
    gsl_odeiv2_system ode_sys = {system_of_ODEs_small, jac, 1, (void *)&ode_param};

    /* ICs for the runge-kutta integration of the cosmic time only */
    x1 = pow(10., log_amin - 2.);  /*  initial value of scale factor a */
    y[0] = 2. / 3. * pow(x1, 1.5); /*  initial value of t(a)*Hubble0   */
    hh = x1 / 10.;                 /*  initial guess of time-step */

    /* this function will be integrated within the loop */
    Function_cosmo.function = &IntegrandComovingDistance;

    /***********************************************/
    /* ODE integration of time-dependent functions */
    /***********************************************/
    for (i = 0, Today = 0; i < NBINS; i++)
    {
      x2 = pow(10., log_amin + i * dloga);
      if (fabs(log_amin + i * dloga) < dloga / 10.)
        x2 = 1.0;

      /* integration of ODE system */
      while (x1 < x2 && status == GSL_SUCCESS)
      {
        status = gsl_odeiv2_evolve_apply(ode_e, ode_c, ode_s, &ode_sys, &x1, x2, &hh, y);
        if (status != GSL_SUCCESS)
        {
          printf("ERROR on task %d: integration of cosmological quantities failed\n", ThisTask);
          fflush(stdout);
          goto fail;
        }
      }

      scalef[i] = x2;
      cosmtime[i] = log10(y[0] * HUBBLETIME_GYR / params.Hubble100);
      if (!Today && x2 >= 1.)
        Today = i;

      /* calculation of comoving distance (Mpc) for a generic cosmology (i.e. flat, open or closed)
         and generic equation of state for the DE component  */
      if (i < NBINS - NBB)
      {
        gsl_integration_qags(&Function_cosmo, 0.0, 1. / x2 - 1., 0.0, TOLERANCE, NWINT, workspace, &result, &error);
        comvdist[i] = SPEEDOFLIGHT * result;
        if (fabs(OmegaK) < 1.e-4)
          diamdist[i] = x2 * comvdist[i];
        else if (OmegaK < 0)
        {
          R0 = SPEEDOFLIGHT / params.Hubble100 / 100. / SqrtOmegaK;
          diamdist[i] = x2 * R0 * sin(comvdist[i] / R0);
        }
        else
        {
          R0 = SPEEDOFLIGHT / params.Hubble100 / 100. / SqrtOmegaK;
          diamdist[i] = x2 * R0 * sinh(comvdist[i] / R0);
        }
      }

      /* closing the loop on integrations */
      x1 = x2;
    }

    /* the growth rates are set by reading the CAMB power spectra */
    if (read_Pk_table_from_CAMB(scalef, grow1, grow2, grow31, grow32, fomega1, fomega2, fomega31, fomega32))
      goto fail;

    gsl_odeiv2_evolve_free(ode_e);
    gsl_odeiv2_control_free(ode_c);
    gsl_odeiv2_step_free(ode_s);
    ode_e = NULL;
    ode_c = NULL;
    ode_s = NULL;
  }
#endif

  /* these quantities will be interpolated logarithmically */
  for (i = 0; i < NBINS; i++)
    scalef[i] = log10(scalef[i]);
  for (j = 0; j < NkBINS; j++)
    for (i = 0; i < NBINS; i++)
    {
      grow1[i + j * NBINS] = log10(grow1[i + j * NBINS]);
      grow2[i + j * NBINS] = log10(grow2[i + j * NBINS]);
      grow31[i + j * NBINS] = log10(grow31[i + j * NBINS]);
      grow32[i + j * NBINS] = log10(grow32[i + j * NBINS]);
    }

  /* initialization of spline interpolations of time-dependent quantities */
  if (checked_spline_init(SPLINE[SP_TIME], scalef, cosmtime, NBINS, "SP_TIME") ||
      checked_spline_init(SPLINE[SP_INVTIME], cosmtime, scalef, NBINS, "SP_INVTIME") ||
      checked_spline_init(SPLINE[SP_COMVDIST], scalef, comvdist, NBINS - NBB, "SP_COMVDIST") ||
      checked_spline_init(SPLINE[SP_DIAMDIST], scalef, diamdist, NBINS - NBB, "SP_DIAMDIST") ||
      /* inverse grow is always defined on the first growth rate */
      checked_spline_init(SPLINE[SP_INVGROW], grow1, scalef, NBINS, "SP_INVGROW"))
    goto fail;

  for (j = 0; j < NkBINS; j++)
  {
    char lbl[64];
    snprintf(lbl, sizeof(lbl), "SP_GROW1[%d]", j);
    if (checked_spline_init(SPLINE[SP_GROW1 + j], scalef, grow1 + j * NBINS, NBINS, lbl))
      goto fail;
    snprintf(lbl, sizeof(lbl), "SP_GROW2[%d]", j);
    if (checked_spline_init(SPLINE[SP_GROW2 + j], scalef, grow2 + j * NBINS, NBINS, lbl))
      goto fail;
    snprintf(lbl, sizeof(lbl), "SP_GROW31[%d]", j);
    if (checked_spline_init(SPLINE[SP_GROW31 + j], scalef, grow31 + j * NBINS, NBINS, lbl))
      goto fail;
    snprintf(lbl, sizeof(lbl), "SP_GROW32[%d]", j);
    if (checked_spline_init(SPLINE[SP_GROW32 + j], scalef, grow32 + j * NBINS, NBINS, lbl))
      goto fail;

    snprintf(lbl, sizeof(lbl), "SP_FOMEGA1[%d]", j);
    if (checked_spline_init(SPLINE[SP_FOMEGA1 + j], scalef, fomega1 + j * NBINS, NBINS, lbl))
      goto fail;
    snprintf(lbl, sizeof(lbl), "SP_FOMEGA2[%d]", j);
    if (checked_spline_init(SPLINE[SP_FOMEGA2 + j], scalef, fomega2 + j * NBINS, NBINS, lbl))
      goto fail;
    snprintf(lbl, sizeof(lbl), "SP_FOMEGA31[%d]", j);
    if (checked_spline_init(SPLINE[SP_FOMEGA31 + j], scalef, fomega31 + j * NBINS, NBINS, lbl))
      goto fail;
    snprintf(lbl, sizeof(lbl), "SP_FOMEGA32[%d]", j);
    if (checked_spline_init(SPLINE[SP_FOMEGA32 + j], scalef, fomega32 + j * NBINS, NBINS, lbl))
      goto fail;
  }

  /* deallocation of vectors for interpolation */
  free(fomega32);
  free(fomega31);
  free(fomega2);
  free(fomega1);
  free(grow32);
  free(grow31);
  free(grow2);
  free(grow1);
  free(diamdist);
  free(comvdist);
  free(cosmtime);
  free(scalef);

  /* normalization of power spectrum */
  if (normalize_PowerSpectrum())
    return 1;

  /* initialization of mass variance with Gaussian filter */
  WindowFunctionType = 0;
  if (initialize_MassVariance())
    return 1;

  /* write out cosmological quantities */
  if (!ThisTask)
  {

    strcpy(filename, "pinocchio.");
    strcat(filename, params.RunFlag);
    strcat(filename, ".cosmology.out");

    fd = fopen(filename, "w");

    fprintf(fd, "# Cosmological quantities used in PINOCCHIO (h=%f)\n", params.Hubble100);
    fprintf(fd, "# TIME-DEPENDENT QUANTITIES\n");
    fprintf(fd, "# 1: scale factor\n");
    fprintf(fd, "# 2: cosmic time (Gyr)\n");
    fprintf(fd, "# 3: comoving distance (Mpc)\n");
    fprintf(fd, "# 4: diameter distance (Mpc)\n");
    fprintf(fd, "# 5: Omega matter\n");
    fprintf(fd, "# 6: dark energy EOS\n");
    fprintf(fd, "# 7: linear growth rate\n");
    fprintf(fd, "# 8: 2nd-order growth rate\n");
    fprintf(fd, "# 9: first 3rd-order growth rate\n");
    fprintf(fd, "#10: second 3rd-order growth rate\n");
    fprintf(fd, "#11: linear d ln D/d ln a\n");
    fprintf(fd, "#12: 2nd-order d ln D/d ln a\n");
    fprintf(fd, "#13: first 3rd-order d ln D/d ln a\n");
    fprintf(fd, "#14: second 3rd-order d ln D/d ln a\n");
    fprintf(fd, "# SCALE-DEPENDENT QUANTITIES\n");
    fprintf(fd, "#15: smoothing scale R (Mpc)\n");
    fprintf(fd, "#16: Gaussian-filtered mass variance sigma_G^2(R)\n");
    fprintf(fd, "#17: Gaussian-filtered displacement variance\n");
    fprintf(fd, "#18: d Log sigma_G^2 / d Log R\n");
    fprintf(fd, "# POWER SPECTRUM\n");
    fprintf(fd, "#19: k (true Mpc^-1)\n");
    fprintf(fd, "#20: P(k)\n");
    fprintf(fd, "#\n");

    for (i = 0; i < NBINS; i++)
    {
      k = pow(10., -4.0 + (double)i / (double)NBINS * 6.0);
      fprintf(fd, " %12lg %12lg %12lg %12lg %12lg %12lg %12lg %12lg %12lg %12lg %12lg %12lg %12lg %12lg   %12lg %12lg %12lg %12lg   %12lg %12lg\n",
              pow(10., SPLINE[SP_TIME]->x[i]),
              pow(10., SPLINE[SP_TIME]->y[i]),
              (i < NBINS - NBB ? SPLINE[SP_COMVDIST]->y[i] : 0.0),
              (i < NBINS - NBB ? SPLINE[SP_DIAMDIST]->y[i] : 0.0),
              OmegaMatter(1. / pow(10., SPLINE[SP_TIME]->x[i]) - 1.0),
              (!params.simpleLambda ? -1 : (NtabEoS ? SPLINE[SP_EOS]->y[i] : DE_EquationOfState(pow(10., SPLINE[SP_TIME]->x[i])))),
              pow(10., SPLINE[SP_GROW1]->y[i]),
              pow(10., SPLINE[SP_GROW2]->y[i]),
              pow(10., SPLINE[SP_GROW31]->y[i]),
              pow(10., SPLINE[SP_GROW32]->y[i]),
              SPLINE[SP_FOMEGA1]->y[i],
              SPLINE[SP_FOMEGA2]->y[i],
              SPLINE[SP_FOMEGA31]->y[i],
              SPLINE[SP_FOMEGA32]->y[i],
              pow(10., SPLINE[SP_MASSVAR]->x[i]),
              pow(10., SPLINE[SP_MASSVAR]->y[i]),
              pow(10., SPLINE[SP_DISPVAR]->y[i]),
              SPLINE[SP_DVARDR]->y[i],
              k, PowerSpectrum(k));
    }
    fclose(fd);

#ifdef SCALE_DEPENDENT
    /* writes scale-dependent growth rates on a file */
    strcpy(filename, "pinocchio.");
    strcat(filename, params.RunFlag);
    strcat(filename, ".scaledep.out");

    fd = fopen(filename, "w");

    fprintf(fd, "# Scale-dependent growth rates\n");
    fprintf(fd, "# Scales considered: ");
    for (ik = 0; ik < NkBINS; ik++)
    {
#ifdef MOD_GRAV_FR
      /* with modified gravity the first wavenumber is set to zero */
      if (!ik)
        k = 0.0;
      else
#endif
        k = pow(10., LOGKMIN + ik * DELTALOGK);
      if (ik == NkBINS - 1)
        fprintf(fd, "%d) k=%8.5f\n", ik + 1, k);
      else
        fprintf(fd, "%d) k=%8.5f, ", ik + 1, k);
    }

    fprintf(fd, "# 1: scale factor\n");
    fprintf(fd, "# %d-%d: linear growth rate\n", 2, NkBINS + 1);
    fprintf(fd, "# %d-%d: 2nd-order growth rate\n", NkBINS + 2, 2 * NkBINS + 1);
    fprintf(fd, "# %d-%d: first 3rd-order growth rate\n", 2 * NkBINS + 2, 3 * NkBINS + 1);
    fprintf(fd, "# %d-%d: second 3rd-order growth rate\n", 3 * NkBINS + 2, 4 * NkBINS + 1);
    fprintf(fd, "# %d-%d: linear d ln D/d ln a\n", 4 * NkBINS + 2, 5 * NkBINS + 1);
    fprintf(fd, "# %d-%d: 2nd-order d ln D/d ln a\n", 5 * NkBINS + 2, 6 * NkBINS + 1);
    fprintf(fd, "# %d-%d: first 3rd-order d ln D/d ln a\n", 6 * NkBINS + 2, 7 * NkBINS + 1);
    fprintf(fd, "# %d-%d: second 3rd-order d ln D/d ln a\n", 7 * NkBINS + 2, 8 * NkBINS + 1);
    fprintf(fd, "#\n");

    for (i = 0; i < NBINS; i++)
    {
      fprintf(fd, " %12lg", pow(10., SPLINE[SP_TIME]->x[i]));
      for (ik = 0; ik < NkBINS; ik++)
        fprintf(fd, " %12lg", pow(10., SPLINE[SP_GROW1 + ik]->y[i]));
      fprintf(fd, "   ");
      for (ik = 0; ik < NkBINS; ik++)
        fprintf(fd, " %12lg", pow(10., SPLINE[SP_GROW2 + ik]->y[i]));
      fprintf(fd, "   ");
      for (ik = 0; ik < NkBINS; ik++)
        fprintf(fd, " %12lg", pow(10., SPLINE[SP_GROW31 + ik]->y[i]));
      fprintf(fd, "   ");
      for (ik = 0; ik < NkBINS; ik++)
        fprintf(fd, " %12lg", pow(10., SPLINE[SP_GROW32 + ik]->y[i]));
      fprintf(fd, "   ");
      for (ik = 0; ik < NkBINS; ik++)
        fprintf(fd, " %12lg", SPLINE[SP_FOMEGA1 + ik]->y[i]);
      fprintf(fd, "   ");
      for (ik = 0; ik < NkBINS; ik++)
        fprintf(fd, " %12lg", SPLINE[SP_FOMEGA2 + ik]->y[i]);
      fprintf(fd, "   ");
      for (ik = 0; ik < NkBINS; ik++)
        fprintf(fd, " %12lg", SPLINE[SP_FOMEGA31 + ik]->y[i]);
      fprintf(fd, "   ");
      for (ik = 0; ik < NkBINS; ik++)
        fprintf(fd, " %12lg", SPLINE[SP_FOMEGA32 + ik]->y[i]);
      fprintf(fd, "\n");
    }
    fclose(fd);
#endif
  }

#ifdef RECOMPUTE_DISPLACEMENTS
  /* here the segmentation of the fragmentation process is defined */
  /* for now, segmentation is taken from the outputs file */
  ScaleDep.nseg = outputs.n;
  for (int i = 0; i < outputs.n; i++) // servono i growth rate? a che scala?
  {
    ScaleDep.z[i] = outputs.z[i];
    ScaleDep.D[i] = GrowingMode(ScaleDep.z[i], 0.0);
    ScaleDep.D2[i] = GrowingMode_2LPT(ScaleDep.z[i], 0.0);
    ScaleDep.D31[i] = GrowingMode_3LPT_1(ScaleDep.z[i], 0.0);
    ScaleDep.D32[i] = GrowingMode_3LPT_2(ScaleDep.z[i], 0.0);
  }
#else
  /* in case the displacements are not to be recomputed,
     there is only one segment that gets to the end*/
  ScaleDep.nseg = 1;
  ScaleDep.z[0] = outputs.zlast;
  ScaleDep.D[0] = GrowingMode(ScaleDep.z[0], 0.0);
  ScaleDep.D2[0] = GrowingMode_2LPT(ScaleDep.z[0], 0.0);
  ScaleDep.D31[0] = GrowingMode_3LPT_1(ScaleDep.z[0], 0.0);
  ScaleDep.D32[0] = GrowingMode_3LPT_2(ScaleDep.z[0], 0.0);
#endif

  return 0;

fail:
  if (ode_e)
    gsl_odeiv2_evolve_free(ode_e);
  if (ode_c)
    gsl_odeiv2_control_free(ode_c);
  if (ode_s)
    gsl_odeiv2_step_free(ode_s);
  free(fomega32);
  free(fomega31);
  free(fomega2);
  free(fomega1);
  free(grow32);
  free(grow31);
  free(grow2);
  free(grow1);
  free(diamdist);
  free(comvdist);
  free(cosmtime);
  free(scalef);
  return 1;
}

#ifdef MOD_GRAV_FR
/* scale-dependent functions for 2LPT term */

double mu(double a, double k)
{
  double B1, B2, emme;
  B1 = params.Omega0 / pow(a, 3.) + 4. * params.OmegaLambda;
  B2 = params.Omega0 + 4. * params.OmegaLambda;
  emme = 0.5 * H_over_c * H_over_c * pow(B1, 3.) / (B2 * B2 * FR0);
  return 1. + k * k / 3. / (k * k + a * a * emme);
}

#endif

/* Helpers to compute E(a)^2 and its derivative from H(z) uniformly */
static inline double a_to_z(double a) { return 1.0 / a - 1.0; }
static inline double Ez2_from_a(double a)
{
  double Ezv = Ez(a_to_z(a));
  return Ezv * Ezv;
}
static inline double dlnE2_da(double a)
{
  /* If we have an external E(z)=H(z)/H0 table, use the spline derivative
     of log10 E(a). The derivative is the same as for log10 H(a) up to an
     additive constant in the logarithm. */
#ifdef READ_HUBBLE_TABLE
  if (SPLINE[SP_EXT_HUBBLE])
  {
    double xlna = log10(a);
    /* my_spline_eval_deriv returns d/dx of the spline y(x); here y=log10 E, x=log10 a */
    double dlogE_dlogA = my_spline_eval_deriv(SPLINE[SP_EXT_HUBBLE], xlna, ACCEL[SP_EXT_HUBBLE]);
    /* ln(E^2) = 2 ln E => d/da ln(E^2) = 2 d/da ln E
       and d/da ln E = (1/a) d/d(log a) ln E; since log10, there is a constant factor that cancels
       when converting derivative of log10 to natural log ratio. Using (2/a) * d log10 E / d log10 a suffices. */
    return (2.0 / a) * dlogE_dlogA;
  }
#endif

  /* Otherwise, use analytic expression consistent with Hubble() when not tabulated */
  const double a2 = a * a;
  const double a3 = a2 * a;
  const double a4 = a2 * a2;
  const double a5 = a4 * a;

  double E2 = params.Omega0 / a3 + OmegaK / a2 + OmegaRad / a4;
  double dE2_da = -3.0 * params.Omega0 / a4 - 2.0 * OmegaK / a3 - 4.0 * OmegaRad / a5;

  if (params.simpleLambda)
  {
    E2 += params.OmegaLambda;
    /* derivative of constant term is zero */
  }
  else
  {
    /* Dark-energy term: params.OmegaLambda a^{-3} exp(3 ∫ w(a) d ln a) */
    double de_eos_int = my_spline_eval(SPLINE[SP_INTEOS], log10(a), ACCEL[SP_INTEOS]);
    double w = DE_EquationOfState(a);
    double fac = params.OmegaLambda * exp(3.0 * de_eos_int);
    E2 += fac / a3;
    dE2_da += -3.0 * (1.0 + w) * fac / a4;
  }

  return dE2_da / E2;
}

int system_of_ODEs(double x, const double y[], double *dydx, void *param)
{

  double a1, b1;

  /*
    0 : cosmic time t
    then loop over k bins:
    1 : dD(1) / da
    2 : D(1)
    3 : dD(2) / da
    4 : D(2)
    5 : dD(3a) / da
    6 : D(3a)
    7 : dD(3b) / da
    8 : D(3b)
  */

  /* Build E(a)^2 and its derivative uniformly via H(z) */
  double E2 = Ez2_from_a(x);
  double dlnE2 = dlnE2_da(x);

  /* coefficients in the growth rate equations */
  a1 = -(3. / x + 0.5 * dlnE2);
  b1 = 1.5 * params.Omega0 / (E2 * pow(x, 5.0));

  /* cosmic time */
  dydx[0] = 1.0 / x / sqrt(E2);

#if !defined(MOD_GRAV_FR) && !defined(FOMEGA_GAMMA)

  /* this is the standard integration */
  for (int j = 0; j < NkBINS; j++)
  {
    dydx[1 + j * 8] = a1 * y[1 + j * 8] + b1 * y[2 + j * 8];                                                                                                /* d^2 D1/ dt^2 */
    dydx[2 + j * 8] = y[1 + j * 8];                                                                                                                         /* d D1/ dt */
    dydx[3 + j * 8] = a1 * y[3 + j * 8] + b1 * y[4 + j * 8] - b1 * y[2 + j * 8] * y[2 + j * 8];                                                             /* d^2 D2/ dt^2 */
    dydx[4 + j * 8] = y[3 + j * 8];                                                                                                                         /* d D2/ dt */
    dydx[5 + j * 8] = a1 * y[5 + j * 8] + b1 * y[6 + j * 8] - 2. * b1 * y[2 + j * 8] * y[2 + j * 8] * y[2 + j * 8];                                         /* d^2 D31/ dt^2 */
    dydx[6 + j * 8] = y[5 + j * 8];                                                                                                                         /* d D31/ dt */
    dydx[7 + j * 8] = a1 * y[7 + j * 8] + b1 * y[8 + j * 8] - 2. * b1 * y[2 + j * 8] * y[4 + j * 8] + 2. * b1 * y[2 + j * 8] * y[2 + j * 8] * y[2 + j * 8]; /* d^2 D32/ dt^2 */
    dydx[8 + j * 8] = y[7 + j * 8];                                                                                                                         /* d D32/ dt */
  }
  return GSL_SUCCESS;
#endif

#ifdef FOMEGA_GAMMA

  /* this is a toy model: D(t) is the one obtained by forcing
     the gamma RSD parameter to a certain value.
     Higher orders are obtained with Bouchet's fits */
  /* E2 from the uniform path was computed above; reuse it here. */
  dydx[1] = 0.;
  dydx[2] = (pow(params.Omega0 / pow(x, 3.0) / E2, FOMEGA_GAMMA)) / x * y[1];

  for (i = 3; i < 9; i++)
    dydx[i] = 0.;

  return GSL_SUCCESS;
#endif

#ifdef MOD_GRAV_FR

  /* integration of growth rates in f(R) */
  /* equations for D2(k,a) as in Moretti et al. (2019) */

  double B1, B2, kkk, PI1, PI2, M2;
  B1 = params.Omega0 + 4. * params.OmegaLambda;
  B2 = params.Omega0 / pow(x, 3.) + 4. * params.OmegaLambda;

  for (int ik = 0; ik < NkBINS; ik++)
  {
    if (!ik)
      kkk = 0.0;
    else
      kkk = pow(10., LOGKMIN + ik * DELTALOGK);
    PI1 = kkk * kkk / x / x + 0.5 * H_over_c * H_over_c * pow(B2, 3.) / (B1 * B1 * FR0);
    PI2 = kkk * kkk / x / x / 2. + 0.5 * H_over_c * H_over_c + pow(B2, 3.) / (B1 * B1 * FR0);
    M2 = params.Omega0 * H_over_c * H_over_c * kkk * kkk * (1.5 * H_over_c / FR0) * (1.5 * H_over_c / FR0) * pow(B2, 5.) / pow(B1, 4.) / (9. * pow(x, 5.));

    dydx[1 + 8 * ik] = a1 * y[1 + 8 * ik] + mu(x, kkk) * b1 * y[2 + 8 * ik];
    dydx[2 + 8 * ik] = y[1 + 8 * ik];
    dydx[3 + 8 * ik] = a1 * y[3 + 8 * ik] + mu(x, kkk) * b1 * y[4 + 8 * ik] - (mu(x, kkk) - M2 / PI1 / PI2 / PI2) * b1 * y[2 + 8 * ik] * y[2 + 8 * ik];
    dydx[4 + 8 * ik] = y[3 + 8 * ik];

    /* third-order growth is not used, it is set to the LCDM one */
    dydx[5 + 8 * ik] = a1 * y[5 + 8 * ik] + b1 * y[6 + 8 * ik] - 2. * b1 * y[2 + 8 * ik] * y[2 + 8 * ik] * y[2 + 8 * ik];                                           /* d^2 D31/ dt^2 */
    dydx[6 + 8 * ik] = y[5 + 8 * ik];                                                                                                                               /* d D31/ dt */
    dydx[7 + 8 * ik] = a1 * y[7 + 8 * ik] + b1 * y[8 + 8 * ik] - 2. * b1 * y[2 + 8 * ik] * y[4 + 8 * ik] + 2. * b1 * y[2 + 8 * ik] * y[2 + 8 * ik] * y[2 + 8 * ik]; /* d^2 D32/ dt^2 */
    dydx[8 + 8 * ik] = y[7 + 8 * ik];                                                                                                                               /* d D32/ dt */
  }

  return GSL_SUCCESS;
#endif

  return GSL_FAILURE;
}

int system_of_ODEs_small(double x, const double y[], double *dydx, void *param)
{

  double E2;

  /*
    0 : cosmic time t
  */

  /* cosmic time using uniform H(z) path */
  E2 = Ez2_from_a(x);
  dydx[0] = 1.0 / x / sqrt(E2);

  return GSL_SUCCESS;
}

int jac(double t, const double y[], double *dfdy, double dfdt[], void *params)
{
  printf("This integration method should not call this function.\n");
  return GSL_FAILURE;
}

/* integrand for the computation of comoving distance */
double IntegrandComovingDistance(double z, void *param)
{
  return 1. / Hubble(z);
}

/***********************************/
/* EQUATION OF STATE OF DE SECTION */
/***********************************/

int read_TabulatedEoS(void)
{
  /* Reads the tabulated Equation of State of dark energy from a file */

  int i;
  int err = 0;
  FILE *fd;
  double a, w;
  double *scalef, *EoS;

  if (!ThisTask)
  {
    if (!(fd = fopen(params.TabulatedEoSfile, "r")))
    {
      printf("ERROR on task 0: can't open tabulated EoS in file '%s'\n", params.TabulatedEoSfile);
      fflush(stdout);
      err = 1;
    }

    if (!err)
    {
      NtabEoS = 0;
      while (1)
      {
        if (fscanf(fd, " %lg %lg ", &a, &w) == 2)
          NtabEoS++;
        else
          break;
      }

      fclose(fd);

      if (!NtabEoS)
      {
        printf("ERROR on task 0: can't read tabulated EoS in file '%s'\n", params.TabulatedEoSfile);
        fflush(stdout);
        err = 1;
      }
      else
      {
        printf("Found %d pairs of values in tabulated EoS file\n", NtabEoS);
        fflush(stdout);
      }
    }
  }

  /* all tasks agree on success/failure before proceeding to broadcasts */
  MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (err)
    return 1;

  /* Task 0 communicates the number of lines to all tasks */
  MPI_Bcast(&NtabEoS, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);

  SPLINE[SP_EOS] = gsl_spline_alloc(gsl_interp_cspline, NtabEoS);

  scalef = (double *)malloc(NtabEoS * sizeof(double));
  EoS = (double *)malloc(NtabEoS * sizeof(double));

  if (!ThisTask)
  {
    fd = fopen(params.TabulatedEoSfile, "r");

    for (i = 0; i < NtabEoS; i++)
    {
      if (fscanf(fd, " %lg %lg ", &a, &w) == 2)
      {
        scalef[i] = log10(a);
        EoS[i] = w;
      }
    }

    fclose(fd);
  }

  /* Task 0 broadcasts the Pk and logk vectors to all tasks */
  MPI_Bcast(scalef, NtabEoS * sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(EoS, NtabEoS * sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);

  if (checked_spline_init(SPLINE[SP_EOS], scalef, EoS, NtabEoS, "SP_EOS"))
  {
    free(EoS);
    free(scalef);
    return 1;
  }

  free(EoS);
  free(scalef);

  return 0;
}

/* parametric equation of state for Dark Energy */
double DE_EquationOfState(double a)
{
  if (!NtabEoS)
    return params.DEw0 + (1 - a) * params.DEwa;
  else
    return my_spline_eval(SPLINE[SP_EOS], log10(a), ACCEL[SP_EOS]);
}
/* optional: read external H(z) table */
#ifdef READ_HUBBLE_TABLE
int read_TabulatedHubble(void)
{
  /* Reads a table with redshift z and E(z)=H(z)/H0 (dimensionless)
     and builds a spline over log10(a). */
  if (!strcmp(params.HubbleTableFile, "no") || !strcmp(params.HubbleTableFile, "\0"))
    return 0; /* nothing to do */

  int n = 0;
  int err = 0;
  FILE *fd;
  double z, H;
  double *ax, *Hz;

  if (!ThisTask)
  {
    if (!(fd = fopen(params.HubbleTableFile, "r")))
    {
      printf("ERROR on task 0: can't open Hubble table file '%s' (READ_HUBBLE_TABLE)\n", params.HubbleTableFile);
      fflush(stdout);
      err = 1;
    }

    if (!err)
    {
      while (fscanf(fd, " %lf %lf", &z, &H) == 2)
        n++;
      fclose(fd);

      if (!n)
      {
        printf("ERROR on task 0: Hubble table file '%s' is empty or malformed (READ_HUBBLE_TABLE)\n", params.HubbleTableFile);
        fflush(stdout);
        err = 1;
      }
    }
  }

  /* all tasks agree on success/failure before proceeding to broadcasts */
  MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (err)
    return 1;

  MPI_Bcast(&n, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
  SPLINE[SP_EXT_HUBBLE] = gsl_spline_alloc(gsl_interp_cspline, n);
  ax = (double *)malloc(n * sizeof(double));
  Hz = (double *)malloc(n * sizeof(double));

  if (!ThisTask)
  {
    fd = fopen(params.HubbleTableFile, "r");
    int i = 0;
    while (i < n && fscanf(fd, " %lf %lf", &z, &H) == 2)
    {
      double a = 1.0 / (1.0 + z);
      ax[i] = log10(a);
      Hz[i] = log10(H); /* H stores E(z)=H(z)/H0, dimensionless */
      i++;
    }
    fclose(fd);

    /* GSL requires strictly increasing x values; if the file lists z in
       ascending order, log10(a) ends up in descending order — reverse. */
    if (n > 1 && ax[0] > ax[n - 1])
    {
      for (int j = 0; j < n / 2; j++)
      {
        double tmp;
        tmp = ax[j];
        ax[j] = ax[n - 1 - j];
        ax[n - 1 - j] = tmp;
        tmp = Hz[j];
        Hz[j] = Hz[n - 1 - j];
        Hz[n - 1 - j] = tmp;
      }
    }

    /* Validate absolute normalization at z=0. This matters because the
       code reconstructs H(z)=100*h*E(z) from the tabulated value.
       If z=0 is missing, add E(0)=1 automatically. */
    if (n > 0)
    {
      const double a_tol = 1.e-10;
      const double e_tol = 1.e-4;
      int have_z0 = 0;
      double e0 = 0.0;

      /* search the whole array for z=0, i.e. log10(a)=0 */
      for (int j = 0; j < n; j++)
        if (fabs(ax[j]) < a_tol)
        {
          have_z0 = 1;
          e0 = pow(10., Hz[j]);
          break;
        }

      if (!have_z0)
      {
        printf("WARNING on task 0: Hubble table '%s' does not include z=0. "
               "Adding E(0)=1 automatically.\n",
               params.HubbleTableFile);
        fflush(stdout);
        double *tmp_ax = (double *)realloc(ax, (n + 1) * sizeof(double));
        double *tmp_Hz = (double *)realloc(Hz, (n + 1) * sizeof(double));
        if (!tmp_ax || !tmp_Hz)
        {
          printf("ERROR on task 0: realloc failed while extending Hubble table\n");
          fflush(stdout);
          if (tmp_ax)
            ax = tmp_ax;
          if (tmp_Hz)
            Hz = tmp_Hz;
          err = 1;
        }
        else
        {
          ax = tmp_ax;
          Hz = tmp_Hz;
          /* find insertion point to keep ax[] strictly increasing */
          int ip = n;
          for (int j = 0; j < n; j++)
            if (ax[j] > 0.0)
            {
              ip = j;
              break;
            }
          /* shift elements from ip..n-1 one position to the right */
          memmove(&ax[ip + 1], &ax[ip], (n - ip) * sizeof(double));
          memmove(&Hz[ip + 1], &Hz[ip], (n - ip) * sizeof(double));
          ax[ip] = 0.0; /* log10(a) = log10(1) = 0 */
          Hz[ip] = 0.0; /* log10(E) = log10(1) = 0 */
          n++;
        }
      }
      else if (fabs(e0 - 1.0) > e_tol)
      {
        printf("ERROR on task 0: Hubble table '%s' is not normalized as "
               "E(z)=H(z)/H0 at z=0: found E(0)=%.8g, expected 1 "
               "within %.1e.\n",
               params.HubbleTableFile, e0, e_tol);
        fflush(stdout);
        err = 1;
      }
    }
  }

  MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (err)
  {
    gsl_spline_free(SPLINE[SP_EXT_HUBBLE]);
    SPLINE[SP_EXT_HUBBLE] = NULL;
    free(Hz);
    free(ax);
    return 1;
  }

  /* n may have changed if z=0 was added; re-broadcast and reallocate spline */
  MPI_Bcast(&n, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
  gsl_spline_free(SPLINE[SP_EXT_HUBBLE]);
  SPLINE[SP_EXT_HUBBLE] = gsl_spline_alloc(gsl_interp_cspline, n);

  if (ThisTask)
  {
    double *tmp_ax = (double *)realloc(ax, n * sizeof(double));
    double *tmp_Hz = (double *)realloc(Hz, n * sizeof(double));
    if (!tmp_ax || !tmp_Hz)
    {
      printf("ERROR on task %d: realloc failed while resizing Hubble table\n", ThisTask);
      fflush(stdout);
      if (tmp_ax)
        ax = tmp_ax;
      if (tmp_Hz)
        Hz = tmp_Hz;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    ax = tmp_ax;
    Hz = tmp_Hz;
  }

  MPI_Bcast(ax, n * sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(Hz, n * sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);

  if (checked_spline_init(SPLINE[SP_EXT_HUBBLE], ax, Hz, n, "SP_EXT_HUBBLE"))
  {
    free(Hz);
    free(ax);
    return 1;
  }
  free(Hz);
  free(ax);
  return 0;
}
#endif

double IntegrandForEoS(double a, void *param)
{
  return DE_EquationOfState(a) / a;
}

/**************************/
/* POWER SPECTRUM SECTION */
/**************************/

/* Part of this code is taken from N-GenIC (see GenIC.c) */
double PowerSpec_Tabulated(double);
double PowerSpec_Efstathiou(double);
double PowerSpec_EH(double);
double PowerSpec_PowerLaw(double);
double transf_EH(double);
double T0(double, double, double);

// #define BODE_ET_AL_A9

double PowerSpectrum(double k)
{
  /* input: k (1/Mpc) */
  /* output: P(k) (Mpc^3) */

  double power, alpha, Tf;

  switch (WhichSpectrum)
  {
  case 1:
    power = PowerSpec_EH(k);
    break;

  case 2:
    power = PowerSpec_Tabulated(k);
    break;

  case 3:
    power = PowerSpec_Efstathiou(k);
    break;

  case 4:
    power = PowerSpec_PowerLaw(k);
    break;

  case 5:
    power = PowerSpec_Tabulated(k);
    break;

  default:
    power = 0.0;
    break;
  }

  if (params.WDM_PartMass_in_kev > 0.)
  {
#ifdef BODE_ET_AL_A9
    /* Eqn. (A9) in Bode, Ostriker & Turok (2001), assuming gX=1.5  */
    /* NB: this is a length in Mpc/h*/
    alpha =
        0.048 * pow((params.Omega0 - params.OmegaBaryon) / 0.4, 0.15) * pow(params.Hubble100 / 0.65, 1.3) * pow(1.0 / params.WDM_PartMass_in_kev, 1.15);
    Tf = pow(1 + pow(alpha * k / params.Hubble100 * (3.085678e24 / UnitLength_in_cm), 2 * 1.2), -5.0 / 1.2);
#else
    /* Eqn. just after (A7) in Bode, Ostriker & Turok (2001), assuming gX=1.5  */
    /* NB: this is a length in Mpc/h*/
    alpha =
        0.05 * pow((params.Omega0 - params.OmegaBaryon) / 0.4, 0.15) * pow(params.Hubble100 / 0.65, 1.3) * pow(1.0 / params.WDM_PartMass_in_kev, 1.15);

    Tf = pow(1 + pow(alpha * k / params.Hubble100 * (3.085678e24 / UnitLength_in_cm), 2), -5.0);
#endif
    power *= Tf * Tf;
  }

  return PkNorm * power;
}

int initialize_PowerSpectrum(void)
{
  /* Different options for the power spectrum */

  if (!strcmp(params.FileWithInputSpectrum, "no") || !strcmp(params.FileWithInputSpectrum, "EH"))
  {
    WhichSpectrum = 1;
    if (!ThisTask)
      printf("Power spectrum will be given by the Einsenstein & Hu fit\n");
  }
  else if (!strcmp(params.FileWithInputSpectrum, "Efstathiou"))
  {
    WhichSpectrum = 3;
    if (!ThisTask)
      printf("Power spectrum will be given by the Efstathiou fit with Gamma=%4.2f\n", SHAPE_EFST);
  }
  else if (!strcmp(params.FileWithInputSpectrum, "PowerLaw"))
  {
    WhichSpectrum = 4;
    if (!ThisTask)
      printf("Power spectrum will be a power law with slope %6.3f\n", params.PrimordialIndex);
  }
  else if (!strcmp(params.FileWithInputSpectrum, "CAMBTable"))
  {
#if defined(SCALE_DEPENDENT) && defined(READ_PK_TABLE)
    WhichSpectrum = 5;
    if (!ThisTask)
      printf("Scale-dependent power spectrum will be read from CAMB files\n");
#else
    if (!ThisTask)
      printf("ERROR: to read CAMBTable P(k) use the SCALE_DEPENDENT and READ_PK_TABLE options\n");
    return 1;
#endif
  }
  else
  {
    WhichSpectrum = 2;
    if (read_Pk_from_file())
      return 1;
  }

  if (params.WDM_PartMass_in_kev > 0. && !ThisTask)
  {
    printf("A WDM cut will be applied to the power spectrum following Bode, Ostriker & Turok\n");
  }

  return 0;
}

int normalize_PowerSpectrum(void)
{
  double tmp;

  WindowFunctionType = 2;
  PkNorm = 1.0;
  if (params.Sigma8 != 0.0 && WhichSpectrum != 5)
  {
    tmp = params.Sigma8 * params.Sigma8 / ComputeMassVariance(8.0 / params.Hubble100);
    PkNorm = tmp;
    if (!ThisTask)
    {
      if (WhichSpectrum == 2)
        printf("Warning: you have read a P(k) from file but set its normalization through the parameter file\nThis is fine as long as you know what you are doing, but if you trust the normalization\nof the P(k) you have provided set Sigma8 to 0\n");
      printf("Normalization constant for the power spectrum: %g\n", PkNorm);
    }
  }
  else
  {
    params.Sigma8 = sqrt(ComputeMassVariance(8.0 / params.Hubble100));
    if (!ThisTask)
      printf("Normalization of the provided P(k): Sigma8=%f\n", params.Sigma8);
  }

  return 0;
}

int read_Pk_from_file(void)
{
  /* This is adapted from N-GenIC */

  int i, oldread = 0;
  int err = 0;
  FILE *fd;
  double k, p;
  double *logk, *Pk;

  if (!ThisTask)
  {
    if (!(fd = fopen(params.FileWithInputSpectrum, "r")))
    {
      printf("ERROR on task 0: can't open input spectrum in file '%s' on task %d\n",
             params.FileWithInputSpectrum, ThisTask);
      fflush(stdout);
      err = 1;
    }

    if (!err)
    {
      NPowerTable = 0;
      do
      {
        if (fscanf(fd, " %lg %lg ", &k, &p) == 2)
          NPowerTable++;
        else
          break;
      } while (1);

      fclose(fd);

      if (!NPowerTable)
      {
        printf("ERROR on task 0: can't read data from input spectrum in file '%s'\n",
               params.FileWithInputSpectrum);
        fflush(stdout);
        err = 1;
      }
      else
      {
        printf("Found %d pairs of values in input spectrum table\n", NPowerTable);
        fflush(stdout);
      }
    }
  }

  /* all tasks agree on success/failure before proceeding to broadcasts */
  MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (err)
    return 1;

  /* Task 0 communicates the number of lines to all tasks */
  MPI_Bcast(&NPowerTable, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);

  SPLINE[SP_PK] = gsl_spline_alloc(gsl_interp_cspline, NPowerTable);

  logk = (double *)malloc(NPowerTable * sizeof(double));
  Pk = (double *)malloc(NPowerTable * sizeof(double));

  if (!ThisTask)
  {
    fd = fopen(params.FileWithInputSpectrum, "r");

    for (i = 0; i < NPowerTable; i++)
    {
      if (fscanf(fd, " %lg %lg ", &k, &p) == 2)
      {

        if (!i)
        {
          if (k < 0.0)
          {
            oldread = 1;
            printf("I am assuming that file %s contains Log k - Log k^3 P(k)\n", params.FileWithInputSpectrum);
          }
          else
          {
            oldread = 0;
            printf("I am assuming that file %s contains k - P(k)\n", params.FileWithInputSpectrum);
          }
        }

        if (oldread)
        {
          logk[i] = k;
          Pk[i] = p;
        }
        else
        {
          logk[i] = log10(k);
          Pk[i] = log10(p * k * k * k);
        }

        /* translates into Htrue */
        logk[i] += log10(params.Hubble100);
        /* fixes the units if necessary */
        if (params.InputSpectrum_UnitLength_in_cm != 0.0)
          logk[i] += log10(params.InputSpectrum_UnitLength_in_cm / UnitLength_in_cm);
      }
    }

    fclose(fd);
  }

  /* Task 0 broadcasts the Pk and logk vectors to all tasks */
  MPI_Bcast(Pk, NPowerTable * sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(logk, NPowerTable * sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);

  if (checked_spline_init(SPLINE[SP_PK], logk, Pk, NPowerTable, "SP_PK (from file)"))
  {
    free(Pk);
    free(logk);
    return 1;
  }

  free(Pk);
  free(logk);

  return 0;
}

#ifdef READ_PK_TABLE
int read_Pk_table_from_CAMB(double *scalef, double *grow1, double *grow2, double *grow31, double *grow32,
                            double *fomega1, double *fomega2, double *fomega31, double *fomega32)
{
  /* Reads P_cb(k) at various redshifts from CAMB outputs and fills:
     - SPLINE[SP_PK] with the z=0 spectrum (stored as log10[k^3 P(k)])
     - scale-dependent growth arrays from ratios of k^3 P(k,z) to k^3 P(k,0)

     Expected CAMB files:
       MatterFile_000.dat, ..., MatterFile_(NCAMB-1).dat
       RedshiftsFile (two columns: index  redshift)
     Units expected from CAMB:
       k in h/Mpc, P(k) in (Mpc/h)^3.
  */

  int i, j, dummy, i1, i2, First, Today;
  int err = 0;
  double kappa, myPk, z, Om, slope;
  char filename[LBLENGTH], buffer[LBLENGTH], *ugo;
  FILE *fd;
  double *logk, *Pk, *CAMBScalefac, *lingrow;

  if (!ThisTask)
  {
    /* count the number of CAMB P(k) files and lines */
    params.camb.NCAMB = 0;
    sprintf(filename, "%s_%03d.dat", params.camb.MatterFile, params.camb.NCAMB);
    while ((fd = fopen(filename, "r")) != 0x0)
    {
      if (!params.camb.NCAMB)
      {
        NPowerTable = 0;
        while (!feof(fd))
        {
          ugo = fgets(buffer, LBLENGTH, fd);
          if (ugo && sscanf(buffer, "%lf", &kappa) == 1)
            NPowerTable++;
        }
      }
      fclose(fd);
      params.camb.NCAMB++;
      sprintf(filename, "%s_%03d.dat", params.camb.MatterFile, params.camb.NCAMB);
    }

    if (!params.camb.NCAMB)
    {
      printf("Error on Task 0: CAMB file %s not found\n", filename);
      err = 1;
    }
    else if (!NPowerTable)
    {
      sprintf(filename, "%s_%03d.dat", params.camb.MatterFile, 0);
      printf("Error on Task 0: problem in reading CAMB file %s\n", filename);
      err = 1;
    }

    if (!err)
      printf("Found %d CAMB matter power files with %d lines each\n",
             params.camb.NCAMB, NPowerTable);
  }

  /* all tasks agree on success/failure before proceeding to broadcasts */
  MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (err)
    return 1;

  /* broadcast sizes */
  MPI_Bcast(&params.camb.NCAMB, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&NPowerTable, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);

  /* allocate */
  SPLINE[SP_PK] = gsl_spline_alloc(gsl_interp_cspline, NPowerTable);
  logk = (double *)malloc(NPowerTable * sizeof(double));
  Pk = (double *)malloc(NPowerTable * sizeof(double));
  CAMBScalefac = (double *)malloc(params.camb.NCAMB * sizeof(double));
  lingrow = (double *)malloc(params.camb.NCAMB * NPowerTable * sizeof(double));

  if (!ThisTask)
  {
    /* read redshifts and convert to scale factors */
    if ((fd = fopen(params.camb.RedshiftsFile, "r")) == 0x0)
    {
      printf("Error: Redshift file %s not found\n", params.camb.RedshiftsFile);
      err = 1;
    }

    if (!err)
    {
      for (i = 0; i < params.camb.NCAMB; i++)
        fscanf(fd, "%d %lf", &dummy, CAMBScalefac + i);
      fclose(fd);

      if (CAMBScalefac[params.camb.NCAMB - 1] != 0.0)
      {
        printf("ERROR on Task 0: last CAMB redshift must be 0.0\n");
        err = 1;
      }
    }

    if (!err)
    {
      for (i = 0; i < params.camb.NCAMB; i++)
        CAMBScalefac[i] = 1. / (1. + CAMBScalefac[i]);

      /* loop over files starting from z=0 one (last index) */
      for (i = params.camb.NCAMB - 1; i >= 0; i--)
      {
        sprintf(filename, "%s_%03d.dat", params.camb.MatterFile, i);
        fd = fopen(filename, "r");
        for (j = 0; j < NPowerTable; j++)
        {
          fscanf(fd, "%lf %lf", &kappa, &myPk); /* kappa: k in h/Mpc; myPk: P in (Mpc/h)^3 */

          if (i == params.camb.NCAMB - 1)
          {
            /* store z=0 table as log10[k^3 P(k)] and k (true 1/Mpc) */
            Pk[j] = log10(kappa * kappa * kappa * myPk);
            logk[j] = log10(kappa * params.Hubble100); /* k_true = h * (k in h/Mpc) */
            if (params.InputSpectrum_UnitLength_in_cm != 0.0)
              logk[j] += log10(params.InputSpectrum_UnitLength_in_cm / UnitLength_in_cm);
            lingrow[i + j * params.camb.NCAMB] = 0.0;
          }
          else
          {
            /* growth from half the difference of log10[k^3 P(k,z)] and z=0 */
            lingrow[i + j * params.camb.NCAMB] =
                0.5 * (log10(kappa * kappa * kappa * myPk) - Pk[j]);
          }
        }
        fclose(fd);
      }

      /* consistency of k-range w.r.t. def_splines.h grid */
      if (logk[0] > LOGKMIN || logk[NPowerTable - 1] < LOGKMIN + DELTALOGK * (NkBINS - 1))
      {
        printf("ERROR: CAMB P(k) tables run from k=%10g to k=%10g 1/Mpc\n",
               pow(10., logk[0]), pow(10., logk[NPowerTable - 1]));
        printf("       while the growth rate is requested from k=%10g to k=%10g 1/Mpc\n",
               pow(10., LOGKMIN), pow(10., LOGKMIN + DELTALOGK * (NkBINS - 1)));
        printf("       please extend the k range in CAMB or fix LOGKMIN, DELTALOGK and NkBINS in def_splines.h\n");
        err = 1;
      }
    }
  }

  /* all tasks agree on success/failure before proceeding to data broadcasts */
  MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (err)
    goto fail_arrays;

  /* broadcast data to all tasks */
  MPI_Bcast(Pk, NPowerTable, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(logk, NPowerTable, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(CAMBScalefac, params.camb.NCAMB, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(lingrow, params.camb.NCAMB * NPowerTable, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  /* spline of P(k) at z=0 (stored as log10[k^3 P]) */
  if (checked_spline_init(SPLINE[SP_PK], logk, Pk, NPowerTable, "SP_PK (from CAMB)"))
    goto fail_arrays;

  /* 2D spline for growth(a,k) defined via lingrow(a,k) */
  const gsl_interp2d_type *T = gsl_interp2d_bicubic;
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();
  gsl_spline2d *AnotherSpline = gsl_spline2d_alloc(T, params.camb.NCAMB, NPowerTable);
  if (checked_spline2d_init(AnotherSpline, CAMBScalefac, logk, lingrow, params.camb.NCAMB, NPowerTable, "CAMB growth(a,k)"))
    goto fail_all;

  /* first usable time index within CAMB range */
  for (First = 0; First < NBINS && scalef[First] < CAMBScalefac[0]; First++)
    ;

  /* fill growth tables for a<=1 using the bicubic interpolation */
  for (i = First; i < NBINS && scalef[i] <= 1; i++)
    for (j = 0; j < NkBINS; j++)
    {
      double logk_req = LOGKMIN + j * DELTALOGK;
      z = 1. / scalef[i] - 1.;
      Om = OmegaMatter(z);
      grow1[i + j * NBINS] = pow(10., gsl_spline2d_eval(AnotherSpline, 1. / (1. + z), logk_req, xacc, yacc));
      grow2[i + j * NBINS] = 3. / 7. * pow(grow1[i + j * NBINS], 2.0) * pow(Om, -1. / 143.);
      grow31[i + j * NBINS] = grow1[i + j * NBINS] * grow1[i + j * NBINS] * grow1[i + j * NBINS] * pow(Om, -4. / 275.) / 9.;
      grow32[i + j * NBINS] = grow1[i + j * NBINS] * grow1[i + j * NBINS] * grow1[i + j * NBINS] * pow(Om, -268. / 17875.) * 5. / 42.;
    }

  Today = i - 1;

  /* extrapolate a>1 as power laws */
  for (j = 0; j < NkBINS; j++)
  {
    slope = log10(grow1[Today + j * NBINS] / grow1[Today - 1 + j * NBINS]) /
            log10(scalef[Today] / scalef[Today - 1]);
    for (i = Today + 1; i < NBINS; i++)
    {
      grow1[i + j * NBINS] = grow1[Today + j * NBINS] * pow(scalef[i] / scalef[Today], 1.0 * slope);
      grow2[i + j * NBINS] = grow2[Today + j * NBINS] * pow(scalef[i] / scalef[Today], 2.0 * slope);
      grow31[i + j * NBINS] = grow31[Today + j * NBINS] * pow(scalef[i] / scalef[Today], 3.0 * slope);
      grow32[i + j * NBINS] = grow32[Today + j * NBINS] * pow(scalef[i] / scalef[Today], 3.0 * slope);
    }
  }

  /* scale with a for earlier times (before CAMB first a) */
  for (j = 0; j < NkBINS; j++)
    for (i = 0; i < First; i++)
    {
      grow1[i + j * NBINS] = grow1[First + j * NBINS] * scalef[i] / scalef[First];
      grow2[i + j * NBINS] = grow2[First + j * NBINS] * pow(scalef[i] / scalef[First], 2.);
      grow31[i + j * NBINS] = grow31[First + j * NBINS] * pow(scalef[i] / scalef[First], 3.);
      grow32[i + j * NBINS] = grow32[First + j * NBINS] * pow(scalef[i] / scalef[First], 3.);
    }

  /* f(Omega) arrays */
  for (j = 0; j < NkBINS; j++)
  {
    for (i = 0; i < Today; i++)
    {
      if (i == 0)
      {
        i1 = 0;
        i2 = 2;
      }
      else
      {
        i1 = i - 1;
        i2 = i + 1;
      }
      fomega1[i + j * NBINS] = (grow1[i2 + j * NBINS] - grow1[i1 + j * NBINS]) / (scalef[i2] - scalef[i1]) * scalef[i] / grow1[i + j * NBINS];
      fomega2[i + j * NBINS] = (grow2[i2 + j * NBINS] - grow2[i1 + j * NBINS]) / (scalef[i2] - scalef[i1]) * scalef[i] / grow2[i + j * NBINS];
      fomega31[i + j * NBINS] = (grow31[i2 + j * NBINS] - grow31[i1 + j * NBINS]) / (scalef[i2] - scalef[i1]) * scalef[i] / grow31[i + j * NBINS];
      fomega32[i + j * NBINS] = (grow32[i2 + j * NBINS] - grow32[i1 + j * NBINS]) / (scalef[i2] - scalef[i1]) * scalef[i] / grow32[i + j * NBINS];
    }
    /* extrapolate to today and future with linear slope */
    slope = (fomega1[Today - 1 + j * NBINS] - fomega1[Today - 2 + j * NBINS]) / (scalef[Today - 1] - scalef[Today - 2]);
    for (i = Today; i < NBINS; i++)
      fomega1[i + j * NBINS] = fomega1[Today - 1 + j * NBINS] + slope * (scalef[i] - scalef[Today - 1]);

    slope = (fomega2[Today - 1 + j * NBINS] - fomega2[Today - 2 + j * NBINS]) / (scalef[Today - 1] - scalef[Today - 2]);
    for (i = Today; i < NBINS; i++)
      fomega2[i + j * NBINS] = fomega2[Today - 1 + j * NBINS] + slope * (scalef[i] - scalef[Today - 1]);

    slope = (fomega31[Today - 1 + j * NBINS] - fomega31[Today - 2 + j * NBINS]) / (scalef[Today - 1] - scalef[Today - 2]);
    for (i = Today; i < NBINS; i++)
      fomega31[i + j * NBINS] = fomega31[Today - 1 + j * NBINS] + slope * (scalef[i] - scalef[Today - 1]);

    slope = (fomega32[Today - 1 + j * NBINS] - fomega32[Today - 2 + j * NBINS]) / (scalef[Today - 1] - scalef[Today - 2]);
    for (i = Today; i < NBINS; i++)
      fomega32[i + j * NBINS] = fomega32[Today - 1 + j * NBINS] + slope * (scalef[i] - scalef[Today - 1]);
  }

  gsl_spline2d_free(AnotherSpline);
  gsl_interp_accel_free(yacc);
  gsl_interp_accel_free(xacc);

  free(lingrow);
  free(CAMBScalefac);
  free(Pk);
  free(logk);

  return 0;

fail_all:
  gsl_spline2d_free(AnotherSpline);
  gsl_interp_accel_free(yacc);
  gsl_interp_accel_free(xacc);
fail_arrays:
  free(lingrow);
  free(CAMBScalefac);
  free(Pk);
  free(logk);
  return 1;
}
#endif

double PowerSpec_Tabulated(double k)
{
  return pow(10., my_spline_eval(SPLINE[SP_PK], log10(k), ACCEL[SP_PK])) / k / k / k;
}

double PowerSpec_Efstathiou(double k)
{
  return pow(k, params.PrimordialIndex) / pow(1 + pow(6.4 / SHAPE_EFST * k + pow(3.0 / SHAPE_EFST * k, 1.5) + pow(1.7 / SHAPE_EFST, 2.0) * k * k, 1.13), 2 / 1.13);
}

double PowerSpec_PowerLaw(double k)
{
  return pow(k, params.PrimordialIndex);
}

double PowerSpec_EH(double k)
{
  return pow(k, params.PrimordialIndex) * pow(transf_EH(k), 2.);
}

double transf_EH(double fk)
{
  static double Teta_27 = 1.0104;
  double q, Omegac, Oh2, b1, b2, zd, Rd, zeq, Req, keq, s, ks, alc, bec, f, Tc, beb, bno, kst, ksi, Tb, Tr, y, alb, Ob2, OB;

  /* Eisenstein & Hu fit of the transfer function */

  OB = (params.OmegaBaryon > 1.e-6 ? params.OmegaBaryon : 1.e-6);
  Omegac = params.Omega0 - OB;
  Oh2 = params.Omega0 * params.Hubble100 * params.Hubble100;
  Ob2 = OB * params.Hubble100 * params.Hubble100;
  b1 = 0.313 * pow(Oh2, -0.419) * (1 + 0.607 * pow(Oh2, 0.674));
  b2 = 0.238 * pow(Oh2, 0.223);
  zd = 1291. * pow(Oh2, 0.251) * (1. + b1 * pow(Ob2, b2)) / (1. + 0.659 * pow(Oh2, 0.828));
  Rd = 31.5 * Ob2 / (pow(Teta_27, 4.0) * 0.001 * zd);
  zeq = 2.5e4 * Oh2 / pow(Teta_27, 4.0);
  Req = 31.5 * Ob2 / (pow(Teta_27, 4.0) * 0.001 * zeq);
  keq = 7.46e-2 * Oh2 / Teta_27 / Teta_27;
  s = 1.633 * log((sqrt(1. + Rd) + sqrt(Rd + Req)) / (1 + sqrt(Req))) / (keq * sqrt(Req)); /* 2sqrt(6)/3=1.633 */
  ks = fk * s;
  q = fk * Teta_27 * Teta_27 / Oh2;
  alc = pow(pow(46.9 * Oh2, 0.670) * (1. + pow(32.1 * Oh2, -0.532)), -OB / params.Omega0) *
        pow(pow(12.0 * Oh2, 0.424) * (1. + pow(45.0 * Oh2, -0.582)), -pow(OB / params.Omega0, 3.0));
  bec = 1. / (1. + (0.944 / (1. + pow(458. * Oh2, -0.708))) * (pow(Omegac / params.Omega0, pow(0.395 * Oh2, -0.0266)) - 1.));
  f = 1. / (1 + pow(ks / 5.4, 4.0));
  Tc = f * T0(q, 1., bec) + (1. - f) * T0(q, alc, bec);
  beb = 0.5 + OB / params.Omega0 + (3. - 2. * OB / params.Omega0) * sqrt(pow(17.2 * Oh2, 2.0) + 1.);
  bno = 8.41 * pow(Oh2, 0.435);
  kst = ks / pow(1. + pow(bno / ks, 3.0), 0.3333);
  ksi = 1.6 * pow(Ob2, 0.52) * pow(Oh2, 0.73) * (1. + pow(10.4 * Oh2, -0.95));
  y = (1. + zeq) / (1 + zd);
  alb = 2.07 * keq * s * pow(1.0 + Rd, -0.75) * (y * (-6. * sqrt(1. + y) + (2. + 3. * y) * log((sqrt(1. + y) + 1.) / (sqrt(1. + y) - 1.))));
  Tb = (T0(q, 1., 1.) / (1. + pow(ks / 5.2, 2.0)) + alb / (1. + pow(beb / ks, 3.0)) * exp(-pow(fk / ksi, 1.4))) * sin(kst) / kst;
  Tr = (OB * Tb + Omegac * Tc) / params.Omega0;

  return Tr;
}

double T0(double q, double a, double b)
{
  double ll, C;

  ll = log(exp(1.) + 1.8 * b * q);
  C = 14.2 / a + 386. / (1. + 69.9 * pow(q, 1.08));

  return ll / (ll + C * q * q);
}

/*************************/
/* MASS VARIANCE SECTION */
/*************************/

double IntegrandForMassVariance(double, void *);
double IntegrandForDisplVariance(double, void *);

int initialize_MassVariance(void)
{
  int i;
  double r;
  double rmin = -6.0, dr = 0.04; /* NB: rmax=2.0 */
  double *rv, *massvar, *dmvdr, *massvarneg, *displv;

  rv = (double *)malloc(NBINS * sizeof(double));
  massvar = (double *)malloc(NBINS * sizeof(double));
  dmvdr = (double *)malloc(NBINS * sizeof(double));
  massvarneg = (double *)malloc(NBINS * sizeof(double));
  displv = (double *)malloc(NBINS * sizeof(double));

  for (i = NBINS - 1; i >= 0; i--)
  {
    rv[i] = (rmin + i * dr);
    r = pow(10., rv[i]);
    massvar[i] = log10(ComputeMassVariance(r));
    /* To allow interpolation of reverse function,
 the variance must increase with decreasing radius */
    if (i < NBINS - 1 && massvar[i] - massvar[i + 1] < 1.e-6)
      massvar[i] = massvar[i + 1] + 1.e-6;

    massvarneg[i] = -massvar[i];
    displv[i] = log10(ComputeDisplVariance(r));
  }

  for (i = 0; i < NBINS; i++)
  {
    if (i == 0)
      dmvdr[i] = (massvar[i + 1] - massvar[i]) / (rv[i + 1] - rv[i]);
    else if (i == NBINS - 1)
      dmvdr[i] = (massvar[i] - massvar[i - 1]) / (rv[i] - rv[i - 1]);
    else
      dmvdr[i] = (massvar[i + 1] - massvar[i - 1]) / (rv[i + 1] - rv[i - 1]);
  }

  /* initialization of splines for interpolation */
  if (checked_spline_init(SPLINE[SP_MASSVAR], rv, massvar, NBINS, "SP_MASSVAR") ||
      checked_spline_init(SPLINE[SP_RADIUS], massvarneg, rv, NBINS, "SP_RADIUS") ||
      checked_spline_init(SPLINE[SP_DVARDR], rv, dmvdr, NBINS, "SP_DVARDR") ||
      checked_spline_init(SPLINE[SP_DISPVAR], rv, displv, NBINS, "SP_DISPVAR"))
    return 1;

  free(displv);
  free(massvarneg);
  free(dmvdr);
  free(massvar);
  free(rv);

  return 0;
}

double ComputeMassVariance(double R)
{
  double result, error;
  gsl_function Function;
  double ThisRadius = R;

  Function.function = &IntegrandForMassVariance;
  Function.params = &ThisRadius;
  gsl_integration_qags(&Function, -10., log(500.0 / R), 0.0, TOLERANCE, NWINT, workspace, &result, &error);

  return result;
}

double IntegrandForMassVariance(double logk, void *param)
{
  double w, k, D;
  double ThisRadius = *(double *)param;

  k = exp(logk);
  w = WindowFunction(k * ThisRadius);
  /* This is superfluous in most cases,
     but is important for scale-dependent cases where D(0,k) is not unity at all k */
  D = GrowingMode(0.0, k);
  return PowerSpectrum(k) * w * w * D * D * k * k * k / (2. * PI * PI);
}

double ComputeDisplVariance(double R)
{
  double result, error;
  gsl_function Function;
  double ThisRadius = R;

  Function.function = &IntegrandForDisplVariance;
  Function.params = &ThisRadius;
  gsl_integration_qags(&Function, -10., log(500.0 / R), 0.0, TOLERANCE, NWINT, workspace, &result, &error);

  return result;
}

double IntegrandForDisplVariance(double logk, void *param)
{
  double w, k, D;
  double ThisRadius = *(double *)param;

  k = exp(logk);
  w = WindowFunction(k * ThisRadius);
  /* This is superfluous in most cases,
     but is important for scale-dependent cases where D(0,k) is not unity at all k */
  D = GrowingMode(0.0, k);
  return PowerSpectrum(k) * w * w * D * D * k / (2. * PI * PI);
}

double WindowFunction(double kr)
{
  /* Window function:

     WindowFunctionType = 0: Gaussian smoothing
     WindowFunctionType = 1: SKS smoothing
     WindowFunctionType = 2: top-hat smoothing

     DIMENSIONLESS
  */
  double window, kr2;

  switch (WindowFunctionType)
  {
  case 0: /* Gaussian */
    window = exp(-kr * kr / 2.);
    break;

  case 1: /* sharp k-space */
    window = (kr < 1 ? 1.0 : 0.0);
    break;

  case 2: /* top-hat */
    if (kr < 1.e-5)
      window = 1.0;
    else
    {
      kr2 = kr * kr;
      window = 3. * (sin(kr) / kr2 / kr - cos(kr) / kr2);
    }
    break;

  default:
    window = 1;
    break;
  }

  return window;
}

double MassVariance(double R)
{
  return pow(10., my_spline_eval(SPLINE[SP_MASSVAR], log10(R), ACCEL[SP_MASSVAR]));
}

double dMassVariance_dr(double R)
{
  return my_spline_eval(SPLINE[SP_DVARDR], log10(R), ACCEL[SP_DVARDR]);
}

double DisplVariance(double R)
{
  return pow(10., my_spline_eval(SPLINE[SP_DISPVAR], log10(R), ACCEL[SP_DISPVAR]));
}

double Radius(double Var)
{
  return pow(10., my_spline_eval(SPLINE[SP_RADIUS], -log10(Var), ACCEL[SP_RADIUS]));
}

/**********************************/
/* COSMOLOGICAL FUNCTIONS SECTION */
/**********************************/

double OmegaMatter(double z)
{
  /* Cosmological mass density parameter as a function of redshift
     DIMENSIONLESS */
  double Ezv = Ez(z);
  return params.Omega0 * pow(1. + z, 3.) / (Ezv * Ezv);
}

double OmegaLambda(double z)
{
  /* Cosmological mass density parameter as a function of redshift
     DIMENSIONLESS */
  double Ezv = Ez(z);
  return params.OmegaLambda / (Ezv * Ezv);
}

double Hubble(double z)
{
  /* Hubble parameter as a function of redshift
     DIMENSION: km/s/Mpc
     If an external E(z)=H(z)/H0 table is provided, reconstruct H(z)
     as 100 * Hubble100 * E(z). */
  double Esq, de_eos;
#ifdef READ_HUBBLE_TABLE
  /* If an external E(z) table is provided, prefer it. */
  if (SPLINE[SP_EXT_HUBBLE])
    return 100 * params.Hubble100 * pow(10, my_spline_eval(SPLINE[SP_EXT_HUBBLE], -log10(1. + z), ACCEL[SP_EXT_HUBBLE]));
#endif

  if (params.simpleLambda)
    Esq = OmegaRad * pow(1. + z, 4.) + params.Omega0 * pow(1. + z, 3.) + OmegaK * pow(1. + z, 2.) + params.OmegaLambda;
  else
  {
    de_eos = my_spline_eval(SPLINE[SP_INTEOS], -log10(1. + z), ACCEL[SP_INTEOS]);
    Esq = OmegaRad * pow(1. + z, 4.) + params.Omega0 * pow(1. + z, 3.) + OmegaK * pow(1. + z, 2.) + params.OmegaLambda * pow(1. + z, 3.) * exp(3. * de_eos);
  }

  return 100. * params.Hubble100 * sqrt(Esq);
}

double Ez(double z)
{
  /* Dimensionless Hubble function E(z) = H(z)/H0 */
  double H0 = Hubble(0.0);
  return Hubble(z) / H0;
}

double Hubble_Gyr(double z)
{
  /* Hubble parameter as a function of redshift
     DIMENSION: Gyr^-1 */

  return Hubble(z) / HUBBLETIME_GYR / 100.;
}

double InterpolateGrowth(double z, double k, int pointer)
{
  /* This function interpolates the table for all scale-dependent growth functions */

#ifdef SCALE_DEPENDENT
  int kk;
  double dk;
  /* NB in modified gravity kmin is set to 0, but this makes log interpolation impossible
     so we leave it to kmin */
  if (k < kmin)
    return my_spline_eval(SPLINE[pointer], -log10(1. + z), ACCEL[pointer]);
  else if (k > kmax)
    return my_spline_eval(SPLINE[pointer + NkBINS - 1], -log10(1. + z), ACCEL[pointer + NkBINS - 1]);
  else
  {
    dk = (log10(k) - LOGKMIN) / DELTALOGK;
    kk = (int)dk;
    dk -= kk;

    return dk * my_spline_eval(SPLINE[pointer + kk + 1], -log10(1. + z), ACCEL[pointer + kk + 1]) +
           (1 - dk) * my_spline_eval(SPLINE[pointer + kk], -log10(1. + z), ACCEL[pointer + kk]);
  }
#else
  /* scale-independent case, just return the interpolation */

  return my_spline_eval(SPLINE[pointer], -log10(1. + z), ACCEL[pointer]);
#endif
}

double fomega(double z, double k)
{
  /* Peebles' f(Omega) function, dlogD/dloga
     DIMENSIONLESS */

  return InterpolateGrowth(z, k, SP_FOMEGA1);
}

double fomega_2LPT(double z, double k)
{
  /* second-order f(Omega) function, dlogD2/dloga
     DIMENSIONLESS */

  return InterpolateGrowth(z, k, SP_FOMEGA2);
}

double fomega_3LPT_1(double z, double k)
{
  /* second-order f(Omega) function, dlogD2/dloga
     DIMENSIONLESS */

  return InterpolateGrowth(z, k, SP_FOMEGA31);
}

double fomega_3LPT_2(double z, double k)
{
  /* second-order f(Omega) function, dlogD2/dloga
     DIMENSIONLESS */

  return InterpolateGrowth(z, k, SP_FOMEGA32);
}

double GrowingMode(double z, double k)
{
  /* linear growing mode, interpolation on the grid
     DIMENSIONLESS */

  return pow(10., InterpolateGrowth(z, k, SP_GROW1));
}

double GrowingMode_2LPT(double z, double k)
{
  /* second-order growing mode, interpolation on the grid
     DIMENSIONLESS */

  return pow(10., InterpolateGrowth(z, k, SP_GROW2));
}

double GrowingMode_3LPT_1(double z, double k)
{
  /* second-order growing mode, interpolation on the grid
     DIMENSIONLESS */

  return -pow(10., InterpolateGrowth(z, k, SP_GROW31));
}

double GrowingMode_3LPT_2(double z, double k)
{
  /* second-order growing mode, interpolation on the grid
     DIMENSIONLESS */

  return pow(10., InterpolateGrowth(z, k, SP_GROW32));
}

#ifdef ELL_CLASSIC
double InverseGrowingMode(double D, int ismooth)
{
  /* redshift corresponding to a linear growing mode, interpolation on the grid
     DIMENSIONLESS */

#ifdef SCALE_DEPENDENT
  return 1. / pow(10., my_spline_eval(SPLINE_INVGROW[ismooth], log10(D), ACCEL_INVGROW[ismooth])) - 1.;
#else
  return 1. / pow(10., my_spline_eval(SPLINE[SP_INVGROW], log10(D), ACCEL[SP_INVGROW])) - 1.;
#endif
}
#endif

double CosmicTime(double z)
{
  /* cosmic time, interpolation on the grid
     Gyr */

  return pow(10., my_spline_eval(SPLINE[SP_TIME], -log10(1. + z), ACCEL[SP_TIME]));
}

double InverseCosmicTime(double t)
{
  /* scale factor corresponding to a cosmic time, interpolation on the grid
     Gyr */

  return pow(10., my_spline_eval(SPLINE[SP_INVTIME], log10(t), ACCEL[SP_INVTIME]));
}

double ComovingDistance(double z)
{
  /* comoving distance, interpolation on the grid
     Mpc */
  if (z <= 0.0)
    return 0.0; /* enforce exact zero at z=0 */
  double chi = my_spline_eval(SPLINE[SP_COMVDIST], -log10(1. + z), ACCEL[SP_COMVDIST]);
  /* Clamp any negative (spline end overshoot) to 0 for physical consistency */
  if (chi < 0.0)
    chi = 0.0;
  return chi;
}

double DiameterDistance(double z)
{
  /* diameter distance, interpolation on the grid
     Mpc */
  if (z <= 0.0)
    return 0.0; /* enforce exact zero at z=0 */
  double da = my_spline_eval(SPLINE[SP_DIAMDIST], -log10(1. + z), ACCEL[SP_DIAMDIST]);
  if (da < 0.0)
    da = 0.0; /* guard against tiny negative overshoot */
  return da;
}

double InverseComovingDistance(double chi)
{
  /* redshift as a function of comoving distance, using an inverse spline
     built once from the (log10 a) -> chi table. */
  if (chi <= 0.0)
    return 0.0;

  if (!INVCOMVDIST_READY)
  {
#ifdef _OPENMP
#pragma omp critical(invchi_init)
    {
      if (!INVCOMVDIST_READY)
      {
#endif
        /* Build inverse mapping chi -> log10(a) from existing SPLINE[SP_COMVDIST]
           Tailored to the PLC redshift span plus a small buffer to maximize accuracy
           where needed while keeping safe linear extrapolation beyond the ends. */
        gsl_spline *fwd = SPLINE[SP_COMVDIST];
        const int n = fwd->size;

        /* Determine PLC redshift span and buffer */
        double z_lo = params.LastzForPLC;
        double z_hi = params.StartingzForPLC;
        if (z_hi < z_lo)
        {
          /* swap if provided in reverse */
          double tmp = z_hi;
          z_hi = z_lo;
          z_lo = tmp;
        }
        /* buffer in redshift: 10% of span, at least 0.02 (and do not go negative) */
        double dz_span = (z_hi > z_lo ? (z_hi - z_lo) : 0.0);
        double dz_buf = dz_span * 0.10;
        if (dz_buf < 0.02)
          dz_buf = 0.02;
        double z_lo_buf = z_lo - dz_buf;
        if (z_lo_buf < 0.0)
          z_lo_buf = 0.0;
        double z_hi_buf = z_hi + dz_buf;

        /* Convert to chi bounds */
        double chi_lo = ComovingDistance(z_lo_buf);
        double chi_hi = ComovingDistance(z_hi_buf);
        if (chi_hi < chi_lo)
        {
          /* numeric safety: ensure chi_lo <= chi_hi */
          double tmp = chi_hi;
          chi_hi = chi_lo;
          chi_lo = tmp;
        }

        /* Count how many forward samples fall within [chi_lo, chi_hi] after reversing */
        int m = 0;
        for (int k = 0; k < n; ++k)
        {
          int kr = n - 1 - k;
          double chi_k = fwd->y[kr];
          if (chi_k + 1e-12 >= chi_lo && chi_k - 1e-12 <= chi_hi)
            ++m;
        }

        /* Ensure a reasonable number of nodes; if too few, fall back to full range */
        int use_full = 0;
        if (m < 8)
          use_full = 1;

        int alloc_n = use_full ? n : m;
        /* allocate temporary arrays */
        double *xchi = (double *)malloc((size_t)alloc_n * sizeof(double));
        double *alog = (double *)malloc((size_t)alloc_n * sizeof(double));
        if (!xchi || !alog)
        {
          if (xchi)
            free(xchi);
          if (alog)
            free(alog);
          /* Mark as unavailable */
          INVCOMVDIST_READY = 2;
        }
        else
        {
          if (use_full)
          {
            /* Reverse entire table so that chi increases (from ~0 to max) */
            for (int k = 0; k < n; ++k)
            {
              int kr = n - 1 - k;
              xchi[k] = fwd->y[kr]; /* chi */
              alog[k] = fwd->x[kr]; /* log10(a) */
            }
          }
          else
          {
            /* Take only the PLC-buffered subrange */
            int t = 0;
            for (int k = 0; k < n; ++k)
            {
              int kr = n - 1 - k;
              double chi_k = fwd->y[kr];
              double alog_k = fwd->x[kr];
              if (chi_k + 1e-12 >= chi_lo && chi_k - 1e-12 <= chi_hi)
              {
                xchi[t] = chi_k;
                alog[t] = alog_k;
                ++t;
              }
            }
            /* Safety: if something went wrong, revert to full */
            if (t < 8)
            {
              /* fill with full */
              for (int k = 0; k < n; ++k)
              {
                int kr = n - 1 - k;
                xchi[k] = fwd->y[kr];
                alog[k] = fwd->x[kr];
              }
              alloc_n = n;
            }
          }

          /* Ensure the first x is non-negative strictly increasing; if tiny negative, clamp */
          if (xchi[0] < 0.0)
            xchi[0] = 0.0;

          /* Optionally inject an exact physical anchor (chi=0 -> a=1 -> log10(a)=0)
             at the start to avoid any low-chi extrapolation and guarantee chi_min=0. */
          int inject_anchor = (xchi[0] > 0.0 + 1e-12);
          double *xchi_use = xchi;
          double *alog_use = alog;
          int n_use = alloc_n;
          if (inject_anchor)
          {
            n_use = alloc_n + 1;
            xchi_use = (double *)malloc((size_t)n_use * sizeof(double));
            alog_use = (double *)malloc((size_t)n_use * sizeof(double));
            if (!xchi_use || !alog_use)
            {
              if (xchi_use)
                free(xchi_use);
              if (alog_use)
                free(alog_use);
              /* fall back to original arrays without anchor */
              xchi_use = xchi;
              alog_use = alog;
              n_use = alloc_n;
              inject_anchor = 0;
            }
            else
            {
              xchi_use[0] = 0.0;
              alog_use[0] = 0.0; /* log10(a=1) */
              for (int k = 0; k < alloc_n; ++k)
              {
                xchi_use[k + 1] = xchi[k];
                alog_use[k + 1] = alog[k];
              }
            }
          }

          /* Allocate and init inverse spline */
          SPLINE_INVCOMVDIST = gsl_spline_alloc(gsl_interp_cspline, n_use);
          ACCEL_INVCOMVDIST = gsl_interp_accel_alloc();
          if (checked_spline_init(SPLINE_INVCOMVDIST, xchi_use, alog_use, n_use, "INVCOMVDIST"))
          {
            /* MPI_Abort is used here because InverseComovingDistance returns double
               and is called from hot loops, so there is no way to propagate an
               int-style error code back to the caller. */
            printf("FATAL on task %d: cannot build inverse comoving distance spline\n", ThisTask);
            fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, 1);
          }

          if (!ThisTask && internal.verbose_level >= VDBG)
          {
            fprintf(stdout,
                    "COSMO: Built %s inverse chi->z over [%.6g, %.6g] Mpc with %d nodes (PLC z [%.6g, %.6g], dz_buf=%.3g)\n",
                    (use_full || alloc_n == n) ? "full" : "tailored",
                    xchi_use[0], xchi_use[n_use - 1], n_use, z_lo, z_hi, dz_buf);
            fflush(stdout);
          }

          if (inject_anchor)
          {
            free(alog_use);
            free(xchi_use);
          }
          free(alog);
          free(xchi);
          INVCOMVDIST_READY = 1;
        }
#ifdef _OPENMP
      }
    }
#endif
  }

  if (INVCOMVDIST_READY != 1)
    return 0.0;

  /* Evaluate log10(a) at chi and convert to z */
  double alog = my_spline_eval(SPLINE_INVCOMVDIST, chi, ACCEL_INVCOMVDIST);
  double a = pow(10., alog);
  if (a <= 0.0)
    return 0.0;
  double z = 1.0 / a - 1.0;
  if (z < 0.0)
    z = 0.0; /* clamp negative round-off */
  return z;
}

double SizeForMass(double m)
{
  /* Radius corresponding to mass m
     m: M_sun

     DIMENSION: Mpc
  */

  switch (WindowFunctionType)
  {
  case 0:
    return pow(m / pow(2.0 * PI, 1.5) / MatterDensity, 0.3333333333333333);
    break;
  case 1:
    return pow(m / (6.0 * PI * PI * MatterDensity), 0.3333333333333333);
    break;
  case 2:
    return pow(m / (4.0 * PI * MatterDensity / 3.0), 0.3333333333333333);
    break;
  default:
    return 0.;
    break;
  }
}

double MassForSize(double size)
{

  switch (WindowFunctionType)
  {
  case 0:
    return MatterDensity * pow(2.0 * PI, 1.5) * pow(size, 3.0);
    break;
  case 1:
    return MatterDensity * 6. * PI * PI * pow(size, 3.0);
    break;
  case 2:
    return MatterDensity * 4. * PI / 3. * pow(size, 3.0);
    break;
  default:
    return 0;
    break;
  }
}

/**********************************/
/* ANALYTIC MASS FUNCTION SECTION */
/**********************************/

#define SQRT2PI ((double)0.39894228)
#define ALPHA ((double)0.569558118758974)

double dOmega_dVariance(double v, double z)
{

  /*
    dOmega/dLambda

    DIMENSIONLESS

    params.AnalyticMassFunction = 0:  Press  & Schechter (1974)
    params.AnalyticMassFunction = 1:  Sheth & Tormen (2001)
    params.AnalyticMassFunction = 2:  Jenkins et al. (2001)
    params.AnalyticMassFunction = 3:  Warren et al. (2006)
    params.AnalyticMassFunction = 4:  Reed et al. (2007)
    params.AnalyticMassFunction = 5:  Crocce et al. (2010)
    params.AnalyticMassFunction = 6:  Tinker et al. (2008)
    params.AnalyticMassFunction = 7:  Courtin et al. (2010)
    params.AnalyticMassFunction = 8:  Angulo et al. (2012)
    params.AnalyticMassFunction = 9:  Watson et al. (2013)
    params.AnalyticMassFunction =10:  Crocce et al. (2010), with forced universality

  */

  double sv, ni, ni2, onepz;

  sv = sqrt(v);
  ni = DELTA_C / sv;

  switch (params.AnalyticMassFunction)
  {
  case 0: // Press & Schechter
    return 2. * exp(-0.5 * ni * ni) * ni * SQRT2PI;
    break;

  case 1: // Sheth & Tormen
    ni2 = sqrt(0.707) * ni;
    return 2. * 0.3222 * SQRT2PI * ni2 * exp(-0.5 * ni2 * ni2) * (1.0 + 1.0 / pow(ni2, 0.6));
    break;

  case 2: // Jenkins et al.
    return 0.315 * exp(-pow(fabs(-log(sv) + 0.61), 3.8));
    break;

  case 3: // Warren et al. (2006)
    return 0.7234 * (pow(sv, -1.625) + 0.2538) * exp(-1.1982 / v);
    break;

  case 4: // Reed et al. (2007)
    ni2 = sqrt(0.707) * ni;
    return 2. * 0.3222 * SQRT2PI * ni2 * exp(-0.54 * ni2 * ni2) *
           (1.0 + 1.0 / pow(ni2, 0.6) + 0.2 * exp(-(pow(-log(sv) - 0.4, 2.0) / 0.72)));
    break;

  case 5: // Crocce et al. (2010)
    onepz = (z < 1 ? 1.0 + z : 2.0);
    return 0.58 * pow(onepz, -0.13) * (pow(sv, -1.37 * pow(onepz, -0.15)) + 0.3 * pow(onepz, -0.084)) * exp(-1.036 * pow(onepz, -0.024) / v);
    break;

  case 6: // Tinker et al. (2010)
    onepz = (z < 2.5 ? 1.0 + z : 3.5);
    return 0.186 * pow(onepz, -0.14) * (pow(2.57 * pow(onepz, -0.569558118758974) / sv, 1.47 * pow(onepz, -0.06)) + 1.) * exp(-1.19 / v);
    break;

  case 7: // Courtin et al. (2010)
    ni2 = sqrt(0.695) * 1.673 / sv;
    return 0.348 * 2. * SQRT2PI * ni2 * (1. + pow(1. / ni2 / ni2, 0.1)) * exp(-ni2 * ni2 / 2.);
    break;

  case 8: // Angulo et al. (2012)
    return 0.201 * (pow(ni * 2.08 / DELTA_C, 1.7) + 1.0) * exp(-1.172 * ni * ni / DELTA_C / DELTA_C);
    break;

  case 9: // Watson et al. (2013)
    return 0.282 * (pow(ni * 1.406 / DELTA_C, 2.163) + 1.0) * exp(-1.210 * ni * ni / DELTA_C / DELTA_C);
    break;

  case 10: // Crocce et al. (2010) universal
    onepz = 1.0;
    return 0.58 * pow(onepz, -0.13) * (pow(sv, -1.37 * pow(onepz, -0.15)) + 0.3 * pow(onepz, -0.084)) * exp(-1.036 * pow(onepz, -0.024) / v);
    break;

  default:
    return 0.0;
    break;
  }
}

double AnalyticMassFunction(double mass, double z)
{
  double r, D;

  r = SizeForMass(mass);
  /* This function still must be adapted to scale-dependent growing mode */
  D = GrowingMode(z, params.k_for_GM);

  return MatterDensity * dOmega_dVariance(MassVariance(r) * D * D, z) * fabs(dMassVariance_dr(r) / 6.0) / mass / mass;
}

double my_spline_eval(gsl_spline *spline, double x, gsl_interp_accel *accel)
{
  /* this function performs linear extrapolation beyond the x-range limits,
     and calls the spline evaluation in between */
  if (x < spline->x[0])
    return spline->y[0] + (x - spline->x[0]) * (spline->y[1] - spline->y[0]) / (spline->x[1] - spline->x[0]);
  else if (x > spline->x[spline->size - 1])
    return spline->y[spline->size - 1] + (x - spline->x[spline->size - 1]) *
                                             (spline->y[spline->size - 1] - spline->y[spline->size - 2]) / (spline->x[spline->size - 1] - spline->x[spline->size - 2]);
  else
    return gsl_spline_eval(spline, x, accel);
}

double my_spline_eval_deriv(gsl_spline *spline, double x, gsl_interp_accel *accel)
{
  /* derivative with linear extrapolation beyond the x-range limits,
     and gsl derivative in between */
  if (x < spline->x[0])
    return (spline->y[1] - spline->y[0]) / (spline->x[1] - spline->x[0]);
  else if (x > spline->x[spline->size - 1])
    return (spline->y[spline->size - 1] - spline->y[spline->size - 2]) / (spline->x[spline->size - 1] - spline->x[spline->size - 2]);
  else
    return gsl_spline_eval_deriv(spline, x, accel);
}
