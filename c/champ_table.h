/*
 * champ_table: C lookup of optimal CHAMP interface coefficients (order-2 / BDF2
 * scheme) from the precomputed table champ_order2.h5. Same algorithm and results
 * as the MATLAB routines get_r_opt.m / get_p_opt.m (see the README).
 *
 *   theta  = K_L/K_R,  beta = D_L/D_R,  lambda = dt*D_L/h^2
 *   interface condition (1-r)*N[T] + r*D[T]; r = 0 Neumann, r = 1 Dirichlet;
 *   Robin weight p = r/(1-r)  (p = INFINITY at r = 1)
 *
 * Requires the HDF5 C library (1.8 or later).
 */
#ifndef CHAMP_TABLE_H
#define CHAMP_TABLE_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct champ_table champ_table;

/* lookup method requested */
enum { CHAMP_LINEAR = 0, CHAMP_NEAREST = 1 };

/* method actually used for a query (champ_result.method) */
enum {
  CHAMP_M_LINEAR = 0,         /* per-family interpolation                         */
  CHAMP_M_NEAREST = 1,        /* nearest table node (verified there)              */
  CHAMP_M_NEARESTEXTRAP = 2,  /* outside the table range: nearest (boundary) node */
  CHAMP_M_NONE = 3,           /* no stable family at the table node: r = NaN       */
  CHAMP_M_INVALID = 4         /* non-positive or non-finite input: r = NaN         */
};

#define CHAMP_NPEAK 4

typedef struct {
  double rL, rR;                  /* coefficients (NaN: none)                          */
  double kappa[CHAMP_NPEAK];      /* tracked peaks of |A(kappa)| of the chosen family,  */
  double A_re[CHAMP_NPEAK];       /*   largest |A| first, from the nearest table node   */
  double A_im[CHAMP_NPEAK];       /*   (NaN = no peak)                                  */
  int    active[CHAMP_NPEAK];     /* peak in the minimax active set                     */
  double Amax;                    /* verified max|A| over all roots and kappa           */
  double Astar;                   /* tracked minimax value A*                           */
  int    family;                  /* chosen family, 1-based (0 none, -1 invalid input)  */
  const char *familyName;         /* its name ("" if none)                             */
  int    method;                  /* CHAMP_M_*                                          */
} champ_result;

/* open / close the table (NULL on failure; a message is printed to stderr) */
champ_table *champ_table_open( const char *path );
void         champ_table_close( champ_table *T );

/* optimal r_L, r_R (and peaks, max|A|, family); method CHAMP_LINEAR or CHAMP_NEAREST.
   Returns out->method. */
int champ_get_r_opt( const champ_table *T, double theta, double beta, double lambda,
                     int method, champ_result *out );

/* optimal Robin weights p = r/(1-r); out (may be NULL) receives the full result */
int champ_get_p_opt( const champ_table *T, double theta, double beta, double lambda,
                     int method, double *pL, double *pR, champ_result *out );

/* name of a CHAMP_M_* code: "linear", "nearest", "nearestExtrap", "none", "invalid" */
const char *champ_method_name( int m );

/* table information */
int         champ_table_nfamilies( const champ_table *T );
const char *champ_table_family_name( const champ_table *T, int f );   /* f 1-based */

#ifdef __cplusplus
}
#endif
#endif
