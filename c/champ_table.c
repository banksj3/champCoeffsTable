/* champ_table.c -- see champ_table.h */
#include "champ_table.h"
#include <hdf5.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXFAM 8

typedef struct {
  char    name[32];
  double *rL, *rR, *Atr, *Amax, *status;
  double *pk, *pre, *pim, *pact;   /* CHAMP_NPEAK per node */
  unsigned char *stable;
} family;

struct champ_table {
  int     n[3];          /* theta, beta, lambda */
  double *g[3];          /* log10 grids */
  double *best;
  int     nF;
  family  F[MAXFAM];
};

/* ------------------------------------------------------------------ reading */
static double *readDataset( hid_t file, const char *name, size_t expect, int *ok )
{
  double *buf = NULL;
  hid_t d = H5Dopen2( file,name,H5P_DEFAULT );
  if( d < 0 ){ fprintf( stderr,"champ_table: cannot open dataset %s\n",name ); *ok = 0; return NULL; }
  hid_t s = H5Dget_space( d );
  hssize_t np = H5Sget_simple_extent_npoints( s );
  if( np <= 0 || (expect && (size_t)np != expect) ){
    fprintf( stderr,"champ_table: dataset %s has %lld values (expected %zu)\n",name,(long long)np,expect );
    *ok = 0;
  } else {
    buf = (double *)malloc( (size_t)np*sizeof(double) );
    if( H5Dread( d,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,buf ) < 0 ){
      fprintf( stderr,"champ_table: cannot read %s\n",name ); free( buf ); buf = NULL; *ok = 0;
    }
  }
  H5Sclose( s ); H5Dclose( d );
  return buf;
}

static size_t datasetSize( hid_t file, const char *name )
{
  hid_t d = H5Dopen2( file,name,H5P_DEFAULT );
  if( d < 0 ) return 0;
  hid_t s = H5Dget_space( d );
  hssize_t np = H5Sget_simple_extent_npoints( s );
  H5Sclose( s ); H5Dclose( d );
  return np > 0 ? (size_t)np : 0;
}

static int readStringAttr( hid_t file, const char *name, char *out, size_t len )
{
  hid_t a = H5Aopen( file,name,H5P_DEFAULT );
  if( a < 0 ) return 0;
  hid_t t = H5Aget_type( a ); int ok = 1;
  if( H5Tis_variable_str(t) > 0 ){
    char *s = NULL; hid_t mt = H5Tcopy( H5T_C_S1 ); H5Tset_size( mt,H5T_VARIABLE );
    if( H5Aread( a,mt,&s ) < 0 || !s ) ok = 0;
    else { strncpy( out,s,len-1 ); out[len-1] = 0; free( s ); }
    H5Tclose( mt );
  } else {
    size_t sz = H5Tget_size( t ); char *s = (char *)calloc( sz+1,1 );
    hid_t mt = H5Tcopy( H5T_C_S1 ); H5Tset_size( mt,sz );
    if( H5Aread( a,mt,s ) < 0 ) ok = 0;
    else { strncpy( out,s,len-1 ); out[len-1] = 0; }
    free( s ); H5Tclose( mt );
  }
  H5Tclose( t ); H5Aclose( a );
  return ok;
}

champ_table *champ_table_open( const char *path )
{
  H5Eset_auto2( H5E_DEFAULT,NULL,NULL );
  hid_t file = H5Fopen( path,H5F_ACC_RDONLY,H5P_DEFAULT );
  if( file < 0 ){ fprintf( stderr,"champ_table: cannot open %s\n",path ); return NULL; }
  champ_table *T = (champ_table *)calloc( 1,sizeof(champ_table) );
  int ok = 1;
  const char *gname[3] = { "/theta","/beta","/lambda" };
  for( int d = 0; d < 3; d++ ){
    T->n[d] = (int)datasetSize( file,gname[d] );
    T->g[d] = readDataset( file,gname[d],0,&ok );
    if( !ok ) break;
    for( int i = 0; i < T->n[d]; i++ ) T->g[d][i] = log10( T->g[d][i] );
  }
  size_t N = ok ? (size_t)T->n[0]*T->n[1]*T->n[2] : 0;
  if( ok ) T->best = readDataset( file,"/best",N,&ok );
  char fams[256] = "";
  if( ok && !readStringAttr( file,"families",fams,sizeof(fams) ) ){
    fprintf( stderr,"champ_table: missing attribute 'families'\n" ); ok = 0;
  }
  if( ok ){
    char *save = NULL;
    for( char *tok = strtok_r( fams,",",&save ); tok && T->nF < MAXFAM; tok = strtok_r( NULL,",",&save ) ){
      family *B = &T->F[T->nF++]; char nm[96];
      strncpy( B->name,tok,sizeof(B->name)-1 );
#define RD(field,ds,cnt) snprintf( nm,sizeof(nm),"/%s/%s",B->name,ds ); B->field = readDataset( file,nm,cnt,&ok ); if( !ok ) break;
      RD( rL,"rL",N ) RD( rR,"rR",N ) RD( Atr,"Atracked",N ) RD( Amax,"Amax",N ) RD( status,"status",N )
      RD( pk,"peakKappa",CHAMP_NPEAK*N ) RD( pre,"peakAreal",CHAMP_NPEAK*N ) RD( pim,"peakAimag",CHAMP_NPEAK*N )
      RD( pact,"peakActive",CHAMP_NPEAK*N )
#undef RD
      B->stable = (unsigned char *)malloc( N );
      for( size_t i = 0; i < N; i++ ) B->stable[i] = B->status[i] > 0 && B->Amax[i] <= 1+1e-6;
    }
  }
  H5Fclose( file );
  if( !ok ){ champ_table_close( T ); return NULL; }
  return T;
}

void champ_table_close( champ_table *T )
{
  if( !T ) return;
  for( int d = 0; d < 3; d++ ) free( T->g[d] );
  free( T->best );
  for( int f = 0; f < T->nF; f++ ){
    family *B = &T->F[f];
    free( B->rL ); free( B->rR ); free( B->Atr ); free( B->Amax ); free( B->status );
    free( B->pk ); free( B->pre ); free( B->pim ); free( B->pact ); free( B->stable );
  }
  free( T );
}

int champ_table_nfamilies( const champ_table *T ){ return T->nF; }
const char *champ_table_family_name( const champ_table *T, int f ){ return (f >= 1 && f <= T->nF) ? T->F[f-1].name : ""; }

const char *champ_method_name( int m )
{
  static const char *nm[] = { "linear","nearest","nearestExtrap","none","invalid" };
  return (m >= 0 && m <= 4) ? nm[m] : "?";
}

/* ------------------------------------------------------------------ lookup */
static double clampd( double v, double a, double b ){ return v < a ? a : (v > b ? b : v); }

int champ_get_r_opt( const champ_table *T, double theta, double beta, double lambda,
                     int method, champ_result *out )
{
  const double X[3] = { theta,beta,lambda };
  int valid = 1;
  for( int d = 0; d < 3; d++ ) valid = valid && isfinite(X[d]) && X[d] > 0;
  double q[3]; int inR = valid;
  for( int d = 0; d < 3; d++ ){
    q[d] = valid ? log10( X[d] ) : 0;
    inR = inR && q[d] >= T->g[d][0]-1e-12 && q[d] <= T->g[d][T->n[d]-1]+1e-12;
  }
  /* nearest node */
  int id[3];
  for( int d = 0; d < 3; d++ ){
    double dg = T->g[d][1]-T->g[d][0];
    double s = round( (q[d]-T->g[d][0])/dg );
    id[d] = (int)clampd( s,0,T->n[d]-1 );
  }
  const size_t n0 = T->n[0], n01 = (size_t)T->n[0]*T->n[1];
  size_t nn = id[0] + n0*id[1] + n01*id[2];
  int bestN = valid ? (int)T->best[nn] : 0;

  double rL = NAN, rR = NAN, Am = NAN, As = NAN; int fam = 0, how = CHAMP_M_NEAREST;

  if( method == CHAMP_LINEAR ){
    int lo[3]; double w[3];
    for( int d = 0; d < 3; d++ ){
      double dg = T->g[d][1]-T->g[d][0], s = (q[d]-T->g[d][0])/dg;
      lo[d] = (int)clampd( floor(s),0,T->n[d]-2 );
      w[d]  = clampd( s-lo[d],0,1 );
    }
    double RL[MAXFAM], RR[MAXFAM], AT[MAXFAM], AM[MAXFAM]; int ok[MAXFAM];
    for( int f = 0; f < T->nF; f++ ){
      const family *B = &T->F[f]; double acc[4] = {0,0,0,0}; int okf = 1;
      for( int c = 0; c < 8; c++ ){
        int b0 = c & 1, b1 = (c>>1) & 1, b2 = (c>>2) & 1;
        double wc = (b0 ? w[0] : 1-w[0])*(b1 ? w[1] : 1-w[1])*(b2 ? w[2] : 1-w[2]);
        if( wc < 1e-12 ) continue;                     /* corners with zero weight do not count */
        size_t li = (lo[0]+b0) + n0*(lo[1]+b1) + n01*(lo[2]+b2);
        okf = okf && B->stable[li];
        if( !B->stable[li] ) continue;
        acc[0] += wc*log( B->rL[li]+0.01 ); acc[1] += wc*log( B->rR[li]+0.01 );
        acc[2] += wc*B->Atr[li];            acc[3] += wc*B->Amax[li];
      }
      ok[f] = okf && inR;
      RL[f] = clampd( exp(acc[0])-0.01,0,1 ); RR[f] = clampd( exp(acc[1])-0.01,0,1 );
      AT[f] = acc[2]; AM[f] = acc[3];
    }
    /* selection rule (smallest A*, cap, collapsed case, then smallest verified max|A|) */
    int qmin = -1; double Amin = INFINITY;
    for( int f = 0; f < T->nF; f++ ) if( ok[f] && AT[f] < Amin ){ Amin = AT[f]; qmin = f; }
    int sel = 0;
    if( qmin >= 0 ){
      double cap = fmax( 1.05*Amin,Amin+0.02 );
      if( Amin < 0.1*AM[qmin] ) cap = fmax( cap,0.1*AM[qmin] );
      double mbest = INFINITY;
      for( int f = 0; f < T->nF; f++ ) if( ok[f] && AT[f] <= cap && AM[f] < mbest ){ mbest = AM[f]; sel = f+1; }
    }
    if( bestN > 0 && ok[bestN-1] ) sel = bestN;          /* the nearest node's family, where valid */
    if( sel > 0 ){
      rL = RL[sel-1]; rR = RR[sel-1]; Am = AM[sel-1]; As = AT[sel-1]; fam = sel; how = CHAMP_M_LINEAR;
    }
  }
  if( isnan(rL) && bestN > 0 ){                          /* nearest-node values */
    const family *B = &T->F[bestN-1];
    rL = B->rL[nn]; rR = B->rR[nn]; Am = B->Amax[nn]; As = B->Atr[nn]; fam = bestN; how = CHAMP_M_NEAREST;
  }
  if( isnan(rL) ){ fam = 0; Am = NAN; As = NAN; }
  if( valid && !inR ) how = CHAMP_M_NEARESTEXTRAP;
  if( valid && fam == 0 ) how = CHAMP_M_NONE;
  if( !valid ){ how = CHAMP_M_INVALID; fam = -1; }

  out->rL = rL; out->rR = rR; out->Amax = Am; out->Astar = As; out->family = fam; out->method = how;
  out->familyName = fam > 0 ? T->F[fam-1].name : "";
  for( int p = 0; p < CHAMP_NPEAK; p++ ){
    if( fam > 0 ){
      const family *B = &T->F[fam-1]; size_t i = (size_t)p + CHAMP_NPEAK*nn;
      out->kappa[p] = B->pk[i]; out->A_re[p] = B->pre[i]; out->A_im[p] = B->pim[i]; out->active[p] = B->pact[i] > 0;
    } else {
      out->kappa[p] = out->A_re[p] = out->A_im[p] = NAN; out->active[p] = 0;
    }
  }
  return how;
}

int champ_get_p_opt( const champ_table *T, double theta, double beta, double lambda,
                     int method, double *pL, double *pR, champ_result *out )
{
  champ_result r; if( !out ) out = &r;
  int m = champ_get_r_opt( T,theta,beta,lambda,method,out );
  *pL = out->rL >= 1 ? INFINITY : out->rL/(1-out->rL);
  *pR = out->rR >= 1 ? INFINITY : out->rR/(1-out->rR);
  return m;
}
