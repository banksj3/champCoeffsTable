/* champ_example: look up optimal CHAMP coefficients.
 *
 *   champ_example [table.h5] theta beta lambda [linear|nearest]
 *   champ_example [table.h5] -          < queries.txt     (lines "theta beta lambda method")
 */
#include "champ_table.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void show( const champ_table *T, double th, double be, double la, int method, int verbose )
{
  champ_result r; double pL, pR;
  champ_get_p_opt( T,th,be,la,method,&pL,&pR,&r );
  if( !verbose ){
    printf( "%.17g %.17g %.17g %.17g %d %d",r.rL,r.rR,r.Amax,r.Astar,r.family,r.method );
    for( int p = 0; p < CHAMP_NPEAK; p++ ) printf( " %.17g %.17g %.17g %d",r.kappa[p],r.A_re[p],r.A_im[p],r.active[p] );
    printf( "\n" );
    return;
  }
  printf( "theta=%g beta=%g lambda=%g (%s)\n",th,be,la,champ_method_name(r.method) );
  printf( "  rL=%.6g rR=%.6g   pL=%.6g pR=%.6g\n",r.rL,r.rR,pL,pR );
  printf( "  family %s, verified max|A| = %.5f, minimax A* = %.5f\n",r.familyName,r.Amax,r.Astar );
  for( int p = 0; p < CHAMP_NPEAK; p++ ){
    if( isnan(r.kappa[p]) ) continue;
    printf( "  peak %d: kappa=%.4g  A=%.5f%+.5fi  |A|=%.5f%s\n",p+1,r.kappa[p],r.A_re[p],r.A_im[p],
            hypot(r.A_re[p],r.A_im[p]),r.active[p] ? "  (active)" : "" );
  }
}

int main( int argc, char **argv )
{
  const char *path = "champ_order2.h5"; int a = 1;
  if( argc > 1 && strstr(argv[1],".h5") ){ path = argv[1]; a = 2; }
  champ_table *T = champ_table_open( path );
  if( !T ) return 1;
  if( argc > a && strcmp(argv[a],"-") == 0 ){
    double th, be, la; char m[32];
    while( scanf( "%lf %lf %lf %31s",&th,&be,&la,m ) == 4 ) show( T,th,be,la,strcmp(m,"nearest") ? CHAMP_LINEAR : CHAMP_NEAREST,0 );
  } else if( argc >= a+3 ){
    int method = (argc > a+3 && strcmp(argv[a+3],"nearest") == 0) ? CHAMP_NEAREST : CHAMP_LINEAR;
    show( T,atof(argv[a]),atof(argv[a+1]),atof(argv[a+2]),method,1 );
  } else {
    fprintf( stderr,"usage: %s [table.h5] theta beta lambda [linear|nearest]\n",argv[0] );
  }
  champ_table_close( T );
  return 0;
}
