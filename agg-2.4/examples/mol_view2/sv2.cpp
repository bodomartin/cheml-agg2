
#include "sv2.h"

#include "sv_math.h"

#include <cmath>

#define  MAXRAD     100.0

// perspective
//TODO RE-ENABLE 2016-05-20 static int pmode = 2;      // 0 : no perspective,  1: pseudo , 2 : true  

static int pmode = 1; // 2016-05-20 disable only for testing  zero will not work!!! TODO

static double dist0 = 15.0; // (initial) distance 

static double dist = 15.0; // distance 

// define a global instance of sv_math (for now)
cheml::gui::sv_math svm;

void atompos( double fac, double p[3], double rad,
              double zp[2], double *zr ) {

  //TEST TEST TEST 2016-05-20
  if (1) {
    // use different projection 
    atompos2( fac, p, rad, zp , zr );
    return;
  }
  
  double y[3], q[3], v1[3], v2[3];
  double xxx, za1, za2, zb1, zb2, a, b;

  if ( pmode == 1 ) {
    // with a factor fac=1, effectivly no perspective .. 
    zp[0] = fac * p[0];
    zp[1] = fac * p[1];
    *zr = fac * rad;
    *zr = MAXRAD;  // ?? 

    if ( dist0 - p[2] > 0 ) { *zr = fac * rad * dist0 / ( dist0 - p[2] ); }

    if ( *zr > MAXRAD ) { *zr = MAXRAD; }

    return;
  }

  svm.vscal( p, 1.0, q );

  q[2] = q[2] - dist;
  svm.vscal( p, 1.0, y );
  xxx = -svm.sp( y, q ) / svm.sp( q, q );
  svm.vsum( y, q, 1.0, xxx, y );

  if ( svm.sp( y, y ) <= 1e-3 ) {
    //bm set axis if y is small 
    y[0] = 1.0; y[1] = 0.0; y[2] = 0.0;
  }

  a = -rad * rad / svm.sp( q, q );

  b = rad * sqrt( ( 1.0 + a ) / svm.sp( y, y ) );

  svm.vsum( q, y, a, b, v1 );
  svm.vsum( q, y, a, -b, v2 );
  svm.vsum( p, v1, 1.0, 1.0, v1 );
  svm.vsum( p, v2, 1.0, 1.0, v2 );

  za1 = fac * v1[0] * dist / ( dist - v1[2] );
  za2 = fac * v1[1] * dist / ( dist - v1[2] );
  zb1 = fac * v2[0] * dist / ( dist - v2[2] );
  zb2 = fac * v2[1] * dist / ( dist - v2[2] );

  // zp: paper xy coordinates 
  zp[0] = 0.5 * ( za1 + zb1 );
  zp[1] = 0.5 * ( za2 + zb2 );

  // zr: radius on paper 
  *zr = ( zb1 - za1 ) * ( zb1 - za1 ) + ( zb2 - za2 ) * ( zb2 - za2 );
  *zr = 0.5 * sqrt( *zr );
}


//----

void atompos2( double fac, double a[3], double rad,
              double zp[2], double *zr ) {


  // https://en.wikipedia.org/wiki/3D_projection

  double c[3] = { 0,0,dist}; // camera position

  double theta[3] = { 0,0,0 }; // orientation of camera ( Tait-Bryan angles)

  // TODO use matrices!
  
  // camera transform
  double d[3];

  // use no camera transform: only a shift
  d[0] = a[0] -c[0];
  d[1] = a[1] -c[1];
  d[2] = a[2] -c[2];

  // viewers position relative to display surface
  // make this the same as the camera ? 
  double e[3];
  e[0] = c[0];
  e[1] = c[1];
  e[2] = c[2];
  
  double b[2];

  b[0] = e[2]/d[2]*d[0]-e[0];
  b[1] = e[2]/d[2]*d[1]-e[1];
  
  
  
  // zp: paper xy coordinates 
  zp[0] = -b[0];
  zp[1] = -b[1];


  *zr = 1.0; // DUMMY VALUE 

}

//----


// SET the rotation matrix 
void rotmat( int ixyz, double rad, double tmat[3][3] ) {

  svm.rotmat( ixyz, rad * 3.141529 / 180.0, tmat );

}

void unit_mat ( double tmat [3][3] ) {

  svm.unit_mat( tmat );

}



