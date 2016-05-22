
#ifndef _SV2_H_
#define _SV2_H_

// find position on paper (x,y coordinates)

// TODO replace with/comment: https://en.wikipedia.org/wiki/3D_projection

void atompos( double fac,  // factor
              double p[3], // 3d (xyz) coordinates 
              double rad,  // radius 
              double zp[2],// xy position on paper 
              double *zr   // returned (effective) radius 
              );


void atompos2( double fac,  // factor
              double p[3], // 3d (xyz) coordinates 
              double rad,  // radius 
              double zp[2],// xy position on paper 
              double *zr   // returned (effective) radius 
              );


// TEST 2016-05-16 rotation 
void rotmat( int ixyz, double rad, double tmat [3][3] );

void unit_mat ( double tmat [3][3] );


#endif
