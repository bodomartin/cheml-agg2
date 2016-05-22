#include <cstdlib>
#include <cmath>
#include <cstdio>


// based on a version ( ~ 1995/1997)
// Written long time ago by Jan Labanowski

// use M_PI and M_RAD ??
//#define RAD 57.29577951
#define  PI 3.14159265358979323846

#define RAD 180.0/PI

#include "sv_math.h"

namespace cheml {
  namespace gui {

    double sv_math::distance( double xyz1[3], double xyz2[3] )
//double xyz1[3], xyz2[3];
    {
      int i;
      double s, x;
      s = 0.0;

      for ( i = 0; i < 3; i++ ) {
        x = xyz1[i] - xyz2[i];
        s = s + x * x;
      }

      return ( ( double )sqrt( ( double )s ) );
    }

    /*===============================================================
       function angle finds angle (in deg) between 3 points:

                 xyz1--------xyz2
                            /
                           /
                          /
                         /
                      xyz3
    =================================================================*/


    double sv_math::angle( double xyz1[3], double xyz2[3], double xyz3[3] )
//double xyz1[3], xyz2[3], xyz3[3];
    {
      double d12, d23, d13, s, c, x;

      d12 = distance( xyz1, xyz2 );
      d23 = distance( xyz2, xyz3 );
      d13 = distance( xyz1, xyz3 );

      x = 2.0 * d12 * d23;

      if ( ( x < 0.00001 ) || ( d13 < 0.00001 ) ) {
        return ( 0.0 );
      }

      c = ( d23 * d23 + d12 * d12 - d13 * d13 ) / x;

      if ( 1.0 - c * c < 0.00001 ) {
        return ( 180.0 );
      }

      return ( ( double )( RAD * acos( ( double )c ) ) );
    }

    /*===================================================================
       function vector_length finds length of a vector
    =====================================================================*/
    double sv_math::vector_length( double vector[3] )
//double vector[3];
    {
      double x, y;
      int i;

      x = 0.0;

      for ( i = 0; i < 3; i++ ) {
        y = vector[i];
        x = x + y * y;
      }

      return ( ( double )sqrt( ( double )x ) );
    }


    /*====================================================================

        calculates unit vector from xyz1 to xyz2 and returns distance
        between xyz1 and xyz2 as function value.
    ======================================================================*/

    double sv_math::unit_vector( double xyz1[3], double xyz2[3], double vector[3] )
//double xyz1[3], xyz2[3], vector[3];
    {
      int k;
      double x;

      x = distance( xyz1, xyz2 );

      if ( x < 0.0001 ) {
        return ( 0.0 );
      }

      for ( k = 0; k < 3; k++ ) {
        vector[k] = ( xyz2[k] - xyz1[k] ) / x;
      }

      return ( x );
    }

//      Normalizes a vector and returns its length
    double sv_math::normalize( double vorig[3], double vnorm[3] ) {

      double x = vector_length( vorig );

      if ( x < 0.0001 ) {
        return ( 0.0 );
      }

      for ( unsigned int k = 0; k < 3; k++ ) {
        vnorm[k] = vorig[k] / x;
      }

      //+bm 2015-11-27 !
      return x;
    }

//   calculated scalar product of v1 and v2
    double sv_math::scalar_product( double v1[3], double v2[3] ) {

      double x = 0.0;

      for ( unsigned int k = 0; k < 3; k++ ) {
        x = x + v1[k] * v2[k];
      }

      return x;
    }

//   calculates vector product of v1 and v2 and stores it in v12.
//   Returns length of v12 as function value

    double sv_math::vector_product( double v1[3], double v2[3], double v12[3] ) {
      v12[0] = v1[1] * v2[2] - v1[2] * v2[1];
      v12[1] = v1[2] * v2[0] - v1[0] * v2[2];
      v12[2] = v1[0] * v2[1] - v1[1] * v2[0];
      return vector_length( v12 );
    }

//   Function returns torsion angle between xyz1--xyz2
//                                                  \
//                                                   xyz3--xyz4
// in degrees!!

    double sv_math::torsion( double xyz1[3], double xyz2[3],
                            double xyz3[3], double xyz4[3] ) {
      int i;
      double  v12[3], v23[3], v34[3], plane123[3],
             plane234[3], v1234[3], c, s, x, y;

      if ( ( unit_vector( xyz1, xyz2, v12 ) <= 0.0001 ) ||
           ( unit_vector( xyz2, xyz3, v23 ) <= 0.0001 ) ||
           ( unit_vector( xyz3, xyz4, v34 ) <= 0.0001 ) ) {
        return 0.0 ;
      }

      if ( ( vector_product( v12, v23, plane123 ) <= 0.0001 ) ||
           ( vector_product( v23, v34, plane234 ) <= 0.0001 ) ) {
        return 0.0;
      }

      x = normalize( plane123, plane123 );
      y = normalize( plane234, plane234 );
      c = scalar_product( plane123, plane234 );

      if ( fabs( c + 1.0 ) < 0.0001 ) {
        return 180.0;
      }
      else if ( fabs( c - 1.0 ) < 0.0001 ) {
        return 0.0;
      }

      s = vector_product( plane123, plane234, v1234 );

      s = normalize( v1234, v1234 );

      s = s * scalar_product( v1234, v23 );

      return ( ( double )( RAD * atan( ( double )s / c ) ) );
    }

    /*
        multiply vector with matrix
        erster Index : Zeile, zweiter Index : Spalte
     */

    void sv_math::vecmat( double v[3], double m[3][3] ) {
      double r[3];

      r[0] = m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2];
      r[1] = m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2];
      r[2] = m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2];

      v[0] = r[0];
      v[1] = r[1];
      v[2] = r[2];
    }

    // multiply matrix with matrix Hainzl S. 195
    void sv_math::matmat( double a[3][3], double b[3][3], double c[3][3] ) {

      for ( unsigned int i = 0; i < 3; i++ ) {
        for ( unsigned int j = 0; j < 3; j++ ) {
          c[i][j] = 0.0;

          for ( unsigned int n = 0; n < 3; n++ ) {
            c[i][j] += a[i][n] * b[n][j];
          }
        }
      }
    }

    // berechnet die determinante einer 3x3 matrix Hainzel S.203
    double sv_math::det( double m[3][3] ) {
      double a, b, c;

      a = m[0][0] * ( m[1][1] * m[2][2] - m[2][1] * m[1][2] );
      b = m[0][1] * ( m[1][0] * m[1][2] - m[2][1] * m[1][2] );
      c = m[0][2] * ( m[1][0] * m[2][1] - m[2][0] * m[1][1] );

      return ( a - b + c );
    }

    /* löscht 3x3 Matrix */
    void sv_math::nullify( double m[3][3] )
//double m[3][3];
    {
      m[0][0] = 0;
      m[0][1] = 0;
      m[0][2] = 0;
      m[1][0] = 0;
      m[1][1] = 0;
      m[1][2] = 0;
      m[2][0] = 0;
      m[2][1] = 0;
      m[2][2] = 0;
    }

    /* berechnet die inverse einer 3x3 Matrix Hainzl S 208*/
    void sv_math::invers( double m[3][3], double n[3][3] )
//double m[3][3], n[3][3];
    {
      double dt, dt2 ;  /* determinanten */
      int a, b, c, d, i, j;

      nullify( n );

      dt = det( m );

      if ( dt == 0 ) { return; }

      for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {

          switch ( j ) {

            case 0:
              a = 1;
              b = 2;
              break;

            case 1:
              a = 0;
              b = 2;
              break;

            case 2:
              a = 0;
              b = 1;
              break;
          }

          switch ( i ) {

            case 0:
              c = 1;
              d = 2;
              break;

            case 1:
              c = 0;
              d = 2;
              break;

            case 2:
              c = 0;
              d = 1;
              break;
          }

          dt2 = m[a][c] * m[b][d] - m[a][d] * m[b][c];

          if ( fabs( dt ) > 0.0001 ) {
            n[i][j] = pow( -1, i + 1 + j + 1 ) * dt2 / dt;
          }

        }
      }
    }


    void sv_math::unit_mat( double m[3][3] )
//double m[3][3];
    {
      m[0][0] = 1;
      m[0][1] = 0;
      m[0][2] = 0;
      m[1][0] = 0;
      m[1][1] = 1;
      m[1][2] = 0;
      m[2][0] = 0;
      m[2][1] = 0;
      m[2][2] = 1;
    }

    /*
         add to rotation matrix tmat alpha in radians!
     */
    void sv_math::rotmat( int ixyz,  // TODO USE ENUM !!
                          double alfa, double tmat[3][3] ) {

      int i, j, k;
      double rot[3][3], w[3][3];

      switch ( ixyz ) {

        case 3:  /* in xy-Ebene rotieren (um die z-Achse) */
          rot[0][0] = cos( alfa );
          rot[0][1] = -sin( alfa );
          rot[0][2] = 0.0;
          rot[1][0] = sin( alfa );
          rot[1][1] = cos( alfa );
          rot[1][2] = 0.0;
          rot[2][0] = 0.0;
          rot[2][1] = 0.0;
          rot[2][2] = 1.0;
          break;

        case 1:  /* in  xz-Ebene rotieren (um die y-Achse) */
          rot[0][0] = cos( alfa );
          rot[0][1] = 0.0;
          rot[0][2] = sin( alfa );
          rot[1][0] = 0.0;
          rot[1][1] = 1.0;
          rot[1][2] = 0.0;
          rot[2][0] = -sin( alfa );
          rot[2][1] = 0.0;
          rot[2][2] = cos( alfa );
          break;

        case 2:  /* in  yz-Ebene rotieren (um die x-Achse) */
          rot[0][0] = 1.0;
          rot[0][1] = 0.0;
          rot[0][2] = 0.0;
          rot[1][0] = 0.0;
          rot[1][1] = cos( alfa );
          rot[1][2] = -sin( alfa );
          rot[2][0] = 0.0;
          rot[2][1] = sin( alfa );
          rot[2][2] = cos( alfa );
          break;

      default: {
        // ERROR (TODO)
        return;

      }
      }

      // matrix multipy (multiply rot with tmat)
      for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
          w[i][j] = 0.0;
          for ( k = 0; k < 3; k++ ) {
            w[i][j] = w[i][j] + rot[i][k] * tmat[k][j];
          }
        }
      }

      // assign result to tmat
      for ( i = 0; i < 3; i++ ) for ( j = 0; j < 3; j++ ) { tmat[i][j] = w[i][j]; }

      return;
    }

    /*
        set up euler angles rotation matrix
    */
    void sv_math::eumat( double alfa, double beta, double gama, double tmat[3][3] )
//double alfa, beta, gama, tmat[3][3];
    {
      int i, j, k;
      double gam[3][3], bet[3][3], alf[3][3], w[3][3];
      /*  printf("Euler angles  %7.1f %7.1f %7.1f\n", alfa,beta,gama);*/

      gam[0][0] = 1.0;
      gam[0][1] = 0.0;
      gam[0][2] = 0.0;
      gam[1][0] = 0.0;
      gam[1][1] = cos( gama );
      gam[1][2] = -sin( gama );
      gam[2][0] = 0.0;
      gam[2][1] = sin( gama );
      gam[2][2] = 1.0;

      bet[0][0] = cos( beta );
      bet[0][1] = 0.0;
      bet[0][2] = sin( beta );
      bet[1][0] = 0.0;
      bet[1][1] = 1.0;
      bet[1][2] = 0.0;
      bet[2][0] = -sin( beta );
      bet[2][1] = 0.0;
      bet[2][2] = cos( beta );

      alf[0][0] = cos( alfa );
      alf[0][1] = -sin( alfa );
      alf[0][2] = 0.0;
      alf[1][0] = sin( alfa );
      alf[1][1] = cos( alfa );
      alf[1][2] = 0.0;
      alf[2][0] = 0.0;
      alf[2][1] = 0.0;
      alf[2][2] = 1.0;

      for ( i = 0; i < 3; i++ ) for ( j = 0; j < 3; j++ ) {
          w[i][j] = 0.0;

          for ( k = 0; k < 3; k++ ) { w[i][j] = w[i][j] + bet[i][k] * alf[k][j]; }
        }

      for ( i = 0; i < 3; i++ ) for ( j = 0; j < 3; j++ ) {
          tmat[i][j] = 0.0;

          for ( k = 0; k < 3; k++ ) { tmat[i][j] = tmat[i][j] + gam[i][k] * w[k][j]; }
        }

      return;
    }

    /* ======================


       return x,y,z position from three starting points and
       angle, distance, dihedral (i.e. z-matrix to xyz
       conversion)

       ======================= */

    void sv_math::zmat2xyz( double p1[3], double p2[3], double p3[3], double dist, double ang, double dihedral, double p4[3] )
//double p1[3], p2[3], p3[3], dist, ang, dihedral, p4[3];
    {
      double a1[3], a2[3], a3[3], x[3], z[3], g, xx[3];
      double tmat[3][3], tinv[3][3], t[3][3];
      int i, j;

      ang = ang * PI / 180.0;
      dihedral = dihedral * PI / 180.0;

      /* translation of a2 to origin */

      for ( i = 0; i < 3 ; i++ ) {
        a1[i] = 0.0;
        a2[i] = p2[i] - p1[i];
        a3[i] = p3[i] - p1[i];
      }

      z[0] = 0;

      z[1] = 0;
      z[2] = 1;
      x[0] = 1;
      x[1] = 0;
      x[2] = 0;

      unit_mat( tmat );

      /* Rotationsmatrix aufbauen um p2 auf die z-Achse zu drehen */
      g = angle( a2, a1, x );
      rotmat( 3, -g * PI / 180.0, tmat );

      g = angle( a2, a1, z );
      rotmat( 1, -g * PI / 180.0, tmat );

      vecmat( a2, tmat );
      vecmat( a3, tmat );

      /* a3 nach an x-Achse ausrichten */

      g = angle( a3, a1, x );
      rotmat( 1, -g * PI / 180.0, tmat );
      vecmat( a3, tmat );

      /* neuer Punkt: Polarkoordinaten -> kathesische Koordinaten */
      p4[0] = dist * cos( ang ) * cos( dihedral );
      p4[1] = dist * cos( ang ) * sin( dihedral );
      p4[2] = dist * sin( dihedral );

      /* und ins globale Koordinaten zurücktransformieren */

      invers( tmat, tinv );
      vecmat( p4, tinv );

      /* und translatieren zum Ansetzpunkt */
      p4[0] += p1[0];
      p4[1] += p1[1];
      p4[2] += p1[2];


      matmat( tmat, tinv, t );

    }


    /* vector product of a and b ,returned in c */
    void sv_math::cross( double a[3], double b[3], double c[3] ) {
      c[0] = a[1] * b[2] - a[2] * b[1];
      c[1] = a[2] * b[0] - a[0] * b[2];
      c[2] = a[0] * b[1] - a[1] * b[0];
      return;
    }

    /* dot product of vectors a and b */
    double sv_math::sp( double a[3], double b[3] ) {
      double sp1;
      sp1 = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
      return sp1;
    }

    /* scale vector a by ca */
    void sv_math::vscal( double a[3], double ca, double v[3] ) {
      v[0] = a[0] * ca;
      v[1] = a[1] * ca;
      v[2] = a[2] * ca;
      return;
    }

    /* add vectors a and b scaled by ca and cb */
    void sv_math::vsum( double a[3], double b[3], double ca, double cb, double v[3] ) {
      v[0] = ca * a[0] + cb * b[0];
      v[1] = ca * a[1] + cb * b[1];
      v[2] = ca * a[2] + cb * b[2];
      return;
    }


  }
}






