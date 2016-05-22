
// simple viewer math

namespace cheml {

  namespace gui {

    // sv = simple viewer

    class sv_math {

      public:

        sv_math() { };

        double distance( double xyz1[3], double xyz2[3] );
        double angle( double xyz1[3], double xyz2[3], double xyz3[3] );
        double vector_length( double vector[3] );
        double unit_vector( double xyz1[3], double xyz2[3], double vector[3] );
        double normalize( double vorig[3], double vnorm[3] );
        double scalar_product( double v1[3], double v2[3] );
        double vector_product( double v1[3], double v2[3], double v12[3] );
        double torsion( double xyz1[3], double xyz2[3], double xyz3[3], double xyz4[3] );
        void vecmat( double v[3], double m[3][3] );
        void matmat( double a[3][3], double b[3][3], double c[3][3] );
        double det( double m[3][3] );
        void nullify( double m[3][3] );
        void invers( double m[3][3], double n[3][3] );
        void unit_mat( double m[3][3] );
        void rotmat( int ixyz, double alfa, double tmat[3][3] );
        void eumat( double alfa, double beta, double gama, double tmat[3][3] );
        void zmat2xyz( double p1[3], double p2[3], double p3[3], double dist, double ang, double dihedral, double p4[3] );

        void cross( double a[3], double b[3], double c[3] );
        double sp( double a[3], double b[3] );
        void vscal( double a[3], double ca, double v[3] );
        void vsum( double a[3], double b[3], double ca, double cb, double v[3] );



    }; // sv_math

  }

}
