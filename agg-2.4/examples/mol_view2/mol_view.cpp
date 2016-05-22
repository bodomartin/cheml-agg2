#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "agg_rendering_buffer.h"

// anti-aliasing rasterizer 
#include "agg_rasterizer_scanline_aa.h"


#include "agg_scanline_p.h"
#include "agg_renderer_scanline.h"
#include "agg_conv_transform.h"
#include "agg_conv_stroke.h"
#include "agg_bspline.h"
#include "agg_ellipse.h"
#include "agg_gsv_text.h"

#include "ctrl/agg_slider_ctrl.h"
#include "ctrl/agg_scale_ctrl.h"

#include "platform/agg_platform_support.h"

// OK 2016-05-18 what does the define do ??
// this is probably at the wrong place now, should be in 'molecule.h' ???  TODO!

/// 2016-05-19 this define is needed for the 'pixel_formats.h' below !!!

#define AGG_BGR24

//#define AGG_RGB24
//#define AGG_BGRA32 
//#define AGG_RGBA32 
//#define AGG_ARGB32 
//#define AGG_ABGR32
//#define AGG_BGR96
//#define AGG_BGRA128
//#define AGG_RGB565
//#define AGG_RGB555

#include "pixel_formats.h"

#include "molecule.h"

#include "sv2.h"

#include "bond_vertex_generator.h"

// DEBUG
#include <iostream>

#include <sstream>

enum default_num_points_e { default_num_points = 20000 };

enum start_size_e
{
    start_width  = 400,
    start_height = 400
};

      


// for object id to color mapping we want the _exact_ value, i.e. we need the rgb8 and _not_ the srgb8 color  
agg::rgba8 object_id_to_rgba8(unsigned long id) {

  unsigned int b = id >> 16;
  unsigned int g = (id >> 8) - (b << 8);
  unsigned int r = id - (b << 16 ) - ( g << 8 );

  //  std::cout << " object_id_to_srgba8 id = " << id << " r = " << r << " g = " << g << " b = " << b << std::endl; 
  
  return agg::rgba8(r,g,b);  // 10: transparency
    
}


// sort atoms back to front, return index vector 
std::vector<unsigned int> sort_atoms(molecule & m) { 

  std::vector<unsigned int> idx;
  idx.resize(m.num_atoms());

  std::vector<bool> f;
  f.resize(m.num_atoms());

  // probably not needed 
  for ( unsigned int k = 0; k < m.num_atoms(); k++ ) {
    f[k] = false;
  }
  
  for ( unsigned int n = 0; n < m.num_atoms(); n++ ) {

    double bot = 1.0e10;
    unsigned int ibot = 0;

    for ( unsigned int k = 0; k < m.num_atoms() ; k++ )
      if ( !f[k] && ( m.atom(k).z < bot ) ) {
        bot = m.atom(k).z ; ibot = k;
      }

    idx[n] = ibot;

    f[ibot] = true; // processed atom 
  }

  return idx;
}



//--------------------------------------------------------------------------------

class the_application : public agg::platform_support
{   
    molecule*                    m_molecules;
    unsigned                     m_num_molecules;
    unsigned                     m_cur_molecule;
    agg::slider_ctrl<color_type> m_thickness;
    agg::slider_ctrl<color_type> m_text_size;
    double     m_pdx;
    double     m_pdy;
    double     m_center_x;
    double     m_center_y;
    double     m_scale;
    double     m_prev_scale;
    double     m_angle;
    double     m_prev_angle;
    bool       m_mouse_move;
  agg::srgba8 m_atom_colors[molecule::end_of_atom_colors];


  double m_current_affine_x;
  double m_current_affine_y;

  int m_current_raw_x;
  int m_current_raw_y;

  bool m_picking_mode;

  // TEST a position (set by mouse)
  double mousepos[3] = {0,0,0};   

  bool m_recenter;
  
public:
    virtual ~the_application()
    {
        delete [] m_molecules;
    }

    the_application(agg::pix_format_e format, bool flip_y, const char* fname) :
        agg::platform_support(format, flip_y),
        m_molecules(new molecule [100]),  // TODO change this would allocate space for up to 100 molecules 
        m_num_molecules(0),
        m_cur_molecule(0),
        m_thickness(5, 5,  start_width-5, 12),
        m_text_size(5, 20, start_width-5, 27),
        m_pdx(0.0),
        m_pdy(0.0),
        m_center_x(start_width / 2),
        m_center_y(start_height / 2),
        m_scale(1.0),
        m_prev_scale(1.0),
        m_angle(0.0),
        m_prev_angle(0.0),
        m_mouse_move(false),

        m_current_affine_x(0),
        m_current_affine_y(0),
        m_current_raw_x(0),
        m_current_raw_y(0),

        // TEST ... 

        m_picking_mode(false)

        
  {
        m_thickness.label("Thickness=%3.2f");
        m_text_size.label("Label Size=%3.2f");

        add_ctrl(m_thickness);
        add_ctrl(m_text_size);

        FILE* fd = fopen(full_file_name(fname), "r");
        if(fd)
        {
            unsigned i;
            for(i = 0; i < 100; i++)
            {
                if(!m_molecules[m_num_molecules].read(fd)) break;
                m_num_molecules++;
            }
            fclose(fd);
        }
        else
        {
            char buf[256];
            sprintf(buf, "File not found: '%s'.\n",
                    fname, fname);
            message(buf);
        }

        // TODO replace 
        memset(m_atom_colors, 0, sizeof(agg::srgba8) * molecule::end_of_atom_colors);
        m_atom_colors[molecule::atom_color_general] = agg::srgba8(0,0,0);
        m_atom_colors[molecule::atom_color_N]       = agg::srgba8(0,0,120);
        m_atom_colors[molecule::atom_color_O]       = agg::srgba8(200,0,0);
        m_atom_colors[molecule::atom_color_S]       = agg::srgba8(120,120,0);
        m_atom_colors[molecule::atom_color_P]       = agg::srgba8(80,50,0);
        m_atom_colors[molecule::atom_color_halogen] = agg::srgba8(0,200,0);
    }


    virtual void on_init()
    {
    }

  // overloaded function
  virtual void on_draw() {

    on_draw_bm();  // call code from original demo 
      }
      
  //----------------------------------------


  //----------------------------------------

  
  // 3d transform 
  void tmat_transform( const double tmat[3][3], // double[3][3]
                       molecule & m
                       ) {

    for ( unsigned int n = 0; n < m.num_atoms(); n++ ) {

      double p[3];      

      for ( unsigned int i = 0; i < 3; i++ ) {
        
        p[i]
          = tmat[i][0] * m.atom(n).x
          + tmat[i][1] * m.atom(n).y
          + tmat[i][2] * m.atom(n).z;
      }
      // overwrite atom positions
      m.set_atom_xyz(n,p[0],p[1],p[2]);
    }
  }

  void tmat_transform_pos( const double tmat[3][3], 
                           double in[3] , double out[3] // result 
                           ) {

      for ( unsigned int i = 0; i < 3; i++ ) {
        
        out[i]
          = tmat[i][0] * in[0]
          + tmat[i][1] * in[1]
          + tmat[i][2] * in[2];
      }

  }

// find molecule / atom nearest to the 2x paper position, return z coordinate
double z_of_nearest_atom(unsigned int mol_idx, double pos[2] ) {

  double dist = 1e50;

  // get the current mol
  molecule & m = m_molecules[mol_idx];

  unsigned int a = 0;
  
  for ( unsigned int i = 0; i < m.num_atoms(); i++ ) {

    double d = sqrt( (m.atom(i).px - pos[0]) * (m.atom(i).px - pos[0]) +
                     (m.atom(i).py - pos[1]) * (m.atom(i).py - pos[1])
                     ); 
    
    if ( d < dist) {
      dist = d;
      a = i;
    }   
  }

  return m.atom(a).z;  // returns the 3D coordinate of nearest atom 
}

  
  
  //----------------------------------------
  
    virtual void on_draw_bm()
    {
        double width = initial_width();
        double height = initial_height();

        
        // antialiased renderer
        agg::rasterizer_scanline_aa<> ras;

        agg::scanline_p8 sl;

        // pixfmt is defined in 'pixel_formats.h'
        typedef agg::renderer_base<pixfmt> renderer_base;

        // anti-aliased solid rendering 
        typedef agg::renderer_scanline_aa_solid<renderer_base> renderer_solid;

        //bm this should not be antialiazed ?? 
        typedef agg::renderer_scanline_bin_solid<renderer_base> renderer_draft_solid;

        
        pixfmt pixf(rbuf_window());

        renderer_base rb(pixf);

        // draw antialised 
        //2016-05-21 use this unless too slow .. (different scanlines ??)
        renderer_solid rs(rb);   

        // use when picking 
        //renderer_draft_solid rs(rb);   
        
        //----------------------------------------
        // flags 
        
        bool draw_labels = true;        

        // draw in false colors for picking 
        //bool draw_picking = false;

        bool draw_picking = m_picking_mode;

        
        bool draw_number_labels = false;   // draw atom labels as numbers 
        bool draw_carbon_labels = true;
        bool draw_hydrogen_labels = true;

        //----------------------------------------
        
        // TODO for picking mode the clip_box can be MUCH SMALLER!!        
        // !!!!

        if (!draw_picking) {
          ras.clip_box(0, 0, rb.width(), rb.height());

        } else {

          // restrict drawing to exactly the pixel we are looking at!
          // (deactivate to see the false colors)
          int box_size_half = 1.0;
          ras.clip_box(m_current_raw_x-box_size_half,
                       m_current_raw_y-box_size_half,
                       m_current_raw_x+box_size_half,
                       m_current_raw_y+box_size_half);
          
          
        }
          
        // clear entire area with white TODO: background color
        
        if (draw_picking) {
          
          rb.clear(object_id_to_rgba8(256*256*256-1) ); // corresponds highest number TEST --- to highest number = while 

        } else {
          // background, this is white; should be configurable!
          rb.clear(agg::rgba(1,1,1));
        }
          
        // get the molecule 
        // NOTE 2016-05-17 we are NOT making a copy, i.e. we are changing atom xyz coordinates
        molecule & mol = m_molecules[m_cur_molecule];

          // print coords
          if (0) { std::cout << " COORDS \n";
            for(unsigned int i = 0; i< mol.num_atoms(); ++i) {
              std::cout << " i : " << i << " : "
                        << mol.atom(i).x << " , "  
                        << mol.atom(i).y << " , "  
                        << mol.atom(i).z << " oid = " << mol.atom(i).oid <<std::endl;  
            }
          }
          
        
        for(unsigned int i = 0; i < mol.num_atoms(); i++) {

          double fac = 1.0;
          double pos[3];
          double radius = 1.0;
          double paper_xy[2];
          double paper_radius;
            
          pos[0] = mol.atom(i).x;
          pos[1] = mol.atom(i).y;
          pos[2] = mol.atom(i).z;
          
          // 3D-> 2D transformation (perspective) TODO should a matrix ???
          atompos( fac, pos, radius, paper_xy, &paper_radius );
          
          // now SET xy positions (note: not setting z coordinate may be a problem)
            //          mol.set_atom_xy(i,paper_xy[0],paper_xy[1]);
          
          mol.set_atom_paper_xy(i,paper_xy[0],paper_xy[1]);
          
          mol.set_atom_paper_radius(i,paper_radius);
          
          //          std::cout << " i : " << i << " pos[0] = " << pos[0] << " paper_x = " << paper_xy[0] << std::endl;
          //          std::cout << " i : " << i << " pos[1] = " << pos[1] << " paper_x = " << paper_xy[1] << std::endl;
          
          
        }  // num_atoms 
        

        
        // this is arbitrary, change 
        double min_x =  1e100;
        double max_x = -1e100;
        double min_y =  1e100;
        double max_y = -1e100;

        for(unsigned int i = 0; i < mol.num_atoms(); i++)
        {
            if(mol.atom(i).px < min_x) min_x = mol.atom(i).px;
            if(mol.atom(i).py < min_y) min_y = mol.atom(i).py;
            if(mol.atom(i).px > max_x) max_x = mol.atom(i).px;
            if(mol.atom(i).py > max_y) max_y = mol.atom(i).py;
        }


        

        // sort atoms back to front for 2d drawing
        auto S = sort_atoms(mol);
        
        // TEST 2016-05-18 do not translate every time by making this static 

        static double initial_x = (max_x + min_x);
        static double initial_y = (max_y + min_y); 

        // setup affine translation/scaling 
        
        agg::trans_affine mtx;

        // affine _translation_ 
        // TODO do not do this every time --> molecule 'pumps' otherwise TODO
        //        mtx *= agg::trans_affine_translation(-(max_x + min_x) * 0.5, -(max_y + min_y) * 0.5);

        mtx *= agg::trans_affine_translation(-(initial_x) * 0.5, -(initial_y) * 0.5);

        
        // TEST 2016-05-18 make this static, effect is that the scale stays constant
        static double scale = width  / (max_x - min_x);

        double t = height / (max_y - min_y);
        if(scale > t) scale = t;
        
        double text_size = mol.average_bond_len() * m_text_size.value() / 4.0;
        double thickness = mol.average_bond_len() / 
                           sqrt(m_scale < 0.0001 ? 0.0001 : m_scale) /
                           8.0;
        
        mtx *= agg::trans_affine_scaling(scale*0.80, scale*0.80);

        // 2016-05-18 disable affine (2d) rotation 
        //---        mtx *= agg::trans_affine_rotation(m_angle);

        // apply manual scale
        mtx *= agg::trans_affine_scaling(m_scale, m_scale);

        // apply manual tranlation in xy plane 
        mtx *= agg::trans_affine_translation(m_center_x, m_center_y);

        // utiliy function in 'platform' : world to screen coordinates 
        mtx *= trans_affine_resizing();

        // done with transformations 
        
        // black - for labels (not used when in picking mode )
        rs.color(agg::rgba(0,0,0));

        //--------------------
        
        // draw bonds 
        for(unsigned int i = 0; i < mol.num_bonds(); i++)
        {
          // 2016-05-18 do not draw hydrogen bonds when hydrogen atom drawing is disabled

          bool draw_bond = true;
          if (!draw_hydrogen_labels) {
            // need to test both atoms
            if (strcmp(mol.atom(mol.bond(i).idx1).label, "H") == 0) {
              draw_bond = false;
            }
            if (strcmp(mol.atom(mol.bond(i).idx2).label, "H") == 0) {
              draw_bond = false;
            }

          }
              
          if (draw_bond) { 
            bond_vertex_generator bond(mol,i,mol.bond(i), 
                                       m_thickness.value() * thickness);
            agg::conv_transform<bond_vertex_generator> tr(bond, mtx);
            ras.add_path(tr);
            agg::render_scanlines(ras, sl, rs);
          }
        }

        //--------------------

        // draw atoms / atom labels 
        if (draw_labels) { 
        
        // draw ellipses for atom labels 
        agg::ellipse ell;
        agg::conv_transform<agg::ellipse> tr(ell, mtx);

        // atom label 
        agg::gsv_text label;
        agg::conv_stroke<agg::gsv_text> ls(label);
        agg::conv_transform<agg::conv_stroke<agg::gsv_text> > lo(ls, mtx);
        ls.line_join(agg::round_join);
        ls.line_cap(agg::round_cap);
        ls.approximation_scale(mtx.scale());


        for(unsigned int i = 0; i < mol.num_atoms(); i++)
        {

          // draw sorted  (low z first)
          
          auto A = mol.atom( S[i] );
          
          bool draw_label = true;

          // TODO pre-calc _once_
          if (!draw_carbon_labels && (strcmp(A.label, "C") == 0))  {
            draw_label = false;
          }
          if (!draw_hydrogen_labels && (strcmp(A.label, "H") == 0))  {
            draw_label = false;
          }

          bool draw_ellipsis = true;
          if (!draw_label && !draw_picking) {
            draw_ellipsis = false;
          }
          
          if (draw_label)   // draw ellipsis and label on top 
            {

              if (draw_ellipsis) { 
              
              bool draw_perspective_radii = true;
              double ell_radius = text_size * 2.5;
              
              if (draw_perspective_radii) {

                // 2016-05-18 
                // this is not yet working as expected ....
                ell_radius += (A.p_atom_radius -1.0); 

                
//TODO                std::cout << " atom i : " << i << " radius    = " << A.p_atom_radius << std::endl; 
//TODO                std::cout << " atom i : " << i << " ellradius = " << ell_radius << std::endl; 
              }

              ell.init(A.px, 
                       A.py, 
                       text_size * 2.5,  // rx
                       text_size * 2.5,  // ry
                       20                // steps 
                       );
              

              // ras = rasterizer
              ras.add_path(tr);
              // white   
              //                rs.color(agg::rgba(1.0, 1.0, 1.0));

              // 'atom' color (the ellipsis behind the label)
              
              if (draw_picking) {
                // draw in a unique color 
                //                std::cout << " atom color id = " << A.oid << std::endl;

                // TEST 2016-05-19 
                //red                rs.color(object_id_to_rgba8( 255 ) ); // FOR TESTING               
                //green                rs.color(object_id_to_rgba8( 256*256 - 255 ) ); // FOR TESTING               
                //blue
                //blue rs.color(object_id_to_srgba8( 256*256*256 - 255*255 - 255 ) ); // FOR TESTING               

                //red
                //                rs.color(object_id_to_rgba8( 254 ) ); // FOR TESTING               

                rs.color(object_id_to_rgba8( A.oid ) ); // FOR TESTING               
                
              } else {

                // color of the ellipses (simplistic 'atom balls') 
                // TODO: we could also draw balls in 'atom color'
                rs.color(agg::rgba(0.8, 0.8, 1.0)); // FOR TESTING

              }
              agg::render_scanlines(ras, sl, rs);

              } // .. draw_ellipsis 
              
              //
              // draw label for atom i 
              // 
              
              bool draw_label = true;
              if (!draw_carbon_labels && (strcmp(A.label, "C") == 0))  {
                draw_label = false;
              }
              if (!draw_hydrogen_labels && (strcmp(A.label, "H") == 0))  {
                draw_label = false;
              }

              // for picking mode we do not want labels in a different color
              // picking by selecting the ellipsis color, i.e. the 'atom', would fail
              if (draw_picking) {
                draw_label = false;
              }
              
              if (draw_label) 
                {
                  ls.width(m_thickness.value() * thickness);
                  
                  if (!draw_number_labels) { 
                    label.text(A.label);
                  } else {
                    std::stringstream num;
                    num << S[i];
                    label.text(num.str().c_str());
                  }
                  
                  label.start_point(A.px - text_size/2*3.0, A.py - text_size/2*3.0);
                  label.size(text_size*3.0);
                  ras.add_path(lo);

                  // labels not drawn in picking mode anyway ... 
                  rs.color(m_atom_colors[A.color_idx]);  // atom label color 
                  agg::render_scanlines(ras, sl, rs);
                }
              
            }
              
        }
        }


        // molecule name (upper left) 
        if (1) { 
          agg::gsv_text label;
          agg::conv_stroke<agg::gsv_text> ls(label);
          agg::conv_transform<agg::conv_stroke<agg::gsv_text> > lo(ls, mtx);
          ls.line_join(agg::round_join);
          ls.line_cap(agg::round_cap);
          ls.approximation_scale(1.0);
          agg::conv_transform<agg::conv_stroke<agg::gsv_text> > name(ls, trans_affine_resizing());
          ls.width(1.5);
          label.text(mol.name());  // molecule name 
          label.size(10.0);        // size 
          label.start_point(10.0, start_height - 20.0);
          ras.reset();
          ras.add_path(name);
          rs.color(agg::rgba(0,0,0)); // back (name of molecule)
          agg::render_scanlines(ras, sl, rs);
        }
        
        // controls at bottom 
        if (1) { 
          agg::render_ctrl(ras, sl, rb, m_thickness);
          agg::render_ctrl(ras, sl, rb, m_text_size);
        }

        if (1) {

          //if (m_mouse_move) { 

          if (1) { 

            //
            // TEST 2016-05-19 draw ellipses at current xy position when mouse button is pressed 
            //
            // start from screen positions , use inverse mtx transform
            // and then apply inverse projection (2d->3d) assuming the z position
            
            agg::ellipse ell;

            // mouse coords 
            double x2 = m_current_raw_x;  // should bre renamed to screen
            double y2 = m_current_raw_y;

            //           std::cout << " center = " << m_center_x << " ,  " << m_center_y << std::endl;
            
            // screen coords --> world coords 
            //--TEST 2016-05-20x OK leave commented           trans_affine_resizing().inverse_transform(&x2, &y2);

            // this _px/y have inverse transform + center correction
//            double x2 = m_pdx;
//            double y2 = m_pdy;
            
            // inverse transform: transform from screen coordinates to world coordinates
            mtx.inverse_transform(&x2, &y2);

            //TODO 2016-05-20 now need to undo affine projection!!

            // OK THIS WORKS ! double wanted_z = 5.0; // TEST 2016-05-21

            // 2016-05-21 this is not perfect --> we should have some sort of graphical feedback
            // TEST 
            double paper_xy[2] = { x2, y2 };           
            double wanted_z = z_of_nearest_atom(m_cur_molecule,paper_xy );  // keep fixed for now

            
            // b: screen coords, 2d
            double b[2];
            b[0] = x2;
            b[1] = y2;
            
            // camera (fixed for now)
            double c[3] = { 0 , 0, 15.0 };
            
            // TESTING inverse projection 
            double a[3];
            //            a[2] = 0.0; // fixed (for now)  --> the WANTED z-component 

            a[2] = wanted_z; // fixed (for now)  --> the WANTED z-component 

            
            a[0] = b[0]/c[2]*(a[2]-c[2]);
            a[1] = b[1]/c[2]*(a[2]-c[2]);
            
            //trans_affine_resizing().inverse_transform(&x2, &y2);
           
            //            std::cout << " x  = " << x  << " y  = " << y << std::endl; 
            std::cout << " x2 = " << x2 << " y2 = " << y2 << std::endl; 
            std::cout << " ax = " << a[0] << " ay = " << a[1] << std::endl; 

            // HMM 2016-05-20 perform mtx on ell
            agg::conv_transform<agg::ellipse> tr(ell, mtx);

            std::cout << " mouse_move: x = " << m_current_affine_x << " y = " << m_current_affine_y << std::endl;

            if (m_mouse_move) {  // TEST only draw this ellipsis on mouse_move
            
            //
            // __local__ coordinates 
            // 

              // draw ellipsis at 2d mouse position, but use local coordinates 
            ell.init(

                     x2,y2,  // local coordinates! : should be {x2,y2,0}

                     0.2, //text_size * 2.5,  // rx
                     0.2, // text_size * 2.5,  // ry
                     20                // steps 
                     );
                        
            // ras = rasterizer
            ras.add_path(tr);

            rs.color(agg::rgba(255, 0, 0)); // FOR TESTING

            agg::render_scanlines(ras, sl, rs);

            }
            
            // 2016-05-20 draw a 3d-dot (given by the backtransform of the 2d dot)
            if (1) {

              // world coordinates, NOT screen coordinates 
              // leave 'z' at zero for now 

              static double mousepos_distance = 16.0;

              static double camera_distance = 15.0; // currently fixed in atompos / sv2.cpp to 15.0

              // only update on mouse_move, but always show
              if (m_mouse_move) { 

                // take paper coordinates and determine atom next to the position (in screen coords)
                
//-2016-05-21                mousepos[0] = x2; //* mousepos_distance / camera_distance;
//-2016-05-21                mousepos[1] = y2; // * mousepos_distance / camera_distance;

                //-                double p[2] = { mousepos[0], mousepos[1] };

                //++ TEST mousepos[2] = z_of_nearest_atom(m_cur_molecule,p );  // keep fixed for now

                // 2016-05-20 works, (depending on the perspective transform done in 'atompos'

                mousepos[0] = -a[0];  // something wrong with the atompos^-1
                mousepos[1] = -a[1];
                mousepos[2] = a[2];

                //mousepos[2] = wanted_z; // 0.0;
                
                // set the 'found' z coordinate for determining the 3d projection / mouse pos coordinates 
//-                mousepos_distance = mousepos[2];
//-                mousepos[0] = x2 / ( mousepos_distance / camera_distance);
//-                mousepos[1] = y2 / ( mousepos_distance / camera_distance);

              }

              // mousepos 3d coordinates
              // 1. project
              // 2. transform with mtx (->done by ell ) 
              // (i.e. do the same thing as donw for the atoms/bonds)
              
              // do 3d-2d projection
              double mp[2], mr;
              atompos( /*fac*/ 1.0, mousepos, /*radius*/ 1.0 , mp, &mr );

              // subtract center in local coords TODO HOW ? 
//---              {
//---                // center 3d point 
//---                // works only because the camera is a { 0,0, z }
//---                // 2016-05-20 we only move center in 3d, the viewing distance is constant 
//---                // BUT we need to transform m_center_x and m_center_y from screen to world coordinates
//---
//---                double cx2 = m_center_x;  // should br renamed to screen
//---                double cy2 = m_center_y;
//---                
//---                // screen coords --> world coords 
//---                //--TEST
//---                //-TEST 2016-05-21 OK ???    trans_affine_resizing().inverse_transform(&cx2, &cy2);
//---
//---                // inverse transform
//---                mtx.inverse_transform(&cx2, &cy2);
//---                
//---                //-TEST
//---                //                double center[3] = { cx2, cy2, 15.0 /* dist */  }; 
//---
//---//TODO?                //+TEST 2016-05-20 the dot drawn at mouse pos should ne local coords {x2,y2,dist}                
//---//TODO?                double center[3] = { x2, y2, 15.0 /* dist */  }; 
//---//TODO?                
//---//TODO?                // do 3d-2d projection
//---//TODO?                double cp[2], cr;
//---//TODO?                atompos( /*fac*/ 1.0, center, /*radius*/ 1.0 , cp, &cr );
//---
//---                /// TEST TEST CONTINUE HERE 2016-05-20 
//---                //mp[0] += cp[0];
//---                //mp[1] += cp[1];
//---                
//---              }

              // draw mouse position with z value
              // the mtx transform is the 'normal' one (without the 3d->2d projection which has to be done first 
              ell.init(
                       mp[0],mp[1],                       
                       
                       0.15, //text_size * 2.5,  // rx
                       0.15, // text_size * 2.5,  // ry
                       20                // steps 
                       );
              
              // ras = rasterizer
              ras.add_path(tr);
              
              rs.color(agg::rgba8(100, 255, 100)); // FOR TESTING
              
              agg::render_scanlines(ras, sl, rs);

              // OK
              if (1) {  // draw center dot (i.e. this one should never move in a rotation)
                
                // center 3d point 
                // works only because the camera is a { 0,0, z }
                //                double center[3] = { m_center_x, m_center_y, camera_distance }; 

                double center[3] = { 0,0,0 }; /// WORLD coordinates (NOT screen coordinates) 

                
                // do 3d-2d projection
                double mp[2], mr;
                atompos( /*fac*/ 1.0, center, /*radius*/ 1.0 , mp, &mr );
                
                ell.init(
                         //                       mousepos[0],mousepos[1],
                         mp[0],mp[1],
                         
                         0.3, //text_size * 2.5,  // rx
                         0.3, // text_size * 2.5,  // ry
                         20                // steps 
                         );
                
                // ras = rasterizer
                ras.add_path(tr);
                
                rs.color(agg::rgba8(100, 40, 100)); // FOR TESTING
                
                agg::render_scanlines(ras, sl, rs);
                
              

              }
              
            }
            
          }
          
        }
          
    }
  
  //----------------------------------------
  
    virtual void on_idle()
    {
      //        m_angle += agg::deg2rad(0.1);

      rotate(m_cur_molecule,1,-2.0);         
      
      double oldp[3] = { mousepos[0],mousepos[1],mousepos[2] };          
      rotate_pos(oldp,1,-2.0, mousepos );

      force_redraw();
    }

    virtual void on_mouse_button_down(int x, int y, unsigned flags)
    {

      std::cout << " mouse button down: x =  " << x << " y = " << y 
                << std::endl;

      // TEST 2016-05-18
      if (1) {
        //      Get the pointer to the beginning of the i-th row (Y-coordinate)
        // and shift it to the i-th position, that is, X-coordinate.
        //---------------

        ///// TODO 2016-05-18 not correct ???!
        ///// TODO draw ellipse/dot at position to check (inverted y ??)

        int yy = y ; // not inverted 
        
        //        unsigned char* ptr = rbuf_window().row_ptr(yy) + x * 3; // i.e. 3 bytes 

        //----color_type col = rb.pixel(x,y); // i.e. 3 bytes 

        
        // assume rgb8 format, one byte per color

        // 2016-05-19 why there is not pix_ptr ? 
        //WORKS        unsigned char* ptr = rbuf_window().row_ptr(y) + x * 3; // i.e. 3 bytes 
        // works with bgr24 format ... TODO make this configurable ??! 
//works        unsigned int b = *ptr++;
//works        unsigned int g = *ptr++;
//works        unsigned int r = *ptr++;

        // this does work too - and we do not have to do pointer arithmetic 
        pixfmt pixf(rbuf_window());
        auto col = pixf.pixel(x,y); // i.e. 3 bytes 

        unsigned int r = col.r;
        unsigned int g = col.g;
        unsigned int b = col.b;
        
        // TODO 2016-05-18 look up the buffer format TODO


        std::cout << " read color = " << " r = " << r << " g = " << g << " b = " << b << std::endl; 
        
      }
      

      
        m_mouse_move = true;


        //bm this should be the transform from screen back to object space coords ...
        
        double x2 = x;
        double y2 = y;
        trans_affine_resizing().inverse_transform(&x2, &y2);

        m_pdx = m_center_x - x2;
        m_pdy = m_center_y - y2;

        std::cout << " pdx = " << m_pdx << " pdy = " << m_pdy << std::endl;
        
        m_prev_scale = m_scale;
        m_prev_angle = m_angle + agg::pi;
        
        force_redraw();
        
          

    }





    virtual void on_mouse_button_up(int x, int y, unsigned flags)
    {
        m_mouse_move = false;
    }

    virtual void on_mouse_move(int x, int y, unsigned flags)
    {
        double x2 = x;
        double y2 = y;
        trans_affine_resizing().inverse_transform(&x2, &y2);

        //2016-05-19 store xy coordinates _after_ transformation 

        m_current_affine_x = x2;
        m_current_affine_y = y2;

        m_current_raw_x = x;
        m_current_raw_y = y;

        
        // BM currently not used, will need this for trackball!
        // (affine rotation in xy plane)
        if(m_mouse_move && (flags & agg::mouse_left) != 0)
        //        if (0) 
        {
            double dx = x2 - m_center_x;
            double dy = y2 - m_center_y;
            m_scale = m_prev_scale * 
                      sqrt(dx * dx + dy * dy) / 
                      sqrt(m_pdx * m_pdx + m_pdy * m_pdy);

            //bm TODO m_angle = m_prev_angle + atan2(dy, dx) - atan2(m_pdy, m_pdx);
            force_redraw();
        }

        // 2016-05-19 this _is_ used for xy translation
        if(m_mouse_move && (flags & agg::mouse_right) != 0)
        {
            m_center_x = x2 + m_pdx;
            m_center_y = y2 + m_pdy;
            force_redraw();
        }

    }


  void rotate( unsigned int mol_idx, int dir, double amount ) {

    int ixyz = dir;
    
    double tmat[3][3];
          
    unit_mat( tmat );
    
    // TEST 2016-05-16 rotation 
    rotmat( ixyz, amount , tmat );

    // get the current mol
    molecule & mol = m_molecules[mol_idx];

    // 3d transform 
    tmat_transform( tmat, mol );

  }

  void rotate_pos( double pos[3], int dir, double amount , double res[3] ) {

    int ixyz = dir;
    
    double tmat[3][3];
          
    unit_mat( tmat );
    
    // TEST 2016-05-16 rotation 
    rotmat( ixyz, amount , tmat );

    // 3d transform 
    tmat_transform_pos( tmat, pos , &res[0] );

  }

  



    virtual void on_key(int, int, unsigned key, unsigned)
    {
        switch(key)
        {
        case agg::key_left:
        case agg::key_up:
        case agg::key_page_up:
            if(m_cur_molecule) --m_cur_molecule;
            force_redraw();
            break;

        case agg::key_right:
        case agg::key_down:
        case agg::key_page_down:
            if(m_cur_molecule < m_num_molecules - 1) ++m_cur_molecule;
            force_redraw();
            break;

        case 'a': {

          std::cout << "a pressed - rotate 'left' \n";

          rotate(m_cur_molecule,1,-2.0);         

          double oldp[3] = { mousepos[0],mousepos[1],mousepos[2] };          
          rotate_pos(oldp,1,-2.0, mousepos );
          
          force_redraw();          

          break;
        }


        case 'd': {

          std::cout << "d pressed - rotate 'right' \n";

          rotate(m_cur_molecule,1,2.0);
          double oldp[3] = { mousepos[0],mousepos[1],mousepos[2] };          
          rotate_pos(oldp,1,2.0, mousepos );

          force_redraw();          

          break;
        }

          
        case 'w': {

          std::cout << "w pressed - rotate 'up' \n";

          rotate(m_cur_molecule,2,-2.0);
          double oldp[3] = { mousepos[0],mousepos[1],mousepos[2] };          
          rotate_pos(oldp,2,-2.0, mousepos );
          force_redraw();          

          break;
        }

        case 's': {

          std::cout << "s pressed - rotate 'down' \n";

          rotate(m_cur_molecule,2,2.0);
          double oldp[3] = { mousepos[0],mousepos[1],mousepos[2] };          
          rotate_pos(oldp,2,2.0, mousepos );
          force_redraw();          

          break;
        }

        case 'q': {

          std::cout << "q pressed - rotate 'in plane left' \n";

          rotate(m_cur_molecule,3,2.0);
          double oldp[3] = { mousepos[0],mousepos[1],mousepos[2] };          
          rotate_pos(oldp,3,2.0, mousepos );
          force_redraw();          

          break;
        }

        case 'e': {

          std::cout << "e pressed - rotate 'in plane right' \n";

          rotate(m_cur_molecule,3,-2.0);
          double oldp[3] = { mousepos[0],mousepos[1],mousepos[2] };          
          rotate_pos(oldp,3,-2.0, mousepos );
          force_redraw();          

          break;
        }

          // TESTING  2016-05-19
        case 'm': {

          std::cout << "mouse down pressed  \n";

          force_redraw();          

          break;
        }

          // TESTING  2016-05-19
        case 'p': {

          std::cout << "p pressed - toggle picking mode \n";

          m_picking_mode = !m_picking_mode;
          
          force_redraw();          

          break;
        }

        case 'r' : {

          // recenter (like space in antique hyperchem)

          m_center_x = start_width / 2;
          m_center_y = start_height / 2;

          
          force_redraw();

          break; 
        }

          
        case 'k' : { // should be - but is intercepted by control; todo remove control
          m_scale = 0.9*m_scale;
          force_redraw();
          break;
        }
        case 'l' : { // should be '+'
          m_scale = 1.1*m_scale;
          force_redraw();
          break;
        }
          
        case ' ':
            wait_mode(!wait_mode());
            break;
        }
    }

};


int agg_main(int argc, char* argv[])
{
    const char* fname = "1.sdf";

    if(argc > 1)
    {
        fname = argv[1];
    }

    enum flip_y_e { flip_y = true };

    the_application app(pix_format, flip_y, fname);
    app.caption("AGG - A Simple SDF Molecular Viewer");

    if(app.init(start_width, start_height, agg::window_resize | agg::window_keep_aspect_ratio ))
    {
        return app.run();
    }

    return 1;
}


