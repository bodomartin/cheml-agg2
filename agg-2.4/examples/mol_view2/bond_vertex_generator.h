
#ifndef _BOND_VERTEX_GENERATOR_H_
#define _BOND_VERTEX_GENERATOR_H_

// TODO not all of these are needed !! 

#include "agg_rendering_buffer.h"
#include "agg_rasterizer_scanline_aa.h"
#include "agg_scanline_p.h"
#include "agg_renderer_scanline.h"
#include "agg_conv_transform.h"
#include "agg_conv_stroke.h"
#include "agg_bspline.h"
#include "agg_ellipse.h"
#include "agg_gsv_text.h"

#include "agg_line.h"

class bond_vertex_generator
{
    enum bond_style_e
    {
        bond_single,
        bond_wedged_solid,
        bond_wedged_dashed,
        bond_double,
        bond_double_left,
        bond_double_right,
        bond_triple
    };

public:
  bond_vertex_generator(molecule & mol,
                        unsigned int bond_index,
                        const molecule::bond_type& bond, double thickness) :
    m_mol(mol),  // 2016-05-17 store the molecule , not the bond ... (need updated atom positions)
    m_bidx(bond_index),
    m_bond(bond),
    m_thickness(thickness),
    m_style(bond_single)
    {
        if(bond.order == 1)
        {
            if(bond.stereo == 1) m_style = bond_wedged_solid;
            if(bond.stereo == 6) m_style = bond_wedged_dashed;
        }
        if(bond.order == 2)
        {
            m_style = bond_double;
            if(bond.topology == 1) m_style = bond_double_left;
            if(bond.topology == 2) m_style = bond_double_right;
        }
        if(bond.order == 3) m_style = bond_triple;
        m_line1.thickness(thickness);
        m_line2.thickness(thickness);
        m_line3.thickness(thickness);
        m_solid_wedge.thickness(thickness);
        m_dashed_wedge.thickness(thickness);
    }

    void rewind(unsigned int )
    {
        double dx, dy, dx1, dy1, dx2, dy2;

        double a1x = m_mol.bond_atom_1_x(m_bidx);
        double a1y = m_mol.bond_atom_1_y(m_bidx);
        double a2x = m_mol.bond_atom_2_x(m_bidx);
        double a2y = m_mol.bond_atom_2_y(m_bidx);
        
        switch(m_style)
        {
        case bond_wedged_solid:
          m_solid_wedge.init(a1x,a1y,a2x,a2y);
            m_solid_wedge.rewind(0);
            break;

        case bond_wedged_dashed:
          m_dashed_wedge.init(a1x,a1y,a2x,a2y);
            m_dashed_wedge.rewind(0);
            break;

        case bond_double:
        case bond_double_left:
        case bond_double_right:
          agg::calc_orthogonal(m_thickness, a1x,a1y,a2x,a2y ,
                                 &dx, &dy);
            dx1 = dy1 = 0;

            // To Do: ring perception and the proper drawing 
            // of the double bonds in the aromatic rings.
            //if(m_style == bond_double)
            {
                dx1 = dx2 = dx;
                dy1 = dy2 = dy;
            }
/*
            else if(m_style == bond_double_left)
            {
                dx2 = dx * 2.0;
                dy2 = dy * 2.0;
            }
            else
            {
                dx2 = -dx * 2.0;
                dy2 = -dy * 2.0;
            }
*/

            m_line1.init(a1x - dx1, 
                         a1y - dy1, 
                         a2x - dx1, 
                         a2y - dy1);
            m_line1.rewind(0);

            m_line2.init(a1x + dx2, 
                         a1y + dy2, 
                         a2x + dx2, 
                         a2y + dy2);
            m_line2.rewind(0);
            m_status = 0;
            break;

        case bond_triple:
            // To Do: triple bonds drawing.

        default:
          m_line1.init(a1x,a1y,a2x,a2y);
            m_line1.rewind(0);
            break;
        }
    }


    unsigned vertex(double* x, double* y)
    {
        unsigned flag = agg::path_cmd_stop;
        switch(m_style)
        {
            case bond_wedged_solid:
                return m_solid_wedge.vertex(x, y);

            case bond_wedged_dashed:
                return m_dashed_wedge.vertex(x, y);

            case bond_double_left:
            case bond_double_right:
            case bond_double:
                if(m_status == 0)
                {
                    flag = m_line1.vertex(x, y);
                    if(flag == agg::path_cmd_stop)
                    {
                        m_status = 1;
                    }
                }

                if(m_status == 1)
                {
                    flag = m_line2.vertex(x, y);
                }
                return flag;

            case bond_triple:
            case bond_single:
                break;
        }
        return m_line1.vertex(x, y);
    }



private:
  
  bond_vertex_generator(molecule & mol,
                        unsigned int bond_index,
                        const bond_vertex_generator&);
  const bond_vertex_generator& operator = (const bond_vertex_generator&);

  /*const*/ molecule & m_mol;
  unsigned int m_bidx; // bond index in m_mol
  
  const molecule::bond_type& m_bond;

  double m_thickness;
  bond_style_e m_style;
  agg::line m_line1;
  agg::line m_line2;
  agg::line m_line3;
  agg::solid_wedge m_solid_wedge;
  agg::dashed_wedge m_dashed_wedge;
  unsigned m_status;

};  // class bond_vertex_generator 


#endif
