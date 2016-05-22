
#ifndef _AGG_LINE_H_
#define _AGG_LINE_H_

// TODO do not use the 'agg' namespace !!

namespace agg
{
    class line
    {
    public:
        line() :
            m_x1(0.0), m_y1(0.0), m_x2(1.0), m_y2(0.0), m_thickness(0.1)
        {
        }

        line(double x1, double y1, double x2, double y2, double thickness) :
            m_x1(x1), m_y1(y1), m_x2(x2), m_y2(y2), m_thickness(thickness)
        {
        }

        void init(double x1, double y1, double x2, double y2)
        {
            m_x1 = x1;
            m_y1 = y1;
            m_x2 = x2;
            m_y2 = y2;
        }

        void thickness(double th)
        {
            m_thickness = th;
        }

        void rewind(unsigned start);
        unsigned vertex(double* x, double* y);

    private:
        double   m_x1;
        double   m_y1;
        double   m_x2;
        double   m_y2;
        double   m_dx;
        double   m_dy;
        double   m_thickness;
        unsigned m_vertex;
    };



    inline void line::rewind(unsigned)
    {
        calc_orthogonal(m_thickness*0.5, m_x1, m_y1, m_x2, m_y2, &m_dx, &m_dy);
        m_vertex = 0;
    }



    inline unsigned line::vertex(double* x, double* y)
    {
        switch(m_vertex)
        {
        case 0:
            *x = m_x1 - m_dx;
            *y = m_y1 - m_dy;
            m_vertex++;
            return path_cmd_move_to;

        case 1:
            *x = m_x2 - m_dx;
            *y = m_y2 - m_dy;
            m_vertex++;
            return path_cmd_line_to;

        case 2:
            *x = m_x2 + m_dx;
            *y = m_y2 + m_dy;
            m_vertex++;
            return path_cmd_line_to;

        case 3:
            *x = m_x1 + m_dx;
            *y = m_y1 + m_dy;
            m_vertex++;
            return path_cmd_line_to;
        }
        return path_cmd_stop;
    }



    class solid_wedge
    {
    public:
        solid_wedge() :
            m_x1(0.0), m_y1(0.0), m_x2(1.0), m_y2(0.0), m_thickness(0.1)
        {
        }

        solid_wedge(double x1, double y1, double x2, double y2, double thickness) :
            m_x1(x1), m_y1(y1), m_x2(x2), m_y2(y2), m_thickness(thickness)
        {
        }

        void init(double x1, double y1, double x2, double y2)
        {
            m_x1 = x1;
            m_y1 = y1;
            m_x2 = x2;
            m_y2 = y2;
        }

        void thickness(double th)
        {
            m_thickness = th;
        }

        void rewind(unsigned start);
        unsigned vertex(double* x, double* y);

    private:
        double   m_x1;
        double   m_y1;
        double   m_x2;
        double   m_y2;
        double   m_dx;
        double   m_dy;
        double   m_thickness;
        unsigned m_vertex;
    };



    inline void solid_wedge::rewind(unsigned)
    {
        calc_orthogonal(m_thickness*2.0, m_x1, m_y1, m_x2, m_y2, &m_dx, &m_dy);
        m_vertex = 0;
    }



    inline unsigned solid_wedge::vertex(double* x, double* y)
    {
        switch(m_vertex)
        {
        case 0:
            *x = m_x1;
            *y = m_y1;
            m_vertex++;
            return path_cmd_move_to;

        case 1:
            *x = m_x2 - m_dx;
            *y = m_y2 - m_dy;
            m_vertex++;
            return path_cmd_line_to;

        case 2:
            *x = m_x2 + m_dx;
            *y = m_y2 + m_dy;
            m_vertex++;
            return path_cmd_line_to;
        }
        return path_cmd_stop;
    }









    class dashed_wedge
    {
    public:
        dashed_wedge() :
            m_x1(0.0), m_y1(0.0), m_x2(1.0), m_y2(0.0), 
            m_thickness(0.1),
            m_num_dashes(8)
        {
        }

        dashed_wedge(double x1, double y1, double x2, double y2, 
                     double thickness, unsigned num_dashes=8) :
            m_x1(x2), m_y1(y2), m_x2(x1), m_y2(y1), 
            m_thickness(thickness),
            m_num_dashes(num_dashes)
        {
        }

        void init(double x1, double y1, double x2, double y2)
        {
            m_x1 = x2;
            m_y1 = y2;
            m_x2 = x1;
            m_y2 = y1;
        }

        void num_dashes(unsigned nd)
        {
            m_num_dashes = nd;
        }

        void thickness(double th)
        {
            m_thickness = th;
        }

        void rewind(unsigned start);
        unsigned vertex(double* x, double* y);

    private:
        double   m_x1;
        double   m_y1;
        double   m_x2;
        double   m_y2;
        double   m_xt2;
        double   m_yt2;
        double   m_xt3;
        double   m_yt3;
        double   m_xd[4];
        double   m_yd[4];
        double   m_thickness;
        unsigned m_num_dashes;
        unsigned m_vertex;
    };



    void dashed_wedge::rewind(unsigned)
    {
        double dx;
        double dy;
        calc_orthogonal(m_thickness*2.0, m_x1, m_y1, m_x2, m_y2, &dx, &dy);
        m_xt2 = m_x2 - dx;
        m_yt2 = m_y2 - dy;
        m_xt3 = m_x2 + dx;
        m_yt3 = m_y2 + dy;
        m_vertex = 0;
    }


    unsigned dashed_wedge::vertex(double* x, double* y)
    {
        if(m_vertex < m_num_dashes * 4)
        {
            if((m_vertex % 4) == 0)
            {
                double k1 = double(m_vertex / 4) / double(m_num_dashes);
                double k2 = k1 + 0.4 / double(m_num_dashes);

                m_xd[0] = m_x1 + (m_xt2 - m_x1) * k1;
                m_yd[0] = m_y1 + (m_yt2 - m_y1) * k1;
                m_xd[1] = m_x1 + (m_xt2 - m_x1) * k2;
                m_yd[1] = m_y1 + (m_yt2 - m_y1) * k2;
                m_xd[2] = m_x1 + (m_xt3 - m_x1) * k2;
                m_yd[2] = m_y1 + (m_yt3 - m_y1) * k2;
                m_xd[3] = m_x1 + (m_xt3 - m_x1) * k1;
                m_yd[3] = m_y1 + (m_yt3 - m_y1) * k1;
                *x = m_xd[0];
                *y = m_yd[0];
                m_vertex++;
                return path_cmd_move_to;
            }
            else
            {
                *x = m_xd[m_vertex % 4];
                *y = m_yd[m_vertex % 4];
                m_vertex++;
                return path_cmd_line_to;
            }
        }
        return path_cmd_stop;
    }

}



#endif 
