
#include "molecule.h"

#include <cmath>
#include <cstring>

#include <iostream>

unsigned long molecule::object_id = 0;  // unique object id 

unsigned trim_cr_lf(char* buf)
{
    unsigned len = strlen(buf);

    // Trim "\n\r" at the beginning 
    unsigned pos = 0;
    while(len && (buf[0] == '\n' || buf[0] == '\r'))
    {
        --len;
        ++pos;
    }
    if(pos) strcpy(buf, buf + pos);

    // Trim "\n\r" at the end 
    while(len && (buf[len-1] == '\n' || buf[len-1] == '\r')) --len;
    buf[len] = 0;
    return len;
}



molecule::~molecule()
{
//    delete [] m_bonds;
//    delete [] m_atoms;
}


molecule::molecule() :

//  m_atoms(0),
//    m_num_atoms(0),
//    m_bonds(0),
//    m_num_bonds(0),

    m_avr_len(0)
{
    m_name[0] = 0;
}


int molecule::get_int(const char* buf, int pos, int len)
{
    char tmp[32];
    return atoi(get_str(tmp, buf, pos, len));
}

double molecule::get_dbl(const char* buf, int pos, int len)
{
    char tmp[32];
    return atof(get_str(tmp, buf, pos, len));
}

char* molecule::get_str(char* dst, const char* buf, int pos, int len)
{
    --pos;
    int buf_len = strlen(buf);
    buf += pos;
    while(pos + len < buf_len && isspace(*buf)) 
    {
        ++buf;
        --len;
    }
    if(len >= 0)
    {
        if(len > 0) memcpy(dst, buf, len);
        dst[len] = 0;
    }
    
    char* ts = dst;
    while(*ts && !isspace(*ts)) ts++;
    *ts = 0;
    return dst;
}


// read in the molecule TODO REMOVE !!

bool molecule::read(FILE* fd)
{
    char buf[512];
    unsigned len;
    if(!fgets(buf, 510, fd)) return false;
    len = trim_cr_lf(buf);
    if(len > 128) len = 128;

    if(len) memcpy(m_name, buf, len);
    m_name[len] = 0;

    if(!fgets(buf, 510, fd)) return false;
    if(!fgets(buf, 510, fd)) return false;
    if(!fgets(buf, 510, fd)) return false;
    trim_cr_lf(buf);

//    m_num_atoms = get_int(buf, 1, 3);
//    m_num_bonds = get_int(buf, 4, 3);

    int m_num_atoms_in = get_int(buf, 1, 3);
    int m_num_bonds_in = get_int(buf, 4, 3);


    std::cout << " number of atoms (in) : " << m_num_atoms_in << std::endl;
    std::cout << " number of bonds (in) : " << m_num_bonds_in << std::endl;
    
    if(m_num_atoms_in == 0 || m_num_bonds_in == 0) return false;

    // so that the elements are allocated 
    m_atoms.resize(m_num_atoms_in);
    m_bonds.resize(m_num_bonds_in);
    
    for(unsigned int i = 0; i < m_num_atoms_in; i++)
    {
        if(!fgets(buf, 510, fd)) return false;
        trim_cr_lf(buf);
/*01234567890123456789012345678901234567890
   -2.6793   -0.2552    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
*/
        m_atoms[i].x = get_dbl(buf, 1,  10);
        m_atoms[i].y = get_dbl(buf, 11, 10);

        //+bm store z coordinate!
        m_atoms[i].z = get_dbl(buf, 21, 10); 
        
        get_str(m_atoms[i].label, buf, 32, 3);

        // get charge from sd file
        m_atoms[i].charge = get_int(buf, 39, 1);
        if(m_atoms[i].charge) m_atoms[i].charge = 4 - m_atoms[i].charge;

        // TODO use real colors 
        if(strcmp(m_atoms[i].label, "N") == 0) m_atoms[i].color_idx = atom_color_N;
        if(strcmp(m_atoms[i].label, "O") == 0) m_atoms[i].color_idx = atom_color_O;
        if(strcmp(m_atoms[i].label, "S") == 0) m_atoms[i].color_idx = atom_color_S;
        if(strcmp(m_atoms[i].label, "P") == 0) m_atoms[i].color_idx = atom_color_P;
        if(strcmp(m_atoms[i].label, "F") == 0 ||
           strcmp(m_atoms[i].label, "Cl") == 0 || 
           strcmp(m_atoms[i].label, "Br") == 0 || 
           strcmp(m_atoms[i].label, "I")  == 0) m_atoms[i].color_idx = atom_color_halogen;
    }

    m_avr_len = 0.0;
    for(unsigned int i = 0; i < m_num_bonds_in; i++)
    {
        if(!fgets(buf, 510, fd)) return false;
        trim_cr_lf(buf);
/*
  1  2  1  0  0  0  0
*/
        // store atom indices 
        m_bonds[i].idx1 = get_int(buf, 1, 3) - 1;
        m_bonds[i].idx2 = get_int(buf, 4, 3) - 1;

        if(m_bonds[i].idx1 >= m_num_atoms_in || m_bonds[i].idx2 >= m_num_bonds_in ) return false;

        m_bonds[i].order    = get_int(buf, 7,  3);
        m_bonds[i].stereo   = get_int(buf, 10, 3);
        m_bonds[i].topology = get_int(buf, 13, 3);

        double a1x = bond_atom_1_x(i);
        double a1y = bond_atom_1_y(i);
        double a2x = bond_atom_2_x(i);
        double a2y = bond_atom_2_y(i);

        // average bond length 
        m_avr_len += sqrt((a1x - a2x) * (a1x - a2x) + 
                          (a1y - a2y) * (a1y - a2y));
        
//        m_avr_len += sqrt((m_bonds[i].x1 - m_bonds[i].x2) * (m_bonds[i].x1 - m_bonds[i].x2) + 
//                          (m_bonds[i].y1 - m_bonds[i].y2) * (m_bonds[i].y1 - m_bonds[i].y2));
    }

    // average bond length; used ? 
    m_avr_len /= double(m_bonds.size());

    while(fgets(buf, 510, fd))
    {
        trim_cr_lf(buf);
        if(buf[0] == '$') return true;
    }

    return false;
}

