
#ifndef _MOLECULE_H_
#define _MOLECULE_H_

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include <vector>

#include "sv2.h"

class molecule {

 private:

  static unsigned long object_id;
  
public:

  
  struct atom_type
  {
    // xyz coordinates; these get changed by 3d rotations 
    double   x;
    double   y;
    double   z;

    // 'paper' coords (i.e. 3d->2d transformed)
    double px;
    double py; 
    // 'paper' atom radius 
    double p_atom_radius; 
    
    char     label[4];
    int      charge;      // currently not used 
    unsigned color_idx;

    unsigned long oid; // object id 

  atom_type() : oid ( ++object_id) {}

  };

  struct bond_type
  {
    unsigned idx1;
    unsigned idx2;
    unsigned order;
    int      stereo;
    int      topology;

    unsigned long oid; // object id 

  bond_type() : oid ( ++object_id) {}
    
  };

  // this is rudimentary at best .. change
enum atom_color_e
{
    atom_color_general = 0,
    atom_color_N       = 1,
    atom_color_O       = 2,
    atom_color_S       = 3,
    atom_color_P       = 4,
    atom_color_halogen = 5,
    end_of_atom_colors
};



    ~molecule();
    molecule();

    bool read(FILE* fd);

    unsigned num_atoms() const { return m_atoms.size(); }
    unsigned num_bonds() const { return m_bonds.size(); }

    const atom_type & atom(unsigned int idx) const { return m_atoms[idx]; }
    const bond_type & bond(unsigned int idx) const { return m_bonds[idx]; }

    // 2016-04-17 DO NOT DO THIS THIS WAY, ONLY FOR TESTING!!!!
 
    double bond_atom_1_x(unsigned int bond_idx)  {

      if (1) { 
      
      // 3D-> 2D transformation
      double pos[3];
      unsigned int idx = m_bonds[bond_idx].idx1;
      double paper_xy[2];  // paper coordinates 
      double paper_radius; 
      
      pos[0] = atom(idx).x;
      pos[1] = atom(idx).y;
      pos[2] = atom(idx).z;        
      atompos( 1.0, pos, 1.0, paper_xy, &paper_radius );

      m_atoms[idx].px = paper_xy[0];
      
      return paper_xy[0];

      }
      
      unsigned int idx = m_bonds[bond_idx].idx1;
      return atom(idx).px;

    }
    double bond_atom_1_y(unsigned int bond_idx)  {

      if (1) { 
      
      // 3D-> 2D transformation
      double pos[3];
      unsigned int idx = m_bonds[bond_idx].idx1;
      double paper_xy[2];  // paper coordinates 
      double paper_radius; 
      
      pos[0] = atom(idx).x;
      pos[1] = atom(idx).y;
      pos[2] = atom(idx).z;        
      atompos( 1.0, pos, 1.0, paper_xy, &paper_radius );

      m_atoms[idx].py = paper_xy[1];
      
      return paper_xy[1];

      }
      
      unsigned int idx = m_bonds[bond_idx].idx1;
      return atom(idx).py;
      
    }
    double bond_atom_2_x(unsigned int bond_idx)  {

      if (1) { 
      // 3D-> 2D transformation
      double pos[3];
      unsigned int idx = m_bonds[bond_idx].idx2;
      double paper_xy[2];  // paper coordinates 
      double paper_radius; 
      
      pos[0] = atom(idx).x;
      pos[1] = atom(idx).y;
      pos[2] = atom(idx).z;        
      atompos( 1.0, pos, 1.0, paper_xy, &paper_radius );

      m_atoms[idx].px = paper_xy[0];
      
      return paper_xy[0];

      }
      
      unsigned int idx = m_bonds[bond_idx].idx2;
      return atom(idx).px;

      
    }
    double bond_atom_2_y(unsigned int bond_idx) {

      if (1) { 
      
      // 3D-> 2D transformation
      double pos[3];
      unsigned int idx = m_bonds[bond_idx].idx2;
      double paper_xy[2];  // paper coordinates 
      double paper_radius; 
      
      pos[0] = atom(idx).x;
      pos[1] = atom(idx).y;
      pos[2] = atom(idx).z;        
      atompos( 1.0, pos, 1.0, paper_xy, &paper_radius );

      m_atoms[idx].py = paper_xy[1];

      return paper_xy[1];      
      
      }
      
      unsigned int idx = m_bonds[bond_idx].idx2;
      return atom(idx).py;
      
    }
    

    void set_atom_paper_xy(unsigned int idx, double x, double y)  {
      m_atoms[idx].px = x;
      m_atoms[idx].py = y;

    }

    void set_atom_paper_radius(unsigned int idx, double r)  {
      m_atoms[idx].p_atom_radius = r;
    }

    
    void set_atom_xyz(unsigned int idx, double x, double y, double z) {
      m_atoms[idx].x = x;
      m_atoms[idx].y = y;
      m_atoms[idx].z = z;
    }

    
    //void set_bond(unsigned idx, double x, double y)  { m_bonds[idx].y = val; }

    
    double average_bond_len() const { return m_avr_len; }

    const char* name() const { return m_name; }

    static int    get_int(const char* buf, int pos, int len);
    static double get_dbl(const char* buf, int pos, int len);
    static char*  get_str(char* dst, const char* buf, int pos, int len);

private:

    std::vector<atom_type> m_atoms;
    std::vector<bond_type> m_bonds;

    char       m_name[128];
    double     m_avr_len;
};


#endif 
