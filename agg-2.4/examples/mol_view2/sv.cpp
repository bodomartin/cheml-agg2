/*
 *  xbs: an X-Window ball&sticks plotting program.
 *  Copyright (C) 1995  Michael Methfessel
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  The author can be contacted as follows:
 *
 *  Michael Methfessel
 *  methfessel@ihp-ffo.de
 *  Institute for Semiconductor Physics, PO Box 409,
 *  D-15204 Frankfurt (Oder), Germany
 *
 *
 *  Additions by Bodo Martin
 *
 *  History:
 *  ========
 *
 *  o  5/99 : additions for Labels: do not display a label (in .bs file) beginning
 *            with an underscore '_'
 *
 *  o  6/99 : callable with option -ps : produced a .ps file from .bs file
 *            without invoking X
 *
 *
 *  o  1/00 : Beginn der Arbeiten an einem Builder (primitiv, aber wirkungsvoll)
 *            Ziel : - Umstellung auf C++
 *                   - Schnittstellen zu semiemp. Programmen etc.
 *                   - *optionale* OpenGL -Darstellung, evtl anderes Toolkit (GTK,TCL?)
 *
 * ----> 2015-11-26 (!) use this as a template for real stuff
 */

#include <cstdio>
#include <csignal>

// 2015-11-26 factor them in a separate file / abstraction
#include <X11/X.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>

#include <cmath>
#include <ctime>
#include <cstring>

#include <cctype> // for isalpha
#include <cstdlib> // exit 


#include <boost/thread.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#define  W_HEIGHT     450      /* startup size */
#define  W_WIDTH      550

#define  NAMAX       2000      /* maximum number of atoms ? */
#define  NBMAX       8000      /* maximum number of bonds ? */
#define  NBTMAX       200
#define  NSPMAX        50

/* already defined under linux */
#define  PI           3.1415926

#define  MAXRAD     100.0
#define  GRAY0      0.5

//#define  PSFAC      1.42

static const float PSFAC = 1.42;

#define  NPOINTS    5
#define  NCOL       31       /* max number of colors to be allocated */
#define  FONT       "8x13"   /* font for info/help */
#define  LABFONT    "6x10"   /* font for labels */
#define  SHADOW     5
#define  BELL_LEVEL 30
#define  SVINP      10

#define  FBMAX    70000   /* max for atoms times frames */
#define  NFRMAX   8000    /* max num of frames */

/* return codes for interpret_keypress routine */
#define K_NOP        0
#define K_QUIT       1
#define K_REPLOT     2
#define K_RESETUP    3
#define K_UPDATE     4
#define K_READ_MORE  5
#define K_READ_DONE  6
#define K_ANIMATE    7

/* codes for gray modes */

#define G_STD    0
#define G_RAMP   1
#define G_LIGHT  2

#include "sv_math.h"
#include "sv_hc.h"

//--------------------------------------------------------------------------------

// define a global instance of sv_math (for now)
cheml::gui::sv_math svm;

cheml::gui::sv_hc svps;

//--------------------------------------------------------------------------------

/* struct for a ball (atom) representation */

struct ballstr {
  float pos[3];
  float rad;
  float gray;
  float r, g, b;
  int col;
  int special;
  char lab[21];
};

/* struct for a stick */

struct stickstr {
  int start;
  int end;
  float rad;
  float gray;
  int col;
};

/* ----- global variables ------- */

/* spec */

struct spec_s {
  char  lab[21];   /* label */
  float rad;
  float r, g, b;
  char  cname[81];
  int   col;
  float gray;
};

spec_s spec[NSPMAX];

/* atom */

struct atom_s {
  char lab[21];   /* label */
  float pos[3];
  float pol[3];
};

atom_s atom[NAMAX];

/* bonds */

struct bond_s {
  char  lab1[21];
  char  lab2[21];
  float min;
  float max;
  float rad;
  float r, g, b;
  char  cname[81];
  int   col;
  float gray;
};

bond_s bonds [NBTMAX];


// TODO use boost!!
void ms_sleep( int t );


int natom, nbond;   /* number of atoms, bonds */

struct ballstr  ball[NAMAX];

struct stickstr stick[NBMAX];

float arc[NPOINTS][2], xbot, xtop, ybot, ytop;

/* X vars */

Display        *dpy;
Window         win;
Drawable       drw;
Pixmap         pixmap, bgrmap;
GC             gc, gcbg, graygc[NCOL], shadowgc, labelgc, labbggc;
unsigned long  fground, bground, gpx[NCOL];
int            screen, depth, ncol;
Screen        *screenptr;
Colormap       cmap;

float p[NAMAX][3];
int nspec, nbas, nbonds, ngray;
int count;
int animate; /* jkl */
int msecs = 200; /* jkl */
int ipr = 10;
int igeo, igx, igy, igw, igh;
float igs;
int midx, midy;
float alat, dist, dist0, amp, dalfa, scale, tmat[3][3], radfac, bndfac;
float taux, tauy, dtaux, dtauy, taux0, tauy0, bg;
float gslope, gz0, light[3];
float center[3];

float   frame[3][FBMAX];
char    frstr[NFRMAX][81];
int     nframe, iframe, saveframe;

int hardcopy, usepixmap, numbers, grayvalues, wrinfo, fstep, wrhelp;
int framephase;
int bline, wire, withbonds, recenter, pmode, gmode, shadow, bondnums;
int color, autocolor, reverse, coords, stippled;
int showaxes;
int replot, resetup, startup, chginfo;
int num_print;
float pr_xoff, pr_yoff;
int lnum, xln;

char inf[81] = "in.bs",         /* input file */
               outf[81] = "Save.bs",      /* output file (default) */
                          prf[81] = "Bs.ps";         /* postscript file name (default) */

char inmv[81], prfsave[81], wname[81], curf[81];
char gmsg[101], emsg[101];
char svinput[SVINP][257];
int  svline, nsvline = 1;


/* (bm) */

int autops = 0;
char pmsg1[81], pmsg2[81];
char *lpos;

/* (bm) end */

/* forward declarations */
void WriteStatus( Drawable draw );
void WriteInfo( Drawable draw );
void WriteHelp();

/*
    procedures for graphics
    =======================

*/

#define MAXRGB 65535

#define HIGHGRAY 65535     /* top of gray scale */
#define LOWGRAY      0     /* bottom of gray scale */

int parse_color( char str[], float * rval, float * gval,
                 float * bval, float * gray ) {

  bool isname;
  int i, rc;
  XColor color, exact;
  float r, g, b;
  char cname[81];

  rc = 1;
  strcpy( cname, str );
  std::size_t j0 = 0;
  isname = false;

  for ( std::size_t j = 0; j < strlen( cname ); j++ ) {
    if ( std::isalpha( cname[j] ) ) { isname = true; }

    if ( cname[j] != ' ' ) { j0 = j; }
  }

  if ( isname ) {

    if ( j0 > 0 ) { cname[j0 + 1] = 0; }

    if ( j0 > 0 ) { str[j0 + 1] = 0; }

    if ( XLookupColor( dpy, cmap, cname, &exact, &color ) ) {
      r = exact.red   / ( float ) MAXRGB;
      g = exact.green / ( float ) MAXRGB;
      b = exact.blue  / ( float ) MAXRGB;
    }
    else {
      /*      printf("parse_color: Cannot look up color <%s>\n", cname); */
      r = g = b = 1.0;
      rc = 0;
    }
  }
  else {
    r = b = g = 999.9;
    sscanf( cname, "%f %f %f", &r, &g, &b );

    if ( r > 999.0 ) { r = 0.0; }

    if ( g > 999.0 ) { g = r; }

    if ( b > 999.0 ) { b = r; }
  }

  if ( reverse ) { r = 1.0 - r; g = 1.0 - g; b = 1.0 - b; }

  *rval = r;

  *gval = g;
  *bval = b;
  *gray = ( r + g + b ) / 3.0;
  /*  printf ("parse_color: rgb %.3f %.3f %.3f  gray %.3f  <%s>\n",
            *rval, *gval, *bval, *gray, cname);  */
  return rc;
}

int GetColorGC( unsigned long gpixel ) {

  int i, j, j0, icol;

  j0 = -1;

  for ( j = 0; j < ncol; j++ ) if ( gpx[j] == gpixel ) { j0 = j; }

  if ( j0 >= 0 ) {
    /*    printf("GetColorGC: already have pixel %d in slot %d\n",gpixel,j0);  */
    return j0;
  }

  if ( ncol >= NCOL ) {
    sprintf( emsg, "no more room for colors.. increase dimension NCOL" );
    return 0;
  }

  /* find a slot and make a graphics context */
  icol = ncol;

  for ( j = 0; j < ncol; j++ )
    if ( gpx[j] == -1 ) { icol = j; }

  /*  if (icol<ncol)
      printf ("Make GC in old slot %d (%d) for pixel %d\n", icol,ncol,gpixel);
    else
      printf ("Make GC in new slot %d (%d) for pixel %d\n", icol,ncol,gpixel); */

  if ( icol == ncol ) { ncol++; }

  graygc[icol] = XCreateGC( dpy, win, 0, 0 );

  XSetForeground( dpy, graygc[icol], gpixel );

  gpx[icol] = gpixel;

  return icol;

}

void FreeColorGC( int icol ) {
  int i, num;

  num = 0;

  for ( i = 0; i < nspec; i++ )
    if ( spec[i].col == icol ) { num++; }

  for ( i = 0; i < nbonds; i++ )
    if ( bonds[i].col == icol ) { num++; }

  if ( num == 0 ) {
    XFreeGC( dpy, graygc[icol] );
    gpx[icol] = -1;
  }

}

void SetColors() {

  int i, j, j0;
  int colptr;
  unsigned long gpixel, bground;
  XColor color, col1;

  bground = WhitePixel( dpy, screen );

  ncol = 0;

  for ( i = 0; i < nspec + nbonds; i++ ) {
    if ( i < nspec ) {
      col1.red   = spec[i].r * MAXRGB;
      col1.green = spec[i].g * MAXRGB;
      col1.blue  = spec[i].b * MAXRGB;
    }
    else {
      col1.red   = bonds[i - nspec].r * MAXRGB;
      col1.green = bonds[i - nspec].g * MAXRGB;
      col1.blue  = bonds[i - nspec].b * MAXRGB;
    }

    color = col1;

    if ( XAllocColor( dpy, cmap, &color ) ) {
      gpixel = color.pixel;
    }
    else {
      printf( "Cannot allocate %6d %6d %6d, use background instead\n",
              col1.red, col1.blue, col1.green );
      gpixel = bground;
    }

    colptr = GetColorGC( gpixel );

    if ( i < nspec ) {
      spec[i].col = colptr;
    }
    else {
      bonds[i - nspec].col = colptr;
    }
  }
}

// forward-declare match
int match( char str[], char pat[] );

int NewSpecColor( char pat[], char cname[], int helpme ) {

  int i, icol, nmatch;
  unsigned long gpixel;
  float f;
  char list[151];
  XColor col1, col2;

  if ( helpme ) {
    sprintf( gmsg, "Usage: color pattern [color]  - query or set color" );
    return 0;
  }

  if ( strlen( pat ) == 0 ) {
    sprintf( emsg, "color: no pattern specified" );
    return 0;
  }

  nmatch = 0;

  strcpy( list, "color" );

  for ( i = 0; i < nspec; i++ ) {
    if ( match( spec[i].lab, pat ) ) {
      nmatch++;
      strcat( list, " " );
      strcat( list, spec[i].lab );

      if ( strlen( cname ) == 0 ) {
        sprintf( gmsg, "Species %s has color <%s> rgb %.2f %.2f %.2f",
                 spec[i].lab, spec[i].cname, spec[i].r, spec[i].g, spec[i].b );
        return 0;
      }

      if ( ! parse_color( cname, &spec[i].r, &spec[i].g,
                          &spec[i].b, &spec[i].gray ) ) {
        sprintf( emsg, "Cannot identify color <%s>", cname );
        return 0;
      }

      strcpy( spec[i].cname, cname );

      if ( color ) {
        col1.red   = spec[i].r * MAXRGB;
        col1.green = spec[i].g * MAXRGB;
        col1.blue  = spec[i].b * MAXRGB;
        col2 = col1;

        if ( XAllocColor( dpy, cmap, &col1 ) ) {
          gpixel = col1.pixel;
        }
        else {
          sprintf( emsg, "Cannot allocate color <%s>", cname );
          return 0;
        }

        icol = spec[i].col;

        spec[i].col = -1;
        FreeColorGC( icol );
        icol = GetColorGC( gpixel );
        spec[i].col = icol;
      }
    }
  }

  if ( nmatch == 0 && ( strlen( cname ) == 0 ) ) {
    sprintf( emsg, "No species matches \'%s\'", pat );
    return 0;
  }

  f = MAXRGB;

  sprintf( gmsg, "%s <%s> rgb %.2f %.2f %.2f shown as %.2f %.2f %.2f",
           list, cname, col2.red / f, col2.green / f, col2.blue / f,
           col1.red / f, col1.green / f, col1.blue / f );

  if ( nmatch == 0 ) { return 0; }

  return 1;
}

void SetSmoothGrays() {

  int i, allocated[NCOL];
  int lastgnum;
  unsigned long gpixel;
  XColor color;
  int g1, g2, gnum;
  float dg;

  ngray = 21;
  g1 = HIGHGRAY;
  g2 = LOWGRAY;
  dg = ( g2 - g1 ) / ( ngray - 1.0 );

  lastgnum = 0;

  for ( i = 0; i < ngray; i++ ) {
    gnum = g1 + i * dg;
    color.red = color.blue = color.green = gnum;
    allocated[i] = XAllocColor( dpy, cmap, &color );

    if ( allocated[i] ) {
      lastgnum = gnum;
    }
    else {
      printf( "%7d: cannot allocate color.. use %d\n", gnum, lastgnum );
      color.red = color.blue = color.green = lastgnum;
      XAllocColor( dpy, cmap, &color );
    }

    graygc[i] = XCreateGC( dpy, win, 0, 0 );

    gpixel = color.pixel;
    XSetForeground( dpy, graygc[i], gpixel );
  }
}

void SetStippled4x4() {
  int i, j, n, depth;
  int gnum1, gnum2;
  unsigned long fg, bg, lg, dg, g1, g2;
  XColor color1, color2;
  Pixmap tile;
  int tw = 4, th = 4;
  static char t0[]  = {( char )0x00, ( char )0x00,
                       ( char )0x00, ( char )0x00
                      };
  static char t2[]  = {( char )0x01, ( char )0x00,
                       ( char )0x04, ( char )0x00
                      };
  static char t4[]  = {( char )0x01, ( char )0x04,
                       ( char )0x01, ( char )0x04
                      };
  static char t6[]  = {( char )0x05, ( char )0x02,
                       ( char )0x05, ( char )0x08
                      };
  static char t8[]  = {( char )0x05, ( char )0x0a,
                       ( char )0x05, ( char )0x0a
                      };
  static char t10[] = {( char )0xfa, ( char )0xfd,
                       ( char )0xfa, ( char )0xf7
                      };
  static char t12[] = {( char )0xfe, ( char )0xfb,
                       ( char )0xfe, ( char )0xfb
                      };
  static char t14[] = {( char )0xfe, ( char )0xff,
                       ( char )0xfb, ( char )0xff
                      };

  depth  = XDefaultDepth( dpy, screen );
  fg = BlackPixel( dpy, screen );
  bg = WhitePixel( dpy, screen );

  gnum1 = 0.33 * HIGHGRAY + 0.67 * LOWGRAY;
  color1.red = color1.blue = color1.green = gnum1;

  if ( XAllocColor( dpy, cmap, &color1 ) ) { dg = color1.pixel; }
  else { printf( "Could not allocate gray value %d\n", gnum1 ); }

  gnum2 = 0.67 * HIGHGRAY + 0.33 * LOWGRAY;

  color2.red = color2.blue = color2.green = gnum2;

  if ( XAllocColor( dpy, cmap, &color2 ) ) { lg = color2.pixel; }
  else { printf( "Could not allocate gray value %d\n", gnum2 ); }

  /*  printf("Gray values: %d %d\n", color1.red, color2.red); */

  n = 0;

  for ( j = 0; j < 3; j++ ) {
    if ( j == 2 ) { g1 = fg; g2 = dg; }

    if ( j == 1 ) { g1 = dg; g2 = lg; }

    if ( j == 0 ) { g1 = lg; g2 = bg; }

    for ( i = 0; i < 8; i++ ) {
      graygc[n] = XCreateGC( dpy, win, 0, 0 );

      if ( i == 0 ) { tile = XCreatePixmapFromBitmapData( dpy, win, t0, tw, th, g1, g2, depth ); }

      if ( i == 1 ) { tile = XCreatePixmapFromBitmapData( dpy, win, t2, tw, th, g1, g2, depth ); }

      if ( i == 2 ) { tile = XCreatePixmapFromBitmapData( dpy, win, t4, tw, th, g1, g2, depth ); }

      if ( i == 3 ) { tile = XCreatePixmapFromBitmapData( dpy, win, t6, tw, th, g1, g2, depth ); }

      if ( i == 4 ) { tile = XCreatePixmapFromBitmapData( dpy, win, t8, tw, th, g1, g2, depth ); }

      if ( i == 5 ) { tile = XCreatePixmapFromBitmapData( dpy, win, t10, tw, th, g1, g2, depth ); }

      if ( i == 6 ) { tile = XCreatePixmapFromBitmapData( dpy, win, t12, tw, th, g1, g2, depth ); }

      if ( i == 7 ) { tile = XCreatePixmapFromBitmapData( dpy, win, t14, tw, th, g1, g2, depth ); }

      XSetTile( dpy, graygc[n], tile );

      XSetFillStyle( dpy, graygc[n], FillTiled );

      n++;
    }
  }

  ngray = n;
}

void SetStippled4x6() {
  int i, depth;
  unsigned long fg, bg, gr;
  XColor color;
  Pixmap tile;
  int tw = 4, th = 6;
  static char t0[]  = {( char )0x00, ( char )0x00, ( char )0x00,
                       ( char )0x00, ( char )0x00, ( char )0x00
                      };
  static char t6[]  = {( char )0x01, ( char )0x04, ( char )0x01,
                       ( char )0x04, ( char )0x01, ( char )0x04
                      };
  static char t12[] = {( char )0x05, ( char )0x0a, ( char )0x05,
                       ( char )0x0a, ( char )0x05, ( char )0x0a
                      };
  static char t18[] = {( char )0xfe, ( char )0xfb, ( char )0xfe,
                       ( char )0xfb, ( char )0xfe, ( char )0xfb
                      };
  static char t24[] = {( char )0xff, ( char )0xff, ( char )0xff,
                       ( char )0xff, ( char )0xff, ( char )0xff
                      };

  depth  = XDefaultDepth( dpy, screen );
  fg = BlackPixel( dpy, screen );
  bg = WhitePixel( dpy, screen );

  if ( XAllocColor( dpy, cmap, &color ) ) {
    gr = color.pixel;
    printf( " gray is %d %d %d\n", color.red, color.blue, color.green );
  }
  else { printf( "Could not allocate gray value\n" ); }

  ngray = 9;

  for ( i = 0; i < ngray; i++ ) {
    graygc[i] = XCreateGC( dpy, win, 0, 0 );

    if ( i == 8 ) {
      tile = XCreatePixmapFromBitmapData( dpy, win, t0, tw, th, gr, bg, depth );
      XSetTile( dpy, graygc[i], tile );
    }

    if ( i == 7 ) {
      tile = XCreatePixmapFromBitmapData( dpy, win, t6, tw, th, gr, bg, depth );
      XSetTile( dpy, graygc[i], tile );
    }

    if ( i == 6 ) {
      tile = XCreatePixmapFromBitmapData( dpy, win, t12, tw, th, gr, bg, depth );
      XSetTile( dpy, graygc[i], tile );
    }

    if ( i == 5 ) {
      tile = XCreatePixmapFromBitmapData( dpy, win, t18, tw, th, gr, bg, depth );
      XSetTile( dpy, graygc[i], tile );
    }

    if ( i == 4 ) {
      tile = XCreatePixmapFromBitmapData( dpy, win, t0, tw, th, fg, gr, depth );
      XSetTile( dpy, graygc[i], tile );
    }

    if ( i == 3 ) {
      tile = XCreatePixmapFromBitmapData( dpy, win, t6, tw, th, fg, gr, depth );
      XSetTile( dpy, graygc[i], tile );
    }

    if ( i == 2 ) {
      tile = XCreatePixmapFromBitmapData( dpy, win, t12, tw, th, fg, gr, depth );
      XSetTile( dpy, graygc[i], tile );
    }

    if ( i == 1 ) {
      tile = XCreatePixmapFromBitmapData( dpy, win, t18, tw, th, fg, gr, depth );
      XSetTile( dpy, graygc[i], tile );
    }

    if ( i == 0 ) {
      tile = XCreatePixmapFromBitmapData( dpy, win, t24, tw, th, fg, gr, depth );
      XSetTile( dpy, graygc[i], tile );
    }

    XSetFillStyle( dpy, graygc[i], FillTiled );
  }
}

GC ChooseColor( float gray, int col0 ) {
  int igray;

  if ( ! grayvalues ) { return gcbg; }

  if ( color ) {
    return graygc[col0];
  }
  else {
    igray = ( 1.0 - gray ) * ngray - 0.49;

    if ( igray > ngray - 1 ) { igray = ngray - 1; }

    if ( igray < 0 ) { igray = 0; }

    return graygc[igray];
  }
}

void showline( Drawable drawable, int x, int y,
               const char s1[], const char s2[], const char s3[] ) {

  char sx[1001];
  strcpy( sx, s1 );
  strcat( sx, s2 );
  strcat( sx, s3 );
  XFillRectangle( dpy, drawable, gcbg, x, igh - y - 9, x + 6 * strlen( sx ), 12 );
  XDrawString( dpy, drawable, labelgc, x, igh - y, sx, strlen( sx ) );

}

void clearline( Drawable drawable, int x, int y ) {
  XFillRectangle( dpy, drawable, gcbg, x, igh - y - 11, igw - x - 20, 14 );
}

void DrawArrow( float x0, float y0, float x1, float y1, float rad,
                const char str[] ) {

  int xx1, yy1, xx0, yy0, xx, yy, dd;

  xx1 = midx + PSFAC * x1;
  yy1 = midy - PSFAC * y1;
  xx0 = midx + PSFAC * x0;
  yy0 = midy - PSFAC * y0;
  xx = 0.9 * xx0 + 0.1 * xx1;
  yy = 0.9 * yy0 + 0.1 * yy1;

  XDrawLine( dpy, drw, shadowgc, xx, yy, xx1, yy1 );
  XDrawLine( dpy, drw, gc, xx0, yy0, xx1, yy1 );

  xx = xx1 + 0.3 * ( xx1 - xx0 );
  yy = yy1 + 0.3 * ( yy1 - yy0 ) + 4;
  dd = strlen( str ) * 6;
  XDrawLine( dpy, drw, labbggc,
             xx + 2 - dd / 2, yy - 4, xx + dd / 2 - 2, yy - 4 );
  XDrawString( dpy, drw, labelgc, xx - dd / 2, yy, str, strlen( str ) );

}

void DrawBall( float gray, float col0, float x, float y, int rad ) {

  int xx, yy, rr, icol;
  GC gcfill;

  rr = PSFAC * rad;
  xx = midx + PSFAC * x;
  yy = midy - PSFAC * y;

  if ( shadow ) {
    XDrawArc( dpy, drw, shadowgc, xx - rr, yy - rr, 2 * rr, 2 * rr, 0, 360 * 64 );
  }

  if ( ! wire ) {
    gcfill = ChooseColor( gray, col0 );
    XFillArc( dpy, drw, gcfill, xx - rr, yy - rr, 2 * rr, 2 * rr, 0, 360 * 64 );
  }

  XDrawArc( dpy, drw, gc, xx - rr, yy - rr, 2 * rr, 2 * rr, 0, 360 * 64 );
}

void DrawStick( float gray, int col0, float m1[6], float m2[6] ) {

  int x1, y1, x2, y2, i, igray;
  float x, y;
  XPoint pp[NPOINTS * 2 + 1];
  GC gcfill;

  if ( bline ) {
    x1 = midx + PSFAC * m1[4];
    y1 = midy - PSFAC * m1[5];
    x2 = midx + PSFAC * m2[4];
    y2 = midy - PSFAC * m2[5];

    if ( shadow ) {
      XDrawLine( dpy, drw, shadowgc, x1, y1, x2, y2 );
    }

    XDrawLine( dpy, drw, gc, x1, y1, x2, y2 );

    return;
  }

  for ( i = 0; i < NPOINTS; i++ ) {
    x = m1[0] * arc[i][0] + m1[2] * arc[i][1] + m1[4];
    y = m1[1] * arc[i][0] + m1[3] * arc[i][1] + m1[5];
    pp[i].x = midx + PSFAC * x;
    pp[i].y = midy - PSFAC * y;
  }

  for ( i = 0; i < NPOINTS; i++ ) {
    x = -m2[0] * arc[i][0] + m2[2] * arc[i][1] + m2[4];
    y = -m2[1] * arc[i][0] + m2[3] * arc[i][1] + m2[5];
    pp[2 * NPOINTS - i - 1].x = midx + PSFAC * x;
    pp[2 * NPOINTS - i - 1].y = midy - PSFAC * y;
  }

  pp[2 * NPOINTS] = pp[0];

  if ( shadow ) {
    XDrawLines( dpy, drw, shadowgc, &pp[0], 2 * NPOINTS + 1, CoordModeOrigin );
  }

  if ( ! wire ) {
    gcfill = ChooseColor( gray, col0 );
    XFillPolygon( dpy, drw, gcfill, &pp[0], 2 * NPOINTS + 1,
                  Nonconvex, CoordModeOrigin );
  }

  XDrawLines( dpy, drw, gc, &pp[0], 2 * NPOINTS + 1, CoordModeOrigin );

}

void LabelBG( float x, float y, float g1, float g2, char str[] ) {
  int xx, yy, dd;

  xx = midx + PSFAC * x;
  yy = midy - PSFAC * y;
  /*  dd=strlen(str)*6+5; */
  /*  XFillRectangle(dpy, drw, gcbg, xx-2, yy-10, dd, 13); */
  dd = strlen( str ) * 6;

  XDrawLine( dpy, drw, labbggc,
             xx + 2 - dd / 2, yy - 4, xx + dd / 2 - 2, yy - 4 );
  XDrawString( dpy, drw, labelgc, xx - dd / 2, yy, str, strlen( str ) );

}




/* print error and exit */
void rx( const char * msg ) {
  printf( "Error: %s\n", msg );
  exit( 1 );
}


/*
   parse strings into individual substrings
   (quoted strings can contain whitespace)    (input parsing)
*/
int parse_args( char str[], char w[8][41] )
//char str[], w[8][41];
{
  int i, num, reading, quoted;
  char *p;

  p = str;
  reading = quoted = 0;
  num = -1;

  while ( *p != 0 ) {
    if ( reading ) {
      if ( ( quoted && ( *p == '\'' ) ) ||
           ( ( !quoted ) && ( *p == ' ' ) ) ) { i++; w[num][i] = 0; reading = 0; }
      else { i++; w[num][i] = *p; }
    }
    else {
      if ( *p != ' ' ) {
        num++;

        if ( *p == '\'' ) { quoted = 1; i = -1; }
        else { i = 0; w[num][i] = *p; }

        reading = 1;
      }
    }

    p++;
  }

  if ( reading ) { i++; w[num][i] = 0; }

  num++;

  return num;
}

/* strip leading and trailing spaces from str , result in str1 */
void strip( char str1[], char str[] )
//char str[], str1[];
{
  int l, i, i1, i2;
  l = strlen( str );

  i1 = 0;

  for ( i = 0; i < l; i++ )
    if ( ( str[i] != ' ' ) && ( str[i] != '\n' ) ) { i1 = i; break; }

  i2 = 0;

  for ( i = l - 1; i >= 0; i-- )
    if ( ( str[i] != ' ' ) && ( str[i] != '\n' ) ) { i2 = i + 1; break; }

  for ( i = i1; i < i2; i++ ) { str1[i - i1] = str[i]; }

  str1[i2 - i1] = 0;

  /*  printf (" l=%d i1=%d i2=%d <%s> <%s>\n", l, i1, i2, str, str1);*/
}

// return 1 if str equals ab. at least nchar chars are compared
int abbrev( const char str[], const char ab[], int nchar ) {
  int i, nc;

  if ( strlen( str ) > strlen( ab ) ) { return 0; }

  nc = strlen( str );

  if ( nc < nchar ) { nc = nchar; }

  for ( i = 0; i < nc; i++ ) if ( str[i] != ab[i] ) { return 0; }

  return 1;
}

/*
    set extension on a file identifier, output=fid1
    force=1 forces change even if fid already has an extension
    force=0 does not change the extension if there already is one
*/
void strext( char fid1[], const char fid[], const char ext[], int force ) {

  int i, l;
  char *p, *q;

  strcpy( fid1, fid );
  l = strlen( fid1 );
  p = fid1;

  for ( i = 0; i < l; i++ )
    if ( fid1[i] == '/' ) { p = fid1 + i; }

  if ( !force ) {
    q = strchr( p, '.' );

    if ( q && ( q != fid1 + strlen( fid1 ) - 1 ) ) { return; }
  }

  if ( !strchr( p, '.' ) ) { strcat( fid1, "." ); }

  q = strchr( p, '.' );

  if ( strlen( ext ) > 0 ) { q++; }

  *q = 0;

  strcat( fid1, ext );

}


/* match  a pattern pat in a str (with wildcards) */
int match( char str[], char pat[] )
//char str[], pat[];
{
  char *p, *s;
  p = pat;
  s = str;

  while ( *p != 0 ) {

    if ( *p == '*' ) {         /* found wildcard '*' in pattern */
      p++;

      while ( *p == '*' ) { p++; }

      if ( *p == 0 ) { return 1; } /* trailing '*' matches all */

      for ( ;; ) {             /* find match to char after '*' */
        if ( *s == 0 ) { return 0; }

        if ( ( *s == *p ) || ( *p == '+' ) )
          if ( match( s + 1, p + 1 ) ) { return 1; }   /* ok if rest matches */

        s++;
      }

      return 0;                /* tried all cases but none worked */
    }

    else {                     /* no wildcard -- char must match */
      if ( *s == 0 ) { return 0; }

      if ( ( *p != *s ) && ( *p != '+' ) ) { return 0; }

      s++;
    }

    p++;
  }

  if ( *s != 0 ) { return 0; }     /* pattern but not string exhausted */

  return 1;
}

/* ----- get_extent ----- */
void get_extent( float * xx1, float * xx2, float * yy1, float * yy2, float * zz1, float * zz2 )
//float *xx1, *xx2, *yy1, *yy2, *zz1, *zz2;
{
  float big, x1, x2, y1, y2, z1, z2;
  int i;

  big = 1000000;

  if ( nbas == 0 ) { big = 0; }

  x1 = y1 = z1 = big;

  x2 = y2 = z2 = -big;

  for ( i = 0; i < nbas; i++ ) {
    if ( p[i][0] < x1 ) { x1 = p[i][0]; }

    if ( p[i][0] > x2 ) { x2 = p[i][0]; }

    if ( p[i][1] < y1 ) { y1 = p[i][1]; }

    if ( p[i][1] > y2 ) { y2 = p[i][1]; }

    if ( p[i][2] < z1 ) { z1 = p[i][2]; }

    if ( p[i][2] > z2 ) { z2 = p[i][2]; }
  }

  *xx1 = x1 + center[0];

  *xx2 = x2 + center[0];
  *yy1 = y1 + center[1];
  *yy2 = y2 + center[1];
  *zz1 = z1 + center[2];
  *zz2 = z2 + center[2];
}

/* ----- atompos: position and radius on paper for an atom --- */

/* zp = position on paper (2 dimensions) : also Projektion in die xy-Ebene */

void atompos( float fac, float p[3], float rad,
              float zp[2], float *zr ) {
  float y[3], q[3], v1[3], v2[3];
  float xxx, za1, za2, zb1, zb2, a, b;

  if ( pmode == 1 ) {
    zp[0] = fac * p[0];
    zp[1] = fac * p[1];
    *zr = fac * rad;
    *zr = MAXRAD;

    if ( dist0 - p[2] > 0 ) { *zr = fac * rad * dist0 / ( dist0 - p[2] ); }

    if ( *zr > MAXRAD ) { *zr = MAXRAD; }

    return;
  }

  svm.vscal( p, 1.0, q );

  q[2] = q[2] - dist;
  svm.vscal( p, 1.0, y );
  xxx = -svm.sp( y, q ) / svm.sp( q, q );
  svm.vsum( y, q, 1.0, xxx, y );

  if ( svm.sp( y, y ) <= 1e-3 ) { y[0] = 1.0; y[1] = 0.0; y[2] = 0.0; }

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
  zp[0] = 0.5 * ( za1 + zb1 );
  zp[1] = 0.5 * ( za2 + zb2 );
  *zr = ( zb1 - za1 ) * ( zb1 - za1 ) + ( zb2 - za2 ) * ( zb2 - za2 );
  *zr = 0.5 * sqrt( *zr );
}

/* ----- readclusterdata ---- */

// forward declare
int readclusterline( char str[], int helpme );

/* read .bs file , i.e. molecule info */

int readclusterdata( char infile[] ) {
  FILE *fp;
  char str[257], token[81];
  char xxx [81];
  char *p;
  int l, nn;

  if ( ( fp = fopen( infile, "r" ) ) == NULL ) { return 0; }

  fgets( str, 257, fp );

  while ( !feof( fp ) ) {
    l = strlen( str );
    str[l + 1] = '\0';
    strcpy( token, "SNOT" );
    sscanf( str, "%s", token );

    if ( ( p = strstr( str, "frame" ) ) ) {
      if ( nframe * nbas > FBMAX ) { rx( "increase internal dimension FBMAX" ); }

      if ( nframe > NFRMAX ) { rx( "increase internal dimension NFRMAX" ); }

      p = p + 6;

      sprintf( frstr[nframe], "%-80s", p );

      for ( std::size_t i = 0; i < strlen( frstr[nframe] ); i++ ) {
        if ( frstr[nframe][i] == '\n' ) { frstr[nframe][i] = '\0'; }
      }

      for ( unsigned int i = 0; i < nbas; i++ ) {
        nn = nframe * nbas + i;
        fscanf( fp, "%f", &frame[0][nn] );
        fscanf( fp, "%f", &frame[1][nn] );
        fscanf( fp, "%f", &frame[2][nn] );
      }

      nframe++;

      if ( nframe % 50 == 0 ) {
        sprintf( xxx, ": frame %d         ", nframe );
        showline( win, 10, 8, "Reading ", inmv, xxx );
        XFlush( dpy );
      }

      fgets( str, 257, fp );
    }
    else {
      readclusterline( str, 0 );
    }

    fgets( str, 257, fp );
  }

  sprintf( curf, "%s ", inf );

  return 1;
}

/*
    readclusterline
    reads in commands from command line: this is where all
    'long' commands are presently defined

*/
int readclusterline( char str[], int helpme ) {
  char token[81], cname[81];
  float atcors[4][3];
  float rval, bval, gval, gray;
  int l, n, i, k, nn, a[4];

  l = strlen( str );

  if ( l < 1 ) { return 0; }

  strcpy( token, "SNOT" );

  sscanf( str, "%s", token );

  if ( !strcmp( token, "SNOT" ) ) { return 0; }  /* empty line */

  if ( token[0] == '*' ) { return 0; }     /* comment -- no error */


  /* print help messages */

  if ( helpme ) {
    if ( abbrev( token, "spec", 4 ) ) {
      sprintf( gmsg, "Usage: spec label radius color  - define species" );
    }
    else if ( abbrev( token, "time", 4 ) ) {
      sprintf( gmsg, "Usage: time msecs - set duration of each frame in animation" );
    }
    else if ( abbrev( token, "tell", 4 ) ) {
      sprintf( gmsg, "Usage: tell i ... - tell xyz, distance, angle, or torsion" );
    }
    else if ( abbrev( token, "atom", 4 ) ) {
      sprintf( gmsg, "Usage: atom label x y z  - place atom at (x,y,z)" );
    }
    else if ( abbrev( token, "bonds", 5 ) ) {
      sprintf( gmsg, "Usage: bonds pat1 pat2 min max rad color  - select bonds" );
    }
    else if ( abbrev( token, "light", 5 ) ) {
      sprintf( gmsg, "Usage: light vx vy vz  - light along vector in bw mode)" );
    }
    else if ( abbrev( token, "inc", 3 ) ) {
      sprintf( gmsg, "Usage: inc degrees  - angle increment for rotation" );
    }
    else if ( abbrev( token, "dist", 5 ) ) {
      sprintf( gmsg, "Usage: dist d  - set distance for perspective" );
    }
    else if ( abbrev( token, "frm", 3 ) ) {
      sprintf( gmsg, "Usage: frm n  - goto frame n" );
    }
    else if ( abbrev( token, "step", 4 ) ) {
      sprintf( gmsg, "Usage: step n  - set step for frames" );
    }
    else if ( abbrev( token, "gramp", 5 ) ) {
      sprintf( gmsg, "Usage: gramp slope [middle]  - set gray ramp in bw mode)" );
    }
    else if ( abbrev( token, "scale", 5 ) ) {
      sprintf( gmsg, "Usage: scale x  - set overall scale factor" );
    }
    else if ( abbrev( token, "rfac", 4 ) ) {
      sprintf( gmsg, "Usage: rfac x  - scale all sphere radii by x" );
    }
    else if ( abbrev( token, "bfac", 4 ) ) {
      sprintf( gmsg, "Usage: bfac x  - scale all bond radii by x" );
    }
    else if ( abbrev( token, "pos", 3 ) ) {
      sprintf( gmsg, "Usage: pos px py  - set position on page" );
    }
    else if ( abbrev( token, "dpos", 4 ) ) {
      sprintf( gmsg, "Usage: dpos x  - set increment for position" );
    }
    else if ( abbrev( token, "time", 4 ) ) {
      sprintf( gmsg, "Usage: time x  - set duration time for frame [msec]" );
    }

    /* additions for builder commands */

    else if ( abbrev( token, "azmat", 5 ) ) {
      sprintf( gmsg, "Usage: add an atom with distance, angle, dihedral " );
    }

    else if ( abbrev( token, "axyz", 4 ) ) {
      sprintf( gmsg, "Usage: add atom  at xyz offset" );
    }

    else {
      sprintf( gmsg, "No help available on %s", token );
    }

    return 0;
  }

  /* the following commands read from command line as well as from infile.bs ! */

  /* Extensions for 'tell' , .i.e. distance, angle and dihedral */

  if ( !strcmp( token, "tell" ) ) {
    a[0] = a[1] = a[2] = a[3] = -1;
    sscanf( str, "%*s %d %d %d %d", &a[0], &a[1], &a[2], &a[3] );
    n = 0;

    for ( i = 0; i < 4; i++ ) {
      if ( a[i] > 0 ) {
        if ( a[i] > nbas ) {
          sprintf( gmsg,
                   "Atom number %d  at position %d is too large", a[i], i + 1 );
          return 0;
        }

        for ( k = 0; k < 3; k++ ) {
          if ( iframe == 0 ) {atcors[i][k] = atom[a[i] - 1].pos[k];}
          else {atcors[i][k] = frame[k][iframe * nbas + a[i] - 1];}
        }

        n++;

      }
      else { break; }
    }

    /* n = number of centers */
    if ( n == 1 ) {
      sprintf( gmsg, "Coordinates X[%d]=%.4f Y[%d]=%.4f Z[%d]=%.4f",
               a[0], atcors[0][0], a[0], atcors[0][1], a[0],
               atcors[0][2] );
    }
    else if ( n == 2 ) {
      sprintf( gmsg, "Distance %d:%d is %.4f A", a[0], a[1],
               svm.distance( atcors[0], atcors[1] ) );
    }
    else if ( n == 3 ) {
      sprintf( gmsg, "Angle %d:%d:%d is %.4f deg", a[0], a[1],  a[2],
               svm.angle( atcors[0], atcors[1], atcors[2] ) );
    }
    else if ( n == 4 ) {
      sprintf( gmsg, "Torsion angle %d:%d:%d:%d is %.4f deg",
               a[0], a[1], a[2], a[3],
               svm.torsion( atcors[0], atcors[1], atcors[2], atcors[3] ) );
    }

    return 0;
  }


  /* read in specification of an atom */

  if ( !strcmp( token, "spec" ) ) {
    sscanf( str, "%*s %s %f %n", spec[nspec].lab, &spec[nspec].rad, &n );
    strip( spec[nspec].cname, str + n );
    nspec++;

    if ( nspec > NSPMAX ) { rx( "increase internal dimension NSPMAX" ); }

    return 2;
  }

  /* read in atom */

  if ( !strcmp( token, "atom" ) ) {
    atom[nbas].pol[0] = 0;
    atom[nbas].pol[1] = 0;
    atom[nbas].pol[2] = 0;
    sscanf( str, "%*s %s %f %f %f %f %f %f", atom[nbas].lab,
            &atom[nbas].pos[0], &atom[nbas].pos[1], &atom[nbas].pos[2],
            &atom[nbas].pol[0], &atom[nbas].pol[1], &atom[nbas].pol[2] );
    nbas++;

    if ( nbas > NAMAX ) { rx( "increase internal dimension NAMAX" ); }

    return 2;
  }

  /* read in bond */

  if ( !strcmp( token, "bonds" ) ) {
    sscanf( str, "%*s %s %s %f %f %f %n", bonds[nbonds].lab1,
            bonds[nbonds].lab2, &bonds[nbonds].min, &bonds[nbonds].max,
            &bonds[nbonds].rad, &n );
    strip( bonds[nbonds].cname, str + n );
    nbonds++;

    if ( nbonds > NBTMAX ) { rx( "increase internal dimension NBTMAX" ); }

    return 2;
  }

  /* read in rotation (?) matrix */

  if ( !strcmp( token, "tmat" ) ) {
    sscanf( str, "%*s %f %f %f %f %f %f %f %f %f",
            &tmat[0][0], &tmat[0][1], &tmat[0][2],
            &tmat[1][0], &tmat[1][1], &tmat[1][2],
            &tmat[2][0], &tmat[2][1], &tmat[2][2] );
    return 1;
  }

  if ( !strcmp( token, "dist" ) || !strcmp( token, "d" ) ) {
    sscanf( str, "%*s %f", &dist0 );
    return 1;
  }

  if ( !strcmp( token, "time" ) ) {
    sscanf( str, "%*s %d", &msecs );
    return 1;
  }

  if ( !strcmp( token, "inc" ) ) {sscanf( str, "%*s %f", &dalfa ); return 3;}

  if ( !strcmp( token, "frm" ) ) {
    sscanf( str, "%*s %d", &nn );

    if ( nn > nframe || nn < 1 ) {
      sprintf( emsg, "No frame %d available", nn );
      return 0;
    }

    if ( iframe == nn - 1 ) { return 0; }

    iframe = nn - 1;

    return 2;
  }

  if ( !strcmp( token, "light" ) ) {
    light[0] = light[1] = light[2] = 0.0;
    sscanf( str, "%*s %f %f %f", &light[0], &light[1], &light[2] );

    if ( light[0]*light[0] + light[1]*light[1] + light[2]*light[2] < .01 ) {
      sprintf( gmsg, "Use standard coloring" );
      gmode = G_STD;;
      return 2;
    }

    if ( color ) {
      sprintf( emsg, "light: only works in b/w mode" );
      return 0;
    }

    gmode = G_LIGHT;

    return 2;
  }

  if ( !strcmp( token, "step" ) )  {sscanf( str, "%*s %d", &fstep ); return 1;}

  if ( !strcmp( token, "scale" ) ) {
    sscanf( str, "%*s %f", &scale );
    scale = scale * igs;
    return 1;
  }

  if ( !strcmp( token, "rfac" ) )  {sscanf( str, "%*s %f", &radfac ); return 1;}

  if ( !strcmp( token, "bfac" ) )  {sscanf( str, "%*s %f", &bndfac ); return 1;}

  if ( !strcmp( token, "amp" ) )   {sscanf( str, "%*s %f", &amp ); return 2;}

  if ( !strcmp( token, "pos" ) )
  {sscanf( str, "%*s %f %f", &taux, &tauy ); return 1;}

  if ( !strcmp( token, "dpos" ) )   {
    sscanf( str, "%*s %f", &dtaux );
    dtauy = dtaux;
    chginfo = 1;
    return 0;
  }

  if ( !strcmp( token, "gramp" ) )   {
    gslope = gz0 = 0;
    sscanf( str, "%*s %f %f", &gslope, &gz0 );

    if ( gslope * gslope < 0.1 ) {
      sprintf( gmsg, "Use standard coloring" );
      gmode = G_STD;
      return 2;
    }

    if ( color ) {
      sprintf( emsg, "gramp: only works in b/w mode" );
      return 0;
    }

    gmode = G_RAMP;

    return 2;
  }

  if ( !strcmp( token, "switches" ) )   {
    sscanf( str, "%*s %d %d %d %d %d %d %d %d %d",
            &usepixmap, &numbers, &grayvalues, &bline, &wire,
            &withbonds, &recenter, &pmode, &shadow );
    return 2;
  }

  sprintf( emsg, "Undefined command: %s", token );

  if ( startup ) { printf( "Cannot understand line: %s\n", str ); }

  return 0;
}

/* ----- writeclusterdata ---- */
/* i.e. write bs file */

void writeclusterdata( char outfile[], int svstep, int svrgb ) {
  FILE *fp;
  char nm[81];
  int i, nn, n, nfrm;
  time_t  ltime;
  char    timestr[41];

  if ( ( fp = fopen( outfile, "w" ) ) == NULL ) {
    sprintf( emsg, "Cannot open file %s\n", outfile );
    return;
  }

  time( &ltime );

  strcpy( timestr,  ctime( &ltime ) );
  timestr[24] = 0;
  fprintf( fp, "* Saved %s from %s\n\n", timestr, inf );

  for ( i = 0; i < nbas; i++ )
    fprintf( fp, "atom %6s  %10.3f %10.3f %10.3f \n",
             atom[i].lab, atom[i].pos[0], atom[i].pos[1],
             atom[i].pos[2] );

  fprintf( fp, "\n" );

  for ( i = 0; i < nspec; i++ )
    if ( svrgb || reverse )
      fprintf( fp, "spec %6s %10.3f   %.2f %.2f %.2f\n",
               spec[i].lab, spec[i].rad, spec[i].r, spec[i].g, spec[i].b );
    else
      fprintf( fp, "spec %6s %10.3f   %s\n",
               spec[i].lab, spec[i].rad, spec[i].cname );

  fprintf( fp, "\n" );

  for ( i = 0; i < nbonds; i++ )
    if ( svrgb || reverse )
      fprintf( fp, "bonds %5s %5s %8.3f %8.3f %8.3f   %.2f %.2f %.2f\n",
               bonds[i].lab1, bonds[i].lab2, bonds[i].min, bonds[i].max,
               bonds[i].rad, bonds[i].r, bonds[i].g, bonds[i].b );
    else
      fprintf( fp, "bonds %5s %5s %8.3f %8.3f %8.3f   %s\n",
               bonds[i].lab1, bonds[i].lab2, bonds[i].min, bonds[i].max,
               bonds[i].rad, bonds[i].cname );

  fprintf( fp, "\ntmat %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
           tmat[0][0], tmat[0][1], tmat[0][2], tmat[1][0], tmat[1][1],
           tmat[1][2], tmat[2][0], tmat[2][1], tmat[2][2] );

  fprintf( fp, "dist  %8.3f\n", dist0 );

  fprintf( fp, "inc   %8.3f\n", dalfa );

  fprintf( fp, "scale %8.3f\nrfac %.2f\nbfac %.2f\n", scale, radfac, bndfac );

  fprintf( fp, "pos %8.3f %8.3f\n", taux, tauy );

  if ( gmode == G_RAMP ) {
    fprintf( fp, "gramp %8.3f %8.3f\n", gslope, gz0 );
  }

  if ( gmode == G_LIGHT ) {
    fprintf( fp, "light %8.3f %8.3f %8.3f\n", light[0], light[1], light[2] );
  }

  fprintf( fp, "switches %d %d %d %d %d %d %d %d %d\n",
           usepixmap, numbers, grayvalues, bline, wire,
           withbonds, recenter, pmode, shadow );

  fclose( fp );

  sprintf( gmsg, "Saved in %s", outfile );

  if ( nframe > 1 ) {
    strext( nm, outfile, "mv", 1 );

    if ( ( fp = fopen( nm, "w" ) ) == NULL ) {
      sprintf( emsg, "Cannot open file %s\n", nm );
      return;
    }

    nfrm = 1;

    for ( i = 1; i < nframe; i = i + svstep ) {
      nfrm++;
      fprintf( fp, "frame %s\n", frstr[i] );

      for ( n = 0; n < nbas; n++ ) {
        nn = i * nbas + n;
        fprintf( fp, "%.3f %.3f %.3f ",
                 frame[0][nn], frame[1][nn], frame[2][nn] );
      }

      fprintf( fp, "\n\n" );
    }

    fclose( fp );

    sprintf( gmsg, "Saved %d frames in %s and %s", nfrm, outfile, nm );

  }

  return;
}

/* ----- parse_all_colors ------- */
void parse_all_colors() {
  int i;

  for ( i = 0; i < nspec; i++ )
    parse_color( spec[i].cname, &spec[i].r, &spec[i].g, &spec[i].b,
                 &spec[i].gray );

  for ( i = 0; i < nbonds; i++ )
    parse_color( bonds[i].cname, &bonds[i].r, &bonds[i].g, &bonds[i].b,
                 &bonds[i].gray );
}

/* ----- set_auto_colors ------ */
void set_auto_colors() {
  int k, i, up1;
  char id[21];
  char *p;

  for ( k = 0; k < nspec; k++ ) {
    /* separate off the identifier and change to upper case */
    up1 = isupper( spec[k].lab[0] );
    p = &spec[k].lab[0];
    i = 0;

    while ( *p != 0 ) {
      id[i] = toupper( *p );
      p++;
      i++;

      if ( !isalpha( *p ) ) { break; }

      if ( up1 && isupper( *p ) ) { break; }
    }

    id[i] = 0;

    /* Set default colors here. Use rgb values or color name */

    strcpy( spec[k].cname, ".68 .85 .90" );

    if ( !strcmp( id, "H" ) ) { strcpy( spec[k].cname, "1.0 0.2 0.2" ); }

    if ( !strcmp( id, "C" ) ) { strcpy( spec[k].cname, "0.65 0.7 0.7" ); }

    /*    if (!strcmp(id,"O"))   strcpy(spec[k].cname, "0.2 0.2 1.0"); */
    if ( !strcmp( id, "O" ) ) { strcpy( spec[k].cname, "blue" ); }

    if ( !strcmp( id, "N" ) ) { strcpy( spec[k].cname, "0.8 0.0 1.0" ); }

    if ( !strcmp( id, "P" ) ) { strcpy( spec[k].cname, "0.0 0.8 0.0" ); }

    if ( !strcmp( id, "CL" ) ) { strcpy( spec[k].cname, "0.9 0.9 0.7" ); }

    if ( !strcmp( id, "TI" ) ) { strcpy( spec[k].cname, "0.2 1.0 1.0" ); }
  }

  for ( k = 0; k < nbonds; k++ ) {
    strcpy( bonds[k].cname, "0.8 0.8 0.8" );
  }
}

/* ----- ball_list ----- */
/* setup list of balls */
int ball_list( struct ballstr ball[], int jpr ) {
  int i, j, k, m;
  float top, bot, sp;

  for ( i = 0; i < nbas; i++ ) {
    k = -1;

    for ( j = 0; j < nspec; j++ ) if ( !strcmp( spec[j].lab, atom[i].lab ) ) { k = j; }

    if ( k == -1 ) {
      if ( !startup ) {
        sprintf( emsg, "Undefined species %s ", atom[i].lab );
      }
      else {
        printf( "Undefined species %s\n", atom[i].lab );
      }

      continue;
    }

    for ( m = 0; m < 3; m++ ) { ball[i].pos[m] = atom[i].pos[m] + amp * atom[i].pol[m]; }

    ball[i].rad = spec[k].rad;

    ball[i].gray = spec[k].gray;

    ball[i].r = spec[k].r;

    ball[i].g = spec[k].g;

    ball[i].b = spec[k].b;

    strcpy( ball[i].lab, spec[k].lab );

    ball[i].special = 0;

    ball[i].col = spec[k].col;

    if ( spec[k].gray < -0.1 ) {ball[i].gray = 1.0; ball[i].special = 1;}
  }

  if ( gmode == G_LIGHT ) {
    top = -1000.0;
    bot = 1000.0;

    for ( i = 0; i < nbas; i++ ) {
      if ( ! ball[i].special ) {
        sp = ball[i].pos[0] * light[0] + ball[i].pos[1] * light[1]
             + ball[i].pos[2] * light[2];

        if ( sp > top ) { top = sp; }

        if ( sp < bot ) { bot = sp; }
      }
    }

    for ( i = 0; i < nbas; i++ ) {
      if ( ! ball[i].special ) {
        sp = ball[i].pos[0] * light[0] + ball[i].pos[1] * light[1]
             + ball[i].pos[2] * light[2];
        ball[i].gray = ( sp - bot ) / ( top - bot );
      }
    }
  }

  return nbas;
}

/* ----- stick_list ------ */
/* set up list of sticks ? */


int stick_list( struct ballstr ball[], struct stickstr stick[] ) {
  int i, j, k, l, m, nbond, kb;
  float dis, dd;

  i = -1;

  for ( k = 0; k < nbas; k++ )
    for ( l = k + 1; l < nbas; l++ ) {
      kb = -1;

      for ( j = 0; j < nbonds; j++ ) {
        if ( match( ball[k].lab, bonds[j].lab1 ) &&
             match( ball[l].lab, bonds[j].lab2 ) ) { kb = j; }

        if ( match( ball[l].lab, bonds[j].lab1 ) &&
             match( ball[k].lab, bonds[j].lab2 ) ) { kb = j; }

        /*        if( (!strcmp(bonds[j].lab1,ball[k].lab)) &&
                    (!strcmp(bonds[j].lab2,ball[l].lab)) ) kb=j;
                if( (!strcmp(bonds[j].lab1,ball[l].lab)) &&
                    (!strcmp(bonds[j].lab2,ball[k].lab)) ) kb=j; */
      }

      if ( kb > -1 ) {
        dis = 0.0;

        for ( m = 0; m < 3; m++ ) {
          dd = ball[k].pos[m] - ball[l].pos[m];
          dis = dis + dd * dd;
        }

        dis = alat * sqrt( dis );

        if ( ( dis >= bonds[kb].min ) && ( dis <= bonds[kb].max ) ) {
          i++;

          if ( i > NBMAX ) { rx( "increase internal dimension NBMAX" ); }

          stick[i].start = k;

          stick[i].end = l;

          stick[i].rad = bonds[kb].rad;

          stick[i].gray = bonds[kb].gray;

          stick[i].col = bonds[kb].col;
        }
      }
    }

  nbond = i + 1;

  return nbond;
}


/* ----- duplicate_atoms ---- */
int duplicate_atoms( float sh[3][6], int helpme )
//float sh[3][6];
//int helpme;
{
  int k, l, iv, nbas1, ndup, fr, nn, nn1;
  float cx, cy, cz;

  if ( helpme ) {
    sprintf( gmsg, "Usage: dup vx vy vz  - duplicate shifted by vector" );
    return 0;
  }

  ndup = 0;

  for ( iv = 0; iv < 6; iv++ ) {
    cx = sh[0][iv];
    cy = sh[1][iv];
    cz = sh[2][iv];

    if ( cx * cx + cy * cy + cz * cz > 0.001 ) { ndup++; }
  }

  if ( ndup == 0 ) {
    sprintf( emsg, "Cannot dup for (0,0,0)" );
    return 0;
  }

  nbas1 = nbas * ( ndup + 1 );

  if ( nframe * nbas1 > FBMAX ) {
    sprintf( emsg, "Cannot dup, internal dimension FBMAX too small" );
    return 0;
  }

  if ( nbas * ( ndup + 1 ) > NAMAX ) {
    sprintf( emsg, "Cannot dup, internal dimension NAMAX too small" );
    return 0;
  }

  for ( iv = 0; iv < ndup; iv++ ) {
    for ( k = 0; k < nbas; k++ ) {
      l = k + nbas * ( iv + 1 );
      strcpy( atom[l].lab, atom[k].lab );
      atom[l].pos[0] = atom[k].pos[0] + sh[0][iv];
      atom[l].pos[1] = atom[k].pos[1] + sh[1][iv];
      atom[l].pos[2] = atom[k].pos[2] + sh[2][iv];
      atom[l].pol[0] = atom[k].pol[0];
      atom[l].pol[1] = atom[k].pol[1];
      atom[l].pol[2] = atom[k].pol[2];
    }
  }

  for ( fr = nframe - 1; fr >= 0; fr-- ) {
    for ( k = 0; k < nbas; k++ ) {
      nn = fr * nbas + k;
      nn1 = fr * nbas1 + k;
      frame[0][nn1] = frame[0][nn];
      frame[1][nn1] = frame[1][nn];
      frame[2][nn1] = frame[2][nn];

      for ( iv = 0; iv < ndup; iv++ ) {
        nn1 = nn1 + nbas;
        frame[0][nn1] = frame[0][nn] + sh[0][iv];
        frame[1][nn1] = frame[1][nn] + sh[1][iv];
        frame[2][nn1] = frame[2][nn] + sh[2][iv];
      }
    }
  }

  sprintf( gmsg, "Increased from %d to %d atoms", nbas, nbas1 );

  nbas = nbas1;
  return 1;
}

/* ----- add atom: add one atom to molecule, dx,dy,dz increments starting from one atom ---- */
int add_atom_xyz( int at, char * label, float dx, float dy, float dz, int helpme )
//float dx, dy, dz;  /* displacement from at */
//int   at;          /* atom name to start from */
//int helpme;
//char *label;
{
  int k, l, iv, nbas1, ndup, fr, nn, nn1;
  float cx, cy, cz;

  /* print help */

  if ( helpme ) {
    sprintf( gmsg, "Usage: axyz at new dx dy dz - add atom new at dx,dy,dz from at" );
    return 0;
  }

  nbas1 = nbas + 1; /* increment number of atoms by one */

  if ( nframe * nbas1 > FBMAX ) {
    sprintf( emsg, "Cannot add_atom_xyz, internal dimension FBMAX too small" );
    return 0;
  }

  if ( nbas1 > NAMAX ) {
    sprintf( emsg, "Cannot add_atom_xyz, internal dimension NAMAX too small" );
    return 0;
  }

  /*
  printf("%d %s %f %f %f ",at,label,dx,dy,dz);
  */

  l = nbas1 - 1;

  atom[l].pos[0] = atom[at].pos[0] + dx;

  atom[l].pos[1] = atom[at].pos[1] + dy;

  atom[l].pos[2] = atom[at].pos[2] + dz;

  strcpy( atom[l].lab, label );

  /* copy atom to all frames too */

  for ( fr = nframe - 1; fr >= 0; fr-- ) {
    nn = fr * nbas + ( at );
    nn1 = fr * nbas1 + l;
    frame[0][nn1] = frame[0][nn] + dx;
    frame[1][nn1] = frame[1][nn] + dy;
    frame[2][nn1] = frame[2][nn] + dz;
  }

  sprintf( gmsg, "Increased from %d to %d atoms", nbas, nbas1 );

  nbas = nbas1;
  return 1;
}

/* ----- add atom: add one atom to molecule, dx,dy,dz increments starting from one atom ---- */
int add_atom_zmat( char * label, int a1, float length, int a2, float angle, int a3, float dihedral, int helpme )
//float length, angle, dihedral;  /* displacement from at */
//int   a1, a2, a3;             /* atom name to start from */
//int helpme;
//char *label;
{
  int k, l, iv, nbas1 = 0, ndup, fr, nn, nn1;
  float p[3];

  /* print help */

  if ( helpme ) {
    sprintf( gmsg, "Usage: azmat at new dist angle dihedral - add atom" );
    return 0;
  }

  if ( nframe * nbas1 > FBMAX ) {
    sprintf( emsg, "Cannot add_atom_xyz, internal dimension FBMAX too small" );
    return 0;
  }

  if ( nbas1 > NAMAX ) {
    sprintf( emsg, "Cannot add_atom_xyz, internal dimension NAMAX too small" );
    return 0;
  }

  /* calculate dx ,dy dz from zmatrix input */

  nbas1 = nbas + 1; /* increment number of atoms by one */

  svm.zmat2xyz( atom[a1].pos, atom[a2].pos, atom[a3].pos, length, angle, dihedral, p );

  l = nbas1 - 1;

  atom[l].pos[0] = p[0];

  atom[l].pos[1] = p[1];

  atom[l].pos[2] = p[2];

  strcpy( atom[l].lab, label );

  /* copy atom to all frames too */

  for ( fr = nframe - 1; fr >= 0; fr-- ) {

    nn1 = fr * nbas1 + l;
    frame[0][nn1] = p[0];
    frame[1][nn1] = p[0];
    frame[2][nn1] = p[0];
  }

  sprintf( gmsg, "Increased from %d to %d atoms", nbas, nbas1 );

  nbas = nbas1;
  return 1;
}

/* ----- cut_atoms ---- */
int cut_atoms( float cut[3], float cut1, float cut2, int helpme )
//float cut[3], cut1, cut2;
//int helpme;
{
  int i, j, nbas1, fr, nn, nn1;
  float c2, s, fuzz;
  int ip[NAMAX];

  if ( helpme ) {
    sprintf( gmsg, "Usage: cut vx vy vz a b  - cut along vector at a and b" );
    return 0;
  }

  fuzz = 0.001;

  c2 = cut[0] * cut[0] + cut[1] * cut[1] + cut[2] * cut[2];

  if ( c2 < 0.01 ) {
    sprintf( emsg, "cut: invalid vector (%.2f,%.2f,%.2f)",
             cut[0], cut[1], cut[2], cut1, cut2 );
    return 0;
  }

  nbas1 = 0;

  for ( i = 0; i < nbas; i++ ) {
    s = atom[i].pos[0] * cut[0] + atom[i].pos[1] * cut[1] + atom[i].pos[2] * cut[2];

    if ( c2 * cut1 - fuzz < s && s < c2 * cut2 + fuzz ) {
      ip[nbas1] = i;
      nbas1 = nbas1 + 1;
    }
  }

  if ( nbas1 == nbas ) {
    sprintf( gmsg, "No atoms were cut" );
    return 0;
  }

  for ( i = 0; i < nbas1; i++ ) {
    j = ip[i];
    atom[i].pos[0] = atom[j].pos[0];
    atom[i].pos[1] = atom[j].pos[1];
    atom[i].pos[2] = atom[j].pos[2];
    atom[i].pol[0] = atom[j].pol[0];
    atom[i].pol[1] = atom[j].pol[1];
    atom[i].pol[2] = atom[j].pol[2];
    strcpy( atom[i].lab, atom[j].lab );
  }

  for ( fr = 0; fr < nframe; fr++ ) {
    for ( i = 0; i < nbas1; i++ ) {
      j = ip[i];
      nn  = fr * nbas + j;
      nn1 = fr * nbas1 + i;
      frame[0][nn1] = frame[0][nn];
      frame[1][nn1] = frame[1][nn];
      frame[2][nn1] = frame[2][nn];
    }
  }

  sprintf( gmsg, "Reduced from %d to %d atoms", nbas, nbas1 );

  nbas = nbas1;
  return 1;
}

/* ----- select bonds to plot ----- */
int selectbonds( int natom, struct ballstr ball[], float blen,
                 float rad, float gray, struct stickstr stick[] ) {
  int i, k, l, m, nbond;
  float dis, dd;

  i = -1;

  for ( k = 0; k < natom; k++ )
    for ( l = k + 1; l < natom; l++ ) {
      dis = 0.0;

      for ( m = 0; m < 3; m++ ) {
        dd = ball[k].pos[m] - ball[l].pos[m];
        dis = dis + dd * dd;
      }

      dis = sqrt( dis );

      if ( dis <= blen ) {
        i++;
        stick[i].start = k;
        stick[i].end = l;
        stick[i].rad = rad;
        stick[i].gray = gray;
      }

    }

  printf( "bond   start  end    radius    gray\n" );

  nbond = i + 1;

  for ( i = 0; i < nbond; i++ )
    printf( "%3d    %3d   %3d    %7.3f  %6.2f\n",
            i + 1, stick[i].start + 1, stick[i].end + 1,
            stick[i].rad, stick[i].gray );

  return nbond;
}

void dbond( float gray, float m1[6], float m2[6] ) {
  float dfac = 0.8;
  float rfac = 2.0;

  float r, ax, ay, bx, by, a1, b1, dx, dy, alf = 0, d1, bb, x;

  m2[4] = dfac * m2[4] + ( 1 - dfac ) * m1[4];
  m2[5] = dfac * m2[5] + ( 1 - dfac ) * m1[5];

  ax = m2[0];
  ay = m2[1];
  bx = m2[2];
  by = m2[3];

  b1 = sqrt( bx * bx + by * by );
  a1 = sqrt( ax * ax + ay * ay );
  r = rfac * b1;
  m2[0] = -r * ax / a1;
  m2[1] = -r * ay / a1;
  m2[2] = r * bx / b1;
  m2[3] = r * by / b1;

  dx = m2[4] - m1[4];
  dy = m2[5] - m1[5];

  if ( m2[0]*dx < 0 || m2[1]*dy < 0 ) {
    m1[0] = -m1[0];
    m1[1] = -m1[1];
    m1[2] = -m1[2];
    m1[3] = -m1[3];
    m2[0] = -m2[0];
    m2[1] = -m2[1];
    m2[2] = -m2[2];
    m2[3] = -m2[3];
  }

  d1 = sqrt( dx * dx + dy * dy );

  bb = sqrt( m1[2] * m1[2] + m1[3] * m1[3] );

  if ( r - bb < d1 ) {
    x = bb * d1 / ( r - bb );
    alf = asin( bb / x ) * 57.3;

    if ( hardcopy ) {
      svps.hardcopy_xdbond( gray, m1, m2, alf );
    }
    else {
      printf( "PSWxdbond.. not changed yet\n" );
    }

    /*      PSWxdbond (gray, m1, m2, alf); */
  }
  else {
    if ( hardcopy ) {
      svps.hardcopy_ydbond( gray, m1, m2, alf );
    }
    else {
      printf( "PSWydbond.. not changed yet\n" );
    }

    /*      PSWydbond (gray, m1, m2, alf); */
  }
}


/* ----- getframe ------ */
void getframe( struct ballstr ball[], int fnum ) {
  register int m, n, nn;
  float sum;

  if ( recenter ) {
    for ( m = 0; m < 3; m++ ) {
      sum = 0.0;

      for ( n = 0; n < nbas; n++ ) {
        nn = fnum * nbas + n;
        sum = sum + frame[m][nn];
      }

      center[m] = sum / nbas;
    }
  }

  for ( n = 0; n < nbas; n++ ) {
    nn = fnum * nbas + n;
    ball[n].pos[0] = frame[0][nn] - center[0];
    ball[n].pos[1] = frame[1][nn] - center[1];
    ball[n].pos[2] = frame[2][nn] - center[2];
  }
}


/* ----- putframe ------ */
void putframe( struct ballstr ball[], int fnum ) {
  int   n, nn;

  for ( n = 0; n < nbas; n++ ) {
    nn = fnum * nbas + n;
    frame[0][nn] = ball[n].pos[0];
    frame[1][nn] = ball[n].pos[1];
    frame[2][nn] = ball[n].pos[2];
  }
}


/* ----- prframes  ------ */
void prframes() {
  printf( "Number of frames: %d\n", nframe );
  printf( "frame %-5d <%s>\n", 1, frstr[0] );

  if ( nframe > 1 ) { printf( "frame %-5d <%s>\n", 2, frstr[1] ); }

  if ( nframe > 2 ) { printf( "frame %-5d <%s>\n", nframe, frstr[nframe - 1] ); }

  return;
}

/* ----- draw_axes ----- */
void draw_axes() {
  float e0[3], e1[3], e2[3], z0[2], z1[2], z2[2];
  float fac, r0, r1, r2, tx, ty, zz[3];
  int m, i, j, i0;

  tx  = ( 70 - midx ) / PSFAC;
  ty  = -( igh - 120  - midy ) / PSFAC;
  fac = 30;

  for ( m = 0; m < 3; m++ ) {
    e0[m] = tmat[m][0];
    e1[m] = tmat[m][1];
    e2[m] = tmat[m][2];
  }

  atompos( fac, e0, 1.0, z0, &r0 );

  atompos( fac, e1, 1.0, z1, &r1 );
  atompos( fac, e2, 1.0, z2, &r2 );

  /* sort vectors (clumsily) back to front */
  zz[0] = e0[2];
  zz[1] = e1[2];
  zz[2] = e2[2];

  for ( i = 0; i < 3; i++ ) {
    i0 = 0;

    if ( zz[1] < zz[i0] ) { i0 = 1; }

    if ( zz[2] < zz[i0] ) { i0 = 2; }

    if ( i0 == 0 ) { DrawArrow( tx, ty, z0[0] + tx, z0[1] + ty, r0, "x 100" ); }

    if ( i0 == 1 ) { DrawArrow( tx, ty, z1[0] + tx, z1[1] + ty, r1, "y 010" ); }

    if ( i0 == 2 ) { DrawArrow( tx, ty, z2[0] + tx, z2[1] + ty, r2, "z 001" ); }

    zz[i0] = 1e10;
  }

}

/* ----- bs_transform ------ */
/* multiplizieren der Atompositionen mit der Rotationsmatrix */

void bs_transform( int natom, struct ballstr ball[] ) {
  register int m, n;

  for ( n = 0; n < natom; n++ ) {
    for ( m = 0; m < 3; m++ )
      p[n][m] = tmat[m][0] * ball[n].pos[0]
                + tmat[m][1] * ball[n].pos[1]
                + tmat[m][2] * ball[n].pos[2];
  }
}

//--------------------------------------------------------------------------------

void bs_kernel( int natom, struct ballstr ball[], int nbond,
                struct stickstr stick[] ) {

  int   flag[NAMAX], ip[NAMAX], j, k, kk, m, n, ibot;
  float br, xx, bx, by, rk, rkk, th1, th2, cth1, cth2, sth1, sth2;
  float w, ww, bb, aa, crit1, crit2, fac, beta, gray;
  char label[81];
  int   ib, note;
  float zp[NAMAX][2], zr[NAMAX];
  float q1[3], q2[3], b[3], d[3], bot, big;
  float fudgefac = 0.6, bmidx, bmidy, dd;
  float m1[6], m2[6];
  int nbx, ibx, jbx, kbx[100], abx[100], pbx[100], fbx[100]; /* faster bond search */
  dist = dist0;

  if ( pmode == 0 || pmode == 1 ) { dist = 10000.0; }

  d[0] = 0.0;

  d[1] = 0.0;

  d[2] = dist;

  fac = scale;

  /* ------- sort atoms back to front ----- */
  for ( k = 0; k < natom; k++ ) { flag[k] = 0; }

  for ( n = 0; n < natom; n++ ) {
    bot = 1.0e10;
    ibot = 0;

    for ( k = 0; k < natom; k++ )
      if ( p[k][2] < bot && !flag[k] ) { bot = p[k][2]; ibot = k; }

    ip[n] = ibot;

    flag[ibot] = 1;
  }

  /* ------- make list of sphere centers and radii ---- */
  big = 1000000;

  if ( nbas == 0 ) { big = 0; }

  xbot = ybot = big;

  xtop = ytop = -big;

  for ( k = 0; k < natom; k++ ) {
    atompos( fac, p[k], ball[k].rad, zp[k], &zr[k] );
    zr[k] = radfac * zr[k];

    if ( zp[k][0] - zr[k] < xbot ) { xbot = zp[k][0] - zr[k]; }

    if ( zp[k][0] + zr[k] > xtop ) { xtop = zp[k][0] + zr[k]; }

    if ( zp[k][1] - zr[k] < ybot ) { ybot = zp[k][1] - zr[k]; }

    if ( zp[k][1] + zr[k] > ytop ) { ytop = zp[k][1] + zr[k]; }
  }

  /*  printf ("bounds x %.3f %.3f  y %.3f %.3f\n", xbot,xtop,ybot,ytop); */


  /* ------- start loop over atoms; plot ball first ----- */
  for ( n = 0; n < natom; n++ ) {
    k = ip[n];
    rk = ball[k].rad;

    if ( !ball[k].special ) {
      beta = exp( gslope * ( p[k][2] - gz0 ) * gslope ) ;

      if ( grayvalues ) {
        if ( gmode == G_RAMP ) {
          gray = beta * ball[k].gray + ( 1 - beta ) * GRAY0;
        }
        else {
          gray = ball[k].gray;
        }
      }
      else { gray = 1.0; }

      if ( hardcopy )
        svps.hardcopy_ball( gray, ball[k].r, ball[k].g, ball[k].b,
                            zp[k][0] + taux, zp[k][1] + tauy, zr[k] , midx, midy );
      else {
        DrawBall( gray, ball[k].col, zp[k][0] + taux, zp[k][1] + tauy, zr[k] );
      }

      if ( numbers || coords ) {
        if ( numbers == 1 ) { sprintf( label, "%d", k + 1 ); }

        if ( numbers == 2 ) { sprintf( label, "%s", ball[k].lab ); }

        if ( coords )
          sprintf( label, "(%.2f,%.2f,%.2f)",
                   ball[k].pos[0] + center[0], ball[k].pos[1] + center[1],
                   ball[k].pos[2] + center[2] );

        if ( hardcopy ) {
          svps.hardcopy_label( zp[k][0] + taux, zp[k][1] + tauy, label , midx, midy );
        }
        else {
          LabelBG( zp[k][0] + taux, zp[k][1] + tauy - 2, 0.0, 1.0, label );
        }
      }
    }

    /*  ------ make list of bonds to this atom ----- */
    if ( !withbonds ) { continue; }

    nbx = 0;

    for ( j = 0; j < nbond; j++ ) {
      if ( k == stick[j].start ) {
        kbx[nbx] = j;
        abx[nbx] = stick[j].end;
        nbx++;
      }
      else if ( k == stick[j].end ) {
        kbx[nbx] = j;
        abx[nbx] = stick[j].start;
        nbx++;
      }
    }

    if ( nbx == 0 ) { continue; }

    for ( m = 0; m < nbx; m++ ) { fbx[m] = 0; }   /* sort mini-list */

    for ( m = 0; m < nbx; m++ ) {
      bot = 1.0e10;
      ibot = 0;

      for ( j = 0; j < nbx; j++ )
        if ( p[abx[j]][2] < bot && !fbx[j] ) { bot = p[abx[j]][2] ; ibot = j; }

      pbx[m] = ibot;

      fbx[ibot] = 1;
    }

    /*  ------ inner loop over bonds ----- */
    for ( ibx = 0; ibx < nbx; ibx++ ) {
      jbx = pbx[ibx];
      kk = abx[jbx];
      ib = kbx[jbx];

      if ( ib < 0 ) { printf( "this cannot happen\n" ); }

      rkk = ball[kk].rad;

      /*  next few lines: the old direct procedure  */
      /*  for (nn=0;nn<natom;nn++) {
            kk = ip[nn];
            rkk = ball[kk].rad;
            ib=-1;
            for (j=0;j<nbond;j++) {
              if(k==stick[j].start && kk==stick[j].end) ib=j;
              if(kk==stick[j].start && k==stick[j].end) ib=j;
            }   */

      if ( ib >= 0 ) {
        br = bndfac * stick[ib].rad;
        bx = zp[kk][0] - zp[k][0];
        by = zp[kk][1] - zp[k][1];
        xx = sqrt( bx * bx + by * by );

        if ( xx * xx < 0.0001 ) { continue; }

        bx = bx / xx;

        by = by / xx;

        svm.vsum( d, p[k],  1.0, -1.0, q1 );
        svm.vsum( d, p[kk], 1.0, -1.0, q2 );
        svm.vsum( p[kk], p[k], 1.0, -1.0, b );

        cth1 =  svm.sp( q1, b ) / sqrt( svm.sp( q1, q1 ) * svm.sp( b, b ) );
        cth2 = -svm.sp( q2, b ) / sqrt( svm.sp( q2, q2 ) * svm.sp( b, b ) );

        th1 = acos( cth1 );
        th2 = acos( cth2 );

        crit1 = asin( br / rk ) * fudgefac;

        if ( crit1 < 0.0 ) { crit1 = 0.0; }

        crit2 = asin( br / rkk ) * fudgefac;

        if ( crit2 < 0.0 ) { crit2 = 0.0; }

        note = 0;

        if ( th2 - 0.5 * PI > crit2 && k < kk ) { note = 1; }

        if ( th1 - 0.5 * PI < crit1 && k > kk ) { note = 2; }

        /*        if(th2-0.5*PI>crit2 && n<nn) note=1;
                if(th1-0.5*PI<crit1 && n>nn) note=2; */

        /* ------- plot a stick ------ */
        if ( note == 1 || note == 2 ) {
          w = sqrt( rk * rk - br * br );
          sth1 = sqrt( 1.0 - cth1 * cth1 );
          ww = w * sth1 * zr[k] / rk;
          bb = br * zr[k] / rk;
          aa = br * cth1 * zr[k] / rk;
          m1[0] = bx * aa;
          m1[1] = by * aa;
          m1[2] = -by * bb;
          m1[3] = bx * bb;
          m1[4] = zp[k][0] + bx * ww + taux;
          m1[5] = zp[k][1] + by * ww + tauy;
          w = sqrt( rkk * rkk - br * br );
          sth2 = sqrt( 1.0 - cth2 * cth2 );
          ww = w * sth2 * zr[kk] / rkk;
          bb = br * zr[kk] / rkk;
          aa = br * cth2 * zr[kk] / rkk;
          m2[0] = bx * aa;
          m2[1] = by * aa;
          m2[2] = -by * bb;
          m2[3] = bx * bb;
          m2[4] = zp[kk][0] - bx * ww + taux;
          m2[5] = zp[kk][1] - by * ww + tauy;

          beta = exp( gslope * ( 0.5 * ( p[k][2] + p[kk][2] ) - gz0 ) * gslope );

          if ( grayvalues ) {
            if ( gmode == G_STD ) {
              gray = stick[ib].gray;
            }
            else if ( gmode == G_LIGHT ) {
              gray = 0.5 * ( ball[k].gray + ball[kk].gray );
            }
            else if ( gmode == G_RAMP ) {
              gray = beta * stick[ib].gray + ( 1 - beta ) * GRAY0;
            }
          }
          else { gray = 1.0; }

          if ( grayvalues && bline && gray > 0.7 ) { gray = 0.7; }

          if ( !grayvalues && bline ) { gray = 0.0; }

          if ( grayvalues && wire ) { gray = 0.0; }

          if ( bline ) { gray = 0.0; }    /* overrides... black if lines */

          if ( ball[k].special ) {
            dbond( 1.0, m2, m1 );
          }
          else if ( ball[kk].special ) {
            dbond( 1.0, m1, m2 );
          }
          else {
            if ( hardcopy ) {
              svps.hardcopy_stick( gray, m1, m2 , midx, midy );
            }
            else
              /*              printf ("draw stick %d  col0=%d\n", ib, stick[ib].col); */
            {
              DrawStick( gray, stick[ib].col, m1, m2 );
            }
          }

          /*  next part writes bond lengths onto the sticks */
          if ( bondnums ) {
            bmidx = 0.5 * ( zp[k][0] + zp[kk][0] ) + taux;
            bmidy = 0.5 * ( zp[k][1] + zp[kk][1] ) + tauy;
            dd = 0;

            for ( m = 0; m < 3; m++ ) {
              dd = dd + pow( ball[k].pos[m] - ball[kk].pos[m], 2 );
            }

            dd = sqrt( dd );

            sprintf( label, "%.2f", dd );

            if ( hardcopy ) {
              svps.hardcopy_label( bmidx, bmidy, label , midx, midy );
            }
            else {
              LabelBG( bmidx, bmidy, 0.0, 1.0, label );
            }
          }

        }
      } /* if (ib!=0) */
    }   /* end loop over nn */
  }     /* end loop over n */

  if ( showaxes ) { draw_axes(); }

}




/* ----- do_ConfigureNotify ---- */
int do_ConfigureNotify( XEvent * eventp ) {

  XConfigureEvent *e = ( XConfigureEvent * ) eventp;
  int  x1, x2, y1, y2, mmx, mmy, h1, w1, xx;

  /*  x1=midx+PSFAC*(taux+xbot);
    x2=midx+PSFAC*(taux+xtop);
    y1=midy+PSFAC*(tauy+ybot);
    y2=midy+PSFAC*(tauy+ytop);
    printf ("old area: X %d %d  Y %d %d\n", x1,x2,y1,y2); */

  mmx = midx;
  mmy = midy;
  w1 = igw;
  h1 = igh;

  if ( ( e->height != igh ) || ( e->width != igw ) ) {
    igw = e->width;
    igh = e->height;
    midx = igw / 2;
    midy = igh / 2 - 20;

    /* make new pixmaps but save contents of bgrmap */
    XCopyArea( dpy, bgrmap, pixmap, gc, 0, 0, w1, h1, 0, 0 );
    XFreePixmap( dpy, bgrmap );
    bgrmap = XCreatePixmap( dpy, win, igw, igh, depth );
    XFillRectangle( dpy, bgrmap, gcbg, 0, 0, igw, igh );
    XCopyArea( dpy, pixmap, bgrmap, gc, 0, 0, w1, h1, 0, 0 );
    XFreePixmap( dpy, pixmap );
    pixmap = XCreatePixmap( dpy, win, igw, igh, depth );

    if ( num_print == 0 ) { /* normally put back to middle on resize */
      taux = taux0 = 0;
      tauy = tauy0 = 0;
    }
    else {          /* but don't move plot about if building a print */
      taux = taux - ( midx - mmx ) / PSFAC;
      taux0 = taux0 - ( midx - mmx ) / PSFAC;
      tauy = tauy - ( midy - mmy - igh + h1 ) / PSFAC;
      tauy0 = tauy0 - ( midy - mmy - igh + h1 ) / PSFAC;
    }

    return 1;
  }
  else {
    return 0;
  }
}


//--------------------------------------------------------------------------------

int close_print( int helpme ) {

  if ( helpme ) {
    sprintf( gmsg, "Usage: close  - close current print file" );
    return 0;
  }

  if ( num_print == 0 ) {
    sprintf( emsg, "close: no print file is open" );
    return 0;
  }

  svps.hardcopy_close();

  if ( autops == 0 ) {
    XFillRectangle( dpy, bgrmap, gcbg, 0, 0, igw, igh );
  }

  num_print = 0;

  return 1; // 'true' ?
}


//--------------------------------------------------------------------------------

int handle_print( char msg1[], char msg2[], int helpme ) {

  char xx[81];
  float x, y, x1, y1, z1, x2, y2, z2;
  int ix, iy, dd;
  Drawable drw1;

  if ( helpme ) {
    sprintf( gmsg, "Usage: print [-T] [-t title] [file]  - print to PS file" );
    return 0;
  }

  sprintf( gmsg, "Print:" );

  if ( num_print > 0 && strcmp( prf, prfsave ) ) {
    close_print( 0 );
    sprintf( xx, " %s closed,", prfsave );
    strcat( gmsg, xx );
    replot = 1;
  }

  if ( num_print == 0 ) {
    svps.hardcopy_init( prf );
    sprintf( xx, " %s opened,", prf );
    strcat( gmsg, xx );
    pr_xoff = 0;
    pr_yoff = ybot - 15;
  }
  else {
    svps.hardcopy_redefine();
  }

  num_print++;

  strcpy( prfsave, prf );
  hardcopy = 1;
  bs_transform( natom, ball );
  bs_kernel( natom, ball, nbond, stick );
  hardcopy = 0;
  drw = bgrmap;

  if ( autops == 0 ) {
    bs_kernel( natom, ball, nbond, stick );
  }

  if ( strlen( msg1 ) > 0 ) {
    x = pr_xoff + taux;
    y = pr_yoff + tauy;
    svps.hardcopy_label( x, y,  msg1 , midx, midy );
    ix = midx + PSFAC * x;
    iy = midy - PSFAC * y;
    dd = strlen( msg1 ) * 6;

    if ( autops == 0 ) {
      XDrawString( dpy, win, labelgc, ix - dd / 2, iy, msg1, strlen( msg1 ) );
      XDrawString( dpy, bgrmap, labelgc, ix - dd / 2, iy, msg1, strlen( msg1 ) );
    }
  }

  if ( strlen( msg2 ) > 0 ) {
    x = pr_xoff + taux;
    y = pr_yoff + tauy - 10;
    svps.hardcopy_label( x, y,  msg2 , midx, midy );
    ix = midx + PSFAC * x;
    iy = midy - PSFAC * y;
    dd = strlen( msg2 ) * 6;

    if ( autops == 0 ) {
      XDrawString( dpy, win, labelgc, ix - dd / 2, iy, msg2, strlen( msg2 ) );
      XDrawString( dpy, bgrmap, labelgc, ix - dd / 2, iy, msg2, strlen( msg2 ) );
    }
  }

  svps.flush();

  sprintf( xx, " write(%d) to %s", num_print, prf );
  strcat( gmsg, xx );


  return 1; // 'true'
}



/* ----- update_from_file---- */
int update_from_file() {
  char msg1[81], msg2[81], str[81], pat[81], w[8][41];
  char *p;
  float sh[3][6], cut[3], cut1, cut2;
  int i, j, n, rc, nw, svstep, svrgb, helpme;

  nbas = 0;
  nbonds = 0;
  nspec = 0;
  nframe = 1;

  if ( !readclusterdata( inf ) ) {
    sprintf( emsg, "Cannot update from file %s", inf );
    return 0;
  }

  sprintf( curf, "%s ", inf );

  clearline( win, 10, 8 );
  showline( win, 10, 8, "Reading ", inmv, " .. please wait" );
  XFlush( dpy );
  sprintf( gmsg, "Updated from %s ", inf );

  if ( readclusterdata( inmv ) ) {
    sprintf( gmsg, "Updated from %s and %s", inf, inmv );
    sprintf( frstr[0], "%s", "start ..\0" );
    strcat( curf, inmv );
  }

  for ( i = 0; i < ncol; i++ ) { FreeColorGC( i ); }

  if ( autocolor ) { set_auto_colors(); }

  parse_all_colors();

  if ( color ) {
    SetColors();
  }
  else {
    if ( stippled ) { SetStippled4x4(); }
    else { SetSmoothGrays(); }
  }

  natom = ball_list( ball, 0 );

  nbond = stick_list( ball, stick );
  putframe( ball, 0 );
  getframe( ball, iframe );
  replot = 1;
  resetup = 1;
  return 1;
}

int interpret_input( char inp[] ) {

  char msg1[81], msg2[81], str[81], pat[81], w[8][41];
  char *p;
  float sh[3][6], cut[3], cut1, cut2;
  int j, n, rc, nw, svstep, svrgb, helpme;

  /* for axyz option */
  int atom_nr;
  float dx, dy, dz;
  char label[10];
  int atom1, atom2, atom3;

  nw = parse_args( inp, w );

  if ( nw == 0 ) {
    sprintf( gmsg, "Empty input" );
    return 0;
  }

  helpme = 0;

  if ( !strcmp( w[1], "?" ) ) { helpme = 1; }

  if ( !strcmp( w[1], "-h" ) ) { helpme = 1; }

  if ( abbrev( w[0], "help", 4 ) ) {
    if ( nw == 1 ) {
      sprintf( gmsg, "Try \'command ?\' or press key h" );
      return 0;
    }
    else {
      helpme = 1;
      sprintf( inp, "%s", w[1] );
      nw = parse_args( inp, w );
    }
  }

  if ( abbrev( w[0], "update", 2 ) ) {
    if ( helpme )
      sprintf( gmsg, "Usage: update [-color] [-rv] [+rv] [-bw] [-st] [-auto] "
               "[file]  - update from file" );
    else {
      unsigned int i = 1;

      while ( i < nw ) {
        if ( abbrev( w[i], "-st", 3 ) )    { stippled = 1; color = 0; }
        else if ( abbrev( w[i], "-bw", 3 ) )    { stippled = 0; color = 0; }
        else if ( abbrev( w[i], "-color", 4 ) ) { color = 1; }
        else if ( abbrev( w[i], "-auto", 5 ) )  { autocolor = 1; color = 1; }
        else if ( abbrev( w[i], "-rv", 3 ) ) { reverse = 1; }
        else if ( abbrev( w[i], "+rv", 3 ) ) { reverse = 0; }
        else if ( w[i][0] == '-' ) {
          sprintf( emsg, "update: unknown flag %s", w[i] );
          return 0;
        }
        else {
          strext( inf, w[i], "bs", 0 );
          strext( inmv, inf, "mv", 1 );
        }

        i++;
      }

      update_from_file();
    }
  }

  else if ( abbrev( w[0], "save", 3 ) ) {
    if ( helpme ) {
      sprintf( gmsg, "Usage: save [-rgb] [-step n] file  - save data" );
    }
    else {
      unsigned int i = 1;
      svstep = 1;
      svrgb = 0;

      while ( i < nw ) {
        if ( abbrev( w[i], "-step", 3 ) )   {i++; svstep = atoi( w[i] ); }
        else if ( abbrev( w[i], "-rgb", 4 ) ) { svrgb = 1; }
        else if ( w[i][0] == '-' ) {
          sprintf( emsg, "save: unknown flag %s", w[i] );
          return 0;
        }
        else { strext( outf, w[i], "bs", 0 ); }

        i++;
      }

      showline( win, 10, 8, "Saving .. ", "" , "" );

      XFlush( dpy );
      writeclusterdata( outf, svstep, svrgb );
    }
  }

  else if ( abbrev( w[0], "print", 2 ) ) {
    if ( helpme ) {
      handle_print( msg1, msg2, helpme );
    }
    else {
      unsigned int i = 1;
      strcpy( msg1, "" );
      strcpy( msg2, "" );

      while ( i < nw ) {
        if ( !strcmp( w[i], "-t" ) ) {
          i++;
          strcpy( msg1, w[i] );
        }
        else if ( !strcmp( w[i], "-T" ) ) {
          sprintf( msg1, "%s frame %d of %d", inf, iframe + 1, nframe );

          if ( strlen( frstr[iframe] ) > 0 ) { strip( msg2, frstr[iframe] ); }
        }
        else if ( w[i][0] == '-' ) {
          sprintf( emsg, "print: unknown flag %s", w[i] );
          return 0;
        }
        else {
          strext( prf, w[i], "ps", 0 );
        }

        i++;
      }

      handle_print( msg1, msg2, helpme );
    }
  }

  else if ( abbrev( w[0], "close", 2 ) ) {
    if ( !helpme )
      sprintf( gmsg, "Print: %s closed after %d write%s",
               prf, num_print, num_print == 1 ? "" : "s" );

    close_print( helpme );

    replot = 1;
  }

  else if ( abbrev( w[0], "dup", 3 ) ) {
    for ( unsigned int i = 0; i < 3; i++ ) for ( j = 0; j < 6; j++ ) { sh[i][j] = 0; }

    sscanf( inp, "%*s %f %f %f %f %f %f %f %f %f %f %f %f",
            &sh[0][0], &sh[1][0], &sh[2][0], &sh[0][1], &sh[1][1], &sh[2][1],
            &sh[0][2], &sh[1][2], &sh[2][2], &sh[0][3], &sh[1][3], &sh[2][3],
            &sh[0][4], &sh[1][4], &sh[2][4], &sh[0][5], &sh[1][5], &sh[2][5] );

    if ( duplicate_atoms( sh, helpme ) ) { resetup = 2; }
  }

  /* add one atom at an dx,dy,dz displacement from an atom */

  else if ( abbrev( w[0], "axyz", 4 ) ) {
    sscanf( inp, "%*s %d %s %f %f %f", &atom_nr, &label[0], &dx, &dy, &dz );

    if ( add_atom_xyz( atom_nr - 1, label, dx, dy, dz, helpme ) ) { resetup = 2; }
  }

  /* add with distance , angele, dihedral */
  else if ( abbrev( w[0], "azmat", 4 ) ) {
    sscanf( inp, "%*s %s %d %f %d %f %d %f", &label[0], &atom1, &dx, &atom2, &dy, &atom3, &dz );

    if ( add_atom_zmat( label, atom1 - 1, dx, atom2 - 1, dy, atom3 - 1, dz, helpme ) ) { resetup = 2; }
  }


  else if ( abbrev( w[0], "cut", 3 ) ) {
    cut[0] = cut[1] = cut[2] = cut1 = cut2 = 0;
    sscanf( inp, "%*s %f %f %f %f %f",
            &cut[0], &cut[1], &cut[2], &cut1, &cut2 );

    if ( cut_atoms( cut, cut1, cut2, helpme ) ) { resetup = 2; }
  }

  else if ( abbrev( w[0], "color", 3 ) ) {
    strcpy( pat, "" );
    sscanf( inp, "%*s %s %n", pat, &n );

    for ( std::size_t i = n; i < strlen( inp ) + 1; i++ ) { str[i - n] = inp[i]; }

    if ( NewSpecColor( pat, str, helpme ) ) { resetup = 1; }
  }

  else {
    rc = readclusterline( inp, helpme );

    if ( rc == 1 ) { replot = 1; }

    if ( rc == 2 ) { resetup = 1; }
  }


  return 0; // 2015-11-27 correct ?
}

//--------------------------------------------------------------------------------
// forward declare
int interpret_keypress( XEvent *ev, int *inpmode, char input[],
                        int *ixyz, float *alfa );

//--------------------------------------------------------------------------------

int main( int argc, char * argv[] ) {

#define bs_icon_width 50
#define bs_icon_height 50

  /* program icon */

  static char bs_icon_bits[] = {
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x7c,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x83, ( char )0x01,
    ( char )0x00, ( char )0x00, ( char )0xf8, ( char )0x03,
    ( char )0x80, ( char )0x00, ( char )0x02, ( char )0x00,
    ( char )0x00, ( char )0xac, ( char )0x06, ( char )0x40,
    ( char )0x00, ( char )0x04, ( char )0x00, ( char )0x00,
    ( char )0x57, ( char )0x1d, ( char )0x40, ( char )0x00,
    ( char )0x04, ( char )0x00, ( char )0x00, ( char )0xab,
    ( char )0x1a, ( char )0x20, ( char )0x00, ( char )0x08,
    ( char )0x00, ( char )0x80, ( char )0x55, ( char )0x35,
    ( char )0x20, ( char )0x00, ( char )0x08, ( char )0x00,
    ( char )0xc0, ( char )0xaa, ( char )0x6a, ( char )0x20,
    ( char )0x06, ( char )0x08, ( char )0x00, ( char )0x40,
    ( char )0x55, ( char )0x55, ( char )0xfe, ( char )0x01,
    ( char )0x08, ( char )0x00, ( char )0xc0, ( char )0xaa,
    ( char )0xea, ( char )0x21, ( char )0x00, ( char )0x08,
    ( char )0x00, ( char )0x40, ( char )0x55, ( char )0x55,
    ( char )0x40, ( char )0x00, ( char )0x04, ( char )0x00,
    ( char )0xc0, ( char )0xaa, ( char )0x6a, ( char )0x40,
    ( char )0x00, ( char )0x04, ( char )0x00, ( char )0x40,
    ( char )0xd5, ( char )0x55, ( char )0x80, ( char )0x00,
    ( char )0x02, ( char )0x00, ( char )0xc0, ( char )0xaa,
    ( char )0x6a, ( char )0x00, ( char )0x83, ( char )0x01,
    ( char )0x00, ( char )0x80, ( char )0xd5, ( char )0x35,
    ( char )0x00, ( char )0x7c, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0xab, ( char )0x1b, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x57, ( char )0x1d, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0xac,
    ( char )0x07, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0xf8, ( char )0x03,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x02, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x02, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x04, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0xfe,
    ( char )0x03, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x01, ( char )0x04,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0xc0, ( char )0x00, ( char )0x18, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x40,
    ( char )0x00, ( char )0x10, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x20, ( char )0x00,
    ( char )0x20, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x10, ( char )0x00, ( char )0x40,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x10, ( char )0x00, ( char )0x40, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x10,
    ( char )0x00, ( char )0x40, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x10, ( char )0x00,
    ( char )0x40, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x10, ( char )0x00, ( char )0x40,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x10, ( char )0x00, ( char )0x40, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x10,
    ( char )0x00, ( char )0x40, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x10, ( char )0x00,
    ( char )0x40, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x10, ( char )0x00, ( char )0x40,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x20, ( char )0x00, ( char )0x20, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x40,
    ( char )0x00, ( char )0x10, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0xc0, ( char )0x00,
    ( char )0x18, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x01, ( char )0x04,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0xfe, ( char )0x03, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00, ( char )0x00, ( char )0x00,
    ( char )0x00, ( char )0x00
  };

  Pixmap        iconmap;

  XEvent        ev, ev1;

  Font          font;
  XSizeHints    hint;
  char          input[257], msg[81];
  int           keytype, inpmode, finished;
  char         *p;

  int   i, j, k, ixyz;
  float alfa, beta, gama, cut[3], cut1, cut2;
  char  str[50];
  char yn[2][5] = { "no", "yes"};

  float phi;

  // ----- set parameters and defaults -----

  /*  bm parameters */
  autops = 0;   /* flag for automatic generation of .ps file from .bs file */

  dist = dist0 = 15.0;
  scale = 15;
  animate = 0;
  msecs = 200;
  framephase = 1;
  igs = 1.0;
  taux0 = taux = tauy0 = tauy = 0;
  dtaux = 10;
  dtauy = 10;
  wire = 0;
  bline = 0;
  radfac = 1.0;
  bndfac = 1.0;
  alfa = 0;
  beta = 0;
  gama = 0;

  // setup euler angles
  svm.eumat( alfa * PI / 180.0, beta * PI / 180.0, gama * PI / 180.0, tmat );

  /* dalpha = step by which molecule is translated with arror keys
              set to 5 for a 'smooth' rotation
              dalfa=90.0;
  */
  dalfa = 5.0;

  nspec = 0;
  nbas = 0;
  nbonds = 0;
  alat = 1.0;
  amp = 0.0;
  hardcopy = 0;
  usepixmap = 0;
  numbers = 0;
  pmode = 0;
  recenter = 1;
  wrhelp = 0;
  shadow = 0;
  grayvalues = 1;
  withbonds = 1;
  bondnums = 0;
  color = 1;
  autocolor = 0;
  reverse = 0;
  fstep = 1;
  wrinfo = 0;
  stippled = 0;
  showaxes = 0;
  saveframe = -99;
  coords = 0;
  gmode = G_STD;
  center[0] = center[1] = center[2] = 0;
  startup = 1;

  hint.x = 0;
  hint.y = 0;
  hint.width =  W_WIDTH;
  hint.height = W_HEIGHT;

  // parse arguments

  // command line arguments

  i = 1;
  j = 0;
  k = -1;
  igeo = 0;

  while ( i < argc ) {
    if ( !strcmp( argv[i], "-geo" ) ) {
      i++;

      //+2015-11-26
      int ui_igx = igx;
      int ui_igy = igy;
      unsigned int ui_igw = igw;
      unsigned int ui_igh = igh;


      igeo = XParseGeometry( argv[i], &ui_igx, &ui_igy, &ui_igw, &ui_igh );

      igx = ui_igx;
      igy = ui_igy;
      igw = ui_igw;
      igh = ui_igh;
    }
    else if ( abbrev( argv[i], "-sc", 3 ) )        {i++; igs = atof( argv[i] ); }
    else if ( abbrev( argv[i], "-t", 2 ) )         {i++; k = i; }
    else if ( abbrev( argv[i], "-st", 3 ) )        { stippled = 1; color = 0; }
    else if ( abbrev( argv[i], "-bw", 3 ) )        { stippled = 0; color = 0; }
    else if ( abbrev( argv[i], "-color", 4 ) ) { color = 1; }
    else if ( abbrev( argv[i], "-autocolor", 5 ) ) { autocolor = 1; color = 1; }
    else if ( abbrev( argv[i], "-rv", 3 ) ) { reverse = 1; }
    else if ( abbrev( argv[i], "-hh", 3 ) )        { WriteHelp(); exit( 0 ); }
    else if ( abbrev( argv[i], "-h", 2 ) ) { j = 1; }

    /* bm automatic postscript option : read spec, print, exit */

    else if ( abbrev( argv[i], "-ps", 3 ) ) {

      autops = 1;

    }

    else if ( argv[i][0] == '-' ) {
      printf( "Unknown flag %s\n", argv[i] ) ;
      exit( 1 );
    }
    else { strext( inf, argv[i], "bs", 0 ); }

    i++;
  }

  p = inf;

  for ( std::size_t ii = 0; ii < strlen( inf ); ii++ )  {
    if ( inf[ii] == '/' ) { p = inf + ii + 1; }
  }

  strcpy( wname, p );

  if ( k > -1 ) { strcpy( wname, argv[k] ); }

  if ( igeo & WidthValue ) { hint.width = igw; }
  else { igw = hint.width; }

  if ( igeo & HeightValue ) { hint.height = igh; }
  else { igh = hint.height; }

  if ( igeo & XValue ) { hint.x = igx; }

  if ( igeo & YValue ) { hint.y = igy; }

  /* need next line to make +0+0 work */
  if ( ( igeo & XValue ) && ( hint.x == 0 ) && ( hint.y == 0 ) ) { hint.x = 1; }

  if ( j == 1 ) {
    printf( "\nUsage:  xbs [flags] id   -- ball-and-sticks plotting program\n"
            "Data is read from files id.bs and id.mv\n"
            "\nFlags:  -geo gg     set window geometry   \n"
            "        -sc x       set scale factor\n"
            "        -t title    set window title\n"
            "        -color      use color\n"
            "        -bw         b/w with smooth grays\n"
            "        -st         b/w with stippled grays\n"
            "        -rv         reverse colors\n"
            "        -autocolor  chose own colors\n"
            "        -ps         non-interactively produce .ps from .bs file\n"
            "        -hh         long help\n"
            "\nHelp: Enter 'xbs -hh' to get an overview.\n"
            "For on-line help, press key 'h' for the overview or use\n"
            "'help <cmd>' or '<cmd> ?' in the input line for information\n"
            "on a specific command (including possible options).\n"
            "Request the input line with key 'i'.\n"
          );
    printf( "\nSettings: geometry %dx%d%+d%+d, scale %.2f,\n"
            "color %s, autocolor %s, stippled %s, reverse %s\n ",
            hint.width, hint.height, hint.x, hint.y, igs,
            yn[color], yn[autocolor], yn[stippled], yn[reverse]

          );
    exit( 0 );
  }

  if ( !strstr( inf, "." ) ) { strcat( inf, ".bs" ); }

  /* ----- read data file ------ */
  sprintf( curf, "%s ", inf );

  nframe = 1;

  strext( inmv, inf, "mv", 1 );

  if ( ! readclusterdata( inf ) ) {
    printf( "Could not open input file %s\n", inf );
    exit( 1 );
  }

  if ( natom <= 80 ) { usepixmap = 1; }

  /* ----- setup ---- */

  if ( autops == 0 ) {

    dpy = XOpenDisplay( "" );

    if ( ! dpy ) {
      printf( "Error: could not open display\n" );
      exit( 1 );
    }

    screen = DefaultScreen( dpy );

    screenptr = DefaultScreenOfDisplay( dpy );
    bground = WhitePixel( dpy, screen );
    fground = BlackPixel( dpy, screen );

    if ( ( igeo & XValue ) && ( igeo & XNegative ) ) {
      hint.x = igx + WidthOfScreen( screenptr ) - hint.width - 18;
    }

    if ( ( igeo & YValue ) && ( igeo & YNegative ) ) {
      hint.y = igy + HeightOfScreen( screenptr ) - hint.height - 34;
    }

    hint.flags = PPosition | PSize;

    win = XCreateSimpleWindow( dpy, DefaultRootWindow( dpy ),
                               hint.x, hint.y, hint.width, hint.height,
                               0, 0, bground );

    midx = igw / 2;

    midy = igh / 2 - 20;

    iconmap = XCreateBitmapFromData( dpy, win,
                                     bs_icon_bits, bs_icon_width, bs_icon_height );

    XSetStandardProperties( dpy, win, wname, wname, iconmap,
                            argv, argc, &hint );

    gc = XCreateGC( dpy, win, 0, 0 );

    XSetBackground( dpy, gc, bground );

    XSetForeground( dpy, gc, fground );

    XSelectInput( dpy, win, ( KeyPressMask | ExposureMask | StructureNotifyMask ) );

    XMapRaised( dpy, win );

    drw = win;

    /* ----- additional setup for graphics  ---- */
    font = XLoadFont( dpy, FONT );

    XSetFont( dpy, gc, font );

    gcbg = XCreateGC( dpy, win, 0, 0 );

    XSetForeground( dpy, gcbg, bground );

    shadowgc = XCreateGC( dpy, win, 0, 0 );

    XSetLineAttributes( dpy, shadowgc, SHADOW,
                        LineSolid, CapRound, JoinBevel );

    XSetForeground( dpy, shadowgc, bground );

    labelgc = XCreateGC( dpy, win, 0, 0 );

    XSetForeground( dpy, labelgc, fground );

    font = XLoadFont( dpy, LABFONT );

    XSetFont( dpy, labelgc, font );

    labbggc = XCreateGC( dpy, win, 0, 0 );

    XSetLineAttributes( dpy, labbggc, 13,
                        LineSolid, CapRound, JoinMiter );

    XSetForeground( dpy, labbggc, bground );

    cmap    = XDefaultColormap( dpy, screen );

  } /* end od if autops == 0 */

  for ( i = 0; i < NPOINTS; i++ ) {
    phi = i * PI / ( NPOINTS - 1.0 );
    arc[i][0] = -sin( phi );
    arc[i][1] =  cos( phi );
  }

  /* ----- allocate colors, set atoms and bonds ------ */
  if ( autocolor ) { set_auto_colors(); }

  parse_all_colors();

  if ( autops == 0 ) {

    if ( color ) {
      SetColors();
    }
    else {
      if ( stippled ) { SetStippled4x4(); }
      else { SetSmoothGrays(); }
    }

  }

  natom = ball_list( ball, 10 );

  nbond = stick_list( ball, stick );
  putframe( ball, 0 );
  getframe( ball, 0 );
  sprintf( frstr[0], "%s", " \0" );



  /* ----- create pixmaps for buffering and background ---- */

  if ( autops == 0 ) {

    depth = XDefaultDepth( dpy, screen );
    pixmap = XCreatePixmap( dpy, win, igw, igh, depth );
    bgrmap = XCreatePixmap( dpy, win, igw, igh, depth );
    XFillRectangle( dpy, bgrmap, gcbg, 0, 0, igw, igh );

  }

  /* ----- start of main loop ---- */
  inpmode = 0;

  finished = 0;

  num_print = 0;

  count = 0;

  iframe = 0;

  keytype = K_NOP;


  /* (bm) hack for printing .ps file non-interactively */

  if ( autops == 1 ) {


    /* remove ending from filename */
    strcpy( pmsg1, wname );
    lpos = strrchr( pmsg1, '.' );
    *lpos = '\0';
    strcat( pmsg1, ".ps\0" );

    /* name of postscript file */
    strcpy( prf, pmsg1 );

    strcpy( pmsg1, "" );  /* title */
    strcpy( pmsg2, "" );

    scale =  60;    /* set scale */
    numbers = 2;   /* print atom numbers (mode 0..2) 2: print atom symbols */

    /* print and close file */
    handle_print( pmsg1, pmsg2, 0 );
    close_print( 0 );

    finished = 1; /* .. and exit */

  }

  /* (bm) end */


  strcpy( prfsave, prf );

  while ( finished == 0 ) {
    /* jkl mess: the single line  XNextEvent (dpy, &ev); was replaced with this */
    if ( ( animate == 0 ) || ( nframe == 1 ) ) { /* jkl */
      XNextEvent( dpy, &ev );
    }
    else {
      while ( animate > 0 ) {
        if ( XPending( dpy ) > 0 ) {
          XNextEvent( dpy, &ev );

          if ( ( ev.type == ConfigureNotify ) ||
               ( ev.type == KeyPress ) ) {
            animate = 0;
            goto ProcessEventNow;
          }
        }

        iframe += fstep * framephase;

        if ( animate == 1 ) {  /* frames shown as: 1, 2, ... N, 1, 2 .. */
          /*  or as N, N-1, N-2, ... 2, 1, N, N-1, N-2 */
          if ( ( iframe >= nframe ) && ( framephase > 0 ) ) {
            iframe = 0;
          }
          else if ( ( iframe < 0 ) && ( framephase < 0 ) ) {
            iframe = nframe - 1;
          }
        }
        else {  /* frames shown as: 1, 2, ... N, N-1, N-2, ... 2, 1, 2, 3... */
          if ( ( iframe >= nframe ) || ( iframe < 0 ) ) {
            framephase *= -1;
            iframe += 2 * framephase * fstep;
          }
        }

        boost::this_thread::sleep( boost::posix_time::milliseconds( msecs ) );
        //ms_sleep( msecs );

        resetup = 1;
        /* [][] */
        natom = ball_list( ball, 0 );
        getframe( ball, iframe );
        nbond = stick_list( ball, stick );
        replot = 1;

        /* ----- redraw the plot -------- */
        svm.rotmat( ixyz, alfa * PI / 180.0, tmat );
        bs_transform( natom, ball );

        if ( usepixmap ) {
          drw = pixmap;
          XCopyArea( dpy, bgrmap, pixmap, gc, 0, 0, igw, igh, 0, 0 );
        }
        else {
          drw = win;

          if ( wrinfo ) { WriteInfo( drw ); }

          if ( wrhelp > 0 ) { WriteHelp(); }

          XCopyArea( dpy, bgrmap, win, gc, 0, 0, igw, igh, 0, 0 );
        }

        showline( win, 10, 8,  "Busy                    ", " ", " " );

        bs_kernel( natom, ball, nbond, stick );

        if ( startup ) { showline( win, 10, 8, "Reading ", inmv, " .." ); }

        WriteStatus( drw );

        if ( wrinfo ) { WriteInfo( drw ); }

        if ( wrhelp > 0 ) { WriteHelp(); }

        if ( usepixmap ) {
          XCopyArea( dpy, pixmap, win, gc, 0, 0, igw, igh, 0, 0 );
        }

        if ( startup ) {
          showline( win, 10, 8, "Reading ", inmv, " .. please wait" );
          XFlush( dpy );
          k = nframe;

          if ( readclusterdata( inmv ) ) {
            sprintf( frstr[0], "%s", "start ..\0" );
            sprintf( gmsg, "%d frames were read from %s", nframe - k, inmv );
            strcat( curf, inmv );
            WriteStatus( win );
          }
          else { sprintf( gmsg, "No frames in %s", inmv ); }

          startup = 0;
        }

        if ( ! usepixmap ) {
          XCopyArea( dpy, win, pixmap, gc, 0, 0, igw, igh, 0, 0 );
        }

        showline( win, 10, 8, "Done                          ", "", " " );

        /* Discard extra keypress events so that rotation stops after
           key is released. Leave one for smooth motion. */
        i = 0;

        while ( XCheckTypedEvent( dpy, KeyPress, &ev1 ) ) {i = 1;}

        if ( i == 1 ) { XPutBackEvent( dpy, &ev1 ); }

        /* [][] */
      }
    }

    /* end of jkl mess for a while */

ProcessEventNow:
    ixyz = 0;

    replot = resetup = chginfo = 0;

    strcpy( gmsg, "" );

    strcpy( emsg, "" );

    switch ( ev.type ) {

      case Expose:

        if ( ev.xexpose.count == 0 ) {
          if ( !startup ) {
            XCopyArea( dpy, pixmap, win, gc, 0, 0, igw, igh, 0, 0 );

            if ( inpmode ) {
              clearline( win, 10, 8 );
              showline( win, 10, 8, "Input: ", input, "_ " );
            }

            if ( wrhelp > 0 ) { WriteHelp(); }
          }
          else {
            replot = 1;
          }
        }

        break;

      case ConfigureNotify:

        if ( do_ConfigureNotify( &ev ) ) { replot = 1; }

        if ( startup ) { replot = 0; }

        break;

      /* ----- process keypress events ----- */

      case KeyPress:
        strcpy( gmsg, "" );

        strcpy( emsg, "" );

        keytype = interpret_keypress( &ev, &inpmode, input, &ixyz, &alfa );

        switch ( keytype ) {

          case K_QUIT:
            finished = 1;
            break;

          case K_REPLOT:
            replot = 1;
            break;

          case K_RESETUP:
            resetup = 1;
            break;

          case K_READ_MORE:
            clearline( win, 10, 8 );
            showline( win, 10, 8, "Input: ", input, "_           " );
            break;

          case K_READ_DONE:
            clearline( win, 10, 8 );

            if ( strlen( input ) > 0 ) {
              for ( k = SVINP - 1; k > 1; k-- ) { strcpy( svinput[k], svinput[k - 1] ); }

              strcpy( svinput[1], input );

              nsvline++;
            }

            interpret_input( input );

            break;

          case K_UPDATE:
            update_from_file();
            break;
        }

        break;
    }  /* switch (ev.type) */

    /* ----- repeat the setup steps ------- */
    if ( resetup ) {
      natom = ball_list( ball, 0 );
      getframe( ball, iframe );
      nbond = stick_list( ball, stick );
      replot = 1;
    }

    /* ----- redraw the plot -------- */
    if ( replot ) {
      svm.rotmat( ixyz, alfa * PI / 180.0, tmat );
      bs_transform( natom, ball );

      if ( usepixmap ) {
        drw = pixmap;
        XCopyArea( dpy, bgrmap, pixmap, gc, 0, 0, igw, igh, 0, 0 );
      }
      else {
        drw = win;

        if ( wrinfo ) { WriteInfo( drw ); }

        if ( wrhelp > 0 ) { WriteHelp(); }

        XCopyArea( dpy, bgrmap, win, gc, 0, 0, igw, igh, 0, 0 );
      }

      showline( win, 10, 8,  "Busy                    ", " ", " " );

      bs_kernel( natom, ball, nbond, stick );

      if ( startup ) { showline( win, 10, 8, "Reading ", inmv, " .." ); }

      WriteStatus( drw );

      if ( wrinfo ) { WriteInfo( drw ); }

      if ( wrhelp > 0 ) { WriteHelp(); }

      if ( usepixmap ) {
        XCopyArea( dpy, pixmap, win, gc, 0, 0, igw, igh, 0, 0 );
      }

      if ( startup ) {
        showline( win, 10, 8, "Reading ", inmv, " .. please wait" );
        XFlush( dpy );
        k = nframe;

        if ( readclusterdata( inmv ) ) {
          sprintf( frstr[0], "%s", "start ..\0" );
          sprintf( gmsg, "%d frames were read from %s", nframe - k, inmv );
          strcat( curf, inmv );
          WriteStatus( win );
        }
        else { sprintf( gmsg, "No frames in %s", inmv ); }

        startup = 0;
      }

      if ( ! usepixmap ) {
        XCopyArea( dpy, win, pixmap, gc, 0, 0, igw, igh, 0, 0 );
      }

      showline( win, 10, 8, "Done                          ", "", " " );

      /* Discard extra keypress events so that rotation stops after
         key is released. Leave one for smooth motion. */
      i = 0;

      while ( XCheckTypedEvent( dpy, KeyPress, &ev1 ) ) {i = 1;}

      if ( i == 1 ) { XPutBackEvent( dpy, &ev1 ); }
    }

    /* ----- handle messages, update info if needed -------- */
    if ( strlen( emsg ) > 0 ) {
      clearline( win, 10, 8 );
      showline( win, 10, 8, "+++ ", emsg, "" );
      showline( pixmap, 10, 8, "+++ ", emsg, "" );
      XBell( dpy, BELL_LEVEL );
    }
    else if ( strlen( gmsg ) > 0 ) {
      clearline( win, 10, 8 );
      showline( win, 10, 8, gmsg, "", "" );
      showline( pixmap, 10, 8, gmsg, "", "" );
    }

    if ( ( !replot ) && wrinfo && chginfo ) { WriteInfo( win ); }

  }    /* while */


  /* ----- clean up and exit -------- */
  if ( num_print > 0 ) { svps.hardcopy_close(); }

  if ( autops == 0 ) {

    XFreeGC( dpy, gc );
    XDestroyWindow( dpy, win );
    XCloseDisplay( dpy );

  }
}

/* ------- intepret_keypress subroutine ------ */
int interpret_keypress( XEvent *ev, int *inpmode, char input[],
                        int *ixyz, float *alfa ) {
  int count, l;
  char buff[8], msg[81];
  unsigned int state;
  KeySym key;

  state = ev->xkey.state;
  count = XLookupString( &ev->xkey, buff, 8, &key, 0 );

  if ( state != 0 ) { /* jkl some translations for keypad keys */
    if ( key == XK_Up )           {key = XK_KP_8;}
    else if ( key == XK_Down )    {key = XK_KP_2;}
    else if ( key == XK_Right )   {key = XK_KP_6;}
    else if ( key == XK_Left )    {key = XK_KP_4;}
  }

  if ( key == XK_R7 )      {key = XK_KP_7;}

  if ( key == XK_Home )    {key = XK_KP_7;}

  if ( key == XK_R6 )      {key = XK_KP_Multiply;}

  if ( key == XK_KP_Subtract )   {key = XK_minus;}

  if ( key == XK_KP_Add )  {key = XK_plus;}

  if ( key == XK_F26 )     {key = XK_KP_Multiply;}

  if ( key == XK_F28 )     {key = XK_KP_8;}

  if ( key == XK_F27 )     {key = XK_KP_7;}

  if ( key == XK_F30 )     {key = XK_KP_4;}

  if ( key == XK_F34 )     {key = XK_KP_2;}

  if ( key == XK_F32 )     {key = XK_KP_6;}

  if ( *inpmode ) {
    l = strlen( input );

    if ( key == XK_Up ) {
      if ( svline == 0 ) { strcpy( svinput[0], input ); }

      if ( svline < nsvline - 1 ) { svline++; }

      strcpy( input, svinput[svline] );

      return K_READ_MORE;
    }

    if ( key == XK_Down ) {
      if ( svline > 0 ) { svline--; }

      strcpy( input, svinput[svline] );

      return K_READ_MORE;
    }

    if ( key == XK_Down ) {strcpy( input, "" ); return K_READ_MORE; }

    if ( key == XK_Return ) { *inpmode = 0; return K_READ_DONE; }

    if ( key == XK_BackSpace || key == XK_Left ) {
      if ( l > 0 ) { input[l - 1] = '\0'; }

      return K_READ_MORE;
    }

    if ( count > 0 ) {input[l] = buff[0]; input[l + 1] = '\0';}

    return K_READ_MORE;
  }

  if ( key == XK_bracketright ) {
    if ( iframe == nframe - 1 ) { sprintf( gmsg, "Last frame" ); return K_NOP; }

    iframe = iframe + fstep;

    if ( iframe >= nframe ) { iframe = nframe - 1; }

    return K_RESETUP;
  }

  /* begin jkl mess */
  if ( key == XK_A ) { /* start animation 1,2,...N,1,2  */
    iframe = 0;

    if ( animate > 0 ) {
      animate = 0;
    }
    else {
      animate = 1;
    }

    return K_RESETUP;
  }

  if ( key == XK_C ) { /* start animation 1,2, ..., N,N-1,N-2, .. 3,1,2,3...*/
    iframe = 0;

    if ( animate > 0 ) {
      animate = 0;
    }
    else {
      animate = 2;
    }

    return K_RESETUP;
  }

  if ( key == XK_less ) {
    sprintf( gmsg, "Animate counterclockwise" );
    framephase = -1;
    return K_NOP;
  }

  if ( key == XK_greater ) {
    sprintf( gmsg, "Animate clockwise" );
    framephase = 1;
    return K_NOP;
  }

  /* end jkl mess */

  if ( key == XK_bracketleft ) {
    if ( iframe == 0 ) {sprintf( gmsg, "First frame" ); return K_NOP; }

    iframe = iframe - fstep;

    if ( iframe <= 0 ) { iframe = 0; }

    return K_RESETUP;
  }

  if ( key == XK_backslash ) {
    saveframe = iframe;
    sprintf( gmsg, "Marked frame %d", saveframe + 1 );
    return K_NOP;
  }

  if ( key == XK_bar ) {
    if ( saveframe == -99 ) {sprintf( emsg, "No frame was marked" ); return K_NOP; }

    if ( saveframe < 0 || saveframe >= nframe ) {
      sprintf( emsg, "Marked frame %d does not exist", saveframe + 1 );
      return K_NOP;
    }

    if ( iframe == saveframe ) {sprintf( gmsg, "Already at frame %d", saveframe + 1 ); return K_NOP; }

    iframe = saveframe;

    sprintf( gmsg, "Go back to frame %d", saveframe + 1 );
    return K_RESETUP;
  }

  if ( key == XK_braceleft ) {iframe = 0; return K_RESETUP; }

  if ( key == XK_braceright ) {iframe = nframe - 1; return K_RESETUP; }

  if ( key == XK_Q ) { return K_QUIT; }

  if ( key == XK_U ) { return K_UPDATE; }

  if ( key == XK_h ) { wrhelp = 1 - wrhelp; return K_REPLOT;}

  if ( key == XK_i ) {strcpy( input, "" ); *inpmode = 1; svline = 0; return K_READ_MORE; }

  *alfa = 0;

  *ixyz = 0;

  if ( key == XK_Right )  {*alfa =  dalfa; *ixyz = 1; return K_REPLOT;}

  if ( key == XK_Left )   {*alfa = -dalfa; *ixyz = 1; return K_REPLOT;}

  if ( key == XK_Up )     {*alfa = -dalfa; *ixyz = 2; return K_REPLOT;}

  if ( key == XK_Down )   {*alfa =  dalfa; *ixyz = 2; return K_REPLOT;}

  if ( key == XK_comma )  {*alfa =  dalfa; *ixyz = 3; return K_REPLOT;}

  if ( key == XK_period ) {*alfa = -dalfa; *ixyz = 3; return K_REPLOT;}

  if ( key == XK_p ) {pmode++; if ( pmode == 3 ) { pmode = 0; }  return K_REPLOT;}

  if ( key == XK_P ) {pmode--; if ( pmode == -1 ) { pmode = 2; } return K_REPLOT;}

  if ( key == XK_d )       {dist0 = dist0 * 1.05 ; return K_REPLOT;}

  if ( key == XK_D )       {dist0 = dist0 / 1.05 ; return K_REPLOT;}

  if ( key == XK_r )       {return K_REPLOT;}

  if ( key == XK_plus || key == XK_KP_Add )
  {scale = scale * 1.05; return K_REPLOT;}

  if ( key == XK_minus || key == XK_KP_Subtract )
  {scale = scale / 1.05; return K_REPLOT;}

  if ( key == XK_KP_6 )   {taux = taux + dtaux; return K_REPLOT;}

  if ( key == XK_KP_4 )   {taux = taux - dtaux; return K_REPLOT;}

  if ( key == XK_KP_8 )   {tauy = tauy + dtauy; return K_REPLOT;}

  if ( key == XK_KP_2 )   {tauy = tauy - dtauy; return K_REPLOT;}

  if ( key == XK_KP_7 )   {taux = taux0; tauy = tauy0; return K_REPLOT;}

  if ( key == XK_KP_Multiply ) {
    taux0 = taux;
    tauy0 = tauy;
    sprintf( gmsg, "New home position: %.2f %.2f", taux0, tauy0 );
    chginfo = 1;
    return K_NOP;
  }

  if ( key == XK_l ) {
    bline = 1 - bline;

    if ( withbonds ) { return K_REPLOT; }

    chginfo = 1;

    if ( bline ) { sprintf( gmsg, "Bonds as lines       " ); }
    else { sprintf( gmsg, "Bonds as cylinders   " ); }

    return K_NOP;
  }

  if ( key == XK_s ) { shadow = 1 - shadow; return K_REPLOT; }

  if ( key == XK_a ) { showaxes = 1 - showaxes; return K_REPLOT; }

  if ( key == XK_w ) { wire = 1 - wire; return K_REPLOT; }

  if ( key == XK_x ) {
    chginfo = 1;
    usepixmap = 1 - usepixmap;

    if ( usepixmap ) { sprintf( gmsg, "Draw to pixmap buffer" ); }
    else { sprintf( gmsg, "Draw to screen" ); }

    chginfo = 1;

    return K_NOP;
  }

  if ( key == XK_space ) {wrinfo = 1 - wrinfo; return K_REPLOT; }

  if ( key == XK_b ) { withbonds = 1 - withbonds; return K_REPLOT; }

  if ( key == XK_n ) {
    numbers = numbers + 1;

    if ( numbers == 3 ) { numbers = 0; }

    return K_REPLOT;
  }

  if ( key == XK_c ) { coords = 1 - coords; return K_REPLOT; }

  if ( key == XK_N ) { bondnums = 1 - bondnums; return K_REPLOT;}

  if ( key == XK_g ) { grayvalues = 1 - grayvalues; return K_REPLOT;}

  return K_NOP;
}


/* ------- wln -----------------  */
void wln( const char x[] ) {

  if ( !startup ) {

    int xx = xln;

    if ( xx < 1 ) { xx = 1; }

    XDrawString( dpy, drw, labelgc, xx, 10 * lnum, x, strlen( x ) );

    lnum++;
  }
  else {
    printf( "%s\n", x );
  }

}

/* ------- WriteHelp -----------------  */
void WriteHelp() {
  lnum = 1;

  /* commands */
  xln = igw - 250;
  wln( "COMMANDS for input line:" );
  wln( "  help  cmd         - on-line help" );
  wln( "  inc   degrees     - rotation increment" );
  wln( "  time  milisecs    - frame duration" );
  wln( "  tell  ijkl        - tell geometry" );
  wln( "  pos   x y         - set position" );
  wln( "  dpos  dxy         - shift increment" );
  wln( "  dist  d           - set distance" );
  wln( "  rfac  fac         - scale all radii" );
  wln( "  bfac  fac         - scale all bonds" );
  wln( "  scale fac         - scale plot" );
  wln( "  gramp slope mid   - gray ramp" );
  wln( "  light x y z       - light direction" );
  wln( "  step  n           - set frame step" );
  wln( "  dup   x y z       - duplicate" );
  wln( "  cut   x y z a b   - planar cut" );
  wln( "  frm   n           - go to frame" );
  wln( "  color pat cname   - query or set color" );
  wln( "  save   fname      - save data" );
  wln( "  update fname      - update from file" );
  wln( "  print  fname      - postscript output" );
  wln( "  close             - close print file" );

  /*  input  */
  xln = igw - 230;
  wln( " " );
  wln( "INPUT data format:" );
  wln( "  atom   label x y z " );
  wln( "  spec   label radius color" );
  wln( "  bonds  lab1 lab2 min max radius color" );
  wln( "  frame  str   x1 ..   - add frame" );
  wln( "  tmat   t1 .. t9      - viewpoint" );
  wln( "Also COMMANDS which set parameters" );

  /*  keys  */
  lnum = 1;
  xln = 10;
  wln( "KEYS (KP=Keypad):" );
  wln( "  i        activate input line" );
  wln( "  h        help" );
  wln( "  space    info" );
  wln( "  cursor   rotate R/L, U/D" );
  wln( "  , .      rotate in plane" );
  wln( "  + -      zoom" );
  wln( "  KP 8642  shift plot" );
  wln( "  KP 7     shift home" );
  wln( "  KP *     set home" );
  wln( "  r        redraw" );
  wln( "  n        atom numbers" );
  wln( "  c        coordinates" );
  wln( "  N        bond lengths" );
  wln( "  w        wire mode" );
  wln( "  g        gray/bw" );
  wln( "  b        bonds" );
  wln( "  l        cylinders/lines" );
  wln( "  p        perspective" );
  wln( "  a        show axes" );
  wln( "  s        line shadows" );
  wln( "  x        pixmap buffer" );
  wln( "  d D      distance" );
  wln( "  [ ]      step frames fw/bk" );
  wln( "  A        animate 1,2,..N,1,2" );
  wln( "  C        animate 1,2,..N,N-1.." );
  wln( "  < >      animate back or forward" );
  wln( "  { }      first/last frame" );
  wln( "  \\ |      mark/goto frame" );
  wln( "  U        update from file" );
  wln( "  Q        quit" );


}

/* ------- WriteInfo -----------------  */
void WriteInfo( Drawable draw ) {
  int i;
  float x1, x2, y1, y2, z1, z2;
  char str[201], stt[201];
  char gmd[3][8] = {"std", "ramp", "light"};

  get_extent( &x1, &x2, &y1, &y2, &z1, &z2 );

  if ( natom == 1 ) { sprintf( str, "%d atom, ", natom ); }
  else { sprintf( str, "%d atoms, ", natom ); }

  if ( nbond == 1 ) { sprintf( stt, "%d bond ", nbond ); }
  else { sprintf( stt, "%d bonds", nbond ); }

  showline( draw, 10, igh - 20, "", str, stt );

  sprintf( str, "extent x %.2f to %.2f,  y %.2f to %.2f,  z %.2f to %.2f",
           x1, x2, y1, y2, z1, z2 );

  showline( draw, 10, igh - 30, "", str, "" );

  sprintf( str, "tmat %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f",
           tmat[0][0], tmat[0][1], tmat[0][2], tmat[1][0], tmat[1][1],
           tmat[1][2], tmat[2][0], tmat[2][1], tmat[2][2] );

  showline( draw, 10, igh - 40, "", str, "" );

  sprintf( str, "dist %.2f, scale %.2f, rfac %.2f, bfac %.2f",
           dist0, scale, radfac, bndfac );

  showline( draw, 10, igh - 50, "", str, "" );

  sprintf( str, "gmode %s, ramp,slope %.2f %.2f, light %.2f %.2f %.2f",
           gmd[gmode], gslope, gz0, light[0], light[1], light[2] );

  showline( draw, 10, igh - 60, "", str, "" );

  sprintf( str, "pos (%.1f,%.1f), dpos %.1f, home (%.1f,%.1f)",
           taux, tauy, dtaux, taux0, tauy0 );

  showline( draw, 10, igh - 70, "", str, "" );

  sprintf( str, "color %d, reverse %d, stippled %d, pixmap %d, gray %d, "
           "bonds %d, lines %d",
           color, reverse, stippled, usepixmap, grayvalues,
           withbonds, bline );

  showline( draw, 10, igh - 80, "", str, "" );

  sprintf( str, "input %s %s, save %s, print %s, nprint=%d",
           inf, inmv, outf, prf, num_print );

  showline( draw, 10, igh - 90, "", str, "" );
}


/* ------- WriteStatus -----------------  */
/* at the lower left side of X-Window    */

void WriteStatus( Drawable draw ) {
  char pers[3][8] = {"off", "pseudo", "true"};
  char str[201];

  sprintf( str, " %5.2f %5.2f %5.2f   inc=%.1f"
           "   d=%.2f   p=%s",
           tmat[2][0], tmat[2][1], tmat[2][2], dalfa,
           dist0, pers[pmode] );

  showline( draw, 10, 20, "View:", str, "" );
  sprintf( str, "Frame %3d of %d (step %d)  <%s>",
           iframe + 1, nframe, fstep, frstr[iframe] );
  showline( draw, 10, 44, "Files: ", curf, "" );
  showline( draw, 10, 32, str, "", "" );

}



