#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
using namespace std;

/////////////////////////////////////////////
// These variables are simulation specific
//////////////////////////////////////////////
#include "params.h"
//////////////////////////////////////////////

template < int w, int p >
struct Format
{
  public:
     static ostream& do_it(ostream& os) {
          os<<fixed<<showpoint<<setw(w)<<setprecision(p);
          return os;
     };
};

template < int w,int p>
inline ostream& operator<<(ostream& os, Format<w,p>& f)
{
     return f.do_it(os);
}


int
main (int argc, char **argv)
{
  Format<15,4> dform;
  long nframe=0;
  const double dr = 0.01;
  const double maxr = a * sqrt (3.0) / 2.0;
  const int nbin = 1 + (int) rint (maxr / dr);
  string label[natom];
  const string olabel ("O");
  const string hlabel ("H");
  double r[natom * 3];
  double v[natom * 3];
  ifstream in ("nwchem.xyz");
  ofstream out("rtraj.dat");
  int i, j, ii, jj, kk, nat;
  while (!in.eof())
    {
      double timx=nframe*25.0*2.41888e-5;
      in >> nat;
      out<<dform<<timx<<" "; 
      if (nat != natom)
	{
	  cerr << "NUMBER OF ATOMS IS FARGED! " << nat << endl;
	}
      for (j = 0; j < natom; ++j)
	{
	  in >> label[j];
	  in >> r[3 * j];
	  in >> r[3 * j + 1];
	  in >> r[3 * j + 2];
	  in >> v[3 * j];
	  in >> v[3 * j + 1];
	  in >> v[3 * j + 2];
	}
      double ux = r[0];
      double uy = r[1];
      double uz = r[2];
      int ocnt = 2;
      for (j = 3; j < natom; ++j)
	{
	  if (label[j] == olabel)
	    {
	      ++ocnt;
	      double rux = r[3 * j] - ux;
	      double ruy = r[3 * j + 1] - uy;
	      double ruz = r[3 * j + 2] - uz;
              double rus=sqrt(rux*rux+ruy*ruy+ruz*ruz);
              out<<dform<<rus<<" ";
	    }
	}
      out<<endl;
      if (ocnt != natom_o)
	{
	  cerr << "Error : did not find all the oxygens" << endl;
	  cerr << "Number of oxygens = " << natom_o << endl;
	  cerr << "Found " << ocnt << endl;
	}
      ++nframe;
    }
  in.close ();
  out.close();
  return EXIT_SUCCESS;
}
