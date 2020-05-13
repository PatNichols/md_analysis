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
//////////////////////////////////////////////
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
  long nframe=0;
  const double dr = 0.01;
  const double maxr = a * sqrt (3.0) / 2.0;
  const int nbin = 1 + (int) rint (maxr / dr);
  double obin[nbin + 1];
  double hbin[nbin + 1];
  memset (obin, 0, sizeof (double) * (nbin + 1));
  memset (hbin, 0, sizeof (double) * (nbin + 1));
  string label[natom];
  const string olabel ("O");
  const string hlabel ("H");
  double r[natom * 3];
  double v[natom * 3];
  ifstream in ("nwchem.xyz");
  ofstream out4 ("ouo_asym.dat");
  ofstream out3 ("ouo_sym.dat");
  ofstream out2 ("ouo_ang.dat");
  int i, j, ii, jj, kk, nat;
  const double AUTIME=2.41888E-5;
  Format<15,6> dform;
  double suma=0.0;
  double sums=0.0;
  double angs=0.0;
  while (!in.eof())
    {
      double tim=nframe*25.0*AUTIME;
      in >> nat;
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
////////////////// find the data for o-u-o angle ////
      double uo1x=ux-r[3];
      double uo1y=uy-r[4];
      double uo1z=uz-r[5];
      double uo2x=ux-r[6];
      double uo2y=uy-r[7];
      double uo2z=uz-r[8];
      double xx1=uo1x*uo1x+uo1y*uo1y+uo1z*uo1z;
      double xx2=uo2x*uo2x+uo2y*uo2y+uo2z*uo2z;
      double x12=uo1x*uo2x+uo1y*uo2y+uo1z*uo2z;
      double ouoang=(180.0/M_PI)*acos(x12/sqrt(xx1*xx2));
      out2<<tim<<" "<<dform<<ouoang<<endl;
      angs+=ouoang; 
      double ax=uo1x+uo2x;
      double ay=uo1y+uo2y;
      double az=uo1z+uo2z;
      double asq=sqrt((ax*ax+ay*ay+az*az)*0.25); 
      suma+=asq; 
      out3<<tim<<" "<<dform<<asq<<endl;
      double bx=uo1x-uo2x;
      double by=uo1y-uo2y;
      double bz=uo1z-uo2z;
      double bsq=sqrt((bx*bx+by*by+bz*bz)*0.25); 
      out4<<tim<<" "<<dform<<bsq<<endl;
      sums+=bsq;
      ++nframe;
    }
  out4.close();
  out3.close();
  out2.close();
  in.close ();
  sums/=nframe;
  suma/=nframe;
  angs/=nframe;
  fprintf(stdout,"Average O-U-O angle = %15.6lf\n",angs);
  fprintf(stdout,"Average O-U-O length symm = %15.6lf\n",sums);
  fprintf(stdout,"Average O-U-O length asym = %15.6lf\n",suma);
  return EXIT_SUCCESS;
}
