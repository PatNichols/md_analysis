#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
using namespace std;

#include "params.h"

template < int w, int p >
struct Format
{
public:
    static ostream& do_it(ostream& os) {
        os<<fixed <<showpoint<<setw(w)<<setprecision(p);
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
  Format<15,6> dform;
  long nframe=0;
  long oframe=0;
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
  ofstream out("nwchem.out.xyz");
  int i, j, ii, jj, kk, nat;
  const double AUTIM=2.418889e-5;
  while (!in.eof())
    {
      double t=AUTIM*25*nframe;
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
      if (t<0.5) {
        ++nframe;
	continue;
      }
      out<<setw(3)<<nat<<endl<<endl;
      for (j=0;j<natom;++j) 
      {
	  out<<label[j]<<" ";
          const double *rc=r+3*j;
          const double *vc=v+3*j;
          out<<dform<<rc[0]<<" ";
          out<<dform<<rc[1]<<" ";
          out<<dform<<rc[2]<<" ";
          out<<dform<<vc[0]<<" ";
          out<<dform<<vc[1]<<" ";
          out<<dform<<vc[2]<<" ";
          out<<endl;
      } 
      ++nframe;
      ++oframe;
    }
  in.close ();
  out.close();
  fprintf (stderr, "Finished reading input!\n");
  fprintf (stderr, "total simulation time = %15.4le picoseconds\n",
	   (nframe * 25. * 2.41889e-5));
  fprintf (stderr, "used  simulation time = %15.4le picoseconds\n",
	   (oframe * 25. * 2.41889e-5));
  return EXIT_SUCCESS;
}
