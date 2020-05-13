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
    const int natom = 195;
    const int natom_o = 66;
    const int natom_h = 128;
    const double a = 12.475329;
//////////////////////////////////////////////
//////////////////////////////////////////////
//////////////////////////////////////////////

template < int w, int p > struct Format
{
public:
  static ostream & do_it (ostream & os)
  {
    os << fixed << showpoint << setw (w) << setprecision (p);
    return os;
  };
};

template < int w, int p > inline ostream &
operator<< (ostream & os, Format < w, p > &f)
{
  return f.do_it (os);
}

inline void oh_bond_dist(const double *__restrict__ ro,
const double *__restrict__ rh,
double& __restrict__ x,
double& __restrict__ y,
double& __restrict__ z,
double& __restrict__ r,
const double& __restrict__ a)
{
register int ix,iy,iz;
register double xa,ya,za,ra;
const double thres=2.0;

   x=rh[0]-ro[0];
   y=rh[1]-ro[1];
   z=rh[2]-ro[2];
   for (ix=-1;ix<=1;++ix)
   {
      for (iy=-1;iy<=1;++iy)
      {
         for (iz=-1;iz<=1;++iz)
         {
            xa=x+ix*a;
            ya=y+iy*a;
            za=z+iz*a;
            ra=xa*xa+ya*ya+za*za;
            if (ra<thres) {
               x=xa;
               y=ya;
               z=za;
               r=sqrt(ra);
               return;
            }
         }
      }
   }
}

inline double h2o_bond_angle(const double *__restrict__ ro,
const double *__restrict__ rh1,const double *__restrict__ rh2,
const double &__restrict__ a)
{
   register double dot12,ang;
   register double x1,y1,z1,r1,x2,y2,z2,r2;

   oh_bond_dist(ro,rh1,x1,y1,z1,r1,a);
   oh_bond_dist(ro,rh2,x2,y2,z2,r2,a);
   dot12=x1*x2+y1*y2+z1*z2;
   ang=acos(dot12/r1/r2)*180.0/M_PI;
   return ang;
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
  ofstream out ("hoh_angle_1Shell.dat");
  ofstream out1 ("hoh_angle_2Shell.dat");
  ofstream out2 ("hoh_angle_tot.dat");
  ofstream out3 ("hoh_angle_Others.dat");
  int i, j, ii, jj, kk, nat;
  const double AUTIME = 2.41888E-5;
  Format < 15, 6 > dform;
  double avg_ang1=0.0;
  double avg_ang2=0.0;
  double avg_ang3=0.0;
  double tang_avg=0.0;
  long ang_knt1=0;
  long ang_knt2=0;
  long ang_knt3=0;
  while (!in.eof())
    {
      double tim = nframe * 25.0 * AUTIME;
      out << dform << tim << " ";
      out1 << dform << tim << " ";
      out2 << dform << tim << " ";
      out3 << dform << tim << " ";
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
      int ocnt = 2;
      for (j = 3; j < natom; ++j)
	{
	  if (label[j] == olabel)
	    {
	      ++ocnt;
///////////////////////////////////////////////////////////////
//////////// caution this only works for no hydrolysis!!!!! ///
	      const double *ro = r + 3 * j;
	      const double *rh1 = r + 3 * j + 3;
	      const double *rh2 = r + 3 * j + 6;
	      double rux = ro[0] - ux;
	      double ruy = ro[1] - uy;
	      double ruz = ro[2] - uz;
	      double rus = rux * rux + ruy * ruy + ruz * ruz;
	      double ang = h2o_bond_angle(ro,rh1,rh2,a);
	      tang_avg += ang / 64.0;
	      if (rus < 9.0)
		{
		  avg_ang1 += ang;
      		  ang_knt1++;
                  out<<dform<<ang<<" ";
		}
	      else
                {
                   if (rus>9.0 && rus<27.36) {
			  avg_ang2 += ang;
			  ang_knt2++;
	                  out1<<dform<<ang<<" ";
                   }else{
			avg_ang3+=ang;
                        ang_knt3++;
			out3<<dform<<ang<<" ";
		   }
                }
	    }
	}
      if (ocnt != natom_o)
	{
	  cerr << "Error : did not find all the oxygens" << endl;
	  cerr << "Number of oxygens = " << natom_o << endl;
	  cerr << "Found " << ocnt << endl;
	}
      out << endl;
      out1 << endl;
      out2 << endl;
      ++nframe;
    }
  out3.close();
  out2.close ();
  out.close ();
  in.close ();
  tang_avg = tang_avg / nframe;
  fprintf (stdout, "Average H-O-H angle = %15.6lf\n", tang_avg);
  double ang1 = avg_ang1 / ang_knt1;
  ang_knt1 = ang_knt1 / nframe / 5;
  fprintf (stdout, "average 1rst shell water angle %15.6lf\n", ang1);
  double ang2 = avg_ang2 / ang_knt2;
  ang_knt2 = ang_knt2 / nframe / 5;
  fprintf (stdout, "average 2nd  shell water angle %15.6lf\n", ang2);
  double ang3 = avg_ang3 / ang_knt3;
  ang_knt3 = ang_knt3 / nframe / 5;
  fprintf (stdout, "average water angle others     %15.6lf\n", ang3);
  fprintf (stderr, "simulation time = %15.4le picoseconds\n",
	   (nframe * 25. * 2.41889e-5));
  cout << ang_knt1<< "  "<<ang_knt2<<"  "<<ang_knt3<<endl;
  return EXIT_SUCCESS;
}
