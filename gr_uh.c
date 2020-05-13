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

int
main (int argc, char **argv)
{
  long nframe=0;
  const double dr = 0.005;
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
  int i, j, ii, jj, kk, nat;
  while (!in.eof())
    {
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
      int ocnt = 0;
      for (j = 0; j < natom; ++j)
	{
	  if (label[j] == olabel)
	    {
	      ++ocnt;
	      double rux = r[3 * j] - ux;
	      double ruy = r[3 * j + 1] - uy;
	      double ruz = r[3 * j + 2] - uz;
	      for (ii = -1; ii < 2; ++ii)
		{
		  for (jj = -1; jj < 2; ++jj)
		    {
		      for (kk = -1; kk < 2; ++kk)
			{
			  double xx = rux + ii * a;
			  double yy = ruy + jj * a;
			  double zz = ruz + kk * a;
			  double rr = sqrt (xx * xx + yy * yy + zz * zz);
			  int k = (int) rint (rr / dr);
			  if (k >= 0 && k <= nbin)
			    obin[k] += 1.0;
			}
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
      int hcnt = 0;
      for (j = 0; j < natom; ++j)
	{
	  if (label[j] == hlabel)
	    {
	      ++hcnt;
	      double rux = r[3 * j] - ux;
	      double ruy = r[3 * j + 1] - uy;
	      double ruz = r[3 * j + 2] - uz;
	      for (ii = -1; ii < 2; ++ii)
		{
		  for (jj = -1; jj < 2; ++jj)
		    {
		      for (kk = -1; kk < 2; ++kk)
			{
			  double xx = rux + ii * a;
			  double yy = ruy + jj * a;
			  double zz = ruz + kk * a;
			  double rr = sqrt (xx * xx + yy * yy + zz * zz);
			  int k = (int) rint (rr / dr);
			  if (k >= 0 && k <= nbin)
			    hbin[k] += 1.0;
			}
		    }
		}
	    }
	}
      if (hcnt != natom_h)
	{
	  cerr << "Error : did not find all the hydrogens" << endl;
	  cerr << "Number of hydrogens = " << natom_h << endl;
	  cerr << "Found " << hcnt << endl;
	}
      ++nframe;
    }
  in.close ();
  fprintf (stderr, "Finished reading input!\n");
  double x = 1.0 / ((double) nframe);
  double sum = 0.0;
  for (i = 0; i <= nbin; ++i)
    {
      sum += obin[i];
      obin[i] *= x;
      hbin[i] *= x;
    }
  fprintf (stderr, "simulation time = %15.4le picoseconds\n",
	   (nframe * 25. * 2.41889e-5));
  fprintf (stderr, "images          = %15.2le\n", (sum / nframe / natom));
  fprintf (stderr, "nframe         = %12d\n", nframe);
  const double pi4 = 4.0 * M_PI;
  const double pi43 = pi4 / 3.0;
  sum = 0.0;
  const double y = a * a * a / natom_o;
  FILE *fp = fopen ("gr_uo.dat", "w");
  for (i = 0; i <= nbin; ++i)
    {
      double r = (i + 0.5) * dr;
      double vbin = (i * i * 3 + i * 3 + 1) * dr * dr * dr * pi43;
      double gr = obin[i] * y / vbin;
      sum += obin[i];
      fprintf (fp, "%15.6le %15.6le %15.6le %15.6le\n",
	       r, gr, sum, (sum - 2.));
    }
  fclose (fp);
  sum = 0.0;
  fp = fopen ("gr_uh.dat", "w");
  const double yh = a * a * a / natom_h;
  for (i = 0; i <= nbin; ++i)
    {
      double r = (i + 0.5) * dr;
      double vbin = (i * i * 3 + i * 3 + 1) * dr * dr * dr * pi43;
      double gr = hbin[i] * yh / vbin;
      sum += hbin[i];
      fprintf (fp, "%15.6le %15.6le %15.6le\n", r, gr, sum);
    }
  fclose (fp);
  sum = 0.;
  double sum1 = 0.;
  double sum2 = 0.;
  double rmax = 0.;
  double ibmax = 0.;
  int shknt = 0;
  double rc[] = { 2.0, 3.2, 5.25 };
  const int shmax = 3;
  for (i = 0; i <= nbin; ++i)
    {
      double r = (i + 0.5) * dr;
      double vbin = (i * i * 3 + i * 3 + 1) * dr * dr * dr * pi43;
      double gr = obin[i] / vbin;
      double vr = pi4 * r * r;
      if (obin[i] > ibmax)
	{
	  ibmax = obin[i];
	  rmax = r;
	}
      sum += gr * vr;
      sum1 += gr * r * vr;
      sum2 += gr * r * r * vr;
      fprintf (fp, "%15.6le %15.6le %15.6le %15.6le\n",
	       r, gr, sum, (sum - 2.));
      if (r > rc[shknt])
	{
/// annouce results
	  sum *= dr;
	  sum1 *= dr;
	  sum2 *= dr;
	  cout << "Results for shell " << shknt << endl;
	  cout << "G(r) integral   = " << sum << endl;
	  cout << "Rmax            = " << rmax << endl;
	  cout << "Ravg            = " << (sum1 / sum) << endl;
	  sum1 /= sum;
	  sum2 /= sum;
	  sum2 = sum2 - sum1 * sum1;
	  cout << "sigma-r         = " << sqrt(sum2) << endl;
	  // reset variables
	  ++shknt;
	  if (shknt > 2)
	    break;
	  rmax = ibmax = 0.0;
	  sum = sum1 = sum2 = 0.0;
	}
    }
  return EXIT_SUCCESS;
}
