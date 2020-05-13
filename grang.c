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


//////////////////////////////////////////////////////////////////////
//  Analylze the GR(r,theta) as GR( R_perp, R_par)
//
//
///////////////////////////////////////////////////////////////////////
int
main (int argc, char **argv)
{
    const int ncut=2;
    long nframe=0;
    const double dr = 0.05;
    const double maxr = a * sqrt (3.0) / 2.0;
    const int nbin = 1 + (int) rint (maxr / dr);
    double obin[nbin][nbin];
    double hbin[nbin][nbin];
    string label[natom];
    const string olabel ("O");
    const string hlabel ("H");
    double r[natom * 3];
    double v[natom * 3];
    ifstream in ("nwchem.xyz");
    int i, j, ii, jj, kk, nat;
    for (int i=0;i<nbin;++i)
    {
        memset (obin[i], 0, sizeof (double) * nbin);
        memset (hbin[i], 0, sizeof (double) * nbin);
    }
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
        double x1 = r[3];
        double y1 = r[4];
        double z1 = r[5];
        double x2 = r[6];
        double y2 = r[7];
        double z2 = r[8];
        double ax=  x1 + x2;
        double ay=  y1 + y1;
        double az=  z1 + z2;
        double bx=  y1*z2-z1*y2;
        double by=  z1*x2-x1*z2;
        double bz=  x1*y2-y1*x2;
        double axx=  ay*bz-az*by;
        double axy=  az*bx-ax*bz;
        double axz=  ax*by-ay*bx;
        double ax2= sqrt(axx*axx+axy*axy+axz*axz);
        axx/=ax2;
        axy/=ax2;
        axz/=ax2;
        int ocnt = 2;
        for (j = 3; j < natom; ++j)
        {
            if (label[j] == olabel)
            {
                ++ocnt;
                double rux = r[3 * j] - ux;
                double ruy = r[3 * j + 1] - uy;
                double ruz = r[3 * j + 2] - uz;
                for (ii = -ncut; ii <= ncut; ++ii)
                {
                    for (jj = -ncut; jj <= ncut; ++jj)
                    {
                        for (kk = -ncut; kk <= ncut; ++kk)
                        {
                            double xx = rux + ii * a;
                            double yy = ruy + jj * a;
                            double zz = ruz + kk * a;
                            double rr = xx * xx + yy * yy + zz * zz;
                            double vpar = fabs(xx*axx+yy*axy + zz*axz);
                            double vper = sqrt(rr-vpar*vpar);
                            int k = (int) rint (vper / dr);
                            int l = (int) rint (vpar / dr);
                            if (k >= 0 && k < nbin) {
                                if (l >= 0 && l < nbin) {
                                    obin[k][l] += 1.0;
                                }
                            }
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
        for (j = 4; j < natom; ++j)
        {
            if (label[j] == hlabel)
            {
                ++hcnt;
                double rux = r[3 * j] - ux;
                double ruy = r[3 * j + 1] - uy;
                double ruz = r[3 * j + 2] - uz;
                for (ii = -ncut; ii <= ncut; ++ii)
                {
                    for (jj = -ncut; jj <= ncut; ++jj)
                    {
                        for (kk = -ncut; kk <= ncut; ++kk)
                        {
                            double xx = rux + ii * a;
                            double yy = ruy + jj * a;
                            double zz = ruz + kk * a;
                            double rr = xx * xx + yy * yy + zz * zz;
                            double vpar = fabs(xx*axx+yy*axy + zz*axz);
                            double vper = sqrt(rr-vpar*vpar);
                            int k = (int) rint (vpar / dr);
                            int l = (int) rint (vper / dr);
                            if (k >= 0 && k < nbin) {
                                if (l >= 0 && l < nbin) {
                                    hbin[k][l] += 1.0;
                                }
                            }
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
    for (i = 0; i < nbin; ++i)
    {
        for (j = 0; j < nbin; ++j)
        {
            sum += obin[i][j];
            obin[i][j] *= x;
            hbin[i][j] *= x;
        }
    }
    fprintf (stderr, "simulation time = %15.4le picoseconds\n",
             (nframe * 25. * 2.41889e-5));
    fprintf (stderr, "images          = %15.2le\n", (sum / nframe / natom));
    fprintf (stderr, "nframe         = %12d\n", nframe);
    const double pi4 = 4.0 * M_PI;
    const double pi43 = pi4 / 3.0;
    sum = 0.0;
    const double y = a * a * a / natom_o;
    const double fct= 1.0/M_PI;
    FILE *fp = fopen ("grphi_uo.dat", "w");
    for (i = 0; i < nbin; ++i)
    {
        for (j = 0; j < nbin ;++ j)
        {
            double vbin=M_PI*2.*(2*j+1)*dr*dr*dr;
            double gr = obin[i][j] * y / vbin;
            fprintf (fp, "%15.6le %15.6le %15.6le\n",
                     i*dr,j*dr,gr);
        }
    }
    fclose (fp);
    sum = 0.0;
    fp = fopen ("grphi_uh.dat", "w");
    const double yh = a * a * a / natom_h;
    for (i = 0; i < nbin; ++i)
    {
        for (j = 0; j < nbin ;++ j)
        {
            double vbin=M_PI*2.*(2*j+1)*dr*dr*dr;
            double gr = hbin[i][j] * yh / vbin;
            fprintf (fp, "%15.6le %15.6le %15.6le\n",
                     i*dr,j*dr,gr);
        }
    }
    fclose (fp);
    return EXIT_SUCCESS;
}
