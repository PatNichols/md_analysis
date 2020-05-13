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
//  Use first shell water plane for analyze
//
///////////////////////////////////////////////////////////////////////
int
main (int argc, char **argv)
{
    const int ncut=1;
    long nframe=0;
    const double dr = 0.1;
    const double maxr = a * sqrt (3.0) / 2.0;
    const int nbin = 1 + (int) rint (maxr/dr);
    const int nbin2 = 1 + (int) rint (2*maxr/dr);
    double obin[nbin2][nbin];
    double hbin[nbin2][nbin];
    string label[natom];
    const string olabel ("O");
    const string hlabel ("H");
    double r[natom * 3];
    double v[natom * 3];
    ifstream in ("nwchem.xyz");
    int i, j, ii, jj, kk, nat, k;
    for (int i=0;i<nbin2;++i)
    {
        memset (obin[i], 0, sizeof (double) * nbin);
        memset (hbin[i], 0, sizeof (double) * nbin);
    }
    while (!in.eof())
    {
        ++nframe;
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
        int ocnt=0;
        double ux = r[0];
        double uy = r[1];
        double uz = r[2];
        double rcut = 3.2;
        double rcut2= rcut*rcut;
//
// Find an average vector normal to the first shell water plane
// 
        double ax,ay,az;
	ax=ay=az=0.0;
        for (j = 1; j< natom; ++j)
        {
            if (label[j] == olabel)
            {
                ++ocnt;
                double rux1 = r[3 * j] - ux;
                double ruy1 = r[3 * j + 1] - uy;
                double ruz1 = r[3 * j + 2] - uz;
                double r12=rux1*rux1+ruy1*ruy1+ruz1*ruz1;
		if (r12>rcut2 || r12<4.0) continue;
                double rn1=1.0/sqrt(r12);
	        for (k = j+1; k< natom; ++k)
        	{
            		if (label[j] == olabel)
            		{
		                double rux2 = r[3 * k] - ux;
                		double ruy2 = r[3 * k + 1] - uy;
                		double ruz2 = r[3 * k + 2] - uz;
                		double r22=rux2*rux2+ruy2*ruy2+ruz2*ruz2;
				if (r22>rcut2 || r22<4.0) continue;
				double fctt=rn1/sqrt(r22);
                                ax+=fctt*(ruy1*ruz2-ruz1*ruy2);
				ay+=fctt*(ruz1*rux2-rux1*ruz2);
				az+=fctt*(rux1*ruy2-ruy1*rux2);
                        }
                }
	    }
	}
	double ar=1.0/sqrt(ax*ax+ay*ay+az*az);
	ax*=ar;
	ay*=ar;
	az*=ar;
        ocnt = 0;
        for (j = 1; j < natom; ++j)
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
                            double vpar = xx*ax + yy*ay + zz*az;
                            double vper = sqrt(rr-vpar*vpar);
                            int k = (int) (vpar / dr);
                            int l = (int) (vper / dr);
                            k+=(nbin-1);
                            if (k >= 0 && k < nbin2) {
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
        for (j = 3; j < natom; ++j)
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
                            double vpar = xx*ax + yy*ay + zz*az;
                            double vper = sqrt(rr-vpar*vpar);
                            int k = (int) (vpar / dr);
                            int l = (int) (vper / dr);
                            k+=(nbin-1);
                            if (k >= 0 && k < nbin2) {
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
    }
    in.close ();
    fprintf (stderr, "Finished reading input!\n");
    double x = 1.0 / ((double) nframe);
    double xo= x / natom_o;
    double xh= x / natom_h;
    double sum = 0.0;
    for (i = 0; i < nbin2; ++i)
    {
        for (j = 0; j < nbin; ++j)
        {
            sum += obin[i][j];
            obin[i][j] *= xo;
            hbin[i][j] *= xh;
        }
    }
    fprintf (stderr, "simulation time = %15.4le picoseconds\n",
             (nframe * 25. * 2.41889e-5));
    fprintf (stderr, "images          = %15.2le\n", (sum / nframe / natom));
    fprintf (stderr, "nframe         = %12d\n", nframe);
    const double pi2 = 2.0 * M_PI;
    sum = 0.0;
    const double y = a * a * a / natom_o;
    FILE *fp = fopen ("grphi_uo2.dat", "w");
    double vfact=a*a*a/(dr*dr*dr*M_PI);
    for (i = 0; i < nbin2; ++i)
    {
        double xi=i*dr-maxr;
        for (j = 0; j < nbin ;++ j)
        {
            double gr = obin[i][j] * vfact/((double)(j+j+1));
            fprintf (fp, "%15.6le %15.6le %15.6le\n",
                     xi,j*dr,gr);
        }
    }
    fclose (fp);
    sum = 0.0;
    fp = fopen ("grphi_uh2.dat", "w"); 
    for (i = 0; i < nbin2; ++i)
    {
        double xi=i*dr-maxr;
        for (j = 0; j < nbin ;++ j)
        {
            double gr = hbin[i][j] * vfact/((double)(j+j+1));
            fprintf (fp, "%15.6le %15.6le %15.6le\n",
                     xi,j*dr,gr);
        }
    }
    fclose (fp);
    fprintf (stdout,"NBIN2 = %12d",nbin2);
    return EXIT_SUCCESS;
}
