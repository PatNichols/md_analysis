#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <cfloat>
using namespace std;

#define ANGLE_CUT 45.
#define DIST_CUT 3.80

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
        os<<fixed <<showpoint<<setw(w)<<setprecision(p);
        return os;
    };
};

template < int w,int p>
inline ostream& operator<<(ostream& os, Format<w,p>& f)
{
    return f.do_it(os);
}

inline double p_dist(const double *ro,
                     const double *rh,
                     const double& a)
{
    register int ix,iy,iz;
    register double xa,ya,za,ra,r;
    double rc=1.e10;
    const double x=rh[0]-ro[0];
    const double y=rh[1]-ro[1];
    const double z=rh[2]-ro[2];
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
                if (ra<rc) rc=ra;
            }
        }
    }
    return sqrt(rc);
}

inline double px_dist(const double *ro,
                      const double *rh,
                      const double& a)
{
    register int ix,iy,iz;
    register double xa,ya,za,ra;
    const double thres=a*a;
    xa=rh[0]-ro[0];
    ya=rh[1]-ro[1];
    za=rh[2]-ro[2];
    ra=xa*xa+ya*ya+za*za;
    return sqrt(ra);
}

inline double h_bond_dist(const double *ro,
                          const double *rh,
                          double *rn,
                          const double& a)
{
    register int ix,iy,iz;
    register double xa,ya,za,ra;
    double rm=1.e10;

    const double x=rh[0]-ro[0];
    const double y=rh[1]-ro[1];
    const double z=rh[2]-ro[2];
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
                if (ra<rm) {
                    rn[0]=xa;
                    rn[1]=ya;
                    rn[2]=za;
                    rm=ra;
                }
            }
        }
    }
    return sqrt(rm);
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
    double grsh[nbin];
    double gr[nbin];
    int shell_ix[natom];
    int sh_knt[4];
    double rh2[3],rh1[3];
    ifstream in ("nwchem.xyz");
    int i, j, ii, jj, kk, k, nat;
    for (i=0;i<nbin;++i) {
        gr[i]=0;
        grsh[i]=0.;
    }
    double sum2=0.;
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
        for (j = 1; j < natom; ++j)
        {
            if (label[j] == olabel)
            {
                ++ocnt;
                double rux = r[3 * j] - ux;
                double ruy = r[3 * j + 1] - uy;
                double ruz = r[3 * j + 2] - uz;
                double rs=sqrt(rux*rux+ruy*ruy+ruz*ruz);
                if (rs>2.0) {
                    if (rs>3.5) {
                        if (rs>5.25) {
                            shell_ix[j]=3;
                        }else{
                            shell_ix[j]=2;
                        }
                    }else{
                        shell_ix[j]=1;
                    }
                }else{
                    shell_ix[j]=0;
                }
            }
        }
        if (ocnt != natom_o)
        {
            cerr << "Error : did not find all the oxygens" << endl;
            cerr << "Number of oxygens = " << natom_o << endl;
            cerr << "Found " << ocnt << endl;
        }
        double sum1=0.;
        for (j=1;j<natom;++j) {
            if (label[j]==olabel) {
                const double *ro1=r+3*j;
                double rmax=px_dist(ro1,r,a);
                long nn=0;
                for (k=1;k<natom;++k) {
                    if (shell_ix[k]==2 && k!=j && label[k]==olabel) {
                        const double *ro2=r+3*k;
                        const double *rh1=ro2+3;
                        const double *rh2=ro2+6;
                        double b1=p_dist(ro2,rh1,a);
                        double b2=p_dist(ro2,rh2,a);
                        double a1=p_dist(ro1,rh1,a);
                        double a2=p_dist(ro1,rh2,a);
                        double cc=p_dist(ro1,ro2,a);
                        double ang=((a1*a1+b1*b1-cc*cc)*0.5)/a1/b1;
                        ang=180.-(acos(ang)*180./M_PI);
                        double rcut= DIST_CUT- 0.00044*ang*ang;
                        if ((ang<=ANGLE_CUT)&&(cc<=rcut)) {
                            ++nn;
                        }
                        ang=((a2*a2+b2*b2-cc*cc)*0.5)/a2/b2;
                        ang=180.-(acos(ang)*180./M_PI);
                        rcut= DIST_CUT- 0.00044*ang*ang;
                        if ((ang<=ANGLE_CUT)&&(cc<=rcut)) {
                            ++nn;
                        }
                    }
                } // end of loop over 2nd oxygens
                sum1+=(double)nn;
                double xx=rmax/dr+0.5;
                k=(int)xx;
                if (k>=0 && k<=nbin) {
                    gr[k]+=nn;
                    grsh[k]+=1;
                }
            }
        } // end of loop over 1rst oxygens
        cout<<" time = "<<(nframe*25.*2.41889e-5)<<" frame "<<setw(6)<<nframe
        <<" ave acceptor "<<dform<<( sum1 )<<endl;
        ++nframe;
        sum2+=sum1;
    }
    in.close ();
    double avg_acc=sum2/nframe;
    cout<<"average acceptors = "<<dform<<avg_acc<<endl;
    FILE *out=fopen("gracc2nd.dat","w");
    for (i=0;i<nbin;++i)
    {
        double r=i*dr;
        fprintf(out,"%15.6le %15.6le\n",r,(gr[i]/(grsh[i]+DBL_EPSILON)));
    }
    fclose(out);
    return EXIT_SUCCESS;
}
