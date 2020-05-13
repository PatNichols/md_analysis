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

inline void cross_prod(const double *a,
const double *b,double * c)
{
	register double cz,cy;
	cz=a[0]*b[1]-a[1]*b[0];
        cy=a[2]*b[0]-a[0]*b[2];
        c[0]=a[1]*b[2]-a[2]*b[1];
        c[1]=cy;
	c[2]=cz;
}

int
main (int argc, char **argv)
{
    Format<15,4> dform;
    long nframe=0;
    const double dtang = 0.5*(M_PI/180.0);
    const double maxr = a * sqrt (3.0) / 2.0;
    const int nbin = 1 + (int) rint (M_PI / dtang);
    string label[natom];
    const string olabel ("O");
    const string hlabel ("H");
    double r[natom * 3];
    double v[natom * 3];
    double rnh1[3],rnh2[3],cprod[3];
    double gtilt[nbin];
    double gtilt1[nbin];
    double gtilt2[nbin];
    double gtilt3[nbin];
    long ntilt1,ntilt2,ntilt3,ntilt,shell_ix;
    const double half_pi=0.5*M_PI;
    ifstream in ("nwchem.xyz");
    int i, j, ii, jj, kk, k, nat;
    for (i=0;i<nbin;++i) {
        gtilt[i]=0;
        gtilt1[i]=0;
        gtilt2[i]=0;
        gtilt3[i]=0;
    }
    ntilt=ntilt1=ntilt2=ntilt3=0;
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
        int ocnt = 2;
        for (j = 3; j < natom; ++j)
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
                            shell_ix=3;
                        }else{
                            shell_ix=2;
                        }
                    }else{
                        shell_ix=1;
                    }
                }else{
                    shell_ix=0;
                }
                double *ro=r+3*j;
                double *rh1=ro+3;
                double *rh2=ro+6;
                double rho1=h_bond_dist(ro,rh1,rnh1,a);
                double rho2=h_bond_dist(ro,rh2,rnh2,a);
                cross_prod(rnh1,rnh2,cprod);
                double rdot=rux*cprod[0]+ruy*cprod[1]+ruz*cprod[2];
                double ram=sqrt(cprod[0]*cprod[0]+cprod[1]*cprod[1]+cprod[2]*cprod[2]);
                double tang=asin(rdot/ram/rs);
                tang=tang+half_pi;
                int kang=static_cast<int>((tang/dtang+0.5));
                gtilt[kang]+=1;
                ++ntilt;
                switch(shell_ix) {
                        case 0:
                                break;
			case 1:
				gtilt1[kang]+=1;
				++ntilt1;
                                break;
			case 2:
				gtilt2[kang]+=1;
				++ntilt2;
                                break;
			case 3:
				gtilt3[kang]+=1;
				++ntilt3;
                                break;
  			default:
				cerr<<"YOU ARE SCREWED!!!!"<<endl;
                                exit(EXIT_FAILURE);
		}    
            }
        }
        if (ocnt != natom_o)
        {
            cerr << "Error : did not find all the oxygens" << endl;
            cerr << "Number of oxygens = " << natom_o << endl;
            cerr << "Found " << ocnt << endl;
        }
        cout<<" time = "<<(nframe*25.*2.41889e-5)<<" frame "<<setw(6)<<nframe<<"\n";
        ++nframe;
    }
    in.close ();
    const double conv=180.0/M_PI;
//  const double offset=half_pi;
    const double offset=0.0;
    double av,av2;
    cout<<"time "<<(nframe*25*2.41889e-5)<<" psec\n";
    FILE *out=fopen("gtiltx.dat","w");
    av=av2=0.0;
    for (i=0;i<nbin;++i)
    {
        double r=(i*dtang-offset)*conv;
        double gn=gtilt[i]/(ntilt+DBL_EPSILON);
        fprintf(out,"%15.6le %15.6le\n",r,gn);
        av+=gn*r;
        av2+=gn*gn*r*r;
    }
    av2=av2-av*av;
    fclose(out);
    fprintf(stdout,"tilt all\n"); 
    fprintf(stdout,"Average Ang= %15.6lf\n",av);
    fprintf(stdout,"Std. Dev   = %15.6lf\n",sqrt(av));
    out=fopen("gtilt1x.dat","w");
    av=av2=0.0;
    for (i=0;i<nbin;++i)
    {
        double r=(i*dtang-offset)*conv;
        double gn=gtilt1[i]/(ntilt1+DBL_EPSILON);
        fprintf(out,"%15.6le %15.6le\n",r,gn);
        av+=gn*r;
        av2+=gn*gn*r*r;
    }
    av2=av2-av*av;
    fclose(out);
    fprintf(stdout,"tilt shell 1\n"); 
    fprintf(stdout,"Average Ang= %15.6lf\n",av);
    fprintf(stdout,"Std. Dev   = %15.6lf\n",sqrt(av));
    out=fopen("gtilt2x.dat","w");
    av=av2=0.0;
    for (i=0;i<nbin;++i)
    {
        double r=(i*dtang-offset)*conv;
        double gn=gtilt2[i]/(ntilt2+DBL_EPSILON);
        fprintf(out,"%15.6le %15.6le\n",r,gn);
        av+=gn*r;
        av2+=gn*gn*r*r;
    }
    av2=av2-av*av;
    fclose(out);
    fprintf(stdout,"tilt shell 2\n"); 
    fprintf(stdout,"Average Ang= %15.6lf\n",av);
    fprintf(stdout,"Std. Dev   = %15.6lf\n",sqrt(av));
    out=fopen("gtilt3x.dat","w");
    av=av2=0.0;
    for (i=0;i<nbin;++i)
    {
        double r=(i*dtang-offset)*conv;
        double gn=gtilt3[i]/(ntilt3+DBL_EPSILON);
        fprintf(out,"%15.6le %15.6le\n",r,gn);
        av+=gn*r;
        av2+=gn*gn*r*r;
    }
    av2=av2-av*av;
    fclose(out);
    fprintf(stdout,"tilt shell 3\n"); 
    fprintf(stdout,"Average Ang= %15.6lf\n",av);
    fprintf(stdout,"Std. Dev   = %15.6lf\n",sqrt(av));
    return EXIT_SUCCESS;
}
