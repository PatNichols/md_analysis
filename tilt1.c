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

inline double pdist(
    const double *__restrict__ r1,
    const double *__restrict__ r2,
    const double& __restrict__ a)
{
    register long i,j,k;
    register double x,y,z,xa,ya,za,ra,rm;

    rm=1.e12;
    x=r1[0]-r2[0];
    y=r1[1]-r2[1];
    z=r1[2]-r2[2];
    for (i=-1;i<=1;++i)
    {
        xa=x+i*a;
        for (j=-1;j<=1;++j)
        {
            ya=y+j*a;
            for (k=-1;k<=1;++k)
            {
                za=z+k*a;
                ra=xa*xa+ya*ya+za*za;
                if (ra<rm) rm=ra;
            }
        }
    }
    return sqrt(rm);
}


inline double pvec_dist(
    const double *__restrict__ r1,
    const double *__restrict__ r2,
    double *__restrict__ rn,
    const double& __restrict__ a)
{
    register long i,j,k;
    register double x,y,z,xa,ya,za,ra,rm;

    rm=1.e12;
    x=r1[0]-r2[0];
    y=r1[1]-r2[1];
    z=r1[2]-r2[2];
    for (i=-1;i<=1;++i)
    {
        xa=x+i*a;
        for (j=-1;j<=1;++j)
        {
            ya=y+j*a;
            for (k=-1;k<=1;++k)
            {
                za=z+k*a;
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

inline void cross_prod(
    const double *__restrict__ a,
    const double *__restrict__ b,
    double *__restrict__ c)
{
    register double cx=a[1]*b[2]-a[2]*b[1];
    register double cy=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
    c[1]=cy;
    c[0]=cx;
}

inline double dot(const double *__restrict__ a,
                  const double *__restrict__ b)
{
    return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}

int
main (int argc, char **argv)
{
    Format<15,4> dform;
    long nframe=0;
    const int natom = 195;
    const int natom_o = 66;
    const int natom_h = 128;
    const double a = 12.475329;
    const double dr = 0.01;
    const double maxr = a * sqrt (3.0) / 2.0;
    string label[natom];
    const string olabel ("O");
    const string hlabel ("H");
    double r[natom * 3];
    double v[natom * 3];
    const double maxth=180.0;
    const double dtheta=maxth/360.;
    const int nbin = 1 + (int) rint (maxth / dtheta);
    double rs3[3],roh1[3],roh2[3],rou[3];
    double gtheta[nbin];
    ifstream in ("nwchem.xyz");
    int i, j, ii, jj, kk, k, nat, m;
    for (i=0;i<nbin;++i) {
        gtheta[i]=0;
    }
    long is_shell[natom];
    double ave,ave2,hcnt;
    ave=ave2=hcnt=0.0;
    double count=0.;
    while (!in.eof()) 
    {
        in >> nat;
        cout<<dform<<(2.41888e-5*25.*nframe);
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
        long oknt=2;
        for (k=3;k<natom;++k)
        {
            if (label[k]==olabel) {
                ++oknt;
                const double *ro=r+3*k;
                const double *rh1=ro+3;
                const double *rh2=ro+6;
                double r0=pvec_dist(ro,r,rou,a);
/////////////////////////// first shell condition ///////////
                if (r0<2.0 || r0>3.5) continue;
/////////////////////////////////////////////////////////////
                double r1=pvec_dist(ro,rh1,roh1,a);
                double r2=pvec_dist(ro,rh2,roh2,a);
                cross_prod(roh1,roh2,rs3);
                double r3=sqrt(dot(rs3,rs3));
                double r03=dot(rs3,rou);
                double theta=acos((r03/r0/r3))*180.0/M_PI;
                double x=theta/dtheta+0.5;
                int m=(int)x;
                if (m>=0 && m<nbin) {
                    gtheta[m]+=1.;
                }
                hcnt+=1.0;
                ave+=theta;
                ave2+=theta*theta;
            }
        }
        cout<<" "<<dform<<hcnt<<" "<<dform<<ave<<endl;
        if (oknt != natom_o)
        {
            cerr << "Error : did not find all the oxygens" << endl;
            cerr << "Number of oxygens = " << natom_o << endl;
            cerr << "Found " << oknt << endl;
        }
	++nframe;
    }
    in.close ();
    ave/=hcnt;
    ave2/=hcnt;
    ave2-=ave*ave;
    cout<<"Sim time        = "<<dform
    <<(nframe*25.0*2.41889e-5)<<" picoseconds\n";
    cout<<"Ave. tilt angle = "<<dform<<ave<<endl;
    cout<<"Ave Dev         = "<<dform<<sqrt(ave2)<<endl;
    cout<<"Ave. Std dev    = "<<dform<<ave2<<endl;
    FILE *out=fopen("tilt_dist1.dat","w");
    const double x=1.0/hcnt;
    for (i=0;i<nbin;++i)
    {
        double r=i*dtheta;
        fprintf(out,"%15.6le %15.6le\n",r,x*gtheta[i]);
    }
    fclose(out);
    return EXIT_SUCCESS;
}
