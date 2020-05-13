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

void itostr(int i,string& s)
{
	char buff[4];
	buff[0]=(char)((i/100)+48);
	buff[1]=(char)(((i/10)%10)+48);
	buff[2]=(char)((i%10)+48);
	buff[3]='\0';
	s=buff;
}

int
main (int argc, char **argv)
{
   int i,j,k;
    Format<15,6> dform;
    if (argc < 2)
    {
        cerr << "PROPER USAGE IS H20_ANAL NFRAMES" << endl;
        exit (EXIT_FAILURE);
    }
    const int nframe = strtol (argv[1], NULL, 10);
    ofstream out;
    double tim[nframe]; 
    double rstor[nframe];
    double rx[natom_o];
    string pref("traj");
    string suff(".dat");
    string nstr;
    string fstr;
    FILE *in=fopen("bintraj.dat","r");
    for (long iat=0;iat<natom_o;++iat)
    {
	itostr(iat,nstr);
	fstr=pref+nstr+suff;
        rewind(in);
	for (j=0;j<nframe;++j) {
		fread(tim+j,sizeof(double),1,in);
                fread(rx,sizeof(double),natom_o,in);
		rstor[j]=rx[iat];
	}
        bool sw=false;
        double rold=rstor[0];
	for (j=1;j<nframe;++j) {
		double rnew=rstor[j];	  
		if (rold<5.35 && rnew>5.35) {
			sw=true;
			break;
                }
		if (rold>5.35 && rnew<5.35) {
			sw=true;
			break;
		}
                rold=rnew;
	}
	if (sw) {
		out.open(fstr.c_str());
		for (j=0;j<nframe;++j)
		{
                        cout<<iat<<" "<<dform<<tim[j]<<endl;
			out<<dform<<tim[j];
			out<<" ";
			out<<dform<<rstor[j];
			out<<endl;
		}
		out.close();
	}
    }
    fclose(in);
}
		
			
		

		
			    