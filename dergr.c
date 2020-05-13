#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

template < int w, int p >
struct Format
{
public:
    static ostream& do_it(ostream& os) {
        os<<scientific <<showpoint<<setw(w)<<setprecision(p);
        return os;
    };
};

template < int w,int p>
inline ostream& operator<<(ostream& os, Format<w,p>& f)
{
    return f.do_it(os);
}

int main()
{
        Format<10,6> dform; 
	double sum1,sum2,gr1,gr0,r1,r0,dgr;
	ifstream in("gr_uo.dat");
	ofstream out("der_gr.dat");
		in>>r0;
		in>>gr0;
		in>>sum1;
		in>>sum2;
        while (!in.eof()) {
		in>>r1;
		in>>gr1;
		in>>sum1;
		in>>sum2;
		dgr=(gr1-gr0)/(r1-r0);
                out<<dform<<r1<<" "<<dform<<dgr<<endl;
		r0=r1;
		gr0=gr1;
	}
	out.close();
        in.close();
}