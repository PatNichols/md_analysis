#include <iostream>
using namespace std;
int main()
{
	char ch;
	for (int i=0;i<127;++i) {
		ch=(char)i;
		cout<<i<<" "<<ch<<"\n";
	}
}