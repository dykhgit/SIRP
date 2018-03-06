#include <iostream>
#include "SubtourModel.h"
#include "BendersModel.h"
using namespace std;


//test some data 
void test1() {
	SubtourModel m("56_6.txt");
	m.solve();
}

void test2() {
	BendersModel m("46_7.txt");
	m.solve();
}

void usage(char *file) {
	cout << "---------usage---------" << endl;
	cout << file << " " << "0" << "filename.... for only lazycallback" << endl;
	cout << file << " " << "1" << "filename.... for usercutcallback and lazycallback" << endl;
	cout << file << " " << "2" << "filename....  for benders decomposation" << endl;
}

void op() {
	int i;
	char a[40] = { '\0' };
	cout << "please input the file :";
	cin >> a;

	cout << "please input solve method"<<endl;
	cout << "0 for lazycallback" << endl;
	cout << "1 for both usercutcallback and lazycallback" << endl;
	cout << "2 for benders decomposation" << endl;
	cin >> i;
	if (i == 0) {
		SubtourModel m(a);
		m.solve();
	}
	else if (i==1){
		SubtourModel m(a,true);
		m.solve();
	}
	else {
		BendersModel m(a);
		m.solve();
	}
	cin >> i;
}

int main(int argc,char *argv[]) {
	int i;
	/*test2();
	cin >> i;
	return 0;*/

	if (argc == 1) {
		op();
	}
	else if( argc==2){
		SubtourModel m(argv[1]);
		m.solve();
	}
	else {
		if (argv[1][0] == '0') {
			SubtourModel m(argv[2]);
			m.solve();
		}
		else if(argv[1][0] == '1'){
			SubtourModel m(argv[2],true);
			m.solve();
		}
		else {
			BendersModel m(argv[2]);
			m.solve();
		}
	}
	return 0;
}