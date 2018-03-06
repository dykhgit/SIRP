#include "IRPData.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include<iomanip>
using namespace std;

IRPData::IRPData():M(1000),epsilon(1e-4)
{
}

IRPData::IRPData(char * filename) : M(1000), epsilon(1e-4)
{
	readData(filename);
}

IRPData::~IRPData()
{
}

//read data from file with the format
// line 1: K N T
// line 2: a  ,K numbers
// line 3: b  ,K numbers
// tau_low£¬   N line £¬each line T numbers
// tau_up,     N line £¬each line T numbers
// d_low,      N line £¬each line T numbers
// d_up       ,N line £¬each line T numbers
// mu ,        N line £¬each line T numbers
// cost,       N+1 line £¬each line N+1 numbers
// B:          1 line ,T numbers
void IRPData::readData(char *filename)
{
	int i, n;
	ifstream ifs(filename,ifstream::in);
	if (!ifs) {
		cout<<"No such a file or have no access to read the file"<<endl;
		return;
	}
	else {
		cout << "file " << filename << " is opened successfully" << endl;
	}

	vector<double> v(100, 0);
	char buf[2001];   
	int buf_size = 2000;
	n = 3;
	getOneLineData(v, 0, n, ifs, buf,buf_size);// K N T
	K = (int)(v[0] + 0.1);
	N = (int)(v[1] + 0.1);
	T = (int)(v[2] + 0.1);

	a.resize(K + 1);
	b.resize(K + 1);
	B.resize(T + 1);
	tau_low.resize(N + 1);
	tau_up.resize(N + 1);
	d_low.resize(N + 1);
	d_up.resize(N + 1);
	mu.resize(N + 1);
	cost.resize(N + 1);
	for (i = 1; i<N + 1; i++) {
		tau_low[i].resize(T + 1);
		tau_up[i].resize(T + 1);
		d_low[i].resize(T + 1);
		d_up[i].resize(T + 1);
		mu[i].resize(T + 1);
	}
	for (i = 0; i<N + 1; i++)  cost[i].resize(N + 1);

	getOneLineData(a, 1, K, ifs, buf, buf_size);  // a
	getOneLineData(b, 1, K, ifs, buf, buf_size);  // b
												  
	//tau_low
	for (i = 1; i<N + 1; i++)  getOneLineData(tau_low[i], 1, T, ifs, buf, buf_size);
	
	//tau_up
	for (i = 1; i<N + 1; i++)  getOneLineData(tau_up[i], 1, T, ifs, buf, buf_size);
	
	//d_low
	for (i = 1; i<N + 1; i++)  getOneLineData(d_low[i], 1, T, ifs, buf, buf_size);
	
	//d_up
	for (i = 1; i<N + 1; i++)  getOneLineData(d_up[i], 1, T, ifs, buf, buf_size);
	
	//mu
	for (i = 1; i<N + 1; i++)  getOneLineData(mu[i], 1, T, ifs, buf, buf_size);
	
	//cost
	for (i = 0; i<N + 1; i++)  getOneLineData(cost[i], 0, N+1, ifs, buf, buf_size);
	
	getOneLineData(B, 1, T, ifs, buf, buf_size); // B
	ifs.close();
}


void IRPData::getOneLineData(vector<double> &to, int beginIdx, int len, ifstream &ifs,char *buf,int buf_size)
{
	readLine(ifs, buf, buf_size);
	std::string stringvalues(buf);
	std::istringstream iss(stringvalues);
	double t;
	for (int i = 0; i < len; i++) {
		iss >> t;
		to[i + beginIdx] = t;
	}
}

void IRPData::readLine(ifstream &ifs,char* chs,int n) {
	int i,j;
	bool hasdata = false;
	while (!hasdata) {
		hasdata = false;
		if (ifs.eof()) break;
		ifs.getline(chs, n);
		j = 0;
		for (i = 0; i < n; i++) {
			if (chs[i] >= '0'&&chs[i] <= '9') {
				j += 1;
			}
			else if (chs[i] == '#' || chs[i] == '\0') {
				break;
			}
		}
		//cout << chs << endl;
		chs[i] = '\0';
		if (j > 0) hasdata = true;
	}
}

void IRPData::print(ostream &f, const vector<double> &d, const char * s,int i) {
	if (i < 0) return;
	f << s << endl;
	for (int j = i; j < d.size(); j++) {
		f << setiosflags(ios::left) << setw(10) << d[j] << " ";
	}
	f << endl;
}

void IRPData::print(ostream &f, const vector<vector<double> >&d, const char * s, int i,int j) {
	if (i < 0 || j < 0) return;
	f << s << endl;
	for (int k1 = i; k1 < d.size(); k1++) {
		for (int k2 = j; k2 < d[k1].size(); k2++) {
			f << setiosflags(ios::left) << setw(10) << d[k1][k2] << " ";
		}
		f << endl;
	}
}
void IRPData::print(ostream &f) {
	f << "# K N T" << endl;
	f << K << " " << N << " " << T << endl;
	print(f, a, "# a", 1);
	print(f, b, "# b", 1);
	print(f, tau_low, "# tau_low", 1, 1);
	print(f, tau_up, "# tau_up", 1, 1);
	print(f, d_low, "# d_low", 1, 1);
	print(f, d_up, "# d_up", 1, 1);
	print(f, mu, "# d_mu", 1, 1);
	print(f, cost, "# cost", 0, 0);
	print(f, B, "# B", 1);
}

//used to debug
void IRPData::print(char *filename)
{
	if (filename == NULL || filename[0] == '\0') {
		filename = "copyInputFile.txt";
	}
	ofstream f(filename);
	int i, j;
	f << "# a" << endl;
	for (i = 1; i<a.size(); i++) f << a[i] << " ";
	f << endl;
	print(f, a, "# a", 1);
	f << "# b" << endl;
	for (i = 1; i<b.size(); i++) f << b[i] << " ";
	f << endl;
	f << "# tau_low" << endl;
	for (i = 1; i<tau_low.size(); i++)
	{
		for (j = 1; j<tau_low[i].size(); j++) f << setiosflags(ios::left) << setw(10)<<tau_low[i][j] << " ";
		f << endl;
	}
	f << "# tau_up" << endl;
	for (i = 1; i<tau_up.size(); i++)
	{
		for (j = 1; j<tau_up[i].size(); j++) f << setiosflags(ios::left) << setw(10) << tau_up[i][j] << " ";
		f << endl;
	}
	f << "# d_low" << endl;
	for (i = 1; i<d_low.size(); i++)
	{
		for (j = 1; j<d_low[i].size(); j++) f << setiosflags(ios::left) << setw(10) << d_low[i][j] << " ";
		f << endl;
	}
	f << "# d_up" << endl;
	for (i = 1; i<d_up.size(); i++)
	{
		for (j = 1; j<d_up[i].size(); j++) f << setiosflags(ios::left) << setw(10) << d_up[i][j] << " ";
		f << endl;
	}
	f << "# d_mu" << endl;
	for (i = 1; i<mu.size(); i++)
	{
		for (j = 1; j<mu[i].size(); j++) f << setiosflags(ios::left) << setw(10) << mu[i][j] << " ";
		f << endl;
	}
	f << "# cost" << endl;
	for (i = 0; i<cost.size(); i++)
	{
		for (j = 0; j<cost[i].size(); j++) f << setiosflags(ios::left) << setw(10) << cost[i][j] << " ";
		f << endl;
	}
	f << "# B" << endl;
	for (i = 1; i<B.size(); i++) f << B[i] << " ";
	f << endl;
	f.close();
}
