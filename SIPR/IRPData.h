#pragma once
#include<vector>
#include<iostream>
using namespace std;

class IRPData
{
public:
	IRPData();
	IRPData(char *filename);
	~IRPData();
	const int M;
	const double epsilon;
	int K, N, T;
	vector<double> a, b, B;
	vector<vector<double> > tau_low, tau_up, d_low, d_up, mu, cost;

	void print(ostream &f);

private:
	void readLine(ifstream &ifs, char* chs, int n);
	void getOneLineData(vector<double> &to, int beginIdx, int len, ifstream &ifs, char *buf, int buf_size);
	void readData(char *filename);  //¶ÁÈ¡Êý¾Ý
	void print(ostream &f, const vector<double> &d, const char * s, int i);
	void print(ostream &f, const vector<vector<double> >&d, const char * s, int i, int j);
	void print(char *filename);
};

