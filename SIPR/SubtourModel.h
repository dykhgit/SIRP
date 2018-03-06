#pragma once
#include "IRPData.h"
#include "type.h"

class SubtourModel
{
public:
	SubtourModel(char *filename, bool isadd=false);
	~SubtourModel();

	char *filename;
	bool isadd;
	IRPData* data;
	double runTime;
	int intCutN;
	int fCutN;

	IloEnv env;
	IloCplex cplex;
	IloModel model;

	IloIntVarMatrix y;   //  (N+1) * T 
	IloNumMatrix yy;     //  (N+1) * T
	IloIntVarMatrix3 z;  // (N+1) * (N+1) * T
	IloNumMatrix3 zz;    // (N+1) * (N+1) * T

	IloNumMatrix e;     // (N+1) * (N+1)
	IloNumVarMatrix alpha, q, s0_up, s0_low; //N * T
	IloNumVarMatrix4 s_up, s_low;  //N * N * T * T
	IloNumVarMatrix5 u_up, u_low, v_up, v_low;                 //K * N * N * T * T 

	void solve();
	void getYandZ();
private:
	void initVaribales();
	void buildModel();
	
	void printStatus(ostream &f);
	void showRoute(ostream &f);
	void printAlpha(ostream &f);
	
	void showResult();
};

