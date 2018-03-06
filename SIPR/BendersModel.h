#pragma once
#include "IRPData.h"
#include "ilcplex/ilocplex.h"
#include "type.h"
#include "Info.h"
using namespace std;

class BendersModel
{
public:
	char *filename;
	IRPData* data;
	double runTime;
	int intCutN;
	int fCutN;
	int bendersOPCutN;
	int bendersInCutN;
	int bendersFloatOpCutN;
	int bendersFloatInCutN;

	IloNumMatrix e;     // (N+1) * (N+1)

	IloEnv masterEnv, workerEnv;
	IloCplex masterCplex, workerCplex;
	IloModel masterModel, workerModel;

	//varibles for master problem
	IloNumVar eta;
	IloIntVarMatrix y;
	IloNumMatrix yy;
	IloIntVarMatrix3 z;
	IloNumMatrix3 zz;

	//variables for dual subproblem
	IloObjective dualObj;
	IloNumVarMatrix theta_up_n1t, theta_low_n1t, phi_nt, delta;
	IloNumVarMatrix3 theta_up_nk2t, theta_low_nk2t;
	IloNumVarMatrix5 theta_up_nik3t, theta_low_nik3t;


	vector<Info *> Infos;
	RangeMatrix alphaRngs;
	RangeMatrix q_nt_Rngs;

	IloExpr workExpr;
	void createMasterILP();
	void initWorkerVariable();  //set variables in dual set problem
	void createWorkerLP();

	void initWorkExpr();
	void rebuildWorkerLP(const IloNumMatrix &yy);     //rebuild the dual subproblem 
	

	IloExpr buildOptimalBendersExpr(); //build benders optimal benders cut
	IloExpr buildInfeasibleBendersExpr(); //build benders infeasiable cut

	void solve();
	void printStatus(ostream &f);
	void getYandZ();
	void printAlpha(ostream &f);
	void showRoute(ostream &f);
	void showResult();

	BendersModel(){}
	BendersModel(char *filename);
	~BendersModel();

private:
	void createMasterILP2();
};

