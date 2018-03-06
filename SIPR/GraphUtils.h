#pragma once
#include <vector>
#include "type.h"
#include "ilcplex\cplex.h"
#include "flow\graph.h"

using namespace std;

class GraphUtils
{
public:
	GraphUtils() {}
	~GraphUtils(){}
	static void getConnectedCompoent(const vector<vector<int> >& graph, vector<vector<int> >& result, bool *buf);

	static IloExpr buildCircleExpr(const vector<int>& Snode, 
		const IloIntVarMatrix3& z, const IloIntVarMatrix& y, const IloNumMatrix& yy, IloInt t,IloEnv &env);

	static vector<IloExpr> buildCircleExpr(const vector<int>& Snode, const IloIntVarMatrix3& z,
		const IloNumMatrix3 &zz, const IloIntVarMatrix& y, const IloNumMatrix& yy, IloInt t, IloEnv &env);

	static void sepNode(Graph<double, double, double> *g,const IloNumMatrix3 &zz, int p, int s, int t, vector<int>&s1, vector<int>&s2);

	static vector<IloExpr> buildUserCutExpr(const IloIntVarMatrix3& z,
		const IloNumMatrix3 &zz, const IloIntVarMatrix& y, const IloNumMatrix& yy, IloEnv &env);

	static void printExpr(const IloExpr &expr);
};

