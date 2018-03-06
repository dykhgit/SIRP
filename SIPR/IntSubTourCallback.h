#pragma once
#include "ilcplex/ilocplex.h"
#include "SubtourModel.h"
#include "GraphUtils.h"
using namespace std;

ILOLAZYCONSTRAINTCALLBACK1(IntSubTourLazyCallback, SubtourModel &, model)
{
	int i, j, t;
	int K = model.data->K, N = model.data->N, T = model.data->T;
	IloIntVarMatrix &y = model.y;
	IloNumMatrix &yy = model.yy;
	IloIntVarMatrix3 &z = model.z;
	IloNumMatrix3 &zz = model.zz;

	for (t = 1; t <= T; t++) {
		for (i = 0; i <= N; i++) {
			if (getValue(y[i][t]) < 0.1)  yy[i][t] = 0;
			else yy[i][t] = 1;
			for (j = i+1; j<=N; j++) {
				zz[j][i][t] = zz[i][j][t] = getValue(z[i][j][t]);
			}
		}
	}
	vector<vector<int> >g(N+1,vector<int>());
	vector<vector<int> > connect;
	bool *buf = new bool[N + 1];
	for (t = 1; t <= T; t++)
	{
		for (i = 0; i <= N; i++) g[i].clear();
		for (i = 0; i <= N; i++) { 
			for (j = 0; j <= N; j++) {
				if (zz[i][j][t] > 0.5) {
					g[i].push_back(j);
					g[j].push_back(i);
				}
			}
		}
		GraphUtils::getConnectedCompoent(g, connect, buf);
		if (connect.size() > 1) {
			for (i = 1; i < connect.size(); i++)
			{
				vector<int> &tmp = connect.at(i);
				vector<IloExpr> cuts = GraphUtils::buildCircleExpr(tmp, z, zz, y, yy, t, getEnv());
				for (j = 0; j < cuts.size(); j++) {
					add(cuts[j] <= 0);
					cuts[j].end();
					model.intCutN += 1;
				}
			}
		}
	}
	delete[] buf;
}
