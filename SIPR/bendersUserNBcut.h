#pragma once
#include "ilcplex/ilocplex.h"
#include "BendersModel.h"
#include "GraphUtils.h"
using namespace std;

ILOUSERCUTCALLBACK1(BendersNBCut, BendersModel &, model)
{
	if (!isAfterCutLoop())
		return;
	if (getNnodes() != 0) return;
	int i, j, t;
	int  N = model.data->N, T = model.data->T;
	IloIntVarMatrix &y = model.y;
	IloNumMatrix &yy = model.yy;
	IloIntVarMatrix3 &z = model.z;
	IloNumMatrix3 &zz = model.zz;

	for (t = 1; t <= T; t++) {
		for (i = 0; i <= N; i++) {
			yy[i][t] = getValue(y[i][t]);
			for (j = i + 1; j <= N; j++) {
				zz[j][i][t] = zz[i][j][t] = getValue(z[i][j][t]);
			}
		}
	}

	vector<IloExpr> cuts = GraphUtils::buildUserCutExpr(z, zz, y, yy, getEnv());
	for (i = 0; i < cuts.size(); i++) {
		add(cuts[i] >= 0);
		cuts[i].end();
		model.fCutN += 1;
	}
}