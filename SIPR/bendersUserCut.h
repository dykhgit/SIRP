#pragma once
#include "ilcplex/ilocplex.h"
#include "BendersModel.h"
#include "GraphUtils.h"
using namespace std;

ILOUSERCUTCALLBACK1(mybensersUserCut, BendersModel &, model) {
	if (!isAfterCutLoop())
		return;
	
	int i, j, t;
	int  N = model.data->N, T = model.data->T;
	IloIntVarMatrix &y = model.y;
	IloNumMatrix &yy = model.yy;
	IloIntVarMatrix3 &z = model.z;
	IloNumMatrix3 &zz = model.zz;

	for (t = 1; t <= T; t++) {
		for (i = 0; i <= N; i++) {
			yy[i][t] = getValue(y[i][t]);
		}
	}

	IloNum masterObjVal = getObjValue();
	model.rebuildWorkerLP(yy);
	model.workerCplex.solve();

	//cout << "worker problem status:" << model.workerCplex.getStatus() << endl;

	if (model.workerCplex.getStatus() == IloAlgorithm::Optimal)
	{
		IloNum workObjVal = model.workerCplex.getObjValue();
		//cout<<"worker problem status:"<<model.workerCplex.getStatus()<<",workerobj="<<workObjVal<<",master Obj="<<masterObjVal<<",diff="<<workObjVal-masterObjVal<<endl;

		if (workObjVal - masterObjVal >= 1e-2) {
			cout<<"user cut callback,try to add benders Optimal cut,"<<endl;
			IloExpr temp = model.buildOptimalBendersExpr();
			add(temp <= model.eta);
			temp.end();
			model.bendersFloatOpCutN += 1;
		}
		else {
			for (t = 1; t <= T; t++) {
				for (i = 0; i <= N; i++) {
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
				cout << "user cut callback,add a subtour cut"<<endl;
			}
		}
	}
	else if (model.workerCplex.getStatus() == IloAlgorithm::Unbounded)
	{
		cout<<"usercut callback In Integer node,try to add benders Infeasible cut"<<endl;
		IloExpr temp = model.buildInfeasibleBendersExpr();
		add(temp <= 0);
		//cout<<"the value of the benders expression is:"<<getValue(temp)<<endl;
		temp.end();
		model.bendersFloatInCutN += 1;
	}
	else {
		cout << "this problem is unbound!!!" << endl;
	}
}