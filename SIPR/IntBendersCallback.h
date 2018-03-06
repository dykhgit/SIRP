#pragma once
#include "ilcplex/ilocplex.h"
#include "BendersModel.h"
#include "GraphUtils.h"
using namespace std;

ILOLAZYCONSTRAINTCALLBACK1(IntBendersLazyCallback, BendersModel &, model)
{
	int i, j, t;
	IRPData *data = model.data;
	IloIntVarMatrix &y = model.y;
	IloNumMatrix &yy = model.yy;
	IloIntVarMatrix3 &z = model.z;
	IloNumMatrix3 &zz = model.zz;
	int K = data->K, N = data->N, T = data->T;

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
		
		if (workObjVal - masterObjVal >= 1e-6){
			//cout<<"lazy constraints callback,try to add benders Optimal cut,";
			IloExpr temp = model.buildOptimalBendersExpr();
			add(temp <= model.eta);
			temp.end();
			model.bendersOPCutN += 1;
		}
		else { 
			for (t = 1; t <= T; t++) {
				for (i = 0; i <= N; i++) {
					for (j = i + 1; j <= N; j++) {
						zz[j][i][t] = zz[i][j][t] = getValue(z[i][j][t]);
					}
				}
			}
			vector<vector<int> >g(N + 1, vector<int>());
			vector<vector<int> > connect;
			bool *buf = new bool[N + 1];
			for (t = 1; t <= T; t++) {
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
					for (i = 1; i < connect.size(); i++) {
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
	}
	else if (model.workerCplex.getStatus() == IloAlgorithm::Unbounded)
	{
		//cout<<"lazy constraints callback In Integer node,try to add benders Infeasible cut"<<endl;
		IloExpr temp = model.buildInfeasibleBendersExpr();
		add(temp <= 0);
		//cout<<"the value of the benders expression is:"<<getValue(temp)<<endl;
		temp.end();
		model.bendersInCutN += 1;
	}
	else {
		cout << "this problem is unbound!!!" << endl;
	}
}
