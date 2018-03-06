#include <math.h>
#include <iostream>
#include <string>
#include "SubtourModel.h"
#include "GraphUtils.h"
#include "IntSubTourCallback.h"
#include "UserCutCallback.h"
using namespace std;


SubtourModel::SubtourModel(char *filename,bool isadd){
	fCutN=intCutN = 0;
	data = new IRPData(filename);
	//cout << "the data you input:" << endl;
	//data->print(cout);

	this->filename = filename;
	model =IloModel(env);
	cplex = IloCplex(env);
	cplex.setParam(IloCplex::TiLim, 7200);
	initVaribales();
	buildModel();
	cplex.extract(model);
	cplex.use(IntSubTourLazyCallback(env,*this));
	if(isadd) cplex.use(myuserCut(env, *this));
	this->isadd = isadd;
}


SubtourModel::~SubtourModel()
{
	delete data;
	model.end();
	cplex.end();
	env.end();
}

void SubtourModel::initVaribales()
{
	int i, j, k, t, tt;
	int K = data->K, N = data->N, T = data->T;
	char name[40];

	e = IloNumMatrix(env, T + 1); //单位矩阵
	for (t = 0; t <= T; t++) {
		e[t] = IloNumArray(env, T + 1);
		for (tt = 0; tt <= T; tt++) {
			if (tt == t) {
				e[t][tt] = 1;
			}
			else {
				e[t][tt] = 0;
			}
		}
	}

	z = IloIntVarMatrix3(env, N + 1);
	zz = IloNumMatrix3(env, N + 1);
	y = IloIntVarMatrix(env, N + 1);
	yy = IloNumMatrix(env, N + 1);
	
	for (i = 0; i <= N; i++) {
		y[i] = IloIntVarArray(env, T + 1, 0, 1);
		yy[i] = IloNumArray(env, T + 1);
		for (t = 1; t <= T; t++)
		{
			sprintf_s(name, "y_%d_%d", i, t);
			y[i][t].setName(name);
		}

		z[i] = IloIntVarMatrix(env, N + 1);
		zz[i] = IloNumMatrix(env, N + 1);
		for (j = 0; j <= N; j++) {
			if (i == 0 || j == 0) {
				z[i][j] = IloIntVarArray(env, T + 1, 0, 2);
			}
			else {
				z[i][j] = IloIntVarArray(env, T + 1, 0, 1);
			}
			zz[i][j] = IloNumArray(env, T + 1);
			for (t = 1; t <= T; t++){
				sprintf_s(name, "z_%d.%d_%d", i, j, t);
				z[i][j][t].setName(name);
				if (i == j) zz[i][j][t] = 0;
			}
		}
	}

	//IloNumVarMatrix alpha, q, s0_up, s0_low; N * T 申请内存(N+1)*(T+1) 使用 1->N * 1->T
	alpha = IloNumVarMatrix(env, N + 1);
	q = IloNumVarMatrix(env, N + 1);
	s0_up = IloNumVarMatrix(env, N + 1);
	s0_low = IloNumVarMatrix(env, N + 1);
	for (i = 1; i <= N; i++) {
		alpha[i] = IloNumVarArray(env, T + 1, data->epsilon, IloInfinity, ILOFLOAT); // 23
		q[i] = IloNumVarArray(env, T + 1, 0, IloInfinity, ILOFLOAT);                 // 16
		s0_up[i] = IloNumVarArray(env, T + 1, -IloInfinity, IloInfinity, ILOFLOAT);
		s0_low[i] = IloNumVarArray(env, T + 1, -IloInfinity, IloInfinity, ILOFLOAT);
	}

	//IloNumVarMatrix4 s_up,s_low,l_up,l_low,h_up,h_low,x,q;  //N * N * T * T
	s_up = IloNumVarMatrix4(env, N + 1);
	s_low = IloNumVarMatrix4(env, N + 1);
	for (i = 1; i <= N; i++) {
		s_up[i] = IloNumVarMatrix3(env, N + 1);
		s_low[i] = IloNumVarMatrix3(env, N + 1);
		for (j = 1; j <= N; j++) {
			s_up[i][j] = IloNumVarMatrix(env, T + 1);
			s_low[i][j] = IloNumVarMatrix(env, T + 1);
			for (t = 1; t <= T; t++) {
				s_up[i][j][t] = IloNumVarArray(env, T + 1, -IloInfinity, IloInfinity, ILOFLOAT);
				s_low[i][j][t] = IloNumVarArray(env, T + 1, -IloInfinity, IloInfinity, ILOFLOAT);
			}
		}
	}

	//IloNumVarMatrix5 u_up,u_low,v_up,v_low;    //K * N * N * T * T 
	u_up = IloNumVarMatrix5(env, K + 1);
	u_low = IloNumVarMatrix5(env, K + 1);
	v_up = IloNumVarMatrix5(env, K + 1);
	v_low = IloNumVarMatrix5(env, K + 1);
	for (k = 1; k <= K; k++){
		u_up[k] = IloNumVarMatrix4(env, N + 1);
		u_low[k] = IloNumVarMatrix4(env, N + 1);
		v_up[k] = IloNumVarMatrix4(env, N + 1);
		v_low[k] = IloNumVarMatrix4(env, N + 1);
		for (i = 1; i <= N; i++){
			u_up[k][i] = IloNumVarMatrix3(env, N + 1);
			u_low[k][i] = IloNumVarMatrix3(env, N + 1);
			v_up[k][i] = IloNumVarMatrix3(env, N + 1);
			v_low[k][i] = IloNumVarMatrix3(env, N + 1);
			for (j = 1; j <= N; j++){
				u_up[k][i][j] = IloNumVarMatrix(env, T + 1);
				u_low[k][i][j] = IloNumVarMatrix(env, T + 1);
				v_up[k][i][j] = IloNumVarMatrix(env, T + 1);
				v_low[k][i][j] = IloNumVarMatrix(env, T + 1);
				for (t = 1; t <= T; t++){
					u_up[k][i][j][t] = IloNumVarArray(env, T + 1, 0.0, IloInfinity, ILOFLOAT);   //9
					u_low[k][i][j][t] = IloNumVarArray(env, T + 1, 0.0, IloInfinity, ILOFLOAT);  //14
					v_up[k][i][j][t] = IloNumVarArray(env, T + 1, 0.0, IloInfinity, ILOFLOAT);   //9
					v_low[k][i][j][t] = IloNumVarArray(env, T + 1, 0.0, IloInfinity, ILOFLOAT);  //14
				}
			}
		}
	}
}

void SubtourModel::buildModel()
{
	int i, j, n, k, t, tt,l;
	int K = data->K, N = data->N, T = data->T;

	for (n = 1; n <= N; n++){
		for (t = 2; t <= T; t++){
			IloExpr temp5(env), temp12(env);
			for (i = 1; i <= N; i++){
				for (tt = 1; tt <= T; tt++){
					temp5 += (data->mu[i][tt] * s_up[n][i][t][tt]);
					temp12 += (data->mu[i][tt] * s_low[n][i][t][tt]);
				}
			}
			model.add((s0_up[n][t] + temp5) >= 0.0);    //5
			model.add((s0_low[n][t] + temp12) >= 0.0);   //10
			temp5.end();
			temp12.end();
		}
	}

	for (k = 1; k <= K; k++){
		for (n = 1; n <= N; n++){
			for (t = 2; t <= T; t++){
				IloExpr temp6(env), temp11(env);
				for (i = 1; i <= N; i++){
					for (tt = 1; tt <= T; tt++){
						temp6 += (data->d_up[i][tt] * u_up[k][n][i][t][tt]-data->d_low[i][tt] * v_up[k][n][i][t][tt]);
						temp11 +=(data->d_up[i][tt] * u_low[k][n][i][t][tt]-data->d_low[i][tt] * v_low[k][n][i][t][tt]);
					}
				}

				temp6 += s0_up[n][t] - data->b[k] * alpha[n][t] - data->a[k] * data->tau_up[n][t];
				temp11 += s0_low[n][t]- data->b[k] * alpha[n][t]+ data->a[k] * data->tau_low[n][t];

				for (l = 1; l <= t - 1; l++) {
					temp6 += data->a[k] * q[n][l];
					temp11-= data->a[k] * q[n][l];
				}
					
				model.add(temp6 == 0);                  //6
				model.add(temp11 == 0);                 //11
				temp6.end();
				temp11.end();
			}
		}
	}
	for (k = 1; k <= K; k++) {
		for (n = 1; n <= N; n++) {
			for (t = 2; t <= T; t++) {
				for (i = 1; i <= N; i++) {
					for (tt = 1; tt <= T; tt++) {
						IloExpr temp7(env), temp12(env);
						temp7 += u_up[k][n][i][t][tt] - v_up[k][n][i][t][tt] - s_up[n][i][t][tt];
						temp12 += (u_low[k][n][i][t][tt] - v_low[k][n][i][t][tt]) - s_low[n][i][t][tt];
						if (i == n) {
							for (l = 1; l <= t - 1; l++) {
								temp7 += data->a[k] * e[l][tt];
								temp12 -= data->a[k] * e[l][tt];
							}
						}
						model.add(temp7 == 0);      //7  8
						model.add(temp12 == 0);     //12 13
						temp7.end();
						temp12.end();
					}
				}
			}
		}
	}

	for (t = 1; t <= T; t++) {
		for (n = 1; n <= N; n++) {
			model.add(q[n][t] <= data->M*y[n][t]); //15
			model.add(y[n][t] <= y[0][t]);        //18
		} 
	}
	// z[i][j][t]  i<j
	for (t = 1; t <= T; t++)
	{
		IloExpr temp17(env);
		for (i = 0; i <= N; i++)
			for (j = i + 1; j <= N; j++)
				temp17 += (data->cost[i][j] * z[i][j][t]);
		model.add(temp17 <= data->B[t]);             //17
		temp17.end();

		for (n = 0; n <= N; n++) {
			IloExpr temp1(env);
			for (i = 0; i <= N; i++) {
				if (n < i) temp1 += z[n][i][t];
				if (n > i) temp1 += z[i][n][t];
			}
			model.add(temp1 == 2 * y[n][t]);   //19
			temp1.end();
		}
	}

	IloExpr obj(env);
	for (n = 1; n <= N; n++) {
		for (t = 2; t <= T; t++) {
			obj += alpha[n][t];
		}
	}
	model.add(IloMinimize(env,obj));          //4
	obj.end();
}

void SubtourModel::solve()
{
	time_t t1 = clock();
	try {
		cplex.solve();
	}
	catch (IloException &e) {
		cerr << "ERROR: " << e.getMessage() << endl;
	}
	catch (...) {
		cerr << "Error" << endl;
	}
	
	time_t t2 = clock();
	runTime = (t2 - t1)*1.0 / 1000;
	if (cplex.getStatus() == IloAlgorithm::Optimal || cplex.getStatus() == IloAlgorithm::Feasible){
		cout<<"ObjValue ="<< cplex.getObjValue()<<endl;
		showResult();
	}
	else if (cplex.getStatus() == IloAlgorithm::Infeasible) {
		cout<<"Infeasible"<<endl;
	}
	else if (cplex.getStatus() == IloAlgorithm::Unbounded) {
		cout<<"Unbounded"<<endl;
	}
	else {
		cout << "state=" << cplex.getStatus() << endl;
	}
}

void SubtourModel::printStatus(ostream &f) {
	int width = 20;
	int precision = 8;
	f << setiosflags(ios::left) << setw(width) << "Solution status:" << cplex.getStatus() << endl;
	f << setiosflags(ios::left) << setw(width) << "Run time:" << runTime << endl;
	f << setiosflags(ios::left) << setw(width) << "int cut Number:" << intCutN << endl;
	f << setiosflags(ios::left) << setw(width) << "float cut Number:" << fCutN << endl;
	f << setiosflags(ios::left) << setw(width) << "Optimal value:" << cplex.getObjValue() << endl;
	f << setiosflags(ios::left) << setw(width) << "MIPRelativeGap:" << cplex.getMIPRelativeGap() << endl;
}

void SubtourModel::getYandZ() {
	int i, j, t;
	int K = data->K, N = data->N, T = data->T;

	for (t = 1; t <= T; t++) {
		for (i = 0; i <= N; i++) {
			yy[i][t] = cplex.getValue(y[i][t]);
			for (j = i + 1; j <= N; j++) {
				zz[j][i][t] = zz[i][j][t] = cplex.getValue(z[i][j][t]);
			}
		}
	}
}

void SubtourModel::printAlpha(ostream &f) {
	int K = data->K, N = data->N, T = data->T;
	int i, t, width=15;
	f << "Optimal Alpha:   " << endl;
	//f << std::fixed << setprecision(6);
	double val;
	for (i = 1; i <= N; i++) {
		for (t = 2; t <= T; t++) {
			val = cplex.getValue(alpha[i][t]);
			f << setiosflags(ios::left) << setw(width)<< val << " ";
		}
		f << endl;
	}
}

void SubtourModel::showRoute(ostream &f) {
	int K = data->K, N = data->N, T = data->T;
	int i, j, t;
	int width = 15;
	double total_cost = 0;
	vector<double> per_cost(T + 1, 0);
	vector<vector<int> > g(N + 1, vector<int>());
	vector<vector<int> > connect;
	bool *buf = new bool[N + 1];
	for (t = 1; t <= T; t++){
		for (i = 0; i <= N; i++) g[i].clear();
		for (i = 0; i <= N; i++) {
			for (j = 0; j <= N; j++) {
				if (zz[i][j][t] > 0.5) {
					g[i].push_back(j);
					g[j].push_back(i);
					per_cost[t] += data->cost[i][j];
				}
			}
		}
		total_cost += per_cost[t];
		f << setiosflags(ios::left) << "period " << t << ":" << endl;
		GraphUtils::getConnectedCompoent(g, connect, buf);
		for (i = 0; i < connect.size(); i++)
		{
			vector<int> &tmp = connect.at(i);
			for (j = 0; j<tmp.size(); j++) f << tmp[j] << "-->";
			f << tmp[0] << endl;
		}
		f << setiosflags(ios::left) << setw(width) << "Budget:" << data->B[t] << endl;
		f << setiosflags(ios::left) << setw(width) << "Cost:" << per_cost[t] / 2 << endl << endl;
	}
	f << "the Total cost is " << total_cost / 2 << endl;
	delete[] buf;
}

void SubtourModel::showResult()
{
	printStatus(cout);
	getYandZ();
	showRoute(cout);
	printAlpha(cout);

	string file(filename);
	if (isadd) {
		file += "_user.result";
	}
	else {
		file += "_lazy.result";
	}
	
	ofstream f(file);
	printStatus(f);
	showRoute(f);
	printAlpha(f);
	f.close();
}

