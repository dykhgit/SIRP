#include "BendersModel.h"
#include "Info.h"
#include "IntBendersCallback.h"
#include "bendersUserCut.h"
#include "bendersUserNBcut.h"

using namespace std;


BendersModel::BendersModel(char * filename)
{
	data = new IRPData(filename);
	this->filename = filename;
	fCutN = intCutN = bendersOPCutN = bendersInCutN = bendersFloatOpCutN = bendersFloatInCutN = 0;

	masterModel = IloModel(masterEnv);
	masterCplex = IloCplex(masterEnv);

	workerModel = IloModel(workerEnv);
	workerCplex = IloCplex(workerEnv);
	
	dualObj = IloObjective(workerEnv);
	dualObj.setSense(IloObjective::Maximize);
	workerModel.add(dualObj);
	workerCplex.setOut(workerEnv.getNullStream());
	workerCplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Algorithm::Primal);

	cout << "begin crate master problem" << endl;
	createMasterILP();
	cout << "finished create master problem" << endl;

	cout << "initialize variable in subdual problem" << endl;
	initWorkerVariable();
	cout << "finished variable in subdual problem" << endl;

	cout << "begin create master problem" << endl;
	createWorkerLP();
	cout << "end create master problem" << endl;

	workExpr = IloExpr(workerEnv);
	initWorkExpr();

	masterCplex.setParam(IloCplex::Param::MIP::Strategy::Search, IloCplex::Traditional);
	masterCplex.setParam(IloCplex::TiLim, 7200);
	masterCplex.use(IntBendersLazyCallback(masterEnv, *this));
	//masterCplex.use(mybensersUserCut(masterEnv, *this));
	//masterCplex.use(BendersNBCut(masterEnv, *this));
}

BendersModel::~BendersModel()
{
}

void BendersModel::createMasterILP() {
	int i, j, n, t;
	int N = data->N, T = data->T;
	char name[40];

	z = IloIntVarMatrix3(masterEnv, N + 1);
	zz = IloNumMatrix3(masterEnv, N + 1);
	y = IloIntVarMatrix(masterEnv, N + 1);
	yy = IloNumMatrix(masterEnv, N + 1);
	eta = IloNumVar(masterEnv);
	eta.setName("eta");
	for (i = 0; i <= N; i++) {
		y[i] = IloIntVarArray(masterEnv, T + 1, 0, 1);
		yy[i] = IloNumArray(masterEnv, T + 1);
		for (t = 1; t <= T; t++){
			sprintf_s(name, "y_%d_%d", i, t);
			y[i][t].setName(name);
		}

		z[i] = IloIntVarMatrix(masterEnv, N + 1);
		zz[i] = IloNumMatrix(masterEnv, N + 1);
		for (j = 0; j <= N; j++) {
			zz[i][j] = IloNumArray(masterEnv, T + 1);
			for (t = 1; t <= T; t++) {
				zz[i][j][t] = 0;
			}
		}

		for (j = i+1; j <=N; j++) {
			if (i == 0) {
				z[i][j] = IloIntVarArray(masterEnv, T + 1, 0, 2);
			}
			else {
				z[i][j] = IloIntVarArray(masterEnv, T + 1, 0, 1);
			}
			for (t = 1; t <= T; t++){
				sprintf_s(name, "z_%d.%d_%d", i, j, t);
				z[i][j][t].setName(name);
			}
		}
	}
	for (t = 1; t <= T; t++) {
		for (n = 1; n <= N; n++) {
			masterModel.add(y[n][t] <= y[0][t]);
		}
	}
	for (t = 1; t <= T; t++)
	{
		IloExpr temp(masterEnv);
		for (i = 0; i <= N; i++) {
			for (j = i + 1; j <= N; j++) {
				temp += data->cost[i][j] * z[i][j][t];
			}	
		}
		masterModel.add(temp <= data->B[t]);
		temp.end();

		for (i = 0; i <= N; i++) {
			IloExpr temp1(masterEnv);
			for (j = 0; j <= N; j++) {
				if (j>i) temp1 += z[i][j][t];
				if (j<i) temp1 += z[j][i][t];
			}
			masterModel.add(temp1 == 2 * y[i][t]);
			temp1.end();
		}
	}
	masterModel.add(eta >= 1e-6);
	masterModel.add(IloMinimize(masterEnv, eta));
	masterCplex.extract(masterModel);
}

// for later code optimize
void BendersModel::createMasterILP2() {
	int i, j, n, t;
	int N = data->N, T = data->T;
	char name[20];
	
	z = IloIntVarMatrix3(masterEnv, T + 1);
	zz = IloNumMatrix3(masterEnv, T + 1);
	y = IloIntVarMatrix(masterEnv, T + 1);
	yy = IloNumMatrix(masterEnv, T + 1);
	eta = IloNumVar(masterEnv);
	eta.setName("eta");

	for (t = 0; t <= T; t++) {
		y[t] = IloIntVarArray(masterEnv, N + 1, 0, 1);
		yy[t] = IloNumArray(masterEnv, N + 1);
		z[i] = IloIntVarMatrix(masterEnv, N + 1);
		zz[i] = IloNumMatrix(masterEnv, N + 1);

		for (i = 0; i <= N; i++)
		{
			sprintf_s(name, "y_%d_%d", t, i);
			y[t][i].setName(name);

			z[t][i] = IloIntVarArray(masterEnv, i);  //lower triangle
			zz[t][i] = IloNumArray(masterEnv, i);
			for (j = 0; j < i; j++) {
				if (j == 0) {
					z[t][i][j].setBounds(0, 2);
				}
				else {
					z[t][i][j].setBounds(0, 1);
				}
				sprintf_s(name, "z_%d.%d_%d", t, i, j);
				z[t][i][j].setName(name);
			}
		}
	}

	for (t = 1; t <= T; t++) {
		for (n = 1; n <= N; n++) {
			masterModel.add(y[t][n] <= y[t][0]);
		}
	}
	// z[i][j][t]  i>j
	for (t = 1; t <= T; t++)
	{
		IloExpr temp(masterEnv);
		for (i = 0; i <= N; i++)
			for (j = 0; j <i; j++)
				temp += (data->cost[i][j] * z[t][i][j]);
		masterModel.add(temp <= data->B[t]);
		temp.end();

		for (i = 0; i <= N; i++) {
			IloExpr temp1(masterEnv);
			for (j = 0; j <= N; j++) {
				if (i>j) temp1 += z[t][i][j];
				if (i<j) temp1 += z[t][j][i];
			}
			masterModel.add(temp1 == 2 * y[t][i]);
			temp1.end();
		}
	}
	masterModel.add(eta >= 1e-6);
	masterModel.add(IloMinimize(masterEnv, eta));
	masterCplex.extract(masterModel);
}

void BendersModel::createWorkerLP() {
	int i, j, k, t, n, m, tt;
	int K=data->K,N = data->N, T = data->T;
	double num = 0;

	try
	{
		for (n = 1; n <= N; n++){
			for (t = 2; t <= T; t++){
				IloExpr temp25(workerEnv), temp26(workerEnv);
				temp25 += theta_up_n1t[n][t];
				temp26 += theta_low_n1t[n][t];
				for (k = 1; k <= K; k++){
					temp25 -= theta_up_nk2t[k][n][t];
					temp26 -= theta_low_nk2t[k][n][t];
				}
				workerModel.add(temp25 == 0);  //  25
				workerModel.add(temp26 == 0);  //  26
				temp25.end();
				temp26.end();
			}
		}

		for (n = 1; n <= N; n++){
			for (i = 1; i <= N; i++){
				for (t = 2; t <= T; t++){
					for (tt = 1; tt <= T; tt++){
						IloExpr temp27(workerEnv), temp28(workerEnv);
						temp27 += theta_up_n1t[n][t] * data->mu[i][tt];
						temp28 += theta_low_n1t[n][t] * data->mu[i][tt];
						for (k = 1; k <= K; k++){
							temp27 -= theta_up_nik3t[k][n][i][t][tt];
							temp28 -= theta_low_nik3t[k][n][i][t][tt];
						}
						workerModel.add(temp27 == 0);   //27
						workerModel.add(temp28 == 0);   //28
						temp27.end();
						temp28.end();
					}
				}
			}
		}

		for (k = 1; k <= K; k++){
			for (n = 1; n <= N; n++){
				for (i = 1; i <= N; i++){
					for (t = 2; t <= T; t++){
						for (tt = 1; tt <= T; tt++){
							workerModel.add(theta_up_nk2t[k][n][t] * data->d_low[i][tt] - theta_up_nik3t[k][n][i][t][tt] <= 0); //29
							workerModel.add(theta_low_nk2t[k][n][t] * data->d_low[i][tt] - theta_low_nik3t[k][n][i][t][tt] <= 0); //30
							workerModel.add(-theta_up_nk2t[k][n][t] * data->d_up[i][tt] + theta_up_nik3t[k][n][i][t][tt] <= 0);   //31
							workerModel.add(-theta_low_nk2t[k][n][t] * data->d_up[i][tt] + theta_low_nik3t[k][n][i][t][tt] <= 0);  //32
						}
					}
				}
			}
		}

		for (n = 1; n <= N; n++){
			for (t = 1; t <= T; t++){
				IloExpr temp33(workerEnv);
				temp33 -= phi_nt[n][t];
				for (k = 1; k <= K; k++){
					for (m = t + 1; m <= T; m++)  temp33 += data->a[k] * (theta_low_nk2t[k][n][m] - theta_up_nk2t[k][n][m]);
				}
				IloRange rng= IloRange(workerEnv, -IloInfinity, temp33, 0);
				q_nt_Rngs[n][t] = rng;
				workerModel.add(rng);  //33
				temp33.end();
			}
		}

		for (n = 1; n <= N; n++){
			for (t = 2; t <= T; t++){
				IloExpr temp34(workerEnv);
				temp34 += delta[n][t];
				for (k = 1; k <= K; k++){
					temp34 += data->b[k] * (theta_up_nk2t[k][n][t] + theta_low_nk2t[k][n][t]);
				}
				IloRange rng= IloRange(workerEnv, 1.0, temp34, 1.0);
				alphaRngs[n][t] = rng;
				workerModel.add(rng);  //  34
				temp34.end();
			}
		}
		workerCplex.extract(workerModel);
	}
	catch (IloException& e){
		std::cout << e << std::endl;
		e.end();
	}catch (...){
		std::cout << "Unknown exception\n";
	}
}

void BendersModel::initWorkExpr() {
	int i, j, t, k, tt, n, m;
	int K = data->K, N = data->N, T = data->T;

	IloNum num = 0;
	for (k = 1; k <= K; k++){
		for (t = 2; t <= T; t++){
			for (n = 1; n <= N; n++){
				workExpr -= data->a[k] * data->tau_up[n][t] * theta_up_nk2t[k][n][t];//1
				workExpr += data->a[k] * data->tau_low[n][t] * theta_low_nk2t[k][n][t];//3
			}
		}
	}
	for (k = 1; k <= K; k++){
		for (t = 2; t <= T; t++){
			for (tt = 1; tt <= T; tt++){
				num = 0;
				for (m = 1; m <= t - 1; m++) num += e[m][tt];
				num *= data->a[k];
				for (n = 1; n <= N; n++){
					workExpr -= num*theta_up_nik3t[k][n][n][t][tt];  //2
					workExpr += num*theta_low_nik3t[k][n][n][t][tt]; //4
				}
			}
		}
	}

	for (n = 1; n <= N; n++){
		for (t = 2; t <= T; t++){
			workExpr += data->epsilon*delta[n][t];      //6
		}
	}
}

void BendersModel::rebuildWorkerLP(const IloNumMatrix &yy) {
	int i, j, t, k, tt, n, m;
	int K = data->K, N = data->N, T = data->T;

	workerModel.remove(dualObj);
	IloExpr objExpr = dualObj.getExpr();
	objExpr.clear();

	objExpr += workExpr;
	for (n = 1; n <= N; n++){
		for (t = 1; t <= T; t++){
			objExpr -= data->M*yy[n][t] * phi_nt[n][t];  //5
		}
	}
	dualObj.setExpr(objExpr);
	workerModel.add(dualObj);
	objExpr.end();
}

void BendersModel::initWorkerVariable() {
	int i, j, n, k, t, tt;
	int K=data->K,N = data->N, T = data->T;
	Info *info;
	char name[40];
	e = IloNumMatrix(workerEnv, T + 1); 
	for (t = 0; t <= T; t++) {
		e[t] = IloNumArray(workerEnv, T + 1);
		for (tt = 0; tt <= T; tt++) {
			if (tt == t) e[t][tt] = 1.0;
			else         e[t][tt] = 0;
		}
	}
	theta_up_n1t = IloNumVarMatrix(workerEnv, N + 1);
	theta_low_n1t = IloNumVarMatrix(workerEnv, N + 1);
	phi_nt = IloNumVarMatrix(workerEnv, N + 1);
	delta = IloNumVarMatrix(workerEnv, N + 1);

	alphaRngs = RangeMatrix(workerEnv, N + 1);
	q_nt_Rngs = RangeMatrix(workerEnv, N + 1);
	for (n = 1; n <= N; n++){
		theta_up_n1t[n] =  IloNumVarArray(workerEnv, T + 1, 0, IloInfinity);  //35
		theta_low_n1t[n] = IloNumVarArray(workerEnv, T + 1, 0, IloInfinity); //36
		phi_nt[n] = IloNumVarArray(workerEnv, T + 1, 0, IloInfinity);  //38
		delta[n] =  IloNumVarArray(workerEnv, T + 1, 0, IloInfinity);   //39

		alphaRngs[n] = IloRangeArray(workerEnv, T + 1);
		q_nt_Rngs[n] = IloRangeArray(workerEnv, T + 1);
		for (t = 1; t <= T; t++){
			info = new Info(Info::phi_nt, n, t);
			Infos.push_back(info);
			phi_nt[n][t].setObject(info);
			sprintf_s(name, "phi_nt..%d_%d", n, t);
			phi_nt[n][t].setName(name);

			if (t >= 2) {
				info = new Info(Info::delta_nt, n, t);
				Infos.push_back(info);
				delta[n][t].setObject(info);
				sprintf_s(name, "delta_nt..%d_%d", n, t);
				delta[n][t].setName(name);

				sprintf_s(name, "theta_up_n1t..%d_%d", n, t);
				theta_up_n1t[n][t].setName(name);

				sprintf_s(name, "theta_low_n1t..%d_%d", n, t);
				theta_low_n1t[n][t].setName(name);
			}
		}
	}

	theta_up_nk2t = IloNumVarMatrix3(workerEnv, K + 1);  
	theta_low_nk2t = IloNumVarMatrix3(workerEnv, K + 1);
	for (k = 1; k <= K; k++){
		theta_up_nk2t[k] = IloNumVarMatrix(workerEnv, N + 1);
		theta_low_nk2t[k] = IloNumVarMatrix(workerEnv, N + 1);
		for (n = 1; n <= N; n++){
			theta_up_nk2t[k][n] = IloNumVarArray(workerEnv, T + 1, 0, IloInfinity);  //36
			theta_low_nk2t[k][n] = IloNumVarArray(workerEnv, T + 1,0, IloInfinity);  //36
			for (t = 2; t <= T; t++){
				info = new Info(Info::theta_up_nk2t, k, n, t);
				Infos.push_back(info);
				theta_up_nk2t[k][n][t].setObject(info);
				sprintf_s(name, "theta_up_nk2t..%d_%d_%d", k, n, t);
				theta_up_nk2t[k][n][t].setName(name);


				info = new Info(Info::theta_low_nk2t, k, n, t);
				Infos.push_back(info);
				theta_low_nk2t[k][n][t].setObject(info);
				sprintf_s(name, "theta_low_nk2t..%d_%d_%d", k, n, t);
				theta_low_nk2t[k][n][t].setName(name);
			}
		}
	}

	theta_up_nik3t = IloNumVarMatrix5(workerEnv, K + 1);  
	theta_low_nik3t = IloNumVarMatrix5(workerEnv, K + 1);
	for (k = 1; k <= K; k++){
		theta_up_nik3t[k] = IloNumVarMatrix4(workerEnv, N + 1);
		theta_low_nik3t[k] = IloNumVarMatrix4(workerEnv, N + 1);
		for (n = 1; n <= N; n++){
			theta_up_nik3t[k][n] = IloNumVarMatrix3(workerEnv, N + 1);
			theta_low_nik3t[k][n] = IloNumVarMatrix3(workerEnv, N + 1);
			for (i = 1; i <= N; i++){
				theta_up_nik3t[k][n][i] = IloNumVarMatrix(workerEnv, T + 1);
				theta_low_nik3t[k][n][i] = IloNumVarMatrix(workerEnv, T + 1);
				for (t = 2; t <= T; t++){
					theta_up_nik3t[k][n][i][t] = IloNumVarArray(workerEnv, T + 1, -IloInfinity, IloInfinity);  //37
					theta_low_nik3t[k][n][i][t] = IloNumVarArray(workerEnv, T + 1, -IloInfinity, IloInfinity); //37

					for (tt = 1; tt <= T; tt++){
						if (n == i){
							info = new Info(Info::theta_up_nik3t, k, n, i, t, tt);
							Infos.push_back(info);
							theta_up_nik3t[k][n][i][t][tt].setObject(info);

							info = new Info(Info::theta_low_nik3t, k, n, i, t, tt);
							Infos.push_back(info);
							theta_low_nik3t[k][n][i][t][tt].setObject(info);
						}

						sprintf_s(name, "theta_up_nik3t..%d_%d_%d_%d_%d", k, n, i, t, tt);
						theta_up_nik3t[k][n][i][t][tt].setName(name);

						sprintf_s(name, "theta_low_nik3t..%d_%d_%d_%d_%d", k, n, i, t, tt);
						theta_low_nik3t[k][n][i][t][tt].setName(name);
					}
				}
			}
		}
	}
}

IloExpr BendersModel::buildOptimalBendersExpr() {
	int n, m, t, k, tt;
	int K = data->K, N = data->N, T = data->T;
	IloNum val = -1, num;
	IloExpr temp(masterEnv);
	IloNum b = 0;
	for (k = 1; k <= K; k++){
		for (t = 2; t <= T; t++){
			for (n = 1; n <= N; n++){
				val = workerCplex.getValue(theta_up_nk2t[k][n][t]);
				val = val*data->a[k] * data->tau_up[n][t];
				b -= val;  //1 

				val = workerCplex.getValue(theta_low_nk2t[k][n][t]);
				val = val*data->a[k] * data->tau_low[n][t];
				b += val;  //3
			}
		}
	}
	for (k = 1; k <= K; k++){
		for (t = 2; t <= T; t++){
			for (tt = 1; tt <= T; tt++){
				num = 0;
				for (m = 1; m <= t - 1; m++) num += e[m][tt];
				num = num*data->a[k];

				val = 0;
				for (n = 1; n <= N; n++)  val += workerCplex.getValue(theta_up_nik3t[k][n][n][t][tt]);
				b -= num*val; //2

				val = 0;
				for (n = 1; n <= N; n++)  val += workerCplex.getValue(theta_low_nik3t[k][n][n][t][tt]);
				b += num*val; //4
			}
		}
	}
	for (n = 1; n <= N; n++){
		for (t = 2; t <= T; t++){
			temp -= data->M*y[n][t] * workerCplex.getValue(phi_nt[n][t]);    //6
			b += data->epsilon*workerCplex.getValue(delta[n][t]);          //7
		}
	}
	temp += b;
	return temp;
}

IloExpr BendersModel::buildInfeasibleBendersExpr() {
	int i, j, m, k, n, t, tt;
	IloNum num, val, b = 0;

	IloNumVarArray vars(workerEnv);
	IloNumArray vals(workerEnv);
	workerCplex.getRay(vals, vars);
	IloExpr temp(masterEnv);

	Info *info;
	void* obj;
	
	for (i = 0; i<vars.getSize(); i++){
		obj = (void*)vars[i].getObject();
		if (obj != NULL){
			info = (Info *)obj;
			k = info->k;
			n = info->n;
			t = info->t;
			tt = info->tt;
			switch (info->type){
				case Info::varType::theta_up_nk2t:{
					b -= vals[i] * data->a[k] * data->tau_up[n][t];     //1
					break;
				}
				case Info::varType::theta_low_nk2t:{
					b += vals[i] * data->a[k] * data->tau_low[n][t];   //3
					break;
				}
				case Info::varType::theta_up_nik3t:{
					num = 0;
					for (m = 1; m <= t - 1; m++) num += e[m][tt];
					b -= num*data->a[k]*vals[i];                         // 2
					break;
				}
				case Info::varType::theta_low_nik3t:{
					num = 0;
					for (m = 1; m <= t - 1; m++) num += e[m][tt];
					b += num*data->a[k]*vals[i];                           // 4
					break;
				}
				case Info::varType::phi_nt:{
					temp -= data->M*y[n][t] * vals[i];       //6
					break;
				}
				case Info::varType::delta_nt:{
					b += data->epsilon*vals[i];  //7
					break;
				}
			}
		}
	}
	temp += b;
	vars.end();
	vals.end();
	return temp;
}

void BendersModel::solve() {
	time_t t1 = clock();
	try {
		masterCplex.solve();
	}
	catch (IloException &e) {
		cerr << "ERROR: " << e.getMessage() << endl;
	}
	catch (...) {
		cerr << "Error" << endl;
	}

	time_t t2 = clock();
	runTime = (t2 - t1)*1.0 / 1000;

	if (masterCplex.getStatus() == IloAlgorithm::Optimal|| masterCplex.getStatus() == IloAlgorithm::Feasible) {
		showResult();
	}
	else if (masterCplex.getStatus() == IloAlgorithm::Infeasible) {
		cout << "Infeasible" << endl;
	}
	else if (masterCplex.getStatus() == IloAlgorithm::Unbounded) {
		cout << "Unbounded" << endl;
	}
	else {
		cout << "state=" << masterCplex.getStatus() << endl;
	}
}

void BendersModel::getYandZ() {
	int i, j, t;
	int K = data->K, N = data->N, T = data->T;

	for (t = 1; t <= T; t++) {
		for (i = 0; i <= N; i++) {
			yy[i][t] = masterCplex.getValue(y[i][t]);
			for (j = i + 1; j <= N; j++) {
				zz[j][i][t] = zz[i][j][t] = masterCplex.getValue(z[i][j][t]);
			}
		}
	}
}

void BendersModel::showRoute(ostream &f) {
	int K = data->K, N = data->N, T = data->T;
	int i, j, t;
	int width = 15;
	double total_cost = 0;
	vector<double> per_cost(T + 1, 0);
	vector<vector<int> > g(N + 1, vector<int>());
	vector<vector<int> > connect;
	bool *buf = new bool[N + 1];
	for (t = 1; t <= T; t++) {
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

void BendersModel::printAlpha(ostream &f) {
	int K = data->K, N = data->N, T = data->T;
	int i, j, t, width = 12;
	double val;
	for (t = 1; t <= T; t++) {
		for (i = 0; i <= N; i++) {
			yy[i][t] = masterCplex.getValue(y[i][t]);
		}
	}
	rebuildWorkerLP(yy);
	workerCplex.solve();
	f << "Optimal alpha:" << endl;
	
	for (i = 1; i <= N; i++){
		for (t = 2; t <= T; t++){
			val = workerCplex.getDual(alphaRngs[i][t]);
			f << setiosflags(ios::left) << setw(width) << val << " ";
		}
		f << endl;
	}
}

void BendersModel::printStatus(ostream &f) {
	int width = 45;
	int precision = 8;
	f << setiosflags(ios::left) << setw(width) << "Solution                        status:" << masterCplex.getStatus() << endl;
	f << setiosflags(ios::left) << setw(width) << "Run                               time:" << runTime << endl;

	f << setiosflags(ios::left) << setw(width) << "subtour                     cut Number:" << intCutN << endl;
	f << setiosflags(ios::left) << setw(width) << "benders         optimal     cut Number:" << bendersOPCutN << endl;
	f << setiosflags(ios::left) << setw(width) << "benders         infeasible  cut Number:" << bendersInCutN << endl;

	f << setiosflags(ios::left) << setw(width) << "fraction         subtour    cut Number:" << fCutN << endl;
	f << setiosflags(ios::left) << setw(width) << "fraction benders optimal    cut Number:" << bendersFloatOpCutN << endl;
	f << setiosflags(ios::left) << setw(width) << "fraction benders infeasible cut Number:" << bendersFloatInCutN << endl;

	f << setiosflags(ios::left) << setw(width) << "Optimal                          value:" << masterCplex.getObjValue() << endl;
	f << setiosflags(ios::left) << setw(width) << "MIP          Relative              Gap:" << masterCplex.getMIPRelativeGap() << endl;
}

void BendersModel::showResult() {
	printStatus(cout);
	getYandZ();
	showRoute(cout);
	printAlpha(cout);

	string file(filename);
	file += "_benders.result";
	ofstream f(file);
	printStatus(f);
	showRoute(f);
	printAlpha(f);
	f.close();
}