#include "GraphUtils.h"
#include "type.h"
#include "flow\graph.h"
#include <vector>
#include <algorithm>

using namespace std;

//graph 邻接表,得到每个连通分支
void GraphUtils::getConnectedCompoent(const vector<vector<int> >& graph, vector<vector<int> >& result,bool *isVisited)
{
	int cur = 0,n=graph.size(),i,j,k;
	fill(isVisited, isVisited + n,false);
	result.clear();
	bool flag;
	for (i = 0; i < n; i++)
	{
		if (!isVisited[i])  //当前节点没被访问过
		{
			cur = i;
			isVisited[cur] = true;
			vector<int> tmp;            // 存储当前 连通分支 节点
			tmp.push_back(cur);
			flag = true;
			while (flag) {
				flag = false;
				for (j = 0; j < graph[cur].size(); j++) {
					k = graph[cur][j];
					if (!isVisited[k]) {
						cur = k;
						tmp.push_back(cur);
						isVisited[cur] = true;
						flag = true;
					}
				}
			}
			if (tmp.size() >= 2) result.push_back(tmp);
		}
	}
}

IloExpr GraphUtils::buildCircleExpr(const vector<int>& Snode, const IloIntVarMatrix3& z, 
	const IloIntVarMatrix& y, const IloNumMatrix& yy, IloInt t,IloEnv &env)
{
	IloExpr temp(env);
	IloInt i, j;
	IloInt max_k = -1;        // 使 y_k^t 最大的k
	IloNum val = -1e8;
	for (i = 0; i < Snode.size(); i++){
		if (yy[Snode[i]][t] > val){
			val = yy[Snode[i]][t];
			max_k = Snode[i];
		}
		for (j = 0; j < Snode.size(); j++){
			if (i != j){
				if (Snode[i]>Snode[j]) temp += z[Snode[j]][Snode[i]][t];   //\sum_{ i,j \in A(S)} z_{ij}^t
				if (Snode[i]<Snode[j]) temp += z[Snode[i]][Snode[j]][t];
			}
		}
	}
	for (i = 0; i < Snode.size(); i++) temp -= 2 * y[Snode[i]][t];     //  \sum_{n \in S} y_n^t
	
	temp += 2 * y[max_k][t];
	return temp;
}

vector<IloExpr> GraphUtils::buildCircleExpr(const vector<int>& Snode, const IloIntVarMatrix3& z,
	const IloNumMatrix3 &zz,const IloIntVarMatrix& y, const IloNumMatrix& yy, IloInt t, IloEnv &env)
{
	int i, j;
	double s = 0;
	vector<IloExpr> result;
	IloExpr expr(env);
	for (i = 0; i < Snode.size(); i++) {
		for (j = 0; j < Snode.size(); j++) {
			if (i != j) {
				if (Snode[i]>Snode[j]) expr += z[Snode[j]][Snode[i]][t];   //\sum_{ i,j \in A(S)} z_{ij}^t
				if (Snode[i]<Snode[j]) expr += z[Snode[i]][Snode[j]][t];
				s += zz[Snode[i]][Snode[j]][t];
			}
		}
	}
	for (i = 0; i < Snode.size(); i++) {
		expr -= 2 * y[Snode[i]][t];
		s -= 2 * yy[Snode[i]][t];
	}
	for (i = 0; i < Snode.size(); i++) {
		if (s + 2 * yy[Snode[i]][t] > 0.01) {
			result.push_back(expr + 2 * y[Snode[i]][t]);
		}
	}
	return result;
}

void GraphUtils::sepNode(Graph<double, double, double> *g,const IloNumMatrix3 &zz,
	int p, int s, int t, vector<int>&s1, vector<int>&s2) 
{
	int N = zz.getSize() - 1,i,j;
	double eps = 1e-3;
	g->reset();
	g->add_node(N +1);
	for (i = 0; i <= N; i++) {
		if (i == s || i == t) continue;
		if (zz[s][i][p] > eps || zz[i][t][p] > eps) {
			g->add_tweights(i, zz[s][i][p], zz[i][t][p]);
		}
	}
	for (i = 0; i <= N; i++) {
		for (j = i+1; j <= N; j++) {
			if (i == s || i == t || j == s || j == t) continue;
			if (zz[i][j][p] > eps) {
				g->add_edge(i, j, zz[i][j][p], zz[i][j][p]);
			}
		}
	}
	g->maxflow();
	s1.clear();
	s2.clear();
	s1.push_back(s);
	s2.push_back(i);
	for (i = 0; i <= N; i++) {
		if (i == s || i == t) continue;
		if (g->what_segment(i) == Graph<double, double, double>::SOURCE) {
			s1.push_back(i);
		}
		else {
			s2.push_back(i);
		}
	}
}


vector<IloExpr> GraphUtils::buildUserCutExpr(const IloIntVarMatrix3 & z, const IloNumMatrix3 & zz, 
	const IloIntVarMatrix & y, const IloNumMatrix & yy, IloEnv & env)
{
	int i, j, k, t,r;
	int N = y.getSize()-1, T = y[0].getSize() - 1;
	double f;
	double eps = 1e-2;
	Graph<double, double, double> *g = new Graph<double, double, double>(N - 1, N*N);
	vector<IloExpr> result;
	vector<int> s1, s2;
	for (t = 1; t <= T; t++) {
		for (i = 1; i <= N; i++) {
			g->reset();
			g->add_node(N - 1);
			for (j = 1; j <i; j++) {
				if (zz[0][j][t] > eps || zz[j][i][t] > eps) {
					g->add_tweights(j - 1, zz[0][j][t], zz[j][i][t]);
				}
			}
			for (j = i + 1; j <= N; j++) {
				if (zz[0][j][t] > eps || zz[j][i][t] > eps) {
					g->add_tweights(j - 2, zz[0][j][t], zz[j][i][t]);
				}
			}

			for (j = 1; j <= N; j++) {
				for (k = j + 1; k <= N; k++) {
					if (j == i || k == i||zz[j][k][t]<eps) continue;
					if (j < i) {
						if (k < i) {
							g->add_edge(j - 1, k - 1, zz[j][k][t], zz[j][k][t]);
						}
						else {
							g->add_edge(j - 1, k - 2, zz[j][k][t], zz[j][k][t]);
						}
					}
					else {
						g->add_edge(j - 2, k - 2, zz[j][k][t], zz[j][k][t]);
					}
				}
			}
			f = g->maxflow();
			s1.clear();
			s2.clear();
			s1.push_back(0);
			s2.push_back(i);
			for (j = 0; j < i - 1;j++) {
				if (g->what_segment(j) == Graph<double, double, double>::SOURCE) {
					s1.push_back(j + 1);
				}
				else {
					s2.push_back(j + 1);
				}
			}
			for (j = i - 1; j <N-1; j++) {
				if (g->what_segment(j) == Graph<double, double, double>::SOURCE) {
					s1.push_back(j + 2);
				}
				else {
					s2.push_back(j + 2);
				}
			}

			//f = 0;
			//IloExpr expr(env);
			//for (j = 0; j < s1.size(); j++) {
			//	for (k = 0; k < s2.size(); k++) {
			//		if (s1[j] < s2[k]) {
			//			expr += z[s1[j]][s2[k]][t];
			//		}
			//		else {
			//			expr += z[s2[k]][s1[j]][t];
			//		}
			//		f += zz[s1[j]][s2[k]][t];
			//	}
			//}
			////add all that don't satisfy
			//for (j = 0; j < s2.size(); j++) {
			//	if (f+0.1< 2*yy[s2[j]][t]) {
			//		result.push_back(expr-2*y[s2[j]][t]);
			//	}
			//}
			
			// only add one constraints
			/*k = 0;
			for (j = 1; j < s2.size(); j++) {
				if (yy[s2[j]][t] > yy[s2[j]][k]) {
					k = j;
				}
			}
			if (f + 0.1< 2 * yy[s2[k]][t]) {
				result.push_back(expr - 2 * y[s2[k]][t]);
			}*/
			
			//expr.end();


			f = 0;
			IloExpr expr(env);
			for (j = 0; j < s2.size(); j++) {
				for (k = 0; k < s2.size(); k++) {
					if (j == k) continue;
					if (s2[j] < s2[k]) {
						expr -= z[s2[j]][s2[k]][t];
					}
					else {
						expr -= z[s2[k]][s2[j]][t];
					}
					f -= zz[s2[j]][s2[k]][t];
				}
				expr += 2*y[s2[j]][t];
				f += 2*yy[s2[j]][t];
			}

			for (j = 0; j < s2.size(); j++) {
				if (f-2*yy[s2[j]][t]+0.01<0) {
					result.push_back(expr - 2*y[s2[j]][t]);
				}
			}
			expr.end();
		}
	}
	delete g;
	return result;
}


void GraphUtils::printExpr(const IloExpr &expr)
{
	typedef IloExpr::LinearIterator ExprItem;
	ExprItem it = expr.getLinearIterator();
	bool first = true;
	while (it.ok()) {
		if (it.getCoef()>0) {
			if (!first) {
				cout << "+";
				first = false;
			}
			cout << it.getCoef() << "*" << it.getVar().getName();
		}
		else {
			cout << it.getCoef() << "*" << it.getVar().getName();
		}
		++it;
	}
}
