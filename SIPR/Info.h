#pragma once
struct Info
{
	enum varType { theta_up_nk2t, theta_low_nk2t, theta_up_nik3t, theta_low_nik3t, phi_nt, delta_nt, unknow };
	varType type;
	int k, n, i, t, tt;
	Info(varType type, int k, int n, int i, int t, int tt) :type(type), k(k), n(n), i(i), t(t), tt(tt) {}

	//theta_up_nik3t, theta_low_nik3t
	Info(varType type, int k, int n, int t, int tt) :type(type), k(k), n(n), t(t), tt(tt), i(n) {}

	//theta_up_nk2t, theta_low_nk2t
	Info(varType type, int k, int n, int t) :type(type), k(k), n(n), t(t), tt(0), i(0) {}

	//phi_nt, delta_nt
	Info(varType type, int n, int t) :type(type), n(n), t(t), k(0), tt(0), i(0) {}
	Info(varType type, int n) :type(type), n(n), k(0), t(0), tt(0), i(0) {}
	Info() :type(varType::unknow), k(0), n(0), t(0), tt(0), i(0) {}
};

