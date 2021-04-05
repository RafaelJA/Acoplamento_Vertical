#include "stdafx.h"
#include "unit_cv.h"

double Pa_Psia(double P_Pa)
{
	return P_Pa*0.0001450377l;
}

double Psia_Pa(double P_Psia)
{
	return P_Psia / 0.0001450377l;
}

double Kel_Far(double T_Kel)
{
	return (T_Kel - 273.15l)*9.l / 5.l + 32.l;
}

double Far_Kel(double T_Far)
{
	return (T_Far - 32.l)* 5.l / 9.l + 273.15l;
}

double Kel_Rank(double T_Kel)
{
	return 1.8l*T_Kel;
}

double m_ft(double L_m)
{
	return 3.280839895l*L_m;
}

double ft_m(double L_ft)
{
	return L_ft / 3.280839895l;
}

double in_ft(double L_in) {
	return L_in * 0.083333333l;
}

double ft_in(double L_ft) {
	return L_ft * 0.083333333l;
}

double ft3_m3(double V_ft3)
{
	return V_ft3 / 35.314666721064;
}

double m3_ft3(double V_m3)
{
	return V_m3 * 35.314666721064;
}

double ft3_bbl(double V_ft3) {
	return V_ft3 / 5.61458l;
}

double bbl_ft3(double V_bbl) {
	return V_bbl * 5.61458l;
}

double ft3s_bbld(double Q_ft3s) {
	return Q_ft3s*15644.97217997475l;
}

double bbld_ft3s(double Q_bbld) {
	return Q_bbld / 15644.97217997475l;
}

double cP_Pas(double mu_cP) {
	return 0.001l*mu_cP;
}

double Pas_cP(double mu_Pas) {
	return mu_Pas / 0.001l;
}

double bbld_m3s(double Q_bbld) {
	return 0.158987l*Q_bbld / 86400.l;
}

double m3s_bbld(double Q_m3s) {
	return 86400.l*Q_m3s / 0.158987l;
}

double lbft3_kgm3(double rho_lbft3) {
	return rho_lbft3*16.018463l;
}

double kgm3_lbft3(double rho_kgm3) {
	return rho_kgm3 / 16.018463l;
}
