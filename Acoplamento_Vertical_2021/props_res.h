#pragma once
//#include <math.h>

//Pressão(psia)(na verdade foda-se)
double calc_phi_p(double phi_ref, double Cr, double p, double pi);

//Pressão(psia)(na verdade foda-se)
double calc_der_phi_p(double phi_ref, double Cr, double p, double pi);

//Pressão(psia), Temperatura(~F}
double calc_water_Rsw_p(double T, double p);

//Pressão(psia), Temperatura(~F}
double calc_der_water_Rsw_p(double T, double p);

//Pressão(psia), Temperatura(~F}
double calc_water_hat_p(double T, double p);

//Pressão(psia), Temperatura(~F}
double calc_der_water_hat_p(double T);

//Pressão(psia), Temperatura(~F}
double calc_water_compr_p(double Tf, double p);

//Pressão(psia), Temperatura(~F}
double calc_der_water_compr_p(double Tf, double p);

//Pressão(psia), Temperatura(~F}, Bw(Rb/STB), rho_ws e saída (g/cm3)
double calc_water_rho_p(double Bwi, double rho_ws, double T, double p, double pi);

//Pressão(psia), Temperatura(~F}, Bw(Rb/STB), rho_ws e saída (g/cm3)
double calc_der_water_rho_p(double Bwi, double rho_ws, double T, double p, double pi);

//Pressão(psia), Temperatura(~F} Saída cp
double calc_water_mi_p(double T, double p);

//Pressão(psia), Temperatura(~F}, Saída cp
double calc_der_water_mi_p(double T, double p);

//Pressão(psia), Temperatura(~F}, Bw(Rb/STB)
double calc_water_Bw_p(double p, double T, double *parametrosBw);

//Pressão(psia), Temperatura(~F}, Bw(Rb/STB)
double calc_der_water_Bw_p(double p, double T, double *parametrosBw);

//
void calc_const_Bw(double *parametrosBw);

//Pressão(psia), Resto adimensionaal
double calc_water_pcw_Sw(double Sw, double Swc, double epsilon, double pcwmin, double pcwmax, double B);

//Pressão(psia), Resto adimensionaal
double calc_der_water_pcw_Sw(double Sw, double Swc, double epsilon, double B);

//Pressão(psia), Resto adimensionaal
double calc_Bpcw(double Swc, double epsilon, double pcwmin, double pcwmax);

//Pressão(psia), Resto adimensionaal
double calc_gas_pcg(double Sw, double So, double Swc, double Sor, double epsilon, double pcgmin, double pcgmax, double B);

//Pressão(psia), Resto adimensionaal
double calc_der_gas_pcg_Swo(double Sw, double So, double Swc, double Sor, double epsilon, double B);

//Pressão(psia), Resto adimensionaal
double calc_Bpcg(double Swc, double Sor, double epsilon, double pcgmin, double pcgmax);

//Pressão(psia), Resto adimensionaal (tá correto)
double calc_water_krw_Sw(double krwmax, double Sw, double Swc, double Sor, double nw);

double calc_water_krw_Sw(double Sw);

//Pressão(psia), Resto adimensionaal (tá certo)
double calc_der_water_krw_Sw(double krwmax, double Sw, double Swc, double Sor, double nw);

double calc_der_water_krw_Sw(double Sw);

//Pressão(psia), Resto adimensionaal (tá certo)
double calc_gas_krg_Swo(double So, double Sw, double Sgr, double Sor, double Swc, double ng);

double calc_gas_krg_Swo(double Sg);

//Pressão(psia), Resto adimensionaal  (tá certo)
double calc_der_gas_krg_Swo(double So, double Sw, double Sgr, double Sor, double Swc, double ng);

double calc_der_gas_krg_Swo(double Sg);

//Pressão(psia), Resto adimensionaal (tá certo)
double calc_water_krow_Sw(double Sw, double Swc, double Sor, double now);

double calc_water_krow_Sw(double Sw);

//Pressão(psia), Resto adimensionaal (tá certo)
double calc_der_water_krow_Sw(double Sw, double Swc, double Sor, double now);

double calc_der_water_krow_Sw(double Sw);

//Pressão(psia), Resto adimensionaal (tá certo)
double calc_gas_krog_Swo(double Sw, double So, double Swc, double Sor, double Sgr, double nog);

double calc_gas_krog_Swo(double Sg);

//Pressão(psia), Resto adimensionaal (tá certo)
double calc_der_gas_krog_Swo(double Sw, double So, double Swc, double Sor, double Sgr, double nog);

double calc_der_gas_krog_Swo(double Sg);

double calc_oil_kro_Swo(double Sw, double Sg);

double calc_der_oil_kro_Sw(double Sw, double Sg);

double calc_der_oil_kro_So(double Sw, double Sg);

double calc_oil_kro_Swo(double krwmax, double Sw, double So, double Sgr, double Swc, double Sor, double nw, double ng, double now, double nog);

double calc_der_oil_kro_Sw(double krwmax, double Sw, double So, double Sgr, double Swc, double Sor, double nw, double ng, double now, double nog);

double calc_der_oil_kro_So(double krwmax, double Sw, double So, double Sgr, double Swc, double Sor, double nw, double ng, double now, double nog);

//Pressão(psia), Temperatura(~F}, Rs(SCF/STB)
double calc_gas_Rs_p(double yg, double p, double pb, double API, double T);

//Pressão(psia), Temperatura(~F}, Rs(SCF/STB)
double calc_gas_Rs_p(double yg, double p, double API, double T);
double calc_der_gas_Rs_p(double yg, double p, double pb, double API, double T);
double calc_der_gas_Rs_pb(double yg, double p, double pb, double API, double T);

//Pressão(psia), Temperatura(~F}, Rs(SCF/STB)
double calc_oil_Bo_p(double yo, double yg, double p, double pb, double API, double T);

double calc_der_oil_Bo_p(double yo, double yg, double p, double API, double T);

double calc_oil_mi_p(double Rs, double T, double API, double p, double pb);

double calc_der_oil_mi_p(double Rs, double dRs, double T, double API, double p, double pb);

double calc_der_oil_mi_pb(double Rs, double dRs, double T, double API, double p, double pb);

//
double calc_oil_rho_p(double Rs, double Bo, double rho_os, double rho_gs);

double calc_der_oil_rho_p(double Rs, double dRs, double Bo, double dBo, double rho_os, double rho_gs);

//
void calc_const_mig(double *parametrosmig, double yg, double yn2, double yco2, double yh2s, double Tf);

//Tpc(R), p(psia)
void calc_pseudocrit(double *ppc, double *Tpc, double yg, double yco2, double yh2s);

//Saída cp
double calc_gas_mig_p(double *parametrosmig, double Tred, double pred);

//Saída cp
double calc_der_gas_mig_p(double *parametrosmig, double Tred, double pred, double ppc, double Tpc);

//Saída sem unidade
double calc_gas_Z_p(double pred, double Tred, double *parametrosZ, double pgr);

//Saída bbl/SCF
double calc_gas_Bg_p(double p);

//Saída bbl/SCF
double calc_der_gas_Bg_p(double p);

double calc_gas_rho_p(double rho_gs, double Bg);

double calc_der_gas_rho_p(double rho_gs, double Bg, double der_Bg);

void calc_props_grid_ant(int N, double *po_0, double *pb_0, double *Sw_0, double *So_0, double *Sg_0, double *X_pb_0,
	double *phi_0, double *Bw_0, double *Rs_0, double *Bo_0, double *Bg_0,
	double *phi_ref, double Cr, double pi, double T, double *parametrosBw,
	double API, double yo, double yg);

void calc_props_grid(int N, double *po, double *pg, double *pw, double *pb, double *Sw, double *So, double *Sg, double *X_pb,
	double *phi, double *der_phi,
	double *rho_w, double *der_rho_w, double *mi_w, double *der_mi_w, double *Bw, double *der_Bw, double *pcw, double *der_pcw, double *krw, double *der_krw,
	double *Rs, double *der_Rs_p, double *der_Rs_pb, double *kro, double *der_kro_Sw, double *der_kro_So, double *Bo, double *der_Bo_p, double *der_Bo_pb,
	double *mi_o, double *der_mi_o_p, double *der_mi_o_pb, double *rho_o, double *der_rho_o_p, double *der_rho_o_pb,
	double *pcg, double *der_pcg_Sw, double *der_pcg_So, double *krg, double *der_krg, double *mi_g, double *der_mi_g, double *Bg, double *der_Bg, double *rho_g, double *der_rho_g,
	double *phi_ref, double Cr, double Bwi, double rho_ws, double pi, double T, double *parametrosBw,
	double Bpcw, double Bpcg, double API, double yo, double yg, double rho_os, double rho_gs,
	double *parametrosmig, double ppc, double Tpc, double pcwmin, double pcwmax,
	double *lambda1_W, double *d_lambda1_po_W, double *d_lambda1_SW_W,
	double *lambda2_W, double *d_lambda2_po_W, double *d_lambda2_SW_W,
	double *lambda1_O, double *d_lambda1_po_O, double *d_lambda1_SW_O, double *d_lambda1_pb_O, double *d_lambda1_SO_O,
	double *lambda2_O, double *d_lambda2_po_O, double *d_lambda2_SW_O, double *d_lambda2_pb_O, double *d_lambda2_SO_O,
	double *lambda1_G, double *d_lambda1_po_G, double *d_lambda1_SW_G, double *d_lambda1_pb_G, double *d_lambda1_SO_G,
	double *lambda2_G, double *d_lambda2_po_G, double *d_lambda2_SW_G, double *d_lambda2_pb_G, double *d_lambda2_SO_G,
	double *lambda3_G, double *d_lambda3_po_G, double *d_lambda3_SW_G, double *d_lambda3_pb_G, double *d_lambda3_SO_G,
	double *lambda4_G, double *d_lambda4_po_G, double *d_lambda4_SW_G, double *d_lambda4_pb_G, double *d_lambda4_SO_G);

double perm_harm(double Kp, double Knm);

//https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
double generateGaussianNoise(double mu, double sigma);

void get_perm(double *k, int dimx, int dimy, int dimz);

void gera_permeabilidades(double *k, int dimx, int dimy, int dimz, double perm, double varperm);

void grava_var(double* p, int fim, int inicio, int res_step, int ite_coup, char directory[], __int8 vartoprint);

void grava_ia(int *p, int fim, int inicio, int j);

void grava_a(double *p, int fim, int inicio, int j);

void grava_b(double *p, int fim, int inicio, int j);

void grava_e(double *p, int fim, int inicio, int j);

void grava_ja(int *p, int fim, int inicio, int j);