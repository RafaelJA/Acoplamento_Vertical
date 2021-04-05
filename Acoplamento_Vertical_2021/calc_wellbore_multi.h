#pragma once

double calc_q_I(double n, double A_I, double U_IS);

void wellbore_multi(int correlation, int N,
	double *PWB, double *Pb, double *H_l,
	double *Q_I_w, double *Q_I_o, double *Q_I_g,
	double BHP, double PWB_0, double *deltaL,
	double D, double D_perf, double perf_dens,
	double API, double yo, double yg,
	double rho_os, double rho_ws, double rho_gs, double T,
	double g, double teta, double e, double gama,
	double *parametrosmug, double *parametrosBw,
	double Tpc, double ppc, double pi, double Bwi);

void ouyang_homogeneous(double *H_l, double *dpdx, double U_sl, double U_sg, double mu_l, double mu_g,
	double rho_l, double rho_g, double  q_I_l, double q_I_g, double teta, double g,
	double D, double deltaL, double e, double p);

void ouyang_multiphase_model_node(double *H_l, double *dpdx, double H_l_b, double rho_g, double rho_l,
	double U_sg, double U_sl, double q_I_l, double q_I_g, double mu_l, double mu_g,
	double D, double g, double teta, double deltax,
	double D_perf, double perf_dens, double e, double gama);

double calc_hlD(double H_l);

double calc_A_G(double hlD, double D);

double calc_d_max_hinze(double U_m, double f_m, double D, double sigma, double rho_l);

double calc_C_d_barry(double mu_g, double mu_l, double rho_g, double rho_l, double V_im, double d);

bool is_dispersed_bubble_flow(double H_l, double mu_g, double mu_l, double rho_l, double rho_g,
	double U_m, double g, double teta, double V_im, double D, double e);

bool is_stratified_flow(double H_l, double D, double rho_l, double rho_g, double g, double teta, double V_im, double U_sg);

int flow_pattern_ouyang(double H_l, double mu_g, double mu_l, double rho_l, double rho_g,
	double U_m, double g, double teta, double V_im, double D, double U_sg, double e);

double calc_Re_w(double q_I_l, double q_I_g, double rho_l, double rho_g, double mu_l, double mu_g);

void bubble_flow_slip_calc(double *H_l, double *dpdx, double U_sl, double U_sg, double mu_l, double mu_g,
	double rho_l, double rho_g, double  q_I_l, double q_I_g, double sigma, double teta, double g,
	double D, double deltaL, double e);

void bubble_flow_noslip_calc(double *H_l, double *dpdx, double U_sl, double U_sg, double mu_l, double mu_g,
	double rho_l, double rho_g, double  q_I_l, double q_I_g, double sigma, double teta, double g,
	double D, double A, double deltaL, double e);

void intermittent_flow_calc(double *H_l, double *dpdx, double U_sl, double U_sg, double mu_l, double mu_g,
	double rho_l, double rho_g, double  q_I_l, double q_I_g, double sigma, double teta, double g,
	double D, double A, double deltaL, double e);

double calc_fi_stratified(double Fr_l, double Re_sl, double rho_l, double rho_g, double g, double D, double U_g);

void stratified_flow_calc(double *H_l, double *dpdx, double U_sl, double U_sg, double mu_l, double mu_g,
	double rho_l, double rho_g, double  q_I_l, double q_I_g, double sigma, double teta, double g,
	double D, double A, double deltaL, double e);

double calc_Fe_annularmist(double mu_l, double U_sg, double rho_g, double rho_l, double sigma, double U_sl);

double calc_fi_annularmist(double f_c, double Re_f, double sigma, double rho_c, double U_c, double D_c);

void annularmist_flow_calc(double *H_l, double *dpdx, double U_sl, double U_sg, double mu_l, double mu_g,
	double rho_l, double rho_g, double  q_I_l, double q_I_g, double sigma, double teta, double g,
	double D, double A, double deltaL, double e);

double riddermethodforlockhartmartinelli(double V_Dmax, double V_Dmin, int nivel, int n_threads, int pattern, double D, double A,
	double U_sl, double U_sg, double mu_l, double mu_g, double rho_l, double rho_g, double Re_sl, double Re_sg, double Re_w,
	double  q_I_l, double q_I_g, double f_sl, double f_sg, double dpdx_sg, double sigma, double g, double X2, double Y, double e);

double calc_G(int pattern, double V_D, double D, double A, double U_sl, double U_sg, double mu_l, double mu_g,
	double rho_l, double rho_g, double Re_sl, double Re_sg, double Re_w, double  q_I_l, double q_I_g,
	double f_sl, double f_sg, double dpdx_sg, double sigma, double g, double X2, double Y, double e);

void beb_calc(double *H_l, double *dpdx, double U_sl, double U_sg, double mu_l, double mu_g,
	double rho_l, double rho_g, double  q_I_l, double q_I_g, double teta, double g,
	double D, double deltaL, double e, double rho_l_2, double rho_g_2);

void beb_mod_calc(double *H_l, double *dpdx, double U_sl, double U_sg, double mu_l, double mu_g,
	double rho_l, double rho_g, double  q_I_l, double q_I_g, double teta, double g,
	double D, double deltaL, double e);

int beb_pattern(double N_fr, double X);

double beb_horizontal(int pattern, double N_fr, double lambda);

double beb_C(int pattern, double N_fr, double N_lv, double lambda, double teta);

double beb_S(double y);

