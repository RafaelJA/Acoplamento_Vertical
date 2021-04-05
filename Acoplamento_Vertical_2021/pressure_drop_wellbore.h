#pragma once

double cal_fric_fact_wellbore(double Re, double Rew, double D, double eps);

double cal_fric_fact_lam_prod_mono(double Re, double Rew);

double cal_fric_fact_lam_inj_mono(double Re, double Rew);

double cal_fric_fact_turb_prod_mono(double Re, double Rew, double f0);

double cal_fric_fact_turb_inj_mono(double Re, double Rew, double f0);

double cal_fric_fact_turb_colle(double Re, double D, double eps);

double cal_Reynolds(double rho, double Q, double mu, double D);

double cal_Reynolds_w(double rho, double Q_I, double mu);

double calc_pres_drop_ouyang_mono(double rho, double mu, double D, double D_perf, double perf_dens, double teta, double g, double e, double gama, double deltaL, double Q_med, double U_med, double Q_I);

double calc_f_0_Chen(double eps_0, double Re);

double calc_pres_drop_zhang(double q_w, double q, double deltax, double S, double phi, double A, double D,
	double rho, double mu, double eps_0);

double calc_pres_drop_siwon(double q_w, double q, double deltax, double S, double phi, double A, double D,
	double rho, double mu, double eps_0);

double calc_pres_drop_asheim(double q_w, double q, double deltax, double S, double phi, double A, double D,
	double rho, double mu, double eps_0, double n);

double calc_pres_drop_yalniz(double q_m, double q_i, double D, double A, double d_i, double deltax, double rho, double mu);

double calc_B_yue(double Re, double teta, double phi);

double calc_pres_drop_yue(double Q, double q, double teta, double phi, double D, double deltax, double rho, double mu);