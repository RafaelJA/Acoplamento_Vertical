#pragma once

void wellbore_mono(int correlation, int N, double *PWB, double *Q_I,
	double BHP, double PWB_0, double *deltaL,
	double rho, double mu,
	double D, double D_perf, double perf_dens,
	double teta, double g, double e, double gama, double teta_yue);