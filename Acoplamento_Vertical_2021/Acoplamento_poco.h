#pragma once

double converge_p_sup(double L_in, int nz_in, double dT_in, double P_in_in, double T_in_in, double OFR_in, double GOR_in, double API_in,
	double SG_G_in, double WC_in, double P_surf_in, double P_initial, double rug_in, int coup_ts, int iter_coup, __int8 correlation, char directory[]);
