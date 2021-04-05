#include "stdafx.h"
#include "calc_wellbore_multi.h"

double calc_q_I(double n, double A_I, double U_IS) {

	//n perforation density
	//Cross-sectional area of each perforation
	//Velocities in the perforation	

	return n * A_I * U_IS;


}

void wellbore_multi(int correlation, int N,
	double* PWB, double* Pb, double* H_l,
	double* Q_I_w, double* Q_I_o, double* Q_I_g,
	double BHP, double PWB_0, double* deltaL,
	double D, double D_perf, double perf_dens,
	double API, double yo, double yg,
	double rho_os, double rho_ws, double rho_gs, double T,
	double g, double teta, double e, double gama,
	double* parametrosmug, double* parametrosBw,
	double Tpc, double ppc, double pi, double Bwi) {


	/*  1- Ouyang Mechanicist
	2- Ouyang Homogeneous
	3- Beggs & Brill Original
	4- Beggs & Brill Modified*/

	double* dpdz = (double*)malloc(N * sizeof(double));
	double* q_I_w = (double*)malloc(N * sizeof(double));
	double* q_I_o = (double*)malloc(N * sizeof(double));
	double* q_I_g = (double*)malloc(N * sizeof(double));
	double error, BHP_CALC;


	if (correlation == 1) {

		double Q_o, Q_g, Q_w; //std
		double Q_o_med, Q_g_med, Q_w_med; //std
		double q_o_med, q_g_med, q_w_med; //WB cond
		double q_I_l, q_l_med; //WB cond

		double U_sl, U_sg;
		double H_l_b, dpdz_ref;

		double Rs, Bo, mu_o, rho_o;
		double Bg, rho_g, mu_g;
		double Bw, rho_w, mu_w;
		double rho_l, mu_l;


		double f1, f0, PWB_1, a;
		error = 10;
		int cont = 0;
		int vira = 1;

		while (error > 1e-7 && cont < 51) {

			PWB[N - 1] = PWB_0;
			Q_w = 0.0l;
			Q_o = 0.0l;
			Q_g = 0.0l;

			for (int i = N - 2; i >= -1; i--) {

				Rs = calc_gas_Rs_p(yg, PWB[i + 1], Pb[i + 1], API, T);
				Bo = calc_oil_Bo_p(yo, yg, PWB[i + 1], Pb[i + 1], API, T);
				mu_o = calc_oil_mi_p(Rs, T, API, PWB[i + 1], Pb[i + 1]);
				rho_o = calc_oil_rho_p(Rs, Bo, rho_os, rho_gs);

				Bg = calc_gas_Bg_p(PWB[i + 1]);
				rho_g = calc_gas_rho_p(rho_gs, Bg);
				mu_g = cP_Pas(calc_gas_mig_p(parametrosmug, T / Tpc, PWB[i + 1] / ppc));

				rho_w = calc_water_rho_p(Bwi, rho_ws, T, PWB[i + 1], pi);
				mu_w = calc_water_mi_p(T, PWB[i + 1]);
				Bw = calc_water_Bw_p(PWB[i + 1], T, parametrosBw);

				q_I_w[i + 1] = ft3_m3(bbld_ft3s(Bw * Q_I_w[i + 1]));
				q_I_o[i + 1] = ft3_m3(bbld_ft3s(Bo * Q_I_o[i + 1]));
				q_I_g[i + 1] = ft3_m3(bbld_ft3s(Bg * (Q_I_g[i + 1] - Rs * Q_I_o[i + 1])));

				Q_w_med = Q_w + 0.5l * Q_I_w[i + 1] * deltaL[i + 1];
				Q_o_med = Q_o + 0.5l * Q_I_o[i + 1] * deltaL[i + 1];
				Q_g_med = Q_g + 0.5l * Q_I_g[i + 1] * deltaL[i + 1];

				q_w_med = ft3_m3(bbld_ft3s(Bw * Q_w_med));
				q_o_med = ft3_m3(bbld_ft3s(Bo * Q_o_med));
				q_g_med = ft3_m3(bbld_ft3s(Bg * (Q_g_med - Rs * Q_o_med)));

				q_I_l = q_I_w[i + 1] + q_I_o[i + 1];
				q_l_med = q_w_med + q_o_med;

				rho_l = rho_w * (q_w_med / q_l_med) + rho_o * (1.l - (q_w_med / q_l_med));
				mu_l = mu_w * (q_w_med / q_l_med) + mu_o * (1.l - (q_w_med / q_l_med));

				U_sl = 4.l * q_l_med / ft_m(ft_m(M_PI * D * D));
				U_sg = 4.l * q_g_med / ft_m(ft_m(M_PI * D * D));

				H_l_b = U_sl / (U_sl + U_sg);
				H_l[i + 1] = H_l_b;
				dpdz_ref = 10000000000000;
				dpdz[i + 1] = 10;

				int it = 0;

				while ((abs(H_l_b - H_l[i + 1]) > 1e-7 || abs(dpdz_ref - dpdz[i + 1]) / abs(dpdz_ref) > 1e-7) && it < 50) {


					if (it > 10) {
						H_l_b = 0.5l * (H_l[i + 1] + H_l_b);
						dpdz_ref = 0.5 * (dpdz[i + 1] + dpdz_ref);
					}
					else {
						H_l_b = H_l[i + 1];
						dpdz_ref = dpdz[i + 1];
					}


					ouyang_multiphase_model_node(&H_l[i + 1], &dpdz[i + 1], H_l_b, lbft3_kgm3(rho_g), lbft3_kgm3(rho_l),
						U_sg, U_sl, q_I_l, q_I_g[i + 1], cP_Pas(mu_l), cP_Pas(mu_g), ft_m(D), ft_m(g), teta, deltaL[i + 1],
						D_perf, perf_dens, e, gama);


					it++;
				}

				if (i == -1) {
					BHP_CALC = PWB[i + 1] + Pa_Psia(dpdz[i + 1]);//checar sinal
				}
				else {
					PWB[i] = PWB[i + 1] + Pa_Psia(dpdz[i + 1]);//checar sinal
				}

				Q_w = Q_w + Q_I_w[i + 1] * deltaL[i + 1];
				Q_o = Q_o + Q_I_o[i + 1] * deltaL[i + 1];
				Q_g = Q_g + Q_I_g[i + 1] * deltaL[i + 1];


			}

			error = abs(BHP_CALC - BHP) / BHP;

			//cout << "cont: " << cont << endl;
			//cout << "error: " << (BHP_CALC - BHP) / BHP << endl;

			if (cont == 0) {
				f1 = BHP_CALC - BHP;
				PWB_1 = PWB_0;
				if (BHP_CALC > BHP) {
					PWB_0 = PWB_0 - 0.5l;
				}
				else {
					PWB_0 = PWB_0 + 0.5l;
				}
			}
			else {
				f0 = BHP_CALC - BHP;
				a = (PWB_1 - PWB_0) / (f1 - f0);
				PWB_1 = PWB_0;
				f1 = f0;
				PWB_0 = PWB_0 - f0 / a;

			}

			cont++;

		}


	}

	if (correlation == 2) {

		double Q_o, Q_g, Q_w; //std
		double Q_o_med, Q_g_med, Q_w_med; //std
		double q_o_med, q_g_med, q_w_med; //WB cond
		double q_I_l, q_l_med; //WB cond

		double U_sl, U_sg;
		double H_l_b, dpdz_ref;

		double Rs, Bo, mu_o, rho_o;
		double Bg, rho_g, mu_g;
		double Bw, rho_w, mu_w;
		double rho_l, mu_l;


		double f1, f0, PWB_1, a;
		error = 10;
		int cont = 0;
		int vira = 1;

		while (error > 1e-7 && cont < 51) {

			PWB[N - 1] = PWB_0;
			Q_w = 0.0l;
			Q_o = 0.0l;
			Q_g = 0.0l;

			for (int i = N - 2; i >= -1; i--) {


				Rs = ft3_bbl(calc_gas_Rs_p(yg, PWB[i + 1], Pb[i + 1], API, T));
				Bo = calc_oil_Bo_p(yo, yg, PWB[i + 1], Pb[i + 1], API, T);
				mu_o = calc_oil_mi_p(Rs, T, API, PWB[i + 1], Pb[i + 1]);
				rho_o = calc_oil_rho_p(Rs, Bo, rho_os, rho_gs);

				Bg = bbl_ft3(calc_gas_Bg_p(PWB[i + 1]));
				rho_g = calc_gas_rho_p(rho_gs, Bg);
				mu_g = cP_Pas(calc_gas_mig_p(parametrosmug, T / Tpc, PWB[i + 1] / ppc));

				rho_w = calc_water_rho_p(Bwi, rho_ws, T, PWB[i + 1], pi);
				mu_w = calc_water_mi_p(T, PWB[i + 1]);
				Bw = calc_water_Bw_p(PWB[i + 1], T, parametrosBw);

				q_I_w[i + 1] = ft3_m3(bbld_ft3s(Bw * Q_I_w[i + 1]));
				q_I_o[i + 1] = ft3_m3(bbld_ft3s(Bo * Q_I_o[i + 1]));
				q_I_g[i + 1] = ft3_m3(bbld_ft3s(Bg * (Q_I_g[i + 1] - Rs * Q_I_o[i + 1])));

				Q_w_med = Q_w + 0.5l * Q_I_w[i + 1] * deltaL[i + 1];
				Q_o_med = Q_o + 0.5l * Q_I_o[i + 1] * deltaL[i + 1];
				Q_g_med = Q_g + 0.5l * Q_I_g[i + 1] * deltaL[i + 1];

				q_w_med = ft3_m3(bbld_ft3s(Bw * Q_w_med));
				q_o_med = ft3_m3(bbld_ft3s(Bo * Q_o_med));
				q_g_med = ft3_m3(bbld_ft3s(Bg * (Q_g_med - Rs * Q_o_med)));

				/*cout << "qg/qo: " << Q_g_med / Q_o_med << endl;
				cout << "rs: " << calc_gas_Rs_p(yg, PWB[i + 1], Pb[i + 1], API, T) << endl;
				cout << "owb: " <<  PWB[i + 1] << endl;
				cout << "pb: " << Pb[i + 1] << endl;
				cout << "yg: " << yg << endl;
				cout << "AO: " << API << endl;
				cout << "T: " << T << endl;*/

				q_I_l = q_I_w[i + 1] + q_I_o[i + 1];
				q_l_med = q_w_med + q_o_med;

				rho_l = rho_w * (q_w_med / q_l_med) + rho_o * (1.l - (q_w_med / q_l_med));
				mu_l = mu_w * (q_w_med / q_l_med) + mu_o * (1.l - (q_w_med / q_l_med));

				U_sl = 4.l * q_l_med / ft_m(ft_m(M_PI * D * D));
				U_sg = 4.l * q_g_med / ft_m(ft_m(M_PI * D * D));

				H_l_b = U_sl / (U_sl + U_sg);
				H_l[i + 1] = H_l_b;
				dpdz_ref = 10000000000000;
				dpdz[i + 1] = 10;

				int it = 0;

				while ((abs(H_l_b - H_l[i + 1]) > 1e-7 || abs(dpdz_ref - dpdz[i + 1]) / abs(dpdz_ref) > 1e-7) && it < 50) {


					if (it > 10) {
						H_l_b = 0.5l * (H_l[i + 1] + H_l_b);
						dpdz_ref = 0.5 * (dpdz[i + 1] + dpdz_ref);
					}
					else {
						H_l_b = H_l[i + 1];
						dpdz_ref = dpdz[i + 1];
					}
					// cout << "gft:" << g << endl;
					ouyang_homogeneous(&H_l[i + 1], &dpdz[i + 1], U_sl, U_sg, cP_Pas(mu_l), cP_Pas(mu_g),
						lbft3_kgm3(rho_l), lbft3_kgm3(rho_g), q_I_l, q_I_g[i + 1], teta, ft_m(g),
						ft_m(D), deltaL[i + 1], e, Psia_Pa(PWB[i + 1]));

					it++;
				}

				if (i == -1) {
					BHP_CALC = PWB[i + 1] + Pa_Psia(dpdz[i + 1]);//checar sinal
				}
				else {
					PWB[i] = PWB[i + 1] + Pa_Psia(dpdz[i + 1]);//checar sinal
					if (isnan(PWB[i])) {
						cout << "PWB[" << i << "]: " << PWB[i] << endl;
						cout << "H_l_b: " << H_l_b << endl;
						system("pause");
					}
				}

				Q_w = Q_w + Q_I_w[i + 1] * deltaL[i + 1];
				Q_o = Q_o + Q_I_o[i + 1] * deltaL[i + 1];
				Q_g = Q_g + Q_I_g[i + 1] * deltaL[i + 1];


			}


			error = abs(BHP_CALC - BHP) / BHP;


			if (cont == 0) {
				f1 = BHP_CALC - BHP;
				PWB_1 = PWB_0;
				if (BHP_CALC > BHP) {
					PWB_0 = PWB_0 - 0.5l;
				}
				else {
					PWB_0 = PWB_0 + 0.5l;
				}
			}
			else {
				f0 = BHP_CALC - BHP;
				a = (PWB_1 - PWB_0) / (f1 - f0);
				PWB_1 = PWB_0;
				f1 = f0;
				PWB_0 = PWB_0 - f0 / a;

			}

			cont++;

		}


	}

	if (correlation == 3) {


		double Q_o, Q_g, Q_w; //std
		double Q_o_med, Q_g_med, Q_w_med; //std
		double q_o_med, q_g_med, q_w_med; //WB cond
		double q_I_l, q_l_med; //WB cond

		double U_sl, U_sg;
		double H_l_b, dpdz_ref;

		double Rs, Bo, mu_o, rho_o;
		double Bg, rho_g, mu_g;
		double Bw, rho_w, mu_w;
		double rho_l, mu_l;

		double Rs_2, Bo_2, mu_o_2, rho_o_2;
		double Bg_2, rho_g_2, mu_g_2;
		double Bw_2, rho_w_2, mu_w_2;
		double rho_l_2, mu_l_2;


		double f1, f0, PWB_1, a;
		error = 10;
		double cont = 0;

		while (error > 1e-7 && cont < 51) {

			PWB[N - 1] = PWB_0;
			Q_w = 0.0l;
			Q_o = 0.0l;
			Q_g = 0.0l;

			for (int i = N - 2; i >= -1; i--) {

				Rs = calc_gas_Rs_p(yg, PWB[i + 1], Pb[i + 1], API, T);
				Bo = calc_oil_Bo_p(yo, yg, PWB[i + 1], Pb[i + 1], API, T);
				mu_o = calc_oil_mi_p(Rs, T, API, PWB[i + 1], Pb[i + 1]);
				rho_o = calc_oil_rho_p(Rs, Bo, rho_os, rho_gs);

				Bg = calc_gas_Bg_p(PWB[i + 1]);
				rho_g = calc_gas_rho_p(rho_gs, Bg);
				mu_g = cP_Pas(calc_gas_mig_p(parametrosmug, T / Tpc, PWB[i + 1] / ppc));

				rho_w = calc_water_rho_p(Bwi, rho_ws, T, PWB[i + 1], pi);
				mu_w = calc_water_mi_p(T, PWB[i + 1]);
				Bw = calc_water_Bw_p(PWB[i + 1], T, parametrosBw);

				q_I_w[i + 1] = ft3_m3(bbld_ft3s(Bw * Q_I_w[i + 1]));
				q_I_o[i + 1] = ft3_m3(bbld_ft3s(Bo * Q_I_o[i + 1]));
				q_I_g[i + 1] = ft3_m3(bbld_ft3s(Bg * (Q_I_g[i + 1] - Rs * Q_I_o[i + 1])));

				Q_w_med = Q_w + 0.5l * Q_I_w[i + 1] * deltaL[i + 1];
				Q_o_med = Q_o + 0.5l * Q_I_o[i + 1] * deltaL[i + 1];
				Q_g_med = Q_g + 0.5l * Q_I_g[i + 1] * deltaL[i + 1];

				q_w_med = ft3_m3(bbld_ft3s(Bw * Q_w_med));
				q_o_med = ft3_m3(bbld_ft3s(Bo * Q_o_med));
				q_g_med = ft3_m3(bbld_ft3s(Bg * (Q_g_med - Rs * Q_o_med)));

				q_I_l = q_I_w[i + 1] + q_I_o[i + 1];
				q_l_med = q_w_med + q_o_med;

				rho_l = rho_w * (q_w_med / q_l_med) + rho_o * (1.l - (q_w_med / q_l_med));
				mu_l = mu_w * (q_w_med / q_l_med) + mu_o * (1.l - (q_w_med / q_l_med));

				U_sl = 4.l * q_l_med / ft_m(ft_m(M_PI * D * D));
				U_sg = 4.l * q_g_med / ft_m(ft_m(M_PI * D * D));

				H_l_b = U_sl / (U_sl + U_sg);


				Rs_2 = calc_gas_Rs_p(yg, PWB[i + 1], Pb[i + 1], API, T);
				Bo_2 = calc_oil_Bo_p(yo, yg, PWB[i + 1], Pb[i + 1], API, T);
				mu_o_2 = calc_oil_mi_p(Rs, T, API, PWB[i + 1], Pb[i + 1]);
				rho_o_2 = calc_oil_rho_p(Rs, Bo, rho_os, rho_gs);

				Bg_2 = calc_gas_Bg_p(PWB[i + 1]);
				rho_g_2 = calc_gas_rho_p(rho_gs, Bg);
				mu_g_2 = cP_Pas(calc_gas_mig_p(parametrosmug, T / Tpc, PWB[i + 1] / ppc));

				rho_w_2 = calc_water_rho_p(Bwi, rho_ws, T, PWB[i + 1], pi);
				mu_w_2 = calc_water_mi_p(T, PWB[i + 1]);
				Bw_2 = calc_water_Bw_p(PWB[i + 1], T, parametrosBw);

				rho_l_2 = rho_w_2 * (q_w_med / q_l_med) + rho_o_2 * (1.l - (q_w_med / q_l_med));
				H_l[i + 1] = H_l_b;
				dpdz_ref = 10000000000000;
				dpdz[i + 1] = 10;

				while (abs(H_l_b - H_l[i + 1]) > 1e-7 || abs(dpdz_ref - dpdz[i + 1]) / abs(dpdz_ref) > 1e-7) {

					//cout << "i:" << i << endl;
					H_l_b = H_l[i + 1];
					dpdz_ref = dpdz[i + 1];

					beb_calc(&H_l[i + 1], &dpdz[i + 1], U_sl, U_sg, cP_Pas(mu_l), cP_Pas(mu_g),
						lbft3_kgm3(rho_l), lbft3_kgm3(rho_g), q_I_l, q_I_g[i + 1], teta, ft_m(g),
						ft_m(D), (deltaL[i + 1]), e, lbft3_kgm3(rho_l_2), lbft3_kgm3(rho_g_2));

					Rs_2 = calc_gas_Rs_p(yg, PWB[i + 1] + Pa_Psia(dpdz[i + 1]), Pb[i + 1], API, T);
					Bo_2 = calc_oil_Bo_p(yo, yg, PWB[i + 1] + Pa_Psia(dpdz[i + 1]), Pb[i + 1], API, T);
					mu_o_2 = calc_oil_mi_p(Rs, T, API, PWB[i + 1] + Pa_Psia(dpdz[i + 1]), Pb[i + 1]);
					rho_o_2 = calc_oil_rho_p(Rs, Bo, rho_os, rho_gs);

					Bg_2 = calc_gas_Bg_p(PWB[i + 1] + Pa_Psia(dpdz[i + 1]));
					rho_g_2 = calc_gas_rho_p(rho_gs, Bg);
					mu_g_2 = cP_Pas(calc_gas_mig_p(parametrosmug, T / Tpc, PWB[i + 1] + Pa_Psia(dpdz[i + 1]) / ppc));

					rho_w_2 = calc_water_rho_p(Bwi, rho_ws, T, PWB[i + 1] + Pa_Psia(dpdz[i + 1]), pi);
					mu_w_2 = calc_water_mi_p(T, PWB[i + 1] + Pa_Psia(dpdz[i + 1]));
					Bw_2 = calc_water_Bw_p(PWB[i + 1] + Pa_Psia(dpdz[i + 1]), T, parametrosBw);

				}

				if (i == -1) {
					BHP_CALC = PWB[i + 1] + Pa_Psia(dpdz[i + 1]);//checar sinal
				}
				else {
					PWB[i] = PWB[i + 1] + Pa_Psia(dpdz[i + 1]);//checar sinal
				}

				Q_w = Q_w + Q_I_w[i + 1] * deltaL[i + 1];
				Q_o = Q_o + Q_I_o[i + 1] * deltaL[i + 1];
				Q_g = Q_g + Q_I_g[i + 1] * deltaL[i + 1];


			}

			//cout << "BHP_CALC:" << BHP_CALC << endl;
			error = abs(BHP_CALC - BHP) / BHP;

			if (cont == 0) {
				f1 = BHP_CALC - BHP;
				PWB_1 = PWB_0;
				if (BHP_CALC > BHP) {
					PWB_0 = PWB_0 - 0.5l;
				}
				else {
					PWB_0 = PWB_0 + 0.5l;
				}
			}
			else {
				f0 = BHP_CALC - BHP;
				a = (PWB_1 - PWB_0) / (f1 - f0);
				PWB_1 = PWB_0;
				f1 = f0;
				PWB_0 = PWB_0 - f0 / a;
			}

			cont++;


		}
	}

	if (correlation == 4) {

		double Q_o, Q_g, Q_w; //std
		double Q_o_med, Q_g_med, Q_w_med; //std
		double q_o_med, q_g_med, q_w_med; //WB cond
		double q_I_l, q_l_med; //WB cond

		double U_sl, U_sg;
		double H_l_b, dpdz_ref;

		double Rs, Bo, mu_o, rho_o;
		double Bg, rho_g, mu_g;
		double Bw, rho_w, mu_w;
		double rho_l, mu_l;


		double f1, f0, PWB_1, a;
		error = 10;
		double cont = 0;

		while (error > 1e-7 && cont < 51) {

			PWB[N - 1] = PWB_0;
			Q_w = 0.0l;
			Q_o = 0.0l;
			Q_g = 0.0l;

			for (int i = N - 2; i >= -1; i--) {

				Rs = ft3_bbl(calc_gas_Rs_p(yg, PWB[i + 1], Pb[i + 1], API, T));
				Bo = calc_oil_Bo_p(yo, yg, PWB[i + 1], Pb[i + 1], API, T);
				mu_o = calc_oil_mi_p(Rs, T, API, PWB[i + 1], Pb[i + 1]);
				rho_o = calc_oil_rho_p(Rs, Bo, rho_os, rho_gs);

				Bg = bbl_ft3(calc_gas_Bg_p(PWB[i + 1]));
				rho_g = calc_gas_rho_p(rho_gs, Bg);
				mu_g = cP_Pas(calc_gas_mig_p(parametrosmug, T / Tpc, PWB[i + 1] / ppc));

				rho_w = calc_water_rho_p(Bwi, rho_ws, T, PWB[i + 1], pi);
				mu_w = calc_water_mi_p(T, PWB[i + 1]);
				Bw = calc_water_Bw_p(PWB[i + 1], T, parametrosBw);

				q_I_w[i + 1] = ft3_m3(bbld_ft3s(Bw * Q_I_w[i + 1]));
				q_I_o[i + 1] = ft3_m3(bbld_ft3s(Bo * Q_I_o[i + 1]));
				q_I_g[i + 1] = ft3_m3(bbld_ft3s(Bg * (Q_I_g[i + 1] - Rs * Q_I_o[i + 1])));

				Q_w_med = Q_w + 0.5l * Q_I_w[i + 1] * deltaL[i + 1];
				Q_o_med = Q_o + 0.5l * Q_I_o[i + 1] * deltaL[i + 1];
				Q_g_med = Q_g + 0.5l * Q_I_g[i + 1] * deltaL[i + 1];

				q_w_med = ft3_m3(bbld_ft3s(Bw * Q_w_med));
				q_o_med = ft3_m3(bbld_ft3s(Bo * Q_o_med));
				q_g_med = ft3_m3(bbld_ft3s(Bg * (Q_g_med - Rs * Q_o_med)));

				q_I_l = q_I_w[i + 1] + q_I_o[i + 1];
				q_l_med = q_w_med + q_o_med;

				rho_l = rho_w * (q_w_med / q_l_med) + rho_o * (1.l - (q_w_med / q_l_med));
				mu_l = mu_w * (q_w_med / q_l_med) + mu_o * (1.l - (q_w_med / q_l_med));

				U_sl = 4.l * q_l_med / ft_m(ft_m(M_PI * D * D));
				U_sg = 4.l * q_g_med / ft_m(ft_m(M_PI * D * D));

				H_l_b = U_sl / (U_sl + U_sg);
				H_l[i + 1] = H_l_b;
				dpdz_ref = 10000000000000;
				dpdz[i + 1] = 10;

				while (abs(H_l_b - H_l[i + 1]) > 1e-7 || abs(dpdz_ref - dpdz[i + 1]) / abs(dpdz_ref) > 1e-7) {

					//cout << "gft:" << g << endl;
					H_l_b = H_l[i + 1];
					dpdz_ref = dpdz[i + 1];

					beb_mod_calc(&H_l[i + 1], &dpdz[i + 1], U_sl, U_sg, cP_Pas(mu_l), cP_Pas(mu_g),
						lbft3_kgm3(rho_l), lbft3_kgm3(rho_g), q_I_l, q_I_g[i + 1], teta, ft_m(g),
						ft_m(D), (deltaL[i + 1]), e);

				}

				if (i == -1) {
					BHP_CALC = PWB[i + 1] + Pa_Psia(dpdz[i + 1]);//checar sinal
				}
				else {
					PWB[i] = PWB[i + 1] + Pa_Psia(dpdz[i + 1]);//checar sinal
				}

				Q_w = Q_w + Q_I_w[i + 1] * deltaL[i + 1];
				Q_o = Q_o + Q_I_o[i + 1] * deltaL[i + 1];
				Q_g = Q_g + Q_I_g[i + 1] * deltaL[i + 1];


			}

			error = abs(BHP_CALC - BHP) / BHP;

			if (cont == 0) {
				f1 = BHP_CALC - BHP;
				PWB_1 = PWB_0;
				if (BHP_CALC > BHP) {
					PWB_0 = PWB_0 - 0.5l;
				}
				else {
					PWB_0 = PWB_0 + 0.5l;
				}
			}
			else {
				f0 = BHP_CALC - BHP;
				a = (PWB_1 - PWB_0) / (f1 - f0);
				PWB_1 = PWB_0;
				f1 = f0;
				PWB_0 = PWB_0 - f0 / a;
			}

			cont++;

		}


	}

	free(dpdz);
	free(q_I_w);
	free(q_I_o);
	free(q_I_g);


}

void ouyang_homogeneous(double* H_l, double* dpdx, double U_sl, double U_sg, double mu_l, double mu_g,
	double rho_l, double rho_g, double  q_I_l, double q_I_g, double teta, double g,
	double D, double deltaL, double e, double p) {

	//SLIP CONSTANTS
	double sigma = Tensuper(rho_l, rho_g);
	double c_w = 1.2l;
	double c_0 = c_w - 0.2l * sqrt(rho_g / rho_l);
	double U_0 = 1.53l * sqrt(sqrt(g * (rho_l - rho_g) * sigma) / rho_l) * sin(teta);
	//	

	double U_m = U_sg + U_sl;

	*H_l = 1.0l - U_sg / (c_0 * U_m + U_0);

	double mu_tp = mu_l * *H_l + mu_g * (1.0l - *H_l);
	double rho_tp = rho_l * *H_l + rho_g * (1.0l - *H_l);
	double U_tp = (rho_l * U_sl + rho_g * U_sg) / rho_tp;
	double q_I_tp = (rho_l * q_I_l + rho_g * q_I_g) / rho_tp;
	double q_I_m = q_I_l + q_I_g;


	double Re_tp = rho_tp * U_tp * D / mu_tp;
	double Re_w = calc_Re_w(q_I_l, q_I_g, rho_l, rho_g, mu_l, mu_g);

	double f_tp = cal_fric_fact_wellbore(Re_tp, Re_w, D, D * e);
	double tau_w = 0.5l * f_tp * rho_tp * U_tp * U_tp;


	double dpdxf = -(4.0l * tau_w / D);
	double dpdxg = -rho_tp * g * sin(teta);
	double dpdxaw1 = -4.0l * rho_tp * (U_m * q_I_tp + U_tp * q_I_m) / (M_PI * D * D);
	double dpdxaw2 = -8.0l * rho_tp * U_tp * q_I_tp / (M_PI * D * D);

	double omega = 0.8l;
	double dpdxaw = omega * dpdxaw1 + (1.l - omega) * dpdxaw2;

	double Betae = rho_tp * U_m * U_sg / p;
	double dpdxae = (Betae / (1.l - Betae)) * (dpdxf + dpdxg + dpdxaw);

	/*cout << "Atrito: " << dpdxf << endl;
	cout << "Acelera: " << dpdxaw + dpdxae << endl;
	cout << "Grave: " << dpdxg << endl;
	cout << "H_l: " << *H_l << endl;
	cout << "rho_l: " << rho_l << endl;
	cout << "rho_g: " << rho_g << endl;
	cout << "g: " << g << endl;
	system("pause");*/
	//cout << "rho_l: " << rho_l << endl;
	/*cout << "fator: " << (Betae / (1.l - Betae)) << endl;
	system("pause");*/

	* dpdx = (dpdxf + dpdxg + dpdxaw + dpdxae) * deltaL;

}

void ouyang_multiphase_model_node(double* H_l, double* dpdx, double H_l_b, double rho_g, double rho_l,
	double U_sg, double U_sl, double q_I_l, double q_I_g, double mu_l, double mu_g,
	double D, double g, double teta, double deltax,
	double D_perf, double perf_dens, double e, double gama) {


	double V_im = (q_I_l + q_I_g) / (M_PI * D);
	double U_m = U_sg + U_sl;
	double A = 0.25l * M_PI * D * D;
	double sigma = Tensuper(rho_l, rho_g);

	int pattern = flow_pattern_ouyang(H_l_b, mu_g, mu_l, rho_l, rho_g, U_m, g, teta, V_im, D, U_sg, e);

	//SINGLE-PHASE FLOW
	if (pattern == 0) {

		double Q_m = A * U_m;
		*dpdx = calc_pres_drop_ouyang_mono(rho_l, mu_l, D, D_perf, perf_dens, teta, g, e, gama, deltax, Q_m, U_m, q_I_l + q_I_g);
		*H_l = H_l_b;

	}

	//BUBBLE FLOW
	else if (pattern == 1) {

		bubble_flow_slip_calc(H_l, dpdx, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, q_I_l, q_I_g, sigma, teta, g, D, deltax, e);
		//cout << "bolinha" << endl;

	}

	//STRATIFIED FLOW
	else if (pattern == 2) {

		stratified_flow_calc(H_l, dpdx, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, q_I_l, q_I_g, sigma, teta, g, D, A, deltax, e);
		//cout << "extratinho" << endl;
	}

	//ANNULAR-MIST FLOW
	else if (pattern == 3) {

		annularmist_flow_calc(H_l, dpdx, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, q_I_l, q_I_g, sigma, teta, g, D, A, deltax, e);
		//cout << "anelzinho" << endl;

	}

	//INTERMITTENT FLOW
	else {

		intermittent_flow_calc(H_l, dpdx, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, q_I_l, q_I_g, sigma, teta, g, D, A, deltax, e);
		//cout << "slugzinho" << endl;

	}


}

double calc_hlD(double H_l) {

	double h_lD;

	if ((0.0l < H_l && H_l <= 0.164l) || (0.228l < H_l && H_l <= 0.337l)) {
		h_lD = 0.7769409734 * pow(H_l, 0.6875921155);
	}
	else if ((0.164l < H_l && H_l <= 0.228l) || (0.385l < H_l && H_l < 0.500l)) {
		h_lD = -0.6896121749 * H_l * H_l + 1.2645448768 * H_l + 0.0286615347;
	}
	else if (0.337l < H_l && H_l <= 0.385l) {
		h_lD = 0.9441422481 * H_l + 0.0496745234;
	}
	else {
		h_lD = H_l;
	}

	double teta, S_i;
	double H_l_ref = 0.5l;

	teta = 2.0l * acos(1.0l - 2 * h_lD);
	S_i = sqrt(1.l - (2.l * h_lD - 1.l) * (2.l * h_lD - 1.l));
	H_l_ref = (M_PI - acos(2.l * h_lD - 1.l) + (2.l * h_lD - 1.l) * S_i) / M_PI;
	double h_lD_2, H_l_2, S_i_2;
	h_lD_2 = 0.01l * (1.l - h_lD) + h_lD;
	S_i_2 = sqrt(1.l - (2.l * h_lD_2 - 1.l) * (2.l * h_lD_2 - 1.l));
	H_l_2 = (M_PI - acos(2.l * h_lD_2 - 1.l) + (2.l * h_lD_2 - 1.l) * S_i_2) / M_PI;

	while (abs((H_l_ref - H_l) / H_l) > 1e-7) {

		teta = max(0.0l, min(2 * M_PI, teta - 4.l * (H_l_ref - H_l) / (sin(0.5 * teta) * (H_l_ref - H_l_2) / (h_lD - h_lD_2))));
		//teta=2.0l*M_PI*H_l+sin(teta);	
		h_lD = 0.5 * (1.0l - cos(0.5l * teta));
		S_i = sqrt(1.l - (2.l * h_lD - 1.l) * (2.l * h_lD - 1.l));
		H_l_ref = (M_PI - acos(2.l * h_lD - 1.l) + (2.l * h_lD - 1.l) * S_i) / M_PI;

		h_lD_2 = 0.01l * (1.l - h_lD) + h_lD;
		S_i_2 = sqrt(1.l - (2.l * h_lD_2 - 1.l) * (2.l * h_lD_2 - 1.l));
		H_l_2 = (M_PI - acos(2.l * h_lD_2 - 1.l) + (2.l * h_lD_2 - 1.l) * S_i_2) / M_PI;

	}


	return h_lD;
}

double calc_A_G(double h_lD, double D) {
	double S_i = sqrt(1.0l - pow(2.0l * h_lD - 1.0l, 2.0l));
	return 0.25l * D * D * (acos(2.l * h_lD - 1.l) - (2.l * h_lD - 1.l) * S_i);//A-A_L 
}

double calc_d_max_hinze(double U_m, double f_m, double D, double sigma, double rho_l) {

	double k = 1.14l;
	double kappa = 2.0l * f_m * U_m * U_m * U_m / D;

	return k * pow(sigma / rho_l, 0.6l) * pow(kappa, -0.4l);


}

double calc_C_d_barry(double mu_g, double mu_l, double rho_g, double rho_l, double V_im, double d) {

	double A = 1.10535l;
	double alpha = 2.5891;
	double beta = 0.9879;
	double X = mu_g / mu_l;
	double P = rho_g / rho_l;

	double Z = 1.0l + ((alpha - 2.0l) * sqrt(X * P) + (beta - 1.0l) * X * P) / (1 + 2 * sqrt(X * P) + X * P);

	double AZ = A * Z;
	double C_1 = (8.l / 9.l) * (4.l + 3.l * X) / (1.l + X);
	double tau_i = sqrt(C_1 - 1.4 * AZ);//1.4 is OK for the first step, not good, not bad, just OK
	double tau = 100000000000;

	while (abs(tau_i - tau) > 0.000000001l) {

		tau = tau_i;
		tau_i = sqrt(C_1 - tau * AZ);

	}
	tau = tau_i;

	double lambda = 1.l / tau;
	double omega = 1.l / (3.l * (1.l + X));

	double R_0 = rho_l * V_im * d / mu_l;

	return 48.l / (R_0 * (1.l + X)) * (1.l + 1.5l * X) * (sqrt(R_0) + AZ + tau * exp(-lambda * sqrt(R_0))) / (R_0 / (1.l + X) + 3.l * AZ + 3.l * tau * exp(-omega * sqrt(R_0)));

}

bool is_dispersed_bubble_flow(double H_l, double mu_g, double mu_l, double rho_l, double rho_g,
	double U_m, double g, double teta, double V_im, double D, double e) {

	if (H_l < 0.48l) {
		return 0;
	}

	double sigma = Tensuper(rho_l, rho_g);
	double f_m = cal_fric_fact_turb_colle(rho_l * U_m * D / mu_l, D, e);
	double d = calc_d_max_hinze(U_m, f_m, D, sigma, rho_l);

	//Ouyang 2002
	double C_id = 0.8l;
	double lambda = mu_g / mu_l;
	double A = (3.0l / 8.0l) * (rho_l / (rho_l - rho_g)) * U_m * U_m * f_m / (g * cos(teta));
	double B = (3.0l * lambda + 2.0l) / (lambda + 1.0l) * 6.0l * C_id * mu_l * V_im / ((rho_l - rho_g) * g * cos(teta));

	double d_cb = max(0.5l * (A + pow(A * A + 4.0l * B, 0.5l)), 0.5l * (A - pow(A * A + 4.0l * B, 0.5l)));

	if (d > d_cb) {
		return 0;
	}
	//


	////Alves 2020
	////double 
	//C_id=1.0l;
	//
	////double
	//d_cb=(3.0/8.0l)*(rho_l/(rho_l-rho_g))/(g*cos(teta))*(f_m*U_m*U_m+2.0l*C_id*V_im*V_im*calc_C_d_barry(mu_g, mu_l, rho_g, rho_l, V_im, d));/////CHECAR

	//if(d>d_cb){
	//	cout << "PAU" << endl;
	//	return 0;
	//}
	////


	double d_cd = 2.0l * pow((0.4l * sigma) / ((rho_l - rho_g) * g), 0.5l);

	if (d > d_cd) {
		return 0;
	}

	return 1;

}

bool is_stratified_flow(double H_l, double D, double rho_l, double rho_g, double g, double teta, double V_im, double U_sg) {

	double hlD = calc_hlD(H_l);
	double hl = hlD * D;
	double A_g = calc_A_G(hlD, D);
	double dAldhl = 2.0l * sqrt(hl * (D - hl));
	double U_g = U_sg * 0.25l * M_PI * D * D / A_g;

	double U_g_ref = (1.0l - hlD) * sqrt((A_g / dAldhl) * (rho_l * g * cos(teta) / (rho_g)+0.5l * V_im * abs(V_im) / (D - hl)));
	/*cout << "U_g:" << U_g << endl;
	cout << "U_g_ref:" << U_g_ref << endl;
	cout << "(1.0l-hlD):" << (1.0l - hlD) << endl;
	system("pause");*/

	if (U_g > U_g_ref) {
		/*cout << "U_g:" << U_g << endl;
		cout << "U_g_ref:" << U_g_ref << endl;*/
		return 0;
	}

	return 1;


}

int flow_pattern_ouyang(double H_l, double mu_g, double mu_l, double rho_l, double rho_g,
	double U_m, double g, double teta, double V_im, double D, double U_sg, double e) {

	if (H_l > 0.9999999999999l) {
		return 0;//SINGLE-PHASE FLOW	
	}

	if (is_dispersed_bubble_flow(H_l, mu_g, mu_l, rho_l, rho_g, U_m, g, teta, V_im, D, e)) {
		return 1;//BUBBLE FLOW
	}

	if (is_stratified_flow(H_l, D, rho_l, rho_g, g, teta, V_im, U_sg)) {
		return 2;//STRATIFIED FLOW
	}

	if (H_l < 0.24l) {
		return 3;//ANNULAR-MIST FLOW
	}

	return 4;//INTERMITTENT FLOW

}

double calc_Re_w(double q_I_l, double q_I_g, double rho_l, double rho_g, double mu_l, double mu_g) {

	double C_I_l = q_I_l / (q_I_l + q_I_g);
	double rho_im = rho_l * C_I_l + rho_g * (1.0l - C_I_l);
	double mu_im = mu_l * C_I_l + mu_g * (1.0l - C_I_l);

	return rho_im * (q_I_l + q_I_g) / (M_PI * mu_im);


}

void bubble_flow_slip_calc(double* H_l, double* dpdx, double U_sl, double U_sg, double mu_l, double mu_g,
	double rho_l, double rho_g, double  q_I_l, double q_I_g, double sigma, double teta, double g,
	double D, double deltaL, double e) {

	//SLIP CONSTANTS	
	double c_w = 1.2l;
	double c_0 = c_w - 0.2l * sqrt(rho_g / rho_l);
	double U_0 = 1.53l * sqrt(sqrt(g * (rho_l - rho_g) * sigma) / rho_l) * sin(teta);
	//	

	double U_m = U_sg + U_sl;

	*H_l = 1.0l - U_sg / (c_0 * U_m + U_0);

	double mu_tp = mu_l * *H_l + mu_g * (1.0l - *H_l);
	double rho_tp = rho_l * *H_l + rho_g * (1.0l - *H_l);
	double U_tp = (rho_l * U_sl + rho_g * U_sg) / rho_tp;
	double q_I_tp = (rho_l * q_I_l + rho_g * q_I_g) / rho_tp;

	double Re_tp = rho_tp * U_tp * D / mu_tp;
	double Re_w = calc_Re_w(q_I_l, q_I_g, rho_l, rho_g, mu_l, mu_g);

	double f_tp = cal_fric_fact_wellbore(Re_tp, Re_w, D, D * e);

	double tau_w = 0.5l * f_tp * rho_tp * U_tp * U_tp;

	/*cout << "Bolinha" << endl;
	cout << "Atrito: " << -4.0l*tau_w*deltaL / D << endl;
	cout << "Acelera: " << -8.0l*rho_tp*U_tp*q_I_tp*deltaL / (M_PI*D*D) << endl;
	cout << "Grave isso: " << -rho_tp*g*sin(teta)*deltaL << endl;
	system("pause");*/

	*dpdx = -(4.0l * tau_w / D + 8.0l * rho_tp * U_tp * q_I_tp / (M_PI * D * D) + rho_tp * g * sin(teta)) * deltaL;

}

void bubble_flow_noslip_calc(double* H_l, double* dpdx, double U_sl, double U_sg, double mu_l, double mu_g,
	double rho_l, double rho_g, double  q_I_l, double q_I_g, double sigma, double teta, double g,
	double D, double A, double deltaL, double e) {


	double U_m = U_sg + U_sl;

	*H_l = U_sl / U_m;

	double mu_tp = mu_l * *H_l + mu_g * (1.0l - *H_l);
	double rho_tp = rho_l * *H_l + rho_g * (1.0l - *H_l);
	double U_tp = (rho_l * U_sl + rho_g * U_sg) / rho_tp;
	double q_I_tp = (rho_l * q_I_l + rho_g * q_I_g) / rho_tp;

	double Re_w = calc_Re_w(q_I_l, q_I_g, rho_l, rho_g, mu_l, mu_g);
	double Re_tp = rho_tp * U_tp * D / mu_tp;

	double f_tp = cal_fric_fact_wellbore(Re_tp, Re_w, D, D * e);

	double tau_w = 0.5l * f_tp * rho_tp * U_tp * U_tp;

	*dpdx = -(4.0l * tau_w / D + 8.0l * rho_tp * U_tp * q_I_tp / (M_PI * D * D) + rho_tp * g * sin(teta)) * deltaL;

}

void intermittent_flow_calc(double* H_l, double* dpdx, double U_sl, double U_sg, double mu_l, double mu_g,
	double rho_l, double rho_g, double  q_I_l, double q_I_g, double sigma, double teta, double g,
	double D, double A, double deltaL, double e) {

	//SLIP CONSTANTS	
	double c_w = 1.2l;
	double c_0 = c_w - 0.2l * sqrt(rho_g / rho_l);
	//	

	double U_m = U_sg + U_sl;

	//velocity of dispersed bubbles in liquid slug
	double U_b = 1.53l * sqrt(sqrt(g * (rho_l - rho_g) * sigma) / rho_l) * sin(teta);
	double U_gdb = c_0 * U_m + U_b;
	//

	//translational velocity of the elongated bubbles
	double Bo = (rho_l - rho_g) * g * D * D / sigma;//BOND, NUMBER OF BOND
	double U_DH_inf = (0.54l - 1.76l * pow(Bo, -0.56l)) * sqrt(g * D * (rho_l - rho_g) / rho_l);
	double U_DV_inf = 0.345l * (1.l - exp(0.337l - 0.1l * Bo)) * sqrt(g * D * (rho_l - rho_g) / rho_l);
	double Re_inf = 0.5l * rho_l * (U_DH_inf * cos(teta) + U_DV_inf * sin(teta)) * D / mu_l;
	double f_m = min(1.0l, 0.316l * sqrt(Re_inf));
	double U_d = f_m * (U_DH_inf * cos(teta) + U_DV_inf * sin(teta));
	double U_t = c_0 * U_m + U_d;
	//

	//liquid fraction in the liquid slug
	double H_ls = 1.0l / (1.0l + pow(U_m / 8.66l, 1.39l));

	*H_l = H_ls + (U_gdb * (1.0l - H_ls) - U_sg) / U_t;

	double mu_tp = mu_l * *H_l + mu_g * (1.0l - *H_l);
	double rho_tp = rho_l * *H_l + rho_g * (1.0l - *H_l);
	double U_tp = (rho_l * U_sl + rho_g * U_sg) / rho_tp;
	double q_I_tp = (rho_l * q_I_l + rho_g * q_I_g) / rho_tp;

	double Re_tp = rho_tp * U_tp * D / mu_tp;
	double Re_w = calc_Re_w(q_I_l, q_I_g, rho_l, rho_g, mu_l, mu_g);

	double f_tp = cal_fric_fact_wellbore(Re_tp, Re_w, D, D * e);

	//cout << "Intermitente" << endl;
	//cout << "Atrito: " << -2.0l*f_tp*rho_tp*U_tp*U_tp / D << endl;
	//cout << "Acelera: " << -2.0l*rho_tp*U_tp*q_I_tp*deltaL / A << endl;
	//cout << "Grave isso: " << rho_tp*g*sin(teta)*deltaL << endl;
	//system("pause");
	//
	*dpdx = (-2.0l * f_tp * rho_tp * U_tp * U_tp / D - 2.0l * rho_tp * U_tp * q_I_tp / A - rho_tp * g * sin(teta)) * deltaL;


}

double calc_fi_stratified(double Fr_l, double Re_sl, double rho_l, double rho_g, double g, double D, double U_g) {

	return (0.004l + 0.5l * Re_sl * 10e-6) * pow(Fr_l, 1.335l) * (rho_l * g * D / (rho_g * U_g * U_g));

}

void stratified_flow_calc(double* H_l, double* dpdx, double U_sl, double U_sg, double mu_l, double mu_g,
	double rho_l, double rho_g, double  q_I_l, double q_I_g, double sigma, double teta, double g,
	double D, double A, double deltaL, double e) {


	double Re_sl = rho_l * U_sl * D / mu_l;
	double Re_sg = rho_g * U_sg * D / mu_g;
	double Re_w = calc_Re_w(q_I_l, q_I_g, rho_l, rho_g, mu_l, mu_g);

	double f_sl = cal_fric_fact_wellbore(Re_sl, Re_w, D, D * e);
	double f_sg = cal_fric_fact_wellbore(Re_sg, Re_w, D, D * e);

	double dpdx_sg = -2.0l * f_sg * U_sg * U_sg * rho_g / D;

	double X2 = f_sl * U_sl * U_sl * rho_l / (f_sg * U_sg * U_sg * rho_g);
	double Y = (rho_l - rho_g) * g * sin(teta) / dpdx_sg;

	double h_ld = riddermethodforlockhartmartinelli(0.99999l, 0.000001l, 1, 8, 1, D, A, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g,
		Re_sl, Re_sg, Re_w, q_I_l, q_I_g, f_sl, f_sg, dpdx_sg, sigma, g, X2, Y, e);

	double S_i = D * sqrt(1.0l - pow(2.0l * h_ld - 1.0l, 2.0l));
	*H_l = (M_PI - acos(2.l * h_ld - 1.l) + (2.l * h_ld - 1.l) * S_i / D) / M_PI;

	double A_g = calc_A_G(h_ld, D);

	double S_l = D * (M_PI - acos(2.0l * h_ld - 1.0l));
	double S_g = M_PI * D - S_l;

	double D_g = 4.0l * A_g / (S_g + S_i);

	double U_g = U_sg * A / A_g;
	double U_l = U_sl * A / (A - A_g);
	double U_i = U_g - U_l;

	double Re_g = rho_g * U_g * D_g / mu_g;
	double Fr_l = U_l / sqrt(g * D * h_ld);

	double f_wg = cal_fric_fact_wellbore(Re_g, Re_w, D, D * e);
	double tau_wg = 0.5l * f_wg * rho_g * U_g * U_g;

	/*double A_l = 0.25*M_PI*D*D-A_g;
	double D_l = 4.0l*(A_l) / (S_l);
	double Re_l = rho_l*U_l*D_g / mu_g;
	double f_wl = cal_fric_fact_wellbore(Re_l, Re_w, D, D*e);
	double tau_wl = 0.5l*f_wl*rho_l*U_l*U_l;*/


	double f_i = calc_fi_stratified(Fr_l, Re_sl, rho_l, rho_g, g, D, U_g);
	double tau_i = 0.5l * f_i * rho_g * U_i * abs(U_i);

	/*cout << "estratificado" << endl;
	cout << "H_l: " << *H_l << endl;
	cout << "-tau_i*S_i: " << -tau_i*S_i<< endl;
	cout << "-tau_wg*S_g: " << -tau_wg*S_g << endl;
	cout << "Atrito: " << (-tau_i*S_i - tau_wg*S_g)*deltaL / A_g << endl;
	cout << "Acelera: " << -2.0l*rho_g*U_g*q_I_g*deltaL / A_g << endl;
	cout << "Atrito2: " << (tau_i*S_i - tau_wl*S_l)*deltaL / A_l << endl;
	cout << "Acelera2: " << -2.0l*rho_l*U_l*q_I_l*deltaL / A_l << endl;

	system("pause");*/


	* dpdx = ((-tau_i * S_i - tau_wg * S_g - rho_g * A_g * g * sin(teta) - 2.0l * rho_g * U_g * q_I_g) / A_g) * deltaL;

}

double calc_Fe_annularmist(double mu_l, double U_sg, double rho_g, double rho_l, double sigma, double U_sl) {

	double A = 0.735l * pow(mu_l * mu_l * U_sg * U_sg * rho_g / (sigma * sigma * rho_l), 0.074l) * pow(U_sg / U_sl, 0.2l);

	return A / (A + 1.0l);

}

double calc_fi_annularmist(double f_c, double Re_f, double sigma, double rho_c, double U_c, double D_c) {

	return 0.24l * f_c * pow(Re_f, 0.305l) * pow(sigma / (rho_c * U_c * U_c * D_c), 0.085l);

}

void annularmist_flow_calc(double* H_l, double* dpdx, double U_sl, double U_sg, double mu_l, double mu_g,
	double rho_l, double rho_g, double  q_I_l, double q_I_g, double sigma, double teta, double g,
	double D, double A, double deltaL, double e) {


	double Re_sl = rho_l * U_sl * D / mu_l;
	double Re_sg = rho_g * U_sg * D / mu_g;
	double Re_w = calc_Re_w(q_I_l, q_I_g, rho_l, rho_g, mu_l, mu_g);

	double f_sl = cal_fric_fact_wellbore(Re_sl, Re_w, D, D * e);
	double f_sg = cal_fric_fact_wellbore(Re_sg, Re_w, D, D * e);

	double dpdx_sg = -2.0l * f_sg * U_sg * U_sg * rho_g / D;

	double X2 = f_sl * U_sl * U_sl * rho_l / (f_sg * U_sg * U_sg * rho_g);
	double Y = (rho_l - rho_g) * g * sin(teta) / dpdx_sg;

	double delta_ld = riddermethodforlockhartmartinelli(0.99999l, 0.000001l, 1, 8, 2, D, A, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g,
		Re_sl, Re_sg, Re_w, q_I_l, q_I_g, f_sl, f_sg, dpdx_sg, sigma, g, X2, Y, e);

	double A_c = A - M_PI * D * D * delta_ld * (1.0l - delta_ld);

	double Fe = calc_Fe_annularmist(mu_l, U_sg, rho_g, rho_l, sigma, U_sl);

	double H_c = A_c / A;
	double H_lc = Fe * U_sl / (U_sg + Fe * U_sl);

	*H_l = 1.0l - H_c * (1.0l - H_lc);

	double S_i = M_PI * D * (1.0l - 2.0l * delta_ld);

	double D_c = 4.0l * A_c / S_i;
	double D_f = 4.0l * D * delta_ld * (1.0l - delta_ld);

	double U_f = U_sl * (1.0l - Fe) * A / (1 - A_c);
	double U_c = (U_sg + Fe * U_sl) * A / A_c;
	double U_i = U_c - U_f;

	double rho_c = rho_l * H_lc + rho_g * (1.0l - H_lc);
	double mu_c = mu_l * H_lc + mu_g * (1.0l - H_lc);

	double Re_c = rho_c * U_c * D_c / mu_c;
	double Re_f = rho_l * U_f * D_f / mu_l;

	double f_c = cal_fric_fact_turb_colle(Re_c, D, e);
	double f_i = calc_fi_annularmist(f_c, Re_f, sigma, rho_c, U_c, D_c);

	double tau_i = 0.5l * f_i * rho_c * U_i * abs(U_i);

	/*cout << "Anular" << endl;
	cout << "Atrito: " << -tau_i*S_i*deltaL / A_c << endl;
	cout << "Acelera: " << -2.0l*rho_c*U_c*q_I_g / (1 - Fe)*deltaL / A_c << endl;
	cout << "Grave isso: " << -rho_c*A_c*g*sin(teta)*deltaL / A_c << endl;
	system("pause");*/

	*dpdx = (-tau_i * S_i - rho_c * A_c * g * sin(teta) - 2.0l * rho_c * U_c * q_I_g / (1 - Fe)) * deltaL / A_c;

}

double riddermethodforlockhartmartinelli(double V_Dmax, double V_Dmin, int nivel, int n_threads, int pattern, double D, double A,
	double U_sl, double U_sg, double mu_l, double mu_g, double rho_l, double rho_g, double Re_sl, double Re_sg, double Re_w,
	double  q_I_l, double q_I_g, double f_sl, double f_sg, double dpdx_sg, double sigma, double g, double X2, double Y, double e) {

	double* V_D_k = (double*)malloc(n_threads * sizeof(double));
	double* G_k = (double*)malloc(n_threads * sizeof(double));

	for (int k = 0; k < n_threads; k++) {

		V_D_k[k] = k * (V_Dmax - V_Dmin) / (n_threads - 1) + V_Dmin;
		G_k[k] = calc_G(pattern, V_D_k[k], D, A, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, Re_sl, Re_sg, Re_w, q_I_l, q_I_g, f_sl, f_sg, dpdx_sg, sigma, g, X2, Y, e);

	}

	int cont = 0;
	int froot = -1;
	for (int k = 0; k < n_threads - 1; k++) {

		if (G_k[k] == 0.0l && k != n_threads - 1) {

			if (k == 0) return V_D_k[0];

			return riddermethodforlockhartmartinelli(V_D_k[k], V_Dmin, nivel + 1, n_threads, pattern, D, A, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, Re_sl, Re_sg, Re_w, q_I_l, q_I_g, f_sl, f_sg, dpdx_sg, sigma, g, X2, Y, e);

		}

		if ((abs(G_k[k + 1]) / G_k[k + 1]) != (abs(G_k[k]) / G_k[k]) && G_k[k] != 0.0l) {

			if (froot == -1) {
				froot = k;
			}

			cont++;
		}

	}


	if (nivel >= 5 && cont == 0 && G_k[n_threads - 1] == 0.0l) return V_D_k[n_threads - 1];


	if (cont == 1 && nivel < 5)
	{

		if (froot + 1 == n_threads - 1) return riddermethodforlockhartmartinelli(V_D_k[froot + 1], V_Dmin, nivel + 1, 2 * n_threads, pattern, D, A, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, Re_sl, Re_sg, Re_w, q_I_l, q_I_g, f_sl, f_sg, dpdx_sg, sigma, g, X2, Y, e);

		else return riddermethodforlockhartmartinelli(V_D_k[froot + 1], V_Dmin, nivel + 1, n_threads, pattern, D, A, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, Re_sl, Re_sg, Re_w, q_I_l, q_I_g, f_sl, f_sg, dpdx_sg, sigma, g, X2, Y, e);

	}
	else if (cont == 0) {

		return riddermethodforlockhartmartinelli(V_Dmax, V_Dmin, nivel + 1, 2 * n_threads, pattern, D, A, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, Re_sl, Re_sg, Re_w, q_I_l, q_I_g, f_sl, f_sg, dpdx_sg, sigma, g, X2, Y, e);

	}
	else {
		double V_D_l = V_D_k[froot];
		double V_D_h = V_D_k[froot + 1];
		double G_l = G_k[froot];
		double G_h = G_k[froot + 1];
		double inte = 0;

		while (1) {

			double V_D_m = 0.5l * (V_D_h + V_D_l);
			double G_m = calc_G(pattern, V_D_m, D, A, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, Re_sl, Re_sg, Re_w, q_I_l, q_I_g, f_sl, f_sg, dpdx_sg, sigma, g, X2, Y, e);

			if (G_m == 0.0l) return V_D_m;

			double s = sqrt(G_m * G_m - G_l * G_h);
			double V_D_new = V_D_m + (V_D_m - V_D_l) * (G_l - G_h) * G_m / (abs(G_l - G_h) * s);
			double G_new = calc_G(pattern, V_D_new, D, A, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, Re_sl, Re_sg, Re_w, q_I_l, q_I_g, f_sl, f_sg, dpdx_sg, sigma, g, X2, Y, e);

			if (G_new == 0.0l) return V_D_new;

			if (abs(G_new) / G_new != abs(G_m) / G_m) {
				V_D_l = V_D_m;
				G_l = G_m;
				V_D_h = V_D_new;
				G_h = G_new;
			}
			else if (abs(G_new) / G_new != abs(G_l) / G_l) {
				V_D_h = V_D_new;
				G_h = G_new;
			}
			else if (abs(G_new) / G_new != abs(G_h) / G_h) {
				V_D_l = V_D_new;
				G_l = G_new;
			}


			if ((abs(V_D_l - V_D_h)) < 0.00000001 || inte > 100) return 0.5l * (V_D_h + V_D_l);

			inte++;

		}
	}

}

double calc_G(int pattern, double V_D, double D, double A, double U_sl, double U_sg, double mu_l, double mu_g,
	double rho_l, double rho_g, double Re_sl, double Re_sg, double Re_w, double  q_I_l, double q_I_g,
	double f_sl, double f_sg, double dpdx_sg, double sigma, double g, double X2, double Y, double e) {

	//Stratified Flow V_D=h_ld
	if (pattern == 1) {

		double S_i = D * sqrt(1.0l - pow(2.0l * V_D - 1.0l, 2.0l));
		double H_l = (M_PI - acos(2.l * V_D - 1.l) + (2.l * V_D - 1.l) * S_i / D) / M_PI;
		//cout << "H_l: " << H_l << endl;
		double A_g = calc_A_G(V_D, D);
		double A_l = A - A_g;

		double S_l = D * (M_PI - acos(2.0l * V_D - 1.0l));
		double S_g = M_PI * D - S_l;

		double D_l = 4.0l * A_l / S_l;
		double D_g = 4.0l * A_g / (S_g + S_i);

		double U_l = U_sl * A / A_l;
		double U_g = U_sg * A / A_g;
		double U_i = U_g - U_l;

		double Re_l = rho_l * U_l * D_l / mu_l;
		double Re_g = rho_g * U_g * D_g / mu_g;
		double Fr_l = U_l / sqrt(g * D * V_D);

		double f_wl = cal_fric_fact_wellbore(Re_l, Re_w, D, D * e);
		double f_wg = cal_fric_fact_wellbore(Re_g, Re_w, D, D * e);
		double f_i = calc_fi_stratified(Fr_l, Re_sl, rho_l, rho_g, g, D, U_g);

		double I = 2.0l * (rho_l * U_l * q_I_l / A_l - rho_g * U_g * q_I_g / A_g) / (dpdx_sg);
		double F_1 = (f_wg / f_sg) * (U_g / U_sg) * (U_g / U_sg) * (D * S_g / A_g + (f_i / f_wg) * U_i * abs(U_i) * D * S_i / (U_g * U_g) * (1.0l / A_l + 1.0l / A_g));
		double F_2 = (f_wl / f_sl) * (U_l / U_sl) * (U_l / U_sl) * (D * S_l / A_l);

		/*cout << "G: " << X2*F_2 - F_1 - 4.0l*(Y + I) << endl;
		system("pause");*/
		return X2 * F_2 - F_1 - 4.0l * (Y + I);

	}


	//Annular Mist Flow V_D=delta_ld
	if (pattern == 2) {

		double A_f = M_PI * D * D * V_D * (1.0l - V_D);
		double A_c = A - A_f;

		double S_l = M_PI * D;
		double S_i = M_PI * D * (1.0l - 2.0l * V_D);

		double D_f = 4.0l * A_f / S_l;
		double D_c = 4.0l * A_c / S_i;

		double Fe = calc_Fe_annularmist(mu_l, U_sg, rho_g, rho_l, sigma, U_sl);

		double U_f = U_sl * (1.0l - Fe) * A / A_f;
		double U_c = (U_sg + Fe * U_sl) * A / A_c;
		double U_i = U_c - U_f;

		double H_c = A_c / A;
		double H_f = A_f / A;
		double H_lc = Fe * U_sl / (U_sg + Fe * U_sl);

		double rho_c = rho_l * H_lc + rho_g * (1.0l - H_lc);
		double mu_c = mu_l * H_lc + mu_g * (1.0l - H_lc);

		double Re_f = rho_l * U_f * D_f / mu_l;
		double Re_c = rho_c * U_c * D_c / mu_c;

		double f_wl = cal_fric_fact_wellbore(Re_f, Re_w, D, D * e);
		double f_c = cal_fric_fact_turb_colle(Re_c, D, e);

		double f_i = calc_fi_annularmist(f_c, Re_f, sigma, rho_c, U_c, D_c);

		double I = (2.0l * (rho_l * U_f / A_f) * (q_I_l - Fe * q_I_g / (1.0l - Fe)) - 2.0l * rho_c * U_c * q_I_g / (A_c * (1.0l - Fe))) / dpdx_sg;
		double F_1 = (f_i / f_sg) * (rho_c / rho_g) * (U_c / U_sg) * (U_c / U_sg) * D * S_i * (1.0l / A_f + 1.0l / A_c);
		double F_2 = (f_wl / f_sl) * (U_f / U_sl) * (U_f / U_sl) * (D * S_l / A_f);

		return X2 * F_2 - F_1 - 4.0l * (Y + I);
	}


}

void beb_calc(double* H_l, double* dpdx, double U_sl, double U_sg, double mu_l, double mu_g,
	double rho_l, double rho_g, double  q_I_l, double q_I_g, double teta, double g,
	double D, double deltaL, double e, double rho_l_2, double rho_g_2) {


	double U_m = U_sg + U_sl;

	double sigma = Tensuper(rho_l, rho_g);
	double lambda = U_sl / U_m;
	double X = log(lambda);
	double N_fr = U_m * U_m / (g * D);
	double N_lv = U_sl * pow(rho_l / (g * sigma), 0.25l);

	int pattern = beb_pattern(N_fr, X);

	double H_l_0 = max(beb_horizontal(pattern, N_fr, lambda), lambda);

	double C = beb_C(pattern, N_fr, N_lv, lambda, teta);

	*H_l = min(1.0l, H_l_0 * (1.l + C * (sin(1.8l * teta) - pow(sin(1.8l * teta), 3.l) / 3.l)));

	double y = lambda / (*H_l * *H_l);

	double S = beb_S(y);

	double G_m = (rho_l * lambda + rho_g * (1.l - lambda)) * U_m;
	double Re_ns = G_m * D / (mu_l * lambda + mu_g * (1.l - lambda));
	double f_ns = pow(2.l * log10(Re_ns / (4.25223l * log10(Re_ns) - 3.8215l)), -2.l);

	double f_tp = exp(S) * f_ns;
	double rho_tp = rho_l * *H_l + rho_g * (1.0l - *H_l);
	double rho_tp_2 = rho_l_2 * *H_l + rho_g_2 * (1.0l - *H_l);
	double U_m2 = U_m * rho_tp / rho_tp_2;

	/*cout << "Atrito: " << 0.5*f_tp*G_m*U_m*deltaL / D << endl;
	cout << "Acelera: " << rho_tp*U_m*(U_m2 - U_m) << endl;
	cout << "Grave: " << rho_tp*g*sin(teta)*deltaL << endl;
	system("pause");*/

	*dpdx = -(0.5 * f_tp * G_m * U_m / D + rho_tp * U_m * (U_m2 - U_m) / deltaL + rho_tp * g * sin(teta)) * deltaL;

}

void beb_mod_calc(double* H_l, double* dpdx, double U_sl, double U_sg, double mu_l, double mu_g,
	double rho_l, double rho_g, double  q_I_l, double q_I_g, double teta, double g,
	double D, double deltaL, double e) {


	double U_m = U_sg + U_sl;
	double q_I_m = (q_I_l + q_I_g);

	double sigma = Tensuper(rho_l, rho_g);
	double lambda = U_sl / U_m;
	double X = log(lambda);
	double N_fr = U_m * U_m / (g * D);
	double N_lv = U_sl * pow(rho_l / (g * sigma), 0.25l);

	int pattern = beb_pattern(N_fr, X);

	double H_l_0 = max(beb_horizontal(pattern, N_fr, lambda), lambda);

	double C = beb_C(pattern, N_fr, N_lv, lambda, teta);

	*H_l = min(1.0l, H_l_0 * (1.l + C * (sin(1.8l * teta) - pow(sin(1.8l * teta), 3.l) / 3.l)));

	//if (*H_l < lambda) cout << "morte da interacao na internet" << endl;

	double y = lambda / (*H_l * *H_l);
	double S = beb_S(y);

	double mu_tp = mu_l * *H_l + mu_g * (1.0l - *H_l);
	double rho_tp = rho_l * *H_l + rho_g * (1.0l - *H_l);
	double U_tp = (rho_l * U_sl + rho_g * U_sg) / rho_tp;
	double q_I_tp = (rho_l * q_I_l + rho_g * q_I_g) / rho_tp;
	double Re_tp = rho_tp * U_tp * D / mu_tp;
	double Re_w = calc_Re_w(q_I_l, q_I_g, rho_l, rho_g, mu_l, mu_g);

	double f_tp = exp(S) * cal_fric_fact_wellbore(Re_tp, Re_w, D, D * e);
	//double f_tp = cal_fric_fact_wellbore(Re_tp, Re_w, D, D*e);

	double tau_w = 0.5l * f_tp * rho_tp * U_tp * U_tp;

	double dpdxaw1 = 4.0l * rho_tp * (U_m * q_I_tp + U_tp * q_I_m) / (M_PI * D * D);
	double dpdxaw2 = 8.0l * rho_tp * U_tp * q_I_tp / (M_PI * D * D);

	double omega = 0.8l;
	double dpdxaw = omega * dpdxaw1 + (1.l - omega) * dpdxaw2;

	/*double Betae = rho_tp*U_m*U_sg / p;
	double dpdxae = (Betae / (1.l - Betae))*(dpdxf + dpdxg + dpdxaw);*/

	/*cout << "Atrito: " << 4.0l*tau_w*deltaL / D << endl;
	cout << "Acelera: " << dpdxaw << endl;
	cout << "Grave isso: " << rho_tp*g*sin(teta)*deltaL << endl;
	system("pause");*/

	*dpdx = -(4.0l * tau_w / D + dpdxaw + rho_tp * g * sin(teta)) * deltaL;

}

int beb_pattern(double N_fr, double X) {

	double L_1 = exp(-4.62l - 3.757l * X - 0.481l * X * X - 0.0207 * X * X * X);
	double L_2 = exp(1.061l - 4.602l * X - 1.609l * X * X - 0.179 * X * X * X + 0.635e-3 * X * X * X * X * X);

	if (N_fr < L_1) return 1; //segregated

	if (N_fr > L_1 && N_fr < L_2) return 2; //intermittent

	if (N_fr > L_1 && N_fr > L_2) return 3; //distributed	


}

double beb_horizontal(int pattern, double N_fr, double lambda) {


	if (pattern == 1) return 0.98l * pow(lambda, 0.4846l) / pow(N_fr, 0.0868l);

	if (pattern == 2) return 0.845l * pow(lambda, 0.5351l) / pow(N_fr, 0.0173l);

	if (pattern == 3) return 1.065l * pow(lambda, 0.5824l) / pow(N_fr, 0.0609l);

}

double beb_C(int pattern, double N_fr, double N_lv, double lambda, double teta) {

	double D, delta, epsilon, ksi;

	if (teta < 0.0l) {
		D = 4.7l;
		delta = -0.3692l;
		epsilon = -0.5056l;
		ksi = 0.1244l;
	}
	else {
		if (pattern == 1) {
			//(pattern == 1 && teta >= 0.0l)
			D = 0.011l;
			delta = -3.768l;
			epsilon = -1.614;
			ksi = 3.539;
		}
		if (pattern == 2) {
			D = 2.96l;
			delta = 0.305l;
			epsilon = 0.0978l;
			ksi = -0.4473l;
		}
		if (pattern == 3) {
			D = 1.l;
			delta = 0.0l;
			epsilon = 0.0l;
			ksi = 0.0l;
		}
	}



	return (1.l - lambda) * log(D * pow(lambda, delta) * pow(N_fr, epsilon) * pow(N_lv, ksi));

}

double beb_S(double y) {

	if (y < 1.2l && y>1.0l) {
		return log(2.2l * y - 1.2l);
	}
	else {
		double logy = log(y);
		return logy / (-0.0523l + 3.182l * logy - 0.8725 * logy * logy + 0.01853 * logy * logy * logy * logy);
	}

}

