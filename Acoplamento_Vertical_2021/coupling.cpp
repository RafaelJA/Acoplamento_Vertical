#include "stdafx.h"
#include "coupling.h"

void coup_step_4(int coup_ts, int ite_coup,
	double BHP_r, double* BHP,
	double THP, double L, int nz, double* deltaz,
	double* Qw_total, double* Qo_total, double* Qg_total,
	double* PWB_ref, double* PWB, double* PWB_0, double* dPWB_D, double* PWB_coup,
	double* Qw_wellbore, double* Qo_wellbore, double* Qg_wellbore,
	double* Qw_wellbore_ref, double* Qo_wellbore_ref, double* Qg_wellbore_ref,
	int N_wellbore, double D_perf, double perf_dens,
	double teta, double rug, double gama, double teta_yue,
	double* po_coup, double* So_coup, double* Sw_coup, double* pb_coup,
	double* po_res, double* pw_res, double* pg_res,
	double* So_res, double* Sw_res, double* Sg_res, double* pb_res,
	double* po_res_0, double* So_res_0, double* Sw_res_0, double* pb_res_0, double* Sg_res_0, double* X_pb_0,
	double* phi_res_0, double* Bw_res_0, double* Rs_res_0, double* Bo_res_0, double* Bg_res_0,
	double* po_res_lv0, double* pb_res_lv0, double* So_res_lv0, double* Sw_res_lv0,
	double* varpo, double* varpb, double* varSo, double* varSw,
	double Qwsinj, double Qosinj, double Qgsinj,
	double* K_F, double* K_B, double* K_N, double* K_S, double* K_W, double* K_E,
	double* X_pb, bool* zonaprodutora,
	int res_timesteps, bool dontstopmenow,
	int p, int y, int u, int l,
	int centerP, int centerSw, int centerSo,
	int i, int k, int j,
	int o, int m, int n, int dim,
	double* Kr, double* Kz, double pi, double T,
	int* IA, int* JA, double* A, double* B, double* X,
	double controladorP, double controladorSw, double controladorSo,
	double C1, double C2, double C3,
	double* Ageo_R, double* Ageo_PHI, double* Ageo_Z, double* Ageo_t,
	double* Adir_F, double* Adir_B, double* Adir_W, double* Adir_E, double* Adir_N, double* Adir_S,
	double* CHI_F, double* CHI_B, double* CHI_W, double* CHI_E, double* CHI_N, double* CHI_S,
	double* R_P, double* PHI_P, double* Z_P,
	double* R_PB, double* R_PF, double* delta_PHI, double* delta_Z, double deltat_res,
	double g, double Re, double Rw,
	double rho_air, double rho_ws, double API, double yo, double rho_os, double yg, double rho_gs,
	double yn2, double yco2, double yh2s,
	double* lambda1_W, double* d_lambda1_po_W, double* d_lambda1_SW_W,
	double* lambda2_W, double* d_lambda2_po_W, double* d_lambda2_SW_W,
	double* lambda1_O, double* d_lambda1_po_O, double* d_lambda1_SW_O, double* d_lambda1_pb_O, double* d_lambda1_SO_O,
	double* lambda2_O, double* d_lambda2_po_O, double* d_lambda2_SW_O, double* d_lambda2_pb_O, double* d_lambda2_SO_O,
	double* lambda1_G, double* d_lambda1_po_G, double* d_lambda1_SW_G, double* d_lambda1_pb_G, double* d_lambda1_SO_G,
	double* lambda2_G, double* d_lambda2_po_G, double* d_lambda2_SW_G, double* d_lambda2_pb_G, double* d_lambda2_SO_G,
	double* lambda3_G, double* d_lambda3_po_G, double* d_lambda3_SW_G, double* d_lambda3_pb_G, double* d_lambda3_SO_G,
	double* lambda4_G, double* d_lambda4_po_G, double* d_lambda4_SW_G, double* d_lambda4_pb_G, double* d_lambda4_SO_G,
	double krwmax, double nw, double ng, double now, double nog,
	double* phi_ref, double Cr,
	double Swc, double Sor, double Sgr, double epsilon,
	double pcwmin, double pcwmax, double pcgmin, double pcgmax,
	double* parametrosBw, double* parametrosmig, double* parametrosZ,
	double ppc, double Tpc, double Bwi, double Bpcw, double Bpcg,
	int* iparm, double* dparm, int solver, int maxfct, int mnum, int phase,
	int error, int msglvl, int mtype, int nrhs, double ddum, int idum, void* pt,
	double* phi, double* der_phi,
	double* rho_w, double* der_rho_w,
	double* mi_w, double* der_mi_w,
	double* Bw, double* der_Bw,
	double* pcw, double* der_pcw,
	double* krw, double* der_krw,
	double* Rs, double* der_Rs_p, double* der_Rs_pb,
	double* kro, double* der_kro_Sw, double* der_kro_So,
	double* Bo, double* der_Bo_p, double* der_Bo_pb,
	double* mi_o, double* der_mi_o_p, double* der_mi_o_pb,
	double* rho_o, double* der_rho_o_p, double* der_rho_o_pb,
	double* pcg, double* der_pcg_Sw, double* der_pcg_So,
	double* krg, double* der_krg,
	double* mi_g, double* der_mi_g, double* Z,
	double* Bg, double* der_Bg,
	double* rho_g, double* der_rho_g,
	int* unW, int* ufW, int* uwbW, int* ueW, int* ucW, int* uwW, int* uebW, int* ubW, int* usW,
	int* unO, int* ufO, int* uwbO, int* ueO, int* ucO, int* uwO, int* uebO, int* ubO, int* usO,
	int* unG, int* ufG, int* uwbG, int* ueG, int* ucG, int* uwG, int* uebG, int* ubG, int* usG,
	int* neg, bool* NpureW, int* contresint,
	char directory[], __int8 W_corre, __int8 WB_corre) {

	double* Pb_WB = (double*)calloc(N_wellbore, sizeof(double));
	double* H_l = (double*)calloc(N_wellbore, sizeof(double));

	double deltaBHP = 50.0l / log(coup_ts + 2);
	bool controldelta = 1;
	int avoidshit = 0;

	double indice1 = 0.0;
	double BHP1 = 0.0;
	double indice2 = 0.0;
	double BHP2 = 0.0;

	double error_WB_R = 100.l;
	double PWB_1, PWB_2, PWB_resp, f1, f2, awb;
	int count_WB_R = 0;

	ite_coup = 1;
	BHP_r = 0;
	*contresint = 0;

	char subdir[50], subdirname[100];
	memset(&subdir[0], 0, sizeof(subdir));
	memset(&subdirname[0], 0, sizeof(subdirname));

	sprintf_s(subdir, directory);
	sprintf_s(subdirname, "\\TimeStep_%d", coup_ts);

	strcat_s(subdir, subdirname);
	_mkdir(subdir);

	char arquive[120];
	sprintf_s(arquive, subdir);
	strcat_s(arquive, "\\BHP_x_BHP_r.m");

	fstream bhp(arquive, ios::out | ios::app);
	bhp << "BHP_guess[psi];" << "BHP_guess[psi];" << "BHPg/BHPr;" << "Qo[bbl/d];" << endl;

	while (fabs((*BHP - BHP_r) / *BHP) > 0.001l || count_WB_R > 60) {

		calc_guess_BHP(ite_coup, coup_ts,
			BHP_r, BHP, &BHP1, &BHP2, &indice1, &indice2,
			&deltaBHP, &controldelta, &avoidshit, po_coup[0], THP);

		error_WB_R = 100.l;
		count_WB_R = 1;
		PWB_1 = PWB[N_wellbore - 1];

		while (error_WB_R > 0.001l && count_WB_R <= 60) {

			error_WB_R = 0.0l;
			p = 0;
			for (k = 0; k <= (o - 1); k++) {
				for (j = 0; j <= (m - 1); j++) {
					for (i = 0; i <= (n - 1); i++) {
						po_res[p] = po_coup[p];
						Sw_res[p] = Sw_coup[p];
						So_res[p] = So_coup[p];
						pb_res[p] = pb_coup[p];
						Sg_res[p] = 1 - So_coup[p] - Sw_coup[p];
						pw_res[p] = po_coup[p];// +calc_water_pcw_Sw(Sw[p], Swc, epsilon, pcwmin, pcwmax, Bpcw);
						pg_res[p] = po_coup[p];// +calc_gas_pcg(Sw[p], So[p], Swc, Sor, epsilon, pcgmin, pcgmax, Bpcg);
						po_res_0[p] = po_coup[p];
						pb_res_0[p] = pb_coup[p];
						Sw_res_0[p] = Sw_coup[p];
						So_res_0[p] = So_coup[p];
						Sg_res_0[p] = 1 - So_coup[p] - Sw_coup[p];
						p++;
					}
				}
			}

			for (k = 0; k < N_wellbore; k++) {
				PWB_0[k] = PWB_coup[k];
			}

			calc_reser_evol(res_timesteps, p, y, u, l, centerP, centerSw, centerSo,
				IA, JA, A, B, X,
				dontstopmenow, i, k, j, o, m, n, dim,
				controladorP, controladorSw, controladorSo,
				lambda1_W, d_lambda1_po_W, d_lambda1_SW_W,
				lambda2_W, d_lambda2_po_W, d_lambda2_SW_W,
				lambda1_O, d_lambda1_po_O, d_lambda1_SW_O, d_lambda1_pb_O, d_lambda1_SO_O,
				lambda2_O, d_lambda2_po_O, d_lambda2_SW_O, d_lambda2_pb_O, d_lambda2_SO_O,
				lambda1_G, d_lambda1_po_G, d_lambda1_SW_G, d_lambda1_pb_G, d_lambda1_SO_G,
				lambda2_G, d_lambda2_po_G, d_lambda2_SW_G, d_lambda2_pb_G, d_lambda2_SO_G,
				lambda3_G, d_lambda3_po_G, d_lambda3_SW_G, d_lambda3_pb_G, d_lambda3_SO_G,
				lambda4_G, d_lambda4_po_G, d_lambda4_SW_G, d_lambda4_pb_G, d_lambda4_SO_G,
				Kr, Kz, pi, T,
				PWB, PWB_0, Pb_WB, N_wellbore,
				po_res, pb_res, pw_res, pg_res, Sw_res, So_res, Sg_res,
				po_res_0, pb_res_0, Sw_res_0, So_res_0, Sg_res_0, X_pb_0,
				phi_res_0, Bw_res_0, Rs_res_0, Bo_res_0, Bg_res_0,
				po_res_lv0, pb_res_lv0, So_res_lv0, Sw_res_lv0,
				varpo, varpb, varSo, varSw,
				C1, C2, C3,
				Ageo_R, Ageo_PHI, Ageo_Z, Ageo_t,
				Adir_F, Adir_B, Adir_W, Adir_E, Adir_N, Adir_S,
				CHI_F, CHI_B, CHI_W, CHI_E, CHI_N, CHI_S,
				R_P, PHI_P, Z_P,
				R_PB, R_PF, delta_PHI, delta_Z, deltat_res,
				g, Re, Rw,
				rho_air, rho_ws, rho_os, rho_gs,
				API, yg, yo,
				phi_ref, Cr,
				Qwsinj, Qosinj, Qgsinj,
				Qw_wellbore, Qo_wellbore, Qg_wellbore, coup_ts, ite_coup,
				phi, der_phi,
				rho_w, der_rho_w,
				mi_w, der_mi_w,
				Bw, der_Bw,
				pcw, der_pcw,
				krw, der_krw,
				Rs, der_Rs_p, der_Rs_pb,
				kro, der_kro_Sw, der_kro_So,
				Bo, der_Bo_p, der_Bo_pb,
				mi_o, der_mi_o_p, der_mi_o_pb,
				rho_o, der_rho_o_p, der_rho_o_pb,
				pcg, der_pcg_Sw, der_pcg_So,
				krg, der_krg, mi_g, der_mi_g,
				Bg, der_Bg, rho_g, der_rho_g,
				K_F, K_B, K_N, K_S, K_W, K_E,
				X_pb, zonaprodutora, 0,
				unW, ufW, uwbW, ueW, ucW, uwW, uebW, ubW, usW,
				unO, ufO, uwbO, ueO, ucO, uwO, uebO, ubO, usO,
				unG, ufG, uwbG, ueG, ucG, uwG, uebG, ubG, usG,
				iparm, dparm, solver, maxfct, mnum, phase, error,
				msglvl, mtype, nrhs, ddum, idum, pt,
				Swc, Sor, Sgr, epsilon,
				pcwmin, pcwmax, pcgmin, pcgmax,
				krwmax, nw, ng, nog, now,
				parametrosBw, parametrosmig, ppc, Tpc,
				Bwi, Bpcw, Bpcg,
				yn2, yco2, yh2s,
				neg, NpureW, subdir);

			*contresint = *contresint + 1;

			for (k = 0; k < N_wellbore; k++) {
				deltaz[k] = ft_m(deltaz[k]);
				Qw_wellbore[k] = Qw_wellbore[k] / deltaz[k];
				Qo_wellbore[k] = Qo_wellbore[k] / deltaz[k];
				Qg_wellbore[k] = Qg_wellbore[k] / deltaz[k];
			}

			wellbore_multi(WB_corre, N_wellbore,
				PWB, Pb_WB, H_l,
				Qw_wellbore, Qo_wellbore, Qg_wellbore,
				*BHP, PWB_1, deltaz,
				2.l * Rw, D_perf, perf_dens,
				API, yo, yg,
				rho_os, rho_ws, rho_gs, T,
				g, teta, rug, gama,
				parametrosmig, parametrosBw,
				Tpc, ppc, pi, Bwi);

			for (k = 0; k < N_wellbore; k++) {
				Qw_wellbore[k] = Qw_wellbore[k] * deltaz[k];
				Qo_wellbore[k] = Qo_wellbore[k] * deltaz[k];
				Qg_wellbore[k] = Qg_wellbore[k] * deltaz[k];
				deltaz[k] = m_ft(deltaz[k]);
			}

			//MÉTODO NUM�RICO 
			if (count_WB_R == 1) {

				PWB_resp = PWB[N_wellbore - 1];
				PWB_2 = PWB_1;
				f2 = 1.0l - PWB_resp / PWB_1;
				PWB_1 = 0.5l * (PWB_1 + PWB_resp);

			}
			else {
				PWB_resp = PWB[N_wellbore - 1];
				f1 = 1.0l - PWB_resp / PWB_1;
				awb = (f1 - f2) / (PWB_1 - PWB_2);

				if (PWB_1 == PWB_2 || abs(f1 / (awb)) > 500.l || (PWB_1 - f1 / awb) > 1.05l * po_coup[0] || isnan(PWB_1 - f1 / awb)) {
					PWB_1 = 0.5l * (PWB_1 + PWB_resp);
				}
				else {
					PWB_2 = PWB_1;
					f2 = f1;
					PWB_1 = PWB_1 - f1 / awb;
					if (isnan(PWB_1)) {
						cout << "PWB_1: " << PWB_1 << endl;
						system("pause");
					}
				}
			}

			//Extract the shape of PWB CURVE
			for (k = N_wellbore - 2; k >= 0; k--) {

				dPWB_D[k] = (PWB[k] - PWB[k + 1]) / (*BHP - PWB_resp);

			}

			*Qo_total = 0.0l;
			*Qw_total = 0.0l;
			*Qg_total = 0.0l;
			for (k = N_wellbore - 1; k >= 0; k--) {

				if (abs((Qo_wellbore_ref[k] - Qo_wellbore[k]) / Qo_wellbore_ref[k]) > error_WB_R) {

					error_WB_R = abs((Qo_wellbore_ref[k] - Qo_wellbore[k]) / Qo_wellbore_ref[k]);

				}
				if (abs((Qw_wellbore_ref[k] - Qw_wellbore[k]) / Qw_wellbore_ref[k]) > error_WB_R) {

					error_WB_R = abs((Qw_wellbore_ref[k] - Qw_wellbore[k]) / Qw_wellbore_ref[k]);

				}
				if (abs((Qg_wellbore_ref[k] - Qg_wellbore[k]) / Qg_wellbore_ref[k]) > error_WB_R) {

					error_WB_R = abs((Qg_wellbore_ref[k] - Qg_wellbore[k]) / Qg_wellbore_ref[k]);

				}
				if (abs((PWB_ref[k] - PWB[k]) / PWB[k]) > error_WB_R) {

					error_WB_R = abs((PWB_ref[k] - PWB[k]) / PWB[k]);

				}

				if (k == N_wellbore - 1) {

					PWB[k] = PWB_1;

				}
				else {

					PWB[k] = PWB[k + 1] + (*BHP - PWB_1) * dPWB_D[k];

				}

				*Qo_total = *Qo_total + Qo_wellbore[k];
				*Qw_total = *Qw_total + Qw_wellbore[k];
				*Qg_total = *Qg_total + Qg_wellbore[k];

				Qo_wellbore_ref[k] = Qo_wellbore[k];
				Qw_wellbore_ref[k] = Qw_wellbore[k];
				Qg_wellbore_ref[k] = Qg_wellbore[k];
				PWB_ref[k] = PWB[k];

			}

			if (coup_ts % 50 == 0) {
				grava_var(PWB, N_wellbore, 0, count_WB_R, ite_coup, subdir, 5);
				grava_var(Qo_wellbore, N_wellbore, 0, count_WB_R, ite_coup, subdir, 6);
				grava_var(Qg_wellbore, N_wellbore, 0, count_WB_R, ite_coup, subdir, 7);
			}

			count_WB_R++;

		}


		BHP_r = Pa_Psia(converge_p_sup(L, nz, ft_m(2.l * Rw), Psia_Pa(*BHP), Far_Kel(T), (*Qw_total + *Qo_total),
			(*Qg_total / *Qo_total), API, yg, (*Qw_total) / (*Qw_total + *Qo_total), THP, pi, rug, coup_ts, ite_coup, W_corre, directory));


		bhp << *BHP << "\t;" << BHP_r << "\t;" << (BHP_r / *BHP) << "\t;" << *Qo_total << "\t;" << endl;

		ite_coup++;
		avoidshit = avoidshit + 1;

	}

	bhp.close();

	p = 0;
	for (k = 0; k <= (o - 1); k++) {
		for (j = 0; j <= (m - 1); j++) {
			for (i = 0; i <= (n - 1); i++) {
				po_coup[p] = po_res_0[p];
				Sw_coup[p] = Sw_res_0[p];
				So_coup[p] = So_res_0[p];
				pb_coup[p] = pb_res_0[p];
				p++;
			}
		}
	}

	for (k = 0; k < N_wellbore; k++) {
		PWB_0[k] = PWB[k];
		PWB_coup[k] = PWB[k];
	}

	free(Pb_WB);
	free(H_l);

}

void coup_step_1(int coup_ts, int ite_coup,
	double BHP_r, double* BHP,
	double THP, double L, int nz, double* deltaz,
	double* Qw_total, double* Qo_total, double* Qg_total,
	double* PWB_ref, double* PWB, double* PWB_0, double* dPWB_D, double* PWB_coup,
	double* Qw_wellbore, double* Qo_wellbore, double* Qg_wellbore,
	double* Qw_wellbore_ref, double* Qo_wellbore_ref, double* Qg_wellbore_ref,
	int N_wellbore, double D_perf, double perf_dens,
	double teta, double rug, double gama, double teta_yue,
	double* po_coup, double* So_coup, double* Sw_coup, double* pb_coup,
	double* po_res, double* pw_res, double* pg_res,
	double* So_res, double* Sw_res, double* Sg_res, double* pb_res,
	double* po_res_0, double* So_res_0, double* Sw_res_0, double* pb_res_0, double* Sg_res_0, double* X_pb_0,
	double* phi_res_0, double* Bw_res_0, double* Rs_res_0, double* Bo_res_0, double* Bg_res_0,
	double* po_res_lv0, double* pb_res_lv0, double* So_res_lv0, double* Sw_res_lv0,
	double* varpo, double* varpb, double* varSo, double* varSw,
	double Qwsinj, double Qosinj, double Qgsinj,
	double* K_F, double* K_B, double* K_N, double* K_S, double* K_W, double* K_E,
	double* X_pb, bool* zonaprodutora,
	int res_timesteps, bool dontstopmenow,
	int p, int y, int u, int l,
	int centerP, int centerSw, int centerSo,
	int i, int k, int j,
	int o, int m, int n, int dim,
	double* Kr, double* Kz, double pi, double T,
	int* IA, int* JA, double* A, double* B, double* X,
	double controladorP, double controladorSw, double controladorSo,
	double C1, double C2, double C3,
	double* Ageo_R, double* Ageo_PHI, double* Ageo_Z, double* Ageo_t,
	double* Adir_F, double* Adir_B, double* Adir_W, double* Adir_E, double* Adir_N, double* Adir_S,
	double* CHI_F, double* CHI_B, double* CHI_W, double* CHI_E, double* CHI_N, double* CHI_S,
	double* R_P, double* PHI_P, double* Z_P,
	double* R_PB, double* R_PF, double* delta_PHI, double* delta_Z, double deltat_res,
	double g, double Re, double Rw,
	double rho_air, double rho_ws, double API, double yo, double rho_os, double yg, double rho_gs,
	double yn2, double yco2, double yh2s,
	double* lambda1_W, double* d_lambda1_po_W, double* d_lambda1_SW_W,
	double* lambda2_W, double* d_lambda2_po_W, double* d_lambda2_SW_W,
	double* lambda1_O, double* d_lambda1_po_O, double* d_lambda1_SW_O, double* d_lambda1_pb_O, double* d_lambda1_SO_O,
	double* lambda2_O, double* d_lambda2_po_O, double* d_lambda2_SW_O, double* d_lambda2_pb_O, double* d_lambda2_SO_O,
	double* lambda1_G, double* d_lambda1_po_G, double* d_lambda1_SW_G, double* d_lambda1_pb_G, double* d_lambda1_SO_G,
	double* lambda2_G, double* d_lambda2_po_G, double* d_lambda2_SW_G, double* d_lambda2_pb_G, double* d_lambda2_SO_G,
	double* lambda3_G, double* d_lambda3_po_G, double* d_lambda3_SW_G, double* d_lambda3_pb_G, double* d_lambda3_SO_G,
	double* lambda4_G, double* d_lambda4_po_G, double* d_lambda4_SW_G, double* d_lambda4_pb_G, double* d_lambda4_SO_G,
	double krwmax, double nw, double ng, double now, double nog,
	double* phi_ref, double Cr,
	double Swc, double Sor, double Sgr, double epsilon,
	double pcwmin, double pcwmax, double pcgmin, double pcgmax,
	double* parametrosBw, double* parametrosmig, double* parametrosZ,
	double ppc, double Tpc, double Bwi, double Bpcw, double Bpcg,
	int* iparm, double* dparm, int solver, int maxfct, int mnum, int phase,
	int error, int msglvl, int mtype, int nrhs, double ddum, int idum, void* pt,
	double* phi, double* der_phi,
	double* rho_w, double* der_rho_w,
	double* mi_w, double* der_mi_w,
	double* Bw, double* der_Bw,
	double* pcw, double* der_pcw,
	double* krw, double* der_krw,
	double* Rs, double* der_Rs_p, double* der_Rs_pb,
	double* kro, double* der_kro_Sw, double* der_kro_So,
	double* Bo, double* der_Bo_p, double* der_Bo_pb,
	double* mi_o, double* der_mi_o_p, double* der_mi_o_pb,
	double* rho_o, double* der_rho_o_p, double* der_rho_o_pb,
	double* pcg, double* der_pcg_Sw, double* der_pcg_So,
	double* krg, double* der_krg,
	double* mi_g, double* der_mi_g, double* Z,
	double* Bg, double* der_Bg,
	double* rho_g, double* der_rho_g,
	int* unW, int* ufW, int* uwbW, int* ueW, int* ucW, int* uwW, int* uebW, int* ubW, int* usW,
	int* unO, int* ufO, int* uwbO, int* ueO, int* ucO, int* uwO, int* uebO, int* ubO, int* usO,
	int* unG, int* ufG, int* uwbG, int* ueG, int* ucG, int* uwG, int* uebG, int* ubG, int* usG,
	int* neg, bool* NpureW, int* contresint,
	char directory[], __int8 W_corre, __int8 WB_corre) {

	double* Pb_WB = (double*)calloc(N_wellbore, sizeof(double));
	double* H_l = (double*)calloc(N_wellbore, sizeof(double));

	double deltaBHP = 50.0l / log(coup_ts + 2);
	bool controldelta = 1;
	int avoidshit = 0;

	double indice1 = 0.0;
	double BHP1 = 0.0;
	double indice2 = 0.0;
	double BHP2 = 0.0;

	double error_WB_R = 100.l;
	double PWB_1, PWB_2, PWB_resp, f1, f2, awb;
	int count_W_R = 0;

	ite_coup = 1;
	PWB_1 = PWB_coup[N_wellbore - 1];
	*contresint = 0;

	char subdir[50], subdirname[100];
	memset(&subdir[0], 0, sizeof(subdir));
	memset(&subdirname[0], 0, sizeof(subdirname));

	sprintf_s(subdir, directory);
	sprintf_s(subdirname, "\\TimeStep_%d", coup_ts);

	strcat_s(subdir, subdirname);
	_mkdir(subdir);

	char arquive[120];
	sprintf_s(arquive, subdir);
	strcat_s(arquive, "\\BHP_x_BHP_r.m");

	fstream bhp(arquive, ios::out | ios::app);
	bhp << "BHP_guess[psi];" << "BHP_guess[psi];" << "BHPg/BHPr;" << "Qo[bbl/d];" << endl;

	while (error_WB_R > 0.001l) {

		count_W_R = 1;

		BHP_r = 0;
		avoidshit = 0;
		deltaBHP = 50.0l;
		controldelta = 1;

		indice1 = 0.0;
		BHP1 = 0.0;
		indice2 = 0.0;
		BHP2 = 0.0;

		while (fabs((*BHP - BHP_r) / *BHP) > 0.001l) {

			calc_guess_BHP(count_W_R, ite_coup + coup_ts - 1,
				BHP_r, BHP, &BHP1, &BHP2, &indice1, &indice2,
				&deltaBHP, &controldelta, &avoidshit, po_coup[0], THP);

			for (k = N_wellbore - 1; k >= 0; k--) {

				if (k == N_wellbore - 1) {

					PWB[k] = PWB_1;

				}
				else {

					PWB[k] = PWB[k + 1] + (*BHP - PWB_1) * dPWB_D[k];

				}

				PWB_ref[k] = PWB[k];
			}

			p = 0;
			for (k = 0; k <= (o - 1); k++) {
				for (j = 0; j <= (m - 1); j++) {
					for (i = 0; i <= (n - 1); i++) {
						po_res[p] = po_coup[p];
						Sw_res[p] = Sw_coup[p];
						So_res[p] = So_coup[p];
						pb_res[p] = pb_coup[p];
						Sg_res[p] = 1 - So_coup[p] - Sw_coup[p];
						pw_res[p] = po_coup[p];// +calc_water_pcw_Sw(Sw[p], Swc, epsilon, pcwmin, pcwmax, Bpcw);
						pg_res[p] = po_coup[p];// +calc_gas_pcg(Sw[p], So[p], Swc, Sor, epsilon, pcgmin, pcgmax, Bpcg);
						po_res_0[p] = po_coup[p];
						pb_res_0[p] = pb_coup[p];
						Sw_res_0[p] = Sw_coup[p];
						So_res_0[p] = So_coup[p];
						Sg_res_0[p] = 1 - So_coup[p] - Sw_coup[p];
						p++;
					}
				}


			}

			for (k = 0; k < N_wellbore; k++) {
				PWB_0[k] = PWB_coup[k];
			}

			calc_reser_evol(res_timesteps, p, y, u, l, centerP, centerSw, centerSo,
				IA, JA, A, B, X,
				dontstopmenow, i, k, j, o, m, n, dim,
				controladorP, controladorSw, controladorSo,
				lambda1_W, d_lambda1_po_W, d_lambda1_SW_W,
				lambda2_W, d_lambda2_po_W, d_lambda2_SW_W,
				lambda1_O, d_lambda1_po_O, d_lambda1_SW_O, d_lambda1_pb_O, d_lambda1_SO_O,
				lambda2_O, d_lambda2_po_O, d_lambda2_SW_O, d_lambda2_pb_O, d_lambda2_SO_O,
				lambda1_G, d_lambda1_po_G, d_lambda1_SW_G, d_lambda1_pb_G, d_lambda1_SO_G,
				lambda2_G, d_lambda2_po_G, d_lambda2_SW_G, d_lambda2_pb_G, d_lambda2_SO_G,
				lambda3_G, d_lambda3_po_G, d_lambda3_SW_G, d_lambda3_pb_G, d_lambda3_SO_G,
				lambda4_G, d_lambda4_po_G, d_lambda4_SW_G, d_lambda4_pb_G, d_lambda4_SO_G,
				Kr, Kz, pi, T,
				PWB, PWB_0, Pb_WB, N_wellbore,
				po_res, pb_res, pw_res, pg_res, Sw_res, So_res, Sg_res,
				po_res_0, pb_res_0, Sw_res_0, So_res_0, Sg_res_0, X_pb_0,
				phi_res_0, Bw_res_0, Rs_res_0, Bo_res_0, Bg_res_0,
				po_res_lv0, pb_res_lv0, So_res_lv0, Sw_res_lv0,
				varpo, varpb, varSo, varSw,
				C1, C2, C3,
				Ageo_R, Ageo_PHI, Ageo_Z, Ageo_t,
				Adir_F, Adir_B, Adir_W, Adir_E, Adir_N, Adir_S,
				CHI_F, CHI_B, CHI_W, CHI_E, CHI_N, CHI_S,
				R_P, PHI_P, Z_P,
				R_PB, R_PF, delta_PHI, delta_Z, deltat_res,
				g, Re, Rw,
				rho_air, rho_ws, rho_os, rho_gs,
				API, yg, yo,
				phi_ref, Cr,
				Qwsinj, Qosinj, Qgsinj,
				Qw_wellbore, Qo_wellbore, Qg_wellbore, coup_ts, ite_coup,
				phi, der_phi,
				rho_w, der_rho_w,
				mi_w, der_mi_w,
				Bw, der_Bw,
				pcw, der_pcw,
				krw, der_krw,
				Rs, der_Rs_p, der_Rs_pb,
				kro, der_kro_Sw, der_kro_So,
				Bo, der_Bo_p, der_Bo_pb,
				mi_o, der_mi_o_p, der_mi_o_pb,
				rho_o, der_rho_o_p, der_rho_o_pb,
				pcg, der_pcg_Sw, der_pcg_So,
				krg, der_krg, mi_g, der_mi_g,
				Bg, der_Bg, rho_g, der_rho_g,
				K_F, K_B, K_N, K_S, K_W, K_E,
				X_pb, zonaprodutora, 0,
				unW, ufW, uwbW, ueW, ucW, uwW, uebW, ubW, usW,
				unO, ufO, uwbO, ueO, ucO, uwO, uebO, ubO, usO,
				unG, ufG, uwbG, ueG, ucG, uwG, uebG, ubG, usG,
				iparm, dparm, solver, maxfct, mnum, phase, error,
				msglvl, mtype, nrhs, ddum, idum, pt,
				Swc, Sor, Sgr, epsilon,
				pcwmin, pcwmax, pcgmin, pcgmax,
				krwmax, nw, ng, nog, now,
				parametrosBw, parametrosmig, ppc, Tpc,
				Bwi, Bpcw, Bpcg,
				yn2, yco2, yh2s,
				neg, NpureW,directory);

			*contresint = *contresint + 1;

			*Qo_total = 0.0l;
			*Qw_total = 0.0l;
			*Qg_total = 0.0l;
			for (k = 0; k < N_wellbore; k++) {
				*Qo_total = *Qo_total + Qo_wellbore[k];
				*Qw_total = *Qw_total + Qw_wellbore[k];
				*Qg_total = *Qg_total + Qg_wellbore[k];
				//cout << "Qo_wellbore[" << k << "]: " << Qo_wellbore[k] << endl;
			}

			BHP_r = Pa_Psia(converge_p_sup(L, nz, ft_m(2.l * Rw), Psia_Pa(*BHP), Far_Kel(T), (*Qw_total + *Qo_total),
				(*Qg_total / *Qo_total), API, yg, (*Qw_total) / (*Qw_total + *Qo_total), THP, pi, rug, coup_ts, count_W_R, W_corre, directory));


			count_W_R++;
			avoidshit = avoidshit + 1;

			bhp << *BHP << "\t;" << BHP_r << "\t;" << (BHP_r / *BHP) << "\t;" << *Qo_total << "\t;" << endl;

		}

		for (k = 0; k < N_wellbore; k++) {
			deltaz[k] = ft_m(deltaz[k]);
			Qw_wellbore[k] = Qw_wellbore[k] / deltaz[k];
			Qo_wellbore[k] = Qo_wellbore[k] / deltaz[k];
			Qg_wellbore[k] = Qg_wellbore[k] / deltaz[k];
		}

		wellbore_multi(WB_corre, N_wellbore,
			PWB, Pb_WB, H_l,
			Qw_wellbore, Qo_wellbore, Qg_wellbore,
			*BHP, PWB_1, deltaz,
			2.l * Rw, D_perf, perf_dens,
			API, yo, yg,
			rho_os, rho_ws, rho_gs, T,
			g, teta, rug, gama,
			parametrosmig, parametrosBw,
			Tpc, ppc, pi, Bwi);

		for (k = 0; k < N_wellbore; k++) {
			Qw_wellbore[k] = Qw_wellbore[k] * deltaz[k];
			Qo_wellbore[k] = Qo_wellbore[k] * deltaz[k];
			Qg_wellbore[k] = Qg_wellbore[k] * deltaz[k];
			deltaz[k] = m_ft(deltaz[k]);
		}

		//M�TODO NUM�RICO 
		if (ite_coup == 1) {

			PWB_resp = PWB[N_wellbore - 1];
			PWB_2 = PWB_1;
			f2 = 1.0l - PWB_resp / PWB_1;
			PWB_1 = 0.5l * (PWB_1 + PWB_resp);

		}
		else {
			PWB_resp = PWB[N_wellbore - 1];
			f1 = 1.0l - PWB_resp / PWB_1;
			awb = (f1 - f2) / (PWB_1 - PWB_2);

			if (PWB_1 == PWB_2 || abs(f1 / (awb)) > 500.l || (PWB_1 - f1 / awb) > 1.05l * po_coup[0] || isnan(PWB_1 - f1 / awb)) {
				PWB_1 = 0.5l * (PWB_1 + PWB_resp);
			}
			else {
				PWB_2 = PWB_1;
				f2 = f1;
				PWB_1 = PWB_1 - f1 / awb;
				if (isnan(PWB_1)) {
					cout << "PWB_1: " << PWB_1 << endl;
					system("pause");
				}
			}
		}
		//

		//Extract the shape of PWB CURVE
		for (k = N_wellbore - 2; k >= 0; k--) {

			dPWB_D[k] = (PWB[k] - PWB[k + 1]) / (*BHP - PWB_resp);

		}

		error_WB_R = 0.0l;
		*Qo_total = 0.0l;
		*Qw_total = 0.0l;
		*Qg_total = 0.0l;
		for (k = N_wellbore - 1; k >= 0; k--) {


			if (abs((Qo_wellbore_ref[k] - Qo_wellbore[k]) / Qo_wellbore_ref[k]) > error_WB_R) {

				error_WB_R = abs((Qo_wellbore_ref[k] - Qo_wellbore[k]) / Qo_wellbore_ref[k]);

			}
			if (abs((Qw_wellbore_ref[k] - Qw_wellbore[k]) / Qw_wellbore_ref[k]) > error_WB_R) {

				error_WB_R = abs((Qw_wellbore_ref[k] - Qw_wellbore[k]) / Qw_wellbore_ref[k]);

			}
			if (abs((Qg_wellbore_ref[k] - Qg_wellbore[k]) / Qg_wellbore_ref[k]) > error_WB_R) {

				error_WB_R = abs((Qg_wellbore_ref[k] - Qg_wellbore[k]) / Qg_wellbore_ref[k]);

			}
			if (abs((PWB_ref[k] - PWB[k]) / PWB[k]) > error_WB_R) {

				error_WB_R = abs((PWB_ref[k] - PWB[k]) / PWB[k]);

			}

			if (k == N_wellbore - 1) {

				PWB[k] = PWB_1;

			}
			else {

				PWB[k] = PWB[k + 1] + (*BHP - PWB_1) * dPWB_D[k];

			}

			*Qo_total = *Qo_total + Qo_wellbore[k];
			*Qw_total = *Qw_total + Qw_wellbore[k];
			*Qg_total = *Qg_total + Qg_wellbore[k];

			Qo_wellbore_ref[k] = Qo_wellbore[k];
			Qw_wellbore_ref[k] = Qw_wellbore[k];
			Qg_wellbore_ref[k] = Qg_wellbore[k];
			PWB_ref[k] = PWB[k];


		}

		if (coup_ts % 50 == 0) {
			grava_var(PWB, N_wellbore, 0, count_W_R, ite_coup, subdir, 5);
			grava_var(Qo_wellbore, N_wellbore, 0, count_W_R, ite_coup, subdir, 6);
			grava_var(Qg_wellbore, N_wellbore, 0, count_W_R, ite_coup, subdir, 7);
		}

		ite_coup++;

	}

	bhp.close();

	p = 0;
	for (k = 0; k <= (o - 1); k++) {
		if (zonaprodutora[k]) PWB_coup[k] = PWB[k];
		for (j = 0; j <= (m - 1); j++) {
			for (i = 0; i <= (n - 1); i++) {
				po_coup[p] = po_res_0[p];
				Sw_coup[p] = Sw_res_0[p];
				So_coup[p] = So_res_0[p];
				pb_coup[p] = pb_res_0[p];
				p++;
			}
		}
	}

	for (k = 0; k < N_wellbore; k++) {
		PWB_0[k] = PWB_coup[k];
	}


	free(Pb_WB);
	free(H_l);

}


void coup_step_2(int coup_ts, int ite_coup,
	double BHP_r, double* BHP,
	double THP, double L, int nz, double* deltaz,
	double* Qw_total, double* Qo_total, double* Qg_total,
	double* PWB_ref, double* PWB, double* PWB_0, double* dPWB_D, double* PWB_coup,
	double* Qw_wellbore, double* Qo_wellbore, double* Qg_wellbore,
	double* Qw_wellbore_ref, double* Qo_wellbore_ref, double* Qg_wellbore_ref,
	int N_wellbore, double D_perf, double perf_dens,
	double teta, double rug, double gama, double teta_yue,
	double* po_coup, double* So_coup, double* Sw_coup, double* pb_coup,
	double* po_res, double* pw_res, double* pg_res,
	double* So_res, double* Sw_res, double* Sg_res, double* pb_res,
	double* po_res_0, double* So_res_0, double* Sw_res_0, double* pb_res_0, double* Sg_res_0, double* X_pb_0,
	double* phi_res_0, double* Bw_res_0, double* Rs_res_0, double* Bo_res_0, double* Bg_res_0,
	double* po_res_lv0, double* pb_res_lv0, double* So_res_lv0, double* Sw_res_lv0,
	double* varpo, double* varpb, double* varSo, double* varSw,
	double Qwsinj, double Qosinj, double Qgsinj,
	double* K_F, double* K_B, double* K_N, double* K_S, double* K_W, double* K_E,
	double* X_pb, bool* zonaprodutora,
	int res_timesteps, bool dontstopmenow,
	int p, int y, int u, int l,
	int centerP, int centerSw, int centerSo,
	int i, int k, int j,
	int o, int m, int n, int dim,
	double* Kr, double* Kz, double pi, double T,
	int* IA, int* JA, double* A, double* B, double* X,
	double controladorP, double controladorSw, double controladorSo,
	double C1, double C2, double C3,
	double* Ageo_R, double* Ageo_PHI, double* Ageo_Z, double* Ageo_t,
	double* Adir_F, double* Adir_B, double* Adir_W, double* Adir_E, double* Adir_N, double* Adir_S,
	double* CHI_F, double* CHI_B, double* CHI_W, double* CHI_E, double* CHI_N, double* CHI_S,
	double* R_P, double* PHI_P, double* Z_P,
	double* R_PB, double* R_PF, double* delta_PHI, double* delta_Z, double deltat_res,
	double g, double Re, double Rw,
	double rho_air, double rho_ws, double API, double yo, double rho_os, double yg, double rho_gs,
	double yn2, double yco2, double yh2s,
	double* lambda1_W, double* d_lambda1_po_W, double* d_lambda1_SW_W,
	double* lambda2_W, double* d_lambda2_po_W, double* d_lambda2_SW_W,
	double* lambda1_O, double* d_lambda1_po_O, double* d_lambda1_SW_O, double* d_lambda1_pb_O, double* d_lambda1_SO_O,
	double* lambda2_O, double* d_lambda2_po_O, double* d_lambda2_SW_O, double* d_lambda2_pb_O, double* d_lambda2_SO_O,
	double* lambda1_G, double* d_lambda1_po_G, double* d_lambda1_SW_G, double* d_lambda1_pb_G, double* d_lambda1_SO_G,
	double* lambda2_G, double* d_lambda2_po_G, double* d_lambda2_SW_G, double* d_lambda2_pb_G, double* d_lambda2_SO_G,
	double* lambda3_G, double* d_lambda3_po_G, double* d_lambda3_SW_G, double* d_lambda3_pb_G, double* d_lambda3_SO_G,
	double* lambda4_G, double* d_lambda4_po_G, double* d_lambda4_SW_G, double* d_lambda4_pb_G, double* d_lambda4_SO_G,
	double krwmax, double nw, double ng, double now, double nog,
	double* phi_ref, double Cr,
	double Swc, double Sor, double Sgr, double epsilon,
	double pcwmin, double pcwmax, double pcgmin, double pcgmax,
	double* parametrosBw, double* parametrosmig, double* parametrosZ,
	double ppc, double Tpc, double Bwi, double Bpcw, double Bpcg,
	int* iparm, double* dparm, int solver, int maxfct, int mnum, int phase,
	int error, int msglvl, int mtype, int nrhs, double ddum, int idum, void* pt,
	double* phi, double* der_phi,
	double* rho_w, double* der_rho_w,
	double* mi_w, double* der_mi_w,
	double* Bw, double* der_Bw,
	double* pcw, double* der_pcw,
	double* krw, double* der_krw,
	double* Rs, double* der_Rs_p, double* der_Rs_pb,
	double* kro, double* der_kro_Sw, double* der_kro_So,
	double* Bo, double* der_Bo_p, double* der_Bo_pb,
	double* mi_o, double* der_mi_o_p, double* der_mi_o_pb,
	double* rho_o, double* der_rho_o_p, double* der_rho_o_pb,
	double* pcg, double* der_pcg_Sw, double* der_pcg_So,
	double* krg, double* der_krg,
	double* mi_g, double* der_mi_g, double* Z,
	double* Bg, double* der_Bg,
	double* rho_g, double* der_rho_g,
	int* unW, int* ufW, int* uwbW, int* ueW, int* ucW, int* uwW, int* uebW, int* ubW, int* usW,
	int* unO, int* ufO, int* uwbO, int* ueO, int* ucO, int* uwO, int* uebO, int* ubO, int* usO,
	int* unG, int* ufG, int* uwbG, int* ueG, int* ucG, int* uwG, int* uebG, int* ubG, int* usG,
	int* neg, bool* NpureW, int* contresint,
	char directory[], __int8 W_corre, __int8 WB_corre) {

	double* Pb_WB = (double*)calloc(N_wellbore + 1, sizeof(double));
	double* H_l = (double*)calloc(N_wellbore + 1, sizeof(double));

	double deltaBHP = deltaBHP = 50.0l / log(coup_ts + 2);
	bool controldelta = 1;
	int avoidshit = 0;

	double indice1 = 0.0;
	double BHP1 = 0.0;
	double indice2 = 0.0;
	double BHP2 = 0.0;

	double error_WB_R = 100.l;
	double PWB_1, PWB_2, PWB_resp, f1, f2, awb;

	ite_coup = 1;
	BHP_r = 0;
	error_WB_R = 100.l;
	PWB_1 = PWB[N_wellbore - 1];
	*contresint = 0;

	char subdir[50], subdirname[100];
	memset(&subdir[0], 0, sizeof(subdir));
	memset(&subdirname[0], 0, sizeof(subdirname));

	sprintf_s(subdir, directory);
	sprintf_s(subdirname, "\\TimeStep_%d", coup_ts);

	strcat_s(subdir, subdirname);
	_mkdir(subdir);

	char arquive[120];
	sprintf_s(arquive, subdir);
	strcat_s(arquive, "\\BHP_x_BHP_r.m");

	fstream bhp(arquive, ios::out | ios::app);
	bhp << "BHP_guess[psi];" << "BHP_guess[psi];" << "BHPg/BHPr;" << "Qo[bbl/d];" << endl;

	while (fabs((*BHP - BHP_r) / *BHP) > 0.001l || error_WB_R > 0.001l) {

		calc_guess_BHP(ite_coup, coup_ts,
			BHP_r, BHP, &BHP1, &BHP2, &indice1, &indice2,
			&deltaBHP, &controldelta, &avoidshit, po_coup[0], THP);

		for (k = N_wellbore - 1; k >= 0; k--) {

			if (k == N_wellbore - 1) {

				PWB[k] = PWB_1;

			}
			else {

				PWB[k] = PWB[k + 1] + (*BHP - PWB_1) * dPWB_D[k];

			}

			PWB_ref[k] = PWB[k];

			if (PWB[k] > po_coup[0]) {
				//cout << "PWB[" << k << "]:" << PWB[k] << endl;			
			}

		}

		p = 0;
		for (k = 0; k <= (o - 1); k++) {
			for (j = 0; j <= (m - 1); j++) {
				for (i = 0; i <= (n - 1); i++) {
					po_res[p] = po_coup[p];
					Sw_res[p] = Sw_coup[p];
					So_res[p] = So_coup[p];
					pb_res[p] = pb_coup[p];
					Sg_res[p] = 1 - So_coup[p] - Sw_coup[p];
					pw_res[p] = po_coup[p];// +calc_water_pcw_Sw(Sw[p], Swc, epsilon, pcwmin, pcwmax, Bpcw);
					pg_res[p] = po_coup[p];// +calc_gas_pcg(Sw[p], So[p], Swc, Sor, epsilon, pcgmin, pcgmax, Bpcg);
					po_res_0[p] = po_coup[p];
					pb_res_0[p] = pb_coup[p];
					Sw_res_0[p] = Sw_coup[p];
					So_res_0[p] = So_coup[p];
					Sg_res_0[p] = 1 - So_coup[p] - Sw_coup[p];
					p++;
				}
			}
		}

		for (k = 0; k < N_wellbore; k++) {
			PWB_0[k] = PWB_coup[k];
		}

		calc_reser_evol(res_timesteps, p, y, u, l, centerP, centerSw, centerSo,
			IA, JA, A, B, X,
			dontstopmenow, i, k, j, o, m, n, dim,
			controladorP, controladorSw, controladorSo,
			lambda1_W, d_lambda1_po_W, d_lambda1_SW_W,
			lambda2_W, d_lambda2_po_W, d_lambda2_SW_W,
			lambda1_O, d_lambda1_po_O, d_lambda1_SW_O, d_lambda1_pb_O, d_lambda1_SO_O,
			lambda2_O, d_lambda2_po_O, d_lambda2_SW_O, d_lambda2_pb_O, d_lambda2_SO_O,
			lambda1_G, d_lambda1_po_G, d_lambda1_SW_G, d_lambda1_pb_G, d_lambda1_SO_G,
			lambda2_G, d_lambda2_po_G, d_lambda2_SW_G, d_lambda2_pb_G, d_lambda2_SO_G,
			lambda3_G, d_lambda3_po_G, d_lambda3_SW_G, d_lambda3_pb_G, d_lambda3_SO_G,
			lambda4_G, d_lambda4_po_G, d_lambda4_SW_G, d_lambda4_pb_G, d_lambda4_SO_G,
			Kr, Kz, pi, T,
			PWB, PWB_0, Pb_WB, N_wellbore,
			po_res, pb_res, pw_res, pg_res, Sw_res, So_res, Sg_res,
			po_res_0, pb_res_0, Sw_res_0, So_res_0, Sg_res_0, X_pb_0,
			phi_res_0, Bw_res_0, Rs_res_0, Bo_res_0, Bg_res_0,
			po_res_lv0, pb_res_lv0, So_res_lv0, Sw_res_lv0,
			varpo, varpb, varSo, varSw,
			C1, C2, C3,
			Ageo_R, Ageo_PHI, Ageo_Z, Ageo_t,
			Adir_F, Adir_B, Adir_W, Adir_E, Adir_N, Adir_S,
			CHI_F, CHI_B, CHI_W, CHI_E, CHI_N, CHI_S,
			R_P, PHI_P, Z_P,
			R_PB, R_PF, delta_PHI, delta_Z, deltat_res,
			g, Re, Rw,
			rho_air, rho_ws, rho_os, rho_gs,
			API, yg, yo,
			phi_ref, Cr,
			Qwsinj, Qosinj, Qgsinj,
			Qw_wellbore, Qo_wellbore, Qg_wellbore, coup_ts, ite_coup,
			phi, der_phi,
			rho_w, der_rho_w,
			mi_w, der_mi_w,
			Bw, der_Bw,
			pcw, der_pcw,
			krw, der_krw,
			Rs, der_Rs_p, der_Rs_pb,
			kro, der_kro_Sw, der_kro_So,
			Bo, der_Bo_p, der_Bo_pb,
			mi_o, der_mi_o_p, der_mi_o_pb,
			rho_o, der_rho_o_p, der_rho_o_pb,
			pcg, der_pcg_Sw, der_pcg_So,
			krg, der_krg, mi_g, der_mi_g,
			Bg, der_Bg, rho_g, der_rho_g,
			K_F, K_B, K_N, K_S, K_W, K_E,
			X_pb, zonaprodutora, 0,
			unW, ufW, uwbW, ueW, ucW, uwW, uebW, ubW, usW,
			unO, ufO, uwbO, ueO, ucO, uwO, uebO, ubO, usO,
			unG, ufG, uwbG, ueG, ucG, uwG, uebG, ubG, usG,
			iparm, dparm, solver, maxfct, mnum, phase, error,
			msglvl, mtype, nrhs, ddum, idum, pt,
			Swc, Sor, Sgr, epsilon,
			pcwmin, pcwmax, pcgmin, pcgmax,
			krwmax, nw, ng, nog, now,
			parametrosBw, parametrosmig, ppc, Tpc,
			Bwi, Bpcw, Bpcg,
			yn2, yco2, yh2s,
			neg, NpureW,directory);

		*contresint = *contresint + 1;

		*Qo_total = 0.0l;
		*Qw_total = 0.0l;
		*Qg_total = 0.0l;
		for (k = 0; k < N_wellbore; k++) {

			*Qo_total = *Qo_total + Qo_wellbore[k];
			*Qw_total = *Qw_total + Qw_wellbore[k];
			*Qg_total = *Qg_total + Qg_wellbore[k];
			//cout << "Qo_wellbore[" << k << "]: " << Qo_wellbore[k] << endl;
		}

		BHP_r = Pa_Psia(converge_p_sup(L, nz, ft_m(2.l * Rw), Psia_Pa(*BHP), Far_Kel(T), (*Qw_total + *Qo_total),
			(*Qg_total / *Qo_total), API, yg, (*Qw_total) / (*Qw_total + *Qo_total), THP, pi, rug, coup_ts, ite_coup, W_corre, directory));

		for (k = 0; k < N_wellbore; k++) {
			deltaz[k] = ft_m(deltaz[k]);
			Qw_wellbore[k] = Qw_wellbore[k] / deltaz[k];
			Qo_wellbore[k] = Qo_wellbore[k] / deltaz[k];
			Qg_wellbore[k] = Qg_wellbore[k] / deltaz[k];
		}

		wellbore_multi(WB_corre, N_wellbore,
			PWB, Pb_WB, H_l,
			Qw_wellbore, Qo_wellbore, Qg_wellbore,
			BHP_r, PWB_1, deltaz,
			2.l * Rw, D_perf, perf_dens,
			API, yo, yg,
			rho_os, rho_ws, rho_gs, T,
			g, teta, rug, gama,
			parametrosmig, parametrosBw,
			Tpc, ppc, pi, Bwi);

		for (k = 0; k < N_wellbore; k++) {
			Qw_wellbore[k] = Qw_wellbore[k] * deltaz[k];
			Qo_wellbore[k] = Qo_wellbore[k] * deltaz[k];
			Qg_wellbore[k] = Qg_wellbore[k] * deltaz[k];
			deltaz[k] = m_ft(deltaz[k]);
		}

		//M�TODO NUM�RICO 
		if (ite_coup == 1) {

			PWB_resp = PWB[N_wellbore - 1];
			PWB_2 = PWB_1;
			f2 = 1.0l - PWB_resp / PWB_1;
			PWB_1 = 0.5l * (PWB_1 + PWB_resp);

		}
		else {
			PWB_resp = PWB[N_wellbore - 1];
			f1 = 1.0l - PWB_resp / PWB_1;
			awb = (f1 - f2) / (PWB_1 - PWB_2);

			if (PWB_1 == PWB_2 || abs(f1 / (awb)) > 500.l || (PWB_1 - f1 / awb) > 1.05l * po_coup[0] || isnan(PWB_1 - f1 / awb)) {
				PWB_1 = 0.5l * (PWB_1 + PWB_resp);
			}
			else {
				PWB_2 = PWB_1;
				f2 = f1;
				PWB_1 = PWB_1 - f1 / awb;
				if (isnan(PWB_1)) {
					cout << "PWB_1: " << PWB_1 << endl;
					system("pause");
				}
			}
		}
		//

		//Extract the shape of PWB CURVE
		for (k = N_wellbore - 2; k >= 0; k--) {

			dPWB_D[k] = (PWB[k] - PWB[k + 1]) / (BHP_r - PWB_resp);
			if (isnan(dPWB_D[k])) {
				cout << "*BHP - PWB_resp: " << BHP_r - PWB_resp << endl;
				//system("pause");
			}

		}

		*Qo_total = 0.0l;
		*Qw_total = 0.0l;
		*Qg_total = 0.0l;
		error_WB_R = 0.0l;
		for (k = N_wellbore - 1; k >= 0; k--) {

			if (abs((Qo_wellbore_ref[k] - Qo_wellbore[k]) / Qo_wellbore_ref[k]) > error_WB_R) {

				error_WB_R = abs((Qo_wellbore_ref[k] - Qo_wellbore[k]) / Qo_wellbore_ref[k]);

			}
			if (abs((Qw_wellbore_ref[k] - Qw_wellbore[k]) / Qw_wellbore_ref[k]) > error_WB_R) {

				error_WB_R = abs((Qw_wellbore_ref[k] - Qw_wellbore[k]) / Qw_wellbore_ref[k]);

			}
			if (abs((Qg_wellbore_ref[k] - Qg_wellbore[k]) / Qg_wellbore_ref[k]) > error_WB_R) {

				error_WB_R = abs((Qg_wellbore_ref[k] - Qg_wellbore[k]) / Qg_wellbore_ref[k]);

			}
			if (abs((PWB_ref[k] - PWB[k]) / PWB[k]) > error_WB_R) {

				error_WB_R = abs((PWB_ref[k] - PWB[k]) / PWB[k]);

			}


			*Qo_total = *Qo_total + Qo_wellbore[k];
			*Qw_total = *Qw_total + Qw_wellbore[k];
			*Qg_total = *Qg_total + Qg_wellbore[k];

			Qo_wellbore_ref[k] = Qo_wellbore[k];
			Qw_wellbore_ref[k] = Qw_wellbore[k];
			Qg_wellbore_ref[k] = Qg_wellbore[k];
			PWB_ref[k] = PWB[k];

		}

		bhp << *BHP << "\t;" << BHP_r << "\t;" << (BHP_r / *BHP) << "\t;" << *Qo_total << "\t;" << endl;

		if (coup_ts % 50 == 0) {
			grava_var(PWB, N_wellbore, 0, 0, ite_coup, subdir, 5);
			grava_var(Qo_wellbore, N_wellbore, 0, 0, ite_coup, subdir, 6);
			grava_var(Qg_wellbore, N_wellbore, 0, 0, ite_coup, subdir, 7);
		}

		ite_coup++;
		avoidshit = avoidshit + 1;


	}

	bhp.close();

	p = 0;
	for (k = 0; k <= (o - 1); k++) {
		for (j = 0; j <= (m - 1); j++) {
			for (i = 0; i <= (n - 1); i++) {
				po_coup[p] = po_res_0[p];
				Sw_coup[p] = Sw_res_0[p];
				So_coup[p] = So_res_0[p];
				pb_coup[p] = pb_res_0[p];
				p++;
			}
		}
	}

	for (k = 0; k < N_wellbore; k++) {
		PWB_0[k] = PWB[k];
		PWB_coup[k] = PWB[k];
	}

	free(Pb_WB);
	free(H_l);

}

void coup_step_3(int coup_ts, int ite_coup,
	double BHP_r, double* BHP,
	double THP, double L, int nz, double* deltaz,
	double* Qw_total, double* Qo_total, double* Qg_total,
	double* PWB_ref, double* PWB, double* PWB_0, double* dPWB_D, double* PWB_coup,
	double* Qw_wellbore, double* Qo_wellbore, double* Qg_wellbore,
	double* Qw_wellbore_ref, double* Qo_wellbore_ref, double* Qg_wellbore_ref,
	int N_wellbore, double D_perf, double perf_dens,
	double teta, double rug, double gama, double teta_yue,
	double* po_coup, double* So_coup, double* Sw_coup, double* pb_coup,
	double* po_res, double* pw_res, double* pg_res,
	double* So_res, double* Sw_res, double* Sg_res, double* pb_res,
	double* po_res_0, double* So_res_0, double* Sw_res_0, double* pb_res_0, double* Sg_res_0, double* X_pb_0,
	double* phi_res_0, double* Bw_res_0, double* Rs_res_0, double* Bo_res_0, double* Bg_res_0,
	double* po_res_lv0, double* pb_res_lv0, double* So_res_lv0, double* Sw_res_lv0,
	double* varpo, double* varpb, double* varSo, double* varSw,
	double Qwsinj, double Qosinj, double Qgsinj,
	double* K_F, double* K_B, double* K_N, double* K_S, double* K_W, double* K_E,
	double* X_pb, bool* zonaprodutora,
	int res_timesteps, bool dontstopmenow,
	int p, int y, int u, int l,
	int centerP, int centerSw, int centerSo,
	int i, int k, int j,
	int o, int m, int n, int dim,
	double* Kr, double* Kz, double pi, double T,
	int* IA, int* JA, double* A, double* B, double* X,
	double controladorP, double controladorSw, double controladorSo,
	double C1, double C2, double C3,
	double* Ageo_R, double* Ageo_PHI, double* Ageo_Z, double* Ageo_t,
	double* Adir_F, double* Adir_B, double* Adir_W, double* Adir_E, double* Adir_N, double* Adir_S,
	double* CHI_F, double* CHI_B, double* CHI_W, double* CHI_E, double* CHI_N, double* CHI_S,
	double* R_P, double* PHI_P, double* Z_P,
	double* R_PB, double* R_PF, double* delta_PHI, double* delta_Z, double deltat_res,
	double g, double Re, double Rw,
	double rho_air, double rho_ws, double API, double yo, double rho_os, double yg, double rho_gs,
	double yn2, double yco2, double yh2s,
	double* lambda1_W, double* d_lambda1_po_W, double* d_lambda1_SW_W,
	double* lambda2_W, double* d_lambda2_po_W, double* d_lambda2_SW_W,
	double* lambda1_O, double* d_lambda1_po_O, double* d_lambda1_SW_O, double* d_lambda1_pb_O, double* d_lambda1_SO_O,
	double* lambda2_O, double* d_lambda2_po_O, double* d_lambda2_SW_O, double* d_lambda2_pb_O, double* d_lambda2_SO_O,
	double* lambda1_G, double* d_lambda1_po_G, double* d_lambda1_SW_G, double* d_lambda1_pb_G, double* d_lambda1_SO_G,
	double* lambda2_G, double* d_lambda2_po_G, double* d_lambda2_SW_G, double* d_lambda2_pb_G, double* d_lambda2_SO_G,
	double* lambda3_G, double* d_lambda3_po_G, double* d_lambda3_SW_G, double* d_lambda3_pb_G, double* d_lambda3_SO_G,
	double* lambda4_G, double* d_lambda4_po_G, double* d_lambda4_SW_G, double* d_lambda4_pb_G, double* d_lambda4_SO_G,
	double krwmax, double nw, double ng, double now, double nog,
	double* phi_ref, double Cr,
	double Swc, double Sor, double Sgr, double epsilon,
	double pcwmin, double pcwmax, double pcgmin, double pcgmax,
	double* parametrosBw, double* parametrosmig, double* parametrosZ,
	double ppc, double Tpc, double Bwi, double Bpcw, double Bpcg,
	int* iparm, double* dparm, int solver, int maxfct, int mnum, int phase,
	int error, int msglvl, int mtype, int nrhs, double ddum, int idum, void* pt,
	double* phi, double* der_phi,
	double* rho_w, double* der_rho_w,
	double* mi_w, double* der_mi_w,
	double* Bw, double* der_Bw,
	double* pcw, double* der_pcw,
	double* krw, double* der_krw,
	double* Rs, double* der_Rs_p, double* der_Rs_pb,
	double* kro, double* der_kro_Sw, double* der_kro_So,
	double* Bo, double* der_Bo_p, double* der_Bo_pb,
	double* mi_o, double* der_mi_o_p, double* der_mi_o_pb,
	double* rho_o, double* der_rho_o_p, double* der_rho_o_pb,
	double* pcg, double* der_pcg_Sw, double* der_pcg_So,
	double* krg, double* der_krg,
	double* mi_g, double* der_mi_g, double* Z,
	double* Bg, double* der_Bg,
	double* rho_g, double* der_rho_g,
	int* unW, int* ufW, int* uwbW, int* ueW, int* ucW, int* uwW, int* uebW, int* ubW, int* usW,
	int* unO, int* ufO, int* uwbO, int* ueO, int* ucO, int* uwO, int* uebO, int* ubO, int* usO,
	int* unG, int* ufG, int* uwbG, int* ueG, int* ucG, int* uwG, int* uebG, int* ubG, int* usG,
	int* neg, bool* NpureW, int* contresint,
	char directory[], __int8 W_corre, __int8 WB_corre) {

	double* Pb_WB = (double*)calloc(N_wellbore, sizeof(double));
	double* H_l = (double*)calloc(N_wellbore, sizeof(double));

	double deltaBHP = 50.0l / log(coup_ts + 2);
	bool controldelta = 1;
	int avoidshit = 0;

	double indice1 = 0.0;
	double BHP1 = 0.0;
	double indice2 = 0.0;
	double BHP2 = 0.0;

	double error_WB_R = 100.l;
	double PWB_1, PWB_2, PWB_resp, f1, f2, awb;

	ite_coup = 1;
	BHP_r = 0;
	error_WB_R = 100.l;
	PWB_1 = PWB[N_wellbore - 1];
	*contresint = 0;

	char subdir[50], subdirname[100];
	memset(&subdir[0], 0, sizeof(subdir));
	memset(&subdirname[0], 0, sizeof(subdirname));

	sprintf_s(subdir, directory);
	sprintf_s(subdirname, "\\TimeStep_%d", coup_ts);

	strcat_s(subdir, subdirname);
	_mkdir(subdir);

	char arquive[120];
	sprintf_s(arquive, subdir);
	strcat_s(arquive, "\\BHP_x_BHP_r.m");

	fstream bhp(arquive, ios::out | ios::app);
	bhp << "BHP_guess[psi];" << "BHP_guess[psi];" << "BHPg/BHPr;" << "Qo[bbl/d];" << endl;	

	while (fabs((*BHP - BHP_r) / *BHP) > 0.01l || error_WB_R > 0.01l) {

		calc_guess_BHP(ite_coup, coup_ts,
			BHP_r, BHP, &BHP1, &BHP2, &indice1, &indice2,
			&deltaBHP, &controldelta, &avoidshit, po_coup[0], THP);

		for (k = N_wellbore - 1; k >= 0; k--) {

			if (k == N_wellbore - 1) {

				PWB[k] = PWB_1;

			}
			else {

				PWB[k] = PWB[k + 1] + (*BHP - PWB_1) * dPWB_D[k];

			}

			PWB_ref[k] = PWB[k];
			if (PWB[k] > po_coup[0]) {
				cout << "PWB[" << k << "]:" << PWB[k] << endl;
			}

		}

		p = 0;
		for (k = 0; k <= (o - 1); k++) {
			for (j = 0; j <= (m - 1); j++) {
				for (i = 0; i <= (n - 1); i++) {
					po_res[p] = po_coup[p];
					Sw_res[p] = Sw_coup[p];
					So_res[p] = So_coup[p];
					pb_res[p] = pb_coup[p];
					Sg_res[p] = 1 - So_coup[p] - Sw_coup[p];
					pw_res[p] = po_coup[p];// +calc_water_pcw_Sw(Sw[p], Swc, epsilon, pcwmin, pcwmax, Bpcw);
					pg_res[p] = po_coup[p];// +calc_gas_pcg(Sw[p], So[p], Swc, Sor, epsilon, pcgmin, pcgmax, Bpcg);
					po_res_0[p] = po_coup[p];
					pb_res_0[p] = pb_coup[p];
					Sw_res_0[p] = Sw_coup[p];
					So_res_0[p] = So_coup[p];
					Sg_res_0[p] = 1 - So_coup[p] - Sw_coup[p];
					p++;
				}
			}
		}

		for (k = 0; k < N_wellbore; k++) {
			if (zonaprodutora[k]) {
				PWB_0[k] = PWB_coup[k];
			}
		}

		calc_reser_evol(res_timesteps, p, y, u, l, centerP, centerSw, centerSo,
			IA, JA, A, B, X,
			dontstopmenow, i, k, j, o, m, n, dim,
			controladorP, controladorSw, controladorSo,
			lambda1_W, d_lambda1_po_W, d_lambda1_SW_W,
			lambda2_W, d_lambda2_po_W, d_lambda2_SW_W,
			lambda1_O, d_lambda1_po_O, d_lambda1_SW_O, d_lambda1_pb_O, d_lambda1_SO_O,
			lambda2_O, d_lambda2_po_O, d_lambda2_SW_O, d_lambda2_pb_O, d_lambda2_SO_O,
			lambda1_G, d_lambda1_po_G, d_lambda1_SW_G, d_lambda1_pb_G, d_lambda1_SO_G,
			lambda2_G, d_lambda2_po_G, d_lambda2_SW_G, d_lambda2_pb_G, d_lambda2_SO_G,
			lambda3_G, d_lambda3_po_G, d_lambda3_SW_G, d_lambda3_pb_G, d_lambda3_SO_G,
			lambda4_G, d_lambda4_po_G, d_lambda4_SW_G, d_lambda4_pb_G, d_lambda4_SO_G,
			Kr, Kz, pi, T,
			PWB, PWB_0, Pb_WB, N_wellbore,
			po_res, pb_res, pw_res, pg_res, Sw_res, So_res, Sg_res,
			po_res_0, pb_res_0, Sw_res_0, So_res_0, Sg_res_0, X_pb_0,
			phi_res_0, Bw_res_0, Rs_res_0, Bo_res_0, Bg_res_0,
			po_res_lv0, pb_res_lv0, So_res_lv0, Sw_res_lv0,
			varpo, varpb, varSo, varSw,
			C1, C2, C3,
			Ageo_R, Ageo_PHI, Ageo_Z, Ageo_t,
			Adir_F, Adir_B, Adir_W, Adir_E, Adir_N, Adir_S,
			CHI_F, CHI_B, CHI_W, CHI_E, CHI_N, CHI_S,
			R_P, PHI_P, Z_P,
			R_PB, R_PF, delta_PHI, delta_Z, deltat_res,
			g, Re, Rw,
			rho_air, rho_ws, rho_os, rho_gs,
			API, yg, yo,
			phi_ref, Cr,
			Qwsinj, Qosinj, Qgsinj,
			Qw_wellbore, Qo_wellbore, Qg_wellbore, coup_ts, ite_coup,
			phi, der_phi,
			rho_w, der_rho_w,
			mi_w, der_mi_w,
			Bw, der_Bw,
			pcw, der_pcw,
			krw, der_krw,
			Rs, der_Rs_p, der_Rs_pb,
			kro, der_kro_Sw, der_kro_So,
			Bo, der_Bo_p, der_Bo_pb,
			mi_o, der_mi_o_p, der_mi_o_pb,
			rho_o, der_rho_o_p, der_rho_o_pb,
			pcg, der_pcg_Sw, der_pcg_So,
			krg, der_krg, mi_g, der_mi_g,
			Bg, der_Bg, rho_g, der_rho_g,
			K_F, K_B, K_N, K_S, K_W, K_E,
			X_pb, zonaprodutora, 0,
			unW, ufW, uwbW, ueW, ucW, uwW, uebW, ubW, usW,
			unO, ufO, uwbO, ueO, ucO, uwO, uebO, ubO, usO,
			unG, ufG, uwbG, ueG, ucG, uwG, uebG, ubG, usG,
			iparm, dparm, solver, maxfct, mnum, phase, error,
			msglvl, mtype, nrhs, ddum, idum, pt,
			Swc, Sor, Sgr, epsilon,
			pcwmin, pcwmax, pcgmin, pcgmax,
			krwmax, nw, ng, nog, now,
			parametrosBw, parametrosmig, ppc, Tpc,
			Bwi, Bpcw, Bpcg,
			yn2, yco2, yh2s,
			neg, NpureW, directory);

		*contresint = *contresint + 1;

		for (k = 0; k < N_wellbore; k++) {
			deltaz[k] = ft_m(deltaz[k]);
			Qw_wellbore[k] = Qw_wellbore[k] / deltaz[k];
			Qo_wellbore[k] = Qo_wellbore[k] / deltaz[k];
			Qg_wellbore[k] = Qg_wellbore[k] / deltaz[k];
		}

		wellbore_multi(WB_corre, N_wellbore,
			PWB, Pb_WB, H_l,
			Qw_wellbore, Qo_wellbore, Qg_wellbore,
			*BHP, PWB_1, deltaz,
			2.l * Rw, D_perf, perf_dens,
			API, yo, yg,
			rho_os, rho_ws, rho_gs, T,
			g, teta, rug, gama,
			parametrosmig, parametrosBw,
			Tpc, ppc, pi, Bwi);

		for (k = 0; k < N_wellbore; k++) {
			Qw_wellbore[k] = Qw_wellbore[k] * deltaz[k];
			Qo_wellbore[k] = Qo_wellbore[k] * deltaz[k];
			Qg_wellbore[k] = Qg_wellbore[k] * deltaz[k];
			deltaz[k] = m_ft(deltaz[k]);
		}


		//M�TODO NUM�RICO 
		if (ite_coup == 1) {

			PWB_resp = PWB[N_wellbore - 1];
			PWB_2 = PWB_1;
			f2 = 1.0l - PWB_resp / PWB_1;
			PWB_1 = 0.5l * (PWB_1 + PWB_resp);

		}
		else {
			PWB_resp = PWB[N_wellbore - 1];
			f1 = 1.0l - PWB_resp / PWB_1;
			awb = (f1 - f2) / (PWB_1 - PWB_2);

			if (PWB_1 == PWB_2 || abs(f1 / (awb)) > 500.l || (PWB_1 - f1 / awb) > 1.05l * po_coup[0] || isnan(PWB_1 - f1 / awb)) {
				PWB_1 = 0.5l * (PWB_1 + PWB_resp);
			}
			else {
				PWB_2 = PWB_1;
				f2 = f1;
				PWB_1 = PWB_1 - f1 / awb;
				if (isnan(PWB_1)) {
					cout << "PWB_1: " << PWB_1 << endl;
					//system("pause");
				}
			}
		}
		//


		//Extract the shape of PWB CURVE
		for (k = N_wellbore - 2; k >= 0; k--) {

			dPWB_D[k] = (PWB[k] - PWB[k + 1]) / (*BHP - PWB_resp);
			if (isnan(dPWB_D[k])) {
				cout << "*BHP - PWB_resp: " << *BHP - PWB_resp << endl;
				/*system("pause");*/
			}

		}

		*Qo_total = 0.0l;
		*Qw_total = 0.0l;
		*Qg_total = 0.0l;
		error_WB_R = 0.0l;
		for (k = N_wellbore - 1; k >= 0; k--) {

			if (abs((Qo_wellbore_ref[k] - Qo_wellbore[k]) / Qo_wellbore_ref[k]) > error_WB_R) {

				error_WB_R = abs((Qo_wellbore_ref[k] - Qo_wellbore[k]) / Qo_wellbore_ref[k]);


			}
			if (abs((Qw_wellbore_ref[k] - Qw_wellbore[k]) / Qw_wellbore_ref[k]) > error_WB_R) {

				error_WB_R = abs((Qw_wellbore_ref[k] - Qw_wellbore[k]) / Qw_wellbore_ref[k]);

			}
			if (abs((Qg_wellbore_ref[k] - Qg_wellbore[k]) / Qg_wellbore_ref[k]) > error_WB_R) {

				error_WB_R = abs((Qg_wellbore_ref[k] - Qg_wellbore[k]) / Qg_wellbore_ref[k]);

			}
			if (abs((PWB_ref[k] - PWB[k]) / PWB[k]) > error_WB_R) {

				error_WB_R = abs((PWB_ref[k] - PWB[k]) / PWB[k]);

			}

			//editei hein
			*Qo_total = *Qo_total + Qo_wellbore[k];
			*Qw_total = *Qw_total + Qw_wellbore[k];
			*Qg_total = *Qg_total + Qg_wellbore[k];

			Qo_wellbore_ref[k] = Qo_wellbore[k];
			Qw_wellbore_ref[k] = Qw_wellbore[k];
			Qg_wellbore_ref[k] = Qg_wellbore[k];
			PWB_ref[k] = PWB[k];

		}

		if (coup_ts % 50 == 0) {
			grava_var(PWB, N_wellbore, 0, 0, ite_coup, subdir, 5);
			grava_var(Qo_wellbore, N_wellbore, 0, 0, ite_coup, subdir, 6);
			grava_var(Qg_wellbore, N_wellbore, 0, 0, ite_coup, subdir, 7);
		}

		BHP_r = Pa_Psia(converge_p_sup(L, nz, ft_m(2.l * Rw), Psia_Pa(*BHP), Far_Kel(T), (*Qw_total + *Qo_total),
			(*Qg_total / *Qo_total), API, yg, (*Qw_total) / (*Qw_total + *Qo_total), THP, pi, rug, coup_ts, ite_coup, W_corre, directory));

		bhp << *BHP << "\t;" << BHP_r << "\t;" << (BHP_r / *BHP) << "\t;" << *Qo_total << "\t;" << endl;

		ite_coup++;
		avoidshit = avoidshit + 1;


	}

	bhp.close();


	p = 0;
	for (k = 0; k <= (o - 1); k++) {
		for (j = 0; j <= (m - 1); j++) {
			for (i = 0; i <= (n - 1); i++) {
				po_coup[p] = po_res_0[p];
				Sw_coup[p] = Sw_res_0[p];
				So_coup[p] = So_res_0[p];
				pb_coup[p] = pb_res_0[p];
				p++;
			}
		}
	}

	for (k = 0; k < N_wellbore; k++) {
		PWB_0[k] = PWB[k];
		PWB_coup[k] = PWB[k];
	}


	free(Pb_WB);
	free(H_l);

}

void calc_guess_BHP(int ite_coup, int coup_ts,
	double BHP_r, double* BHP, double* BHP1, double* BHP2, double* indice1, double* indice2,
	double* deltaBHP, bool* controldelta, int* avoidshit, double P20, double THP)
{

	//determina chute BHP
	if (ite_coup == 1) {

	}
	else {

		if (ite_coup > 1 && ite_coup % 2 == 0) {

			*indice1 = BHP_r / *BHP;
			*BHP1 = *BHP;

		}

		if (ite_coup > 1 && ite_coup % 2 == 1) {

			*indice2 = BHP_r / *BHP;
			*BHP2 = *BHP;

		}

		if ((ite_coup > 2 && coup_ts > 0) || ite_coup > 70) {


			double areta = (*indice1 - *indice2) / (*BHP1 - *BHP2);
			double breta = *indice1 - areta * *BHP1;
			*BHP = (1.0l - breta) / areta;

			if ((*BHP >= P20 || isnan(*BHP) || *BHP < (THP * 0.0001450377))) {
				if (isnan(*BHP) || *BHP < (THP * 0.0001450377)) {
					*BHP = 0.5 * (*BHP1 + *BHP2);
				}
				else {
					*BHP = 0.95l * P20 + rand() % 50 + 1;
					//cout << "P20: " << P20 << endl;
				}
				//system("pause");
			}

		}

		else {

			if (BHP_r > *BHP) {
				if (*controldelta) *deltaBHP = *deltaBHP / 2.0l;
				*BHP = *BHP + *deltaBHP;
				*controldelta = 0;
			}

			else {
				*BHP = *BHP - *deltaBHP;
				*controldelta = 1;
			}

		}

	}


	if (*avoidshit > 30 && coup_ts > 0) {
		*BHP = *BHP + 50.0;
		*avoidshit = 0;
	}

	if (*BHP >= P20 || isnan(*BHP) || *BHP < (THP * 0.0001450377)) {
		*BHP = P20 - 100 + rand() % 50 + 1;
	}

}

void calc_THP(double prod_time, int coup_ts, bool* controlstartup, double* THP, double THP_i, double THP_adj) {

	/*if (*controlstartup){
		*THP = Psia_Pa(780 - THP_adj);
	}
	else{

		if (coup_ts >= 0) *THP = Psia_Pa(max((-0.85979l*(prod_time / 60.l) / (1 + exp(-0.039084l*((prod_time / 60.l) - 120.53l))) + 887.52l), 780.l));

		if (Pa_Psia(*THP)<780.0001l) *controlstartup = 1;

	}*/

	* THP = Psia_Pa(THP_i);

}