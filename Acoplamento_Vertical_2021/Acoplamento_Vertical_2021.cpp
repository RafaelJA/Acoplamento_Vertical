// 3dreser.cpp : Defines the entry point for the console application.// acoplamento.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

int main()
{

	bool start = 1;
	__int8 W_corre = 1;
	__int8 WB_corre = 2;
	double deltat_coup_d = 0.01;//days

	//reser
	double Re_i = 2000;
	double D=3.5l;
	double p_i = 3600;
	double pb_i = 2000;
	double perm = 300.0l;
	double varperm = 0.0l;

	//dados do po�o
	double L = 2000;//Comprimento do poço em metros
	int nz = 501;//Numero de pontos do poço
	double THP_i = 780.l;

	double perf_dens_ft = 10.l;
	double D_perf_in = 0.18;
	double rug = 0.0002;//

	//---------------------------------------------------OK-----------------------------------------------------------------//

	/* Pardiso control parameters. */
	int      iparm[64];
	double   dparm[64];
	int      solver;
	int      maxfct, mnum, phase = 11, error, msglvl;
	int      mtype = 11;        /* Real unsymmetric matrix */
	int      nrhs = 1;          /* Number of right hand sides. */

	int      i = 0;

	double   ddum = 0;              /* Double dummy */
	int      idum = 0;              /* Integer dummy. */

	void* pt[64];

	error = 0;
	solver = 0; /* use sparse direct solver */

	pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);

	if (error != 0)
	{
		if (error == -10)
			printf("No license file found \n");
		if (error == -11)
			printf("License is expired \n");
		if (error == -12)
			printf("Wrong username or hostname \n");
		system("pause");
		return 1;
	}
	else
		printf("[PARDISO]: License check was successful ... \n");

	maxfct = 1;         /* Maximum number of numerical factorizations.  */
	mnum = 1;         /* Which factorization to use. */

	msglvl = 0;         /* Print statistical information  */
	error = 0;         /* Initialize error flag */

	iparm[2] = 1;//PROCESSADORES
	//---------------------------------------------------OK-----------------------------------------------------------------//

	double D_perf = ft_m(in_ft(D_perf_in));
	double perf_dens = m_ft(perf_dens_ft);//32.8 shots/m
	
	system("pause");

	double deltat_res, deltat_coup;

	int coup_ts_in = -1;
	int int_in = 0;
	bool controlstartup = 0;
	double THP_adj = 0;

	int X1 = 0;
	double* tempo_THP = (double*)malloc(20000 * sizeof(double));
	double* THP = (double*)malloc(20000 * sizeof(double));

	//Entrada transiente
	ifstream inFile("tabela_BHP.txt");

	while (!inFile.eof()) {
		inFile >> tempo_THP[X1];
		inFile >> THP[X1++];
	}

	setlocale(LC_ALL, "");

	//dados reservat�rio
	double H, Re, Rw;
	double* deltaz = (double*)malloc(29 * sizeof(double));
	int res_timesteps;

	double pi, pbi, Swi, Soi, T, g, C1, C2, C3;
	double Qwsinj, Qosinj, Qgsinj;
	double rho_air, rho_ws, rho_os, rho_gs, API, yg, yo, yn2, yco2, yh2s;
	double Swc, Sor, Sgr, epsilon, pcwmin, pcwmax, pcgmin, pcgmax;//pressao capilar

	double krwmax, nw, ng, now, nog;//permeabilidade relativa
	double phi_ref_0, Cr;//pedra
	int n, m, o;
	double parametrosBw[9];//fator volume de formacao
	double parametrosmig[25];//viscosidade
	double ppc, Tpc;
	double parametrosZ[7];
	double Bwi, Bpcw, Bpcg;
	bool* zonaprodutora = (bool*)malloc(50 * sizeof(bool));//indica se a zona � produtora	

	inicializa_reser(&res_timesteps,
		&Re, &Rw, deltaz,
		&n, &m, &o,
		&pi, &pbi,
		&T, &g,
		&C1, &C2, &C3,
		&rho_air, &rho_ws,
		&rho_os, &rho_gs,
		&yg, &yo, &API,
		&phi_ref_0, &Cr,
		&Swc, &Sor, &Sgr, &epsilon,
		&pcwmin, &pcwmax, &pcgmin, &pcgmax,
		&krwmax, &nw, &ng, &nog, &now,
		parametrosBw, parametrosmig, &ppc, &Tpc,
		&Bwi, &Bpcw, &Bpcg,
		&yn2, &yco2, &yh2s,
		zonaprodutora,
		Re_i, D, p_i, pb_i);

	//vetores PARDISO
	int dim = n * m * o;
	int* IA = (int*)malloc(3 * (dim + 10) * sizeof(int));
	int* JA = (int*)malloc(63 * (dim + 10) * sizeof(int));
	double* A = (double*)malloc(63 * (dim + 10) * sizeof(double));
	double* X = (double*)calloc(3 * (dim + 10), sizeof(double));
	double* B = (double*)malloc(3 * (dim + 10) * sizeof(double));

	//vetores propiedades(ant= passo anterior)
	double* po_res = (double*)malloc((dim + 10) * sizeof(double));
	double* pw_res = (double*)malloc((dim + 10) * sizeof(double));
	double* pg_res = (double*)malloc((dim + 10) * sizeof(double));
	double* pb_res = (double*)malloc((dim + 10) * sizeof(double));
	double* So_res = (double*)malloc((dim + 10) * sizeof(double));
	double* Sw_res = (double*)malloc((dim + 10) * sizeof(double));
	double* Sg_res = (double*)malloc((dim + 10) * sizeof(double));

	//passo anterior do reservatorio
	double* po_res_0 = (double*)malloc((dim + 10) * sizeof(double));
	double* pb_res_0 = (double*)malloc((dim + 10) * sizeof(double));
	double* So_res_0 = (double*)malloc((dim + 10) * sizeof(double));
	double* Sw_res_0 = (double*)malloc((dim + 10) * sizeof(double));
	double* Sg_res_0 = (double*)malloc((dim + 10) * sizeof(double));
	double* X_pb_0 = (double*)malloc((dim + 10) * sizeof(double));
	double* phi_res_0 = (double*)malloc((dim + 10) * sizeof(double));
	double* Bw_res_0 = (double*)malloc((dim + 10) * sizeof(double));
	double* Rs_res_0 = (double*)malloc((dim + 10) * sizeof(double));
	double* Bo_res_0 = (double*)malloc((dim + 10) * sizeof(double));
	double* Bg_res_0 = (double*)malloc((dim + 10) * sizeof(double));

	//passo anterior do acoplamento
	double* po_coup = (double*)malloc((dim + 10) * sizeof(double));
	double* pb_coup = (double*)malloc((dim + 10) * sizeof(double));
	double* So_coup = (double*)malloc((dim + 10) * sizeof(double));
	double* Sw_coup = (double*)malloc((dim + 10) * sizeof(double));

	//passo anterior do reservatorio no level 0
	double* po_res_lv0 = (double*)malloc((dim + 10) * sizeof(double));
	double* pb_res_lv0 = (double*)malloc((dim + 10) * sizeof(double));
	double* So_res_lv0 = (double*)malloc((dim + 10) * sizeof(double));
	double* Sw_res_lv0 = (double*)malloc((dim + 10) * sizeof(double));

	//varia�ao do passo anterior do reservatorio no level 0
	double* varpo = (double*)calloc((dim + 10), sizeof(double));
	double* varpb = (double*)calloc((dim + 10), sizeof(double));
	double* varSo = (double*)calloc((dim + 10), sizeof(double));
	double* varSw = (double*)calloc((dim + 10), sizeof(double));

	//propiedades calculadas
	double* phi_ref = (double*)malloc((dim + 10) * sizeof(double));
	double* phi = (double*)malloc((dim + 10) * sizeof(double));
	double* der_phi = (double*)malloc((dim + 10) * sizeof(double));
	double* rho_w = (double*)malloc((dim + 10) * sizeof(double));
	double* der_rho_w = (double*)malloc((dim + 10) * sizeof(double));
	double* mi_w = (double*)malloc((dim + 10) * sizeof(double));
	double* der_mi_w = (double*)malloc((dim + 10) * sizeof(double));
	double* Bw = (double*)malloc((dim + 10) * sizeof(double));
	double* der_Bw = (double*)malloc((dim + 10) * sizeof(double));
	double* pcw = (double*)malloc((dim + 10) * sizeof(double));
	double* der_pcw = (double*)malloc((dim + 10) * sizeof(double));
	double* krw = (double*)malloc((dim + 10) * sizeof(double));
	double* der_krw = (double*)malloc((dim + 10) * sizeof(double));
	double* Rs = (double*)malloc((dim + 10) * sizeof(double));
	double* der_Rs_p = (double*)malloc((dim + 10) * sizeof(double));
	double* der_Rs_pb = (double*)malloc((dim + 10) * sizeof(double));
	double* kro = (double*)malloc((dim + 10) * sizeof(double));
	double* der_kro_Sw = (double*)malloc((dim + 10) * sizeof(double));
	double* der_kro_So = (double*)malloc((dim + 10) * sizeof(double));
	double* Bo = (double*)malloc((dim + 10) * sizeof(double));
	double* der_Bo_p = (double*)malloc((dim + 10) * sizeof(double));
	double* der_Bo_pb = (double*)malloc((dim + 10) * sizeof(double));
	double* mi_o = (double*)malloc((dim + 10) * sizeof(double));
	double* der_mi_o_p = (double*)malloc((dim + 10) * sizeof(double));
	double* der_mi_o_pb = (double*)malloc((dim + 10) * sizeof(double));
	double* rho_o = (double*)malloc((dim + 10) * sizeof(double));
	double* der_rho_o_p = (double*)malloc((dim + 10) * sizeof(double));
	double* der_rho_o_pb = (double*)malloc((dim + 10) * sizeof(double));
	double* pcg = (double*)malloc((dim + 10) * sizeof(double));
	double* der_pcg_Sw = (double*)malloc((dim + 10) * sizeof(double));
	double* der_pcg_So = (double*)malloc((dim + 10) * sizeof(double));
	double* krg = (double*)malloc((dim + 10) * sizeof(double));
	double* der_krg = (double*)malloc((dim + 10) * sizeof(double));
	double* mi_g = (double*)malloc((dim + 10) * sizeof(double));
	double* der_mi_g = (double*)malloc((dim + 10) * sizeof(double));
	double* Z = (double*)malloc((dim + 10) * sizeof(double));
	double* Bg = (double*)malloc((dim + 10) * sizeof(double));
	double* der_Bg = (double*)malloc((dim + 10) * sizeof(double));
	double* rho_g = (double*)malloc((dim + 10) * sizeof(double));
	double* der_rho_g = (double*)malloc((dim + 10) * sizeof(double));

	double* Kr = (double*)malloc((dim + 10) * sizeof(double));
	double* Kz = (double*)malloc((dim + 10) * sizeof(double));

	//Como se fosse uma transmissibilidade	
	double* lambda1_W = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda1_po_W = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda1_SW_W = (double*)malloc((dim + 10) * sizeof(double));

	double* lambda2_W = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda2_po_W = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda2_SW_W = (double*)malloc((dim + 10) * sizeof(double));

	double* lambda1_O = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda1_po_O = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda1_SW_O = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda1_pb_O = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda1_SO_O = (double*)malloc((dim + 10) * sizeof(double));

	double* lambda2_O = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda2_po_O = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda2_SW_O = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda2_pb_O = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda2_SO_O = (double*)malloc((dim + 10) * sizeof(double));

	double* lambda1_G = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda1_po_G = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda1_SW_G = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda1_pb_G = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda1_SO_G = (double*)malloc((dim + 10) * sizeof(double));

	double* lambda2_G = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda2_po_G = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda2_SW_G = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda2_pb_G = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda2_SO_G = (double*)malloc((dim + 10) * sizeof(double));

	double* lambda3_G = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda3_po_G = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda3_SW_G = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda3_pb_G = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda3_SO_G = (double*)malloc((dim + 10) * sizeof(double));

	double* lambda4_G = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda4_po_G = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda4_SW_G = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda4_pb_G = (double*)malloc((dim + 10) * sizeof(double));
	double* d_lambda4_SO_G = (double*)malloc((dim + 10) * sizeof(double));

	//saber se o ponto do reservat�rio est� acima do pb
	double* X_pb = (double*)malloc((dim + 10) * sizeof(double));
	int* neg = (int*)malloc((dim + 10) * sizeof(int));
	bool* NpureW = (bool*)malloc((dim + 10) * sizeof(bool));

	//posicionamento
	int* unW = (int*)malloc((dim + 10) * sizeof(int));
	int* ufW = (int*)malloc((dim + 10) * sizeof(int));
	int* uwbW = (int*)malloc((dim + 10) * sizeof(int));
	int* ueW = (int*)malloc((dim + 10) * sizeof(int));
	int* ucW = (int*)malloc((dim + 10) * sizeof(int));
	int* uebW = (int*)malloc((dim + 10) * sizeof(int));
	int* uwW = (int*)malloc((dim + 10) * sizeof(int));
	int* ubW = (int*)malloc((dim + 10) * sizeof(int));
	int* usW = (int*)malloc((dim + 10) * sizeof(int));

	int* unO = (int*)malloc((dim + 10) * sizeof(int));
	int* ufO = (int*)malloc((dim + 10) * sizeof(int));
	int* uwbO = (int*)malloc((dim + 10) * sizeof(int));
	int* ueO = (int*)malloc((dim + 10) * sizeof(int));
	int* ucO = (int*)malloc((dim + 10) * sizeof(int));
	int* uebO = (int*)malloc((dim + 10) * sizeof(int));
	int* uwO = (int*)malloc((dim + 10) * sizeof(int));
	int* ubO = (int*)malloc((dim + 10) * sizeof(int));
	int* usO = (int*)malloc((dim + 10) * sizeof(int));

	int* unG = (int*)malloc((dim + 10) * sizeof(int));
	int* ufG = (int*)malloc((dim + 10) * sizeof(int));
	int* uwbG = (int*)malloc((dim + 10) * sizeof(int));
	int* ueG = (int*)malloc((dim + 10) * sizeof(int));
	int* ucG = (int*)malloc((dim + 10) * sizeof(int));
	int* uebG = (int*)malloc((dim + 10) * sizeof(int));
	int* uwG = (int*)malloc((dim + 10) * sizeof(int));
	int* ubG = (int*)malloc((dim + 10) * sizeof(int));
	int* usG = (int*)malloc((dim + 10) * sizeof(int));

	//Volumes Finitos
	double* Ageo_R = (double*)malloc((dim + 10) * sizeof(double));
	double* Ageo_PHI = (double*)malloc((dim + 10) * sizeof(double));
	double* Ageo_Z = (double*)malloc((dim + 10) * sizeof(double));
	double* Ageo_t = (double*)malloc((dim + 10) * sizeof(double));

	double* Adir_F = (double*)malloc((dim + 10) * sizeof(double));
	double* Adir_B = (double*)malloc((dim + 10) * sizeof(double));
	double* Adir_W = (double*)malloc((dim + 10) * sizeof(double));
	double* Adir_E = (double*)malloc((dim + 10) * sizeof(double));
	double* Adir_N = (double*)malloc((dim + 10) * sizeof(double));
	double* Adir_S = (double*)malloc((dim + 10) * sizeof(double));

	double* CHI_F = (double*)malloc((dim + 10) * sizeof(double));
	double* CHI_B = (double*)malloc((dim + 10) * sizeof(double));
	double* CHI_W = (double*)malloc((dim + 10) * sizeof(double));
	double* CHI_E = (double*)malloc((dim + 10) * sizeof(double));
	double* CHI_N = (double*)malloc((dim + 10) * sizeof(double));
	double* CHI_S = (double*)malloc((dim + 10) * sizeof(double));

	double* K_F = (double*)malloc((dim + 10) * sizeof(double));
	double* K_B = (double*)malloc((dim + 10) * sizeof(double));
	double* K_W = (double*)malloc((dim + 10) * sizeof(double));
	double* K_E = (double*)malloc((dim + 10) * sizeof(double));
	double* K_N = (double*)malloc((dim + 10) * sizeof(double));
	double* K_S = (double*)malloc((dim + 10) * sizeof(double));

	double* R_P = (double*)malloc((dim + 10) * sizeof(double));
	double* PHI_P = (double*)malloc((dim + 10) * sizeof(double));
	double* Z_P = (double*)malloc((dim + 10) * sizeof(double));

	double* R_PB = (double*)malloc((dim + 10) * sizeof(double));
	double* R_PF = (double*)malloc((dim + 10) * sizeof(double));
	double* delta_PHI = (double*)malloc((dim + 10) * sizeof(double));
	double* delta_Z = (double*)malloc((dim + 10) * sizeof(double));

	grid_R(Re, Rw, o, m, n,
		R_P, R_PF, R_PB);

	grid_Z(m_ft(L), deltaz, o, m, n,
		delta_Z, Z_P);

	grid_PHI(o, m, n,
		delta_PHI, PHI_P);

	int p = 0;
	int j = 0;
	int k = 0;

	Qwsinj = 0.0l;//Qwsinj / At1[m - 1];//1/s
	Qosinj = 0.0l;//Qosinj / At1[m - 1];//1/s
	Qgsinj = 0.0l;//Qgsinj / At1[m - 1];//1/s

	
	gera_permeabilidades(Kr, n, m, o, perm, varperm);
	get_perm(Kr, n, m, o);
	gera_permeabilidades(Kz, n, m, o, perm, varperm);
	get_perm(Kz, n, m, o);

	double* Sw_layer = (double*)malloc(dim * sizeof(double));
	double* So_layer = (double*)malloc(dim * sizeof(double));

	if (start) {

		for (k = 0; k <= (o - 1); k++) {

			Sw_layer[k] = 0.2l;
			So_layer[k] = 0.8l;

		}

		ini_fields(m, n, o, yo, rho_ws,
			So_coup, Sw_coup, Sg_res,
			po_coup, pw_res, pg_res, pb_coup, X_pb,
			So_layer, Sw_layer, deltaz, pi, pbi, phi_ref, phi_ref_0);
	}
	else {
		reini_fields(m, n, o, So_coup, Sw_coup, Sg_res,
			po_coup, pw_res, pg_res, pb_coup, X_pb,
			phi_ref, phi_ref_0,
			res_timesteps, coup_ts_in, int_in);
	}

	double Qw_total = 0.0l;
	double Qo_total = 0.0l;
	double Qg_total = 0.0l;

	double* Qw_wellbore = (double*)malloc((o + 5) * sizeof(double));
	double* Qo_wellbore = (double*)malloc((o + 5) * sizeof(double));
	double* Qg_wellbore = (double*)malloc((o + 5) * sizeof(double));

	double* Qw_wellbore_ref = (double*)malloc((o + 5) * sizeof(double));
	double* Qo_wellbore_ref = (double*)malloc((o + 5) * sizeof(double));
	double* Qg_wellbore_ref = (double*)malloc((o + 5) * sizeof(double));

	double* PWB_ref = (double*)malloc((o + 5) * sizeof(double));
	double* PWB = (double*)malloc((o + 5) * sizeof(double));
	double* dPWB_D = (double*)malloc((o + 5) * sizeof(double));
	double* PWB_0 = (double*)malloc((o + 5) * sizeof(double));
	double* PWB_coup = (double*)malloc((o + 5) * sizeof(double));

	char pasta[500];
	char pasta2[500];

	//cria arquivo de producao
	memset(&pasta[0], 0, sizeof(pasta));//preenche o string com 0, this terminates the string with a null
	sprintf_s(pasta, "\D_%.2f_THP_%.2f_shotsft_%.2f", D, THP_i, perf_dens_ft);

	_mkdir(pasta);

	sprintf_s(pasta2, pasta);

	char arqu[100];
	sprintf_s(arqu, "\\vazao_x_tempo.csv");
	strcat_s(pasta, arqu);

	fstream vazao(pasta, ios::out | ios::app);
	vazao.precision(8);

	if (start) vazao << "prod_time[s];" << "Ql_total[bbl/d];" << "Qo_total[bbl/d];" << "Qg_total[scf/d];" << "BHP[psi];" << "THP[psi];" << "Res_Iter;" <<endl;

	double prod_time = 0;
	int ite_coup = 0;
	double BHP_r = pi;
	double BHP = po_coup[0] - 500;

	bool dontstopmenow = 1;
	int y = 0, u = 0, l = 0;
	int centerP = 0;
	int centerSw = 0;
	int centerSo = 0;
	double controladorP = 0.0;
	double controladorSw = 0.0;
	double controladorSo = 0.0;
		
	double gama = M_PI / 4;
	double teta = M_PI / 2;
	double teta_yue = 90;
	int N_wellbore = 0;
	int contresint = 0;

	p = 0;
	if (start) {

		for (k = 0; k <= (o - 1); k++) {


			if (zonaprodutora[k]) {

				N_wellbore = N_wellbore + 1;

				PWB_coup[N_wellbore - 1] = po_coup[p];

				PWB[N_wellbore - 1] = po_coup[p] - 500;

			}


			for (j = 0; j <= (m - 1); j++) {

				for (i = 0; i <= (n - 1); i++) {
					p++;

				}

			}

		}
	}
	else {
		double* PWB_temp = (double*)malloc(2 * m * n * o * sizeof(double));
		char arqin1[50];
		sprintf_s(arqin1, "PWB_%d_%d_%d.m", 0, coup_ts_in, int_in);

		ifstream inFile2(arqin1);
		X1 = 0;
		while (!inFile2.eof()) {
			inFile2 >> PWB_temp[X1++];
		}


		for (k = 0; k <= (o - 1); k++) {


			if (zonaprodutora[k]) {

				N_wellbore = N_wellbore + 1;

				PWB_coup[N_wellbore - 1] = PWB_temp[N_wellbore - 1];
				//cout << "PWB_coup: " << PWB_temp[N_wellbore - 1];

				PWB[N_wellbore - 1] = PWB_temp[N_wellbore - 1];
				//system("pause");

			}


			for (j = 0; j <= (m - 1); j++) {

				for (i = 0; i <= (n - 1); i++) {
					p++;

				}

			}

		}

		free(PWB_temp);
	}

	for (k = N_wellbore - 1; k >= 0; k--) {
		dPWB_D[k] = (PWB[k] - PWB[k + 1]) / (BHP - PWB[N_wellbore - 1]);
	}

	//N_wellbore = N_wellbore + 1;

	for (int coup_ts = coup_ts_in + 1; coup_ts <= 18250; coup_ts++) {


		calc_THP(prod_time, coup_ts, &controlstartup, &THP[coup_ts], THP_i ,THP_adj);
			
		deltat_coup = deltat_coup_d * 24 * 3600;

		deltat_res = deltat_coup / res_timesteps;

		ite_coup = 0;
		coup_step_1(coup_ts, ite_coup,
			BHP_r, &BHP,
			THP[coup_ts], L, nz, deltaz,
			&Qw_total, &Qo_total, &Qg_total,
			PWB_ref, PWB, PWB_0, dPWB_D, PWB_coup,
			Qw_wellbore, Qo_wellbore, Qg_wellbore,
			Qw_wellbore_ref, Qo_wellbore_ref, Qg_wellbore_ref,
			N_wellbore, D_perf, perf_dens,
			teta, rug, gama, teta_yue,
			po_coup, So_coup, Sw_coup, pb_coup,
			po_res, pw_res, pg_res,
			So_res, Sw_res, Sg_res, pb_res,
			po_res_0, So_res_0, Sw_res_0, pb_res_0, Sg_res_0, X_pb_0,
			phi_res_0, Bw_res_0, Rs_res_0, Bo_res_0, Bg_res_0,
			po_res_lv0, pb_res_lv0, So_res_lv0, Sw_res_lv0,
			varpo, varpb, varSo, varSw,
			Qwsinj, Qosinj, Qgsinj,
			K_F, K_B, K_N, K_S, K_W, K_E,
			X_pb, zonaprodutora,
			res_timesteps, dontstopmenow,
			p, y, u, l,
			centerP, centerSw, centerSo,
			i, k, j,
			o, m, n, dim,
			Kr, Kz, pi, T,
			IA, JA, A, B, X,
			controladorP, controladorSw, controladorSo,
			C1, C2, C3,
			Ageo_R, Ageo_PHI, Ageo_Z, Ageo_t,
			Adir_F, Adir_B, Adir_W, Adir_E, Adir_N, Adir_S,
			CHI_F, CHI_B, CHI_W, CHI_E, CHI_N, CHI_S,
			R_P, PHI_P, Z_P,
			R_PB, R_PF, delta_PHI, delta_Z, deltat_res,
			g, Re, Rw,
			rho_air, rho_ws, API, yo, rho_os, yg, rho_gs,
			yn2, yco2, yh2s,
			lambda1_W, d_lambda1_po_W, d_lambda1_SW_W,
			lambda2_W, d_lambda2_po_W, d_lambda2_SW_W,
			lambda1_O, d_lambda1_po_O, d_lambda1_SW_O, d_lambda1_pb_O, d_lambda1_SO_O,
			lambda2_O, d_lambda2_po_O, d_lambda2_SW_O, d_lambda2_pb_O, d_lambda2_SO_O,
			lambda1_G, d_lambda1_po_G, d_lambda1_SW_G, d_lambda1_pb_G, d_lambda1_SO_G,
			lambda2_G, d_lambda2_po_G, d_lambda2_SW_G, d_lambda2_pb_G, d_lambda2_SO_G,
			lambda3_G, d_lambda3_po_G, d_lambda3_SW_G, d_lambda3_pb_G, d_lambda3_SO_G,
			lambda4_G, d_lambda4_po_G, d_lambda4_SW_G, d_lambda4_pb_G, d_lambda4_SO_G,
			krwmax, nw, ng, now, nog,
			phi_ref, Cr,
			Swc, Sor, Sgr, epsilon,
			pcwmin, pcwmax, pcgmin, pcgmax,
			parametrosBw, parametrosmig, parametrosZ,
			ppc, Tpc, Bwi, Bpcw, Bpcg,
			iparm, dparm, solver, maxfct, mnum, phase,
			error, msglvl, mtype, nrhs, ddum, idum, pt,
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
			krg, der_krg,
			mi_g, der_mi_g, Z,
			Bg, der_Bg,
			rho_g, der_rho_g,
			unW, ufW, uwbW, ueW, ucW, uwW, uebW, ubW, usW,
			unO, ufO, uwbO, ueO, ucO, uwO, uebO, ubO, usO,
			unG, ufG, uwbG, ueG, ucG, uwG, uebG, ubG, usG,
			neg, NpureW, &contresint,
			pasta2, W_corre, WB_corre);

		prod_time += deltat_coup;
		vazao << prod_time << "\t;" << Qw_total + Qo_total << "\t;" << Qo_total << "\t;" << 5.614585l*Qg_total << "\t;" << BHP << "\t;" << Pa_Psia(THP[coup_ts]) << "\t;" << contresint << "\t;" << endl;

		if (coup_ts % 20 == 0) {
						
			Sleep(10 * 1000);

			if (coup_ts % 50 == 0) {
				Sleep(60 * 1000);

			}
		}

		/*if (controlstartup && (Qo_total) < 1000 && THP_adj<199) {
			THP_adj = THP_adj + 50.l;
		}*/

	}

	vazao.close();

	return 0;
}