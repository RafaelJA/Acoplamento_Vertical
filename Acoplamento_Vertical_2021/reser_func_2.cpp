#include "stdafx.h"
#include "reser_func_2.h"


void calc_coef_grid(double *Ageo_t, double *Ageo_R, double *Ageo_PHI, double *Ageo_Z,
	double *Adir_F, double *Adir_B, double *Adir_W, double *Adir_E, double *Adir_N, double *Adir_S,
	double *CHI_F, double *CHI_B, double *CHI_W, double *CHI_E, double *CHI_N, double *CHI_S,
	double *R_P, double *PHI_P, double *Z_P,
	double *R_PB, double *R_PF, double *delta_PHI, double *delta_Z, double deltat_res,
	double *K_F, double *K_B, double *K_N, double *K_S, double *K_W, double *K_E,
	double *Kr, double *Kz,
	int o, int m, int n) {

	int p = 0;
	for (int k = 0; k <= (o - 1); k++) {
		for (int j = 0; j <= (m - 1); j++) {
			for (int i = 0; i <= (n - 1); i++) {
				
			Ageo_t[p] = delta_PHI[p]* delta_Z[p] * 0.5l*(R_PB[p]*R_PB[p]-R_PF[p]*R_PF[p]);
			Ageo_R[p] = delta_PHI[p]* delta_Z[p]  * deltat_res;
			Ageo_PHI[p] = log(R_PB[p]/R_PF[p])* delta_Z[p]  * deltat_res;
			Ageo_Z[p] = delta_PHI[p] * 0.5l*(R_PB[p]*R_PB[p]-R_PF[p]*R_PF[p]) * deltat_res;


			if (j != 0) K_F[p] = perm_harm(Kr[p], Kr[p - n]);//Forward
			if (j == 0) K_F[p] = Kr[p];//Forward
			if (j != (m - 1)) K_B[p] = perm_harm(Kr[p], Kr[p + n]);;//Backward
			if (j == (m - 1)) K_B[p] = Kr[p];//Backward

			if (n != 1 && i == (n - 1)) K_W[p] = 0.0l;//West(boundary)	
			if (n != 1 && i != (n - 1)) K_W[p] = 0.0l;//West
			if (n != 1 && i == 0) K_E[p] = 0.0l;//East(boundary)
			if (n != 1 && i != (0)) K_E[p] = 0.0l;//East

			if (k != 0) K_N[p] = 0.5*(delta_Z[p] + delta_Z[p - n*m])*perm_harm(Kz[p] /  delta_Z[p], Kz[p - n*m] / delta_Z[p - n*m]);//North
			if (k != (o - 1)) K_S[p] = 0.5*(delta_Z[p] + delta_Z[p + n*m])*perm_harm(Kz[p] / delta_Z[p], Kz[p + n*m] / delta_Z[p + n*m]);//South

			if (j != 0) Adir_F[p] = R_PF[p] / (R_P[p] - R_P[p - n]);//Forward
			if (j == 0) Adir_F[p] = R_PF[p] / (R_P[p] - R_PF[p]);//Forward
			if (j != (m - 1)) Adir_B[p] = R_PB[p] / (R_P[p + n] - R_P[p]);//Backward

			if (n != 1 && i == (n - 1)) Adir_W[p] = 1.0l / abs(PHI_P[p - (n - 1)] - PHI_P[p]);//West(boundary)	
			if (n != 1 && i != (n - 1)) Adir_W[p] = 1.0l / (PHI_P[p + 1] - PHI_P[p]);//West
			if (n != 1 && i == 0) Adir_E[p] = 1.0l / abs(PHI_P[p] - PHI_P[p + (n - 1)]);//East(boundary)
			if (n != 1 && i != (0)) Adir_E[p] = 1.0l / (PHI_P[p] - PHI_P[p - 1]);//East

			//the minus signal is because Z_P denotes depth
			if (k != 0) Adir_N[p] = -1.0l / (Z_P[p - n*m] - Z_P[p]);//North
			if (k != (o - 1)) Adir_S[p] = -1.0l / (Z_P[p] - Z_P[p + n*m]);//South
			
			if (j != 0) CHI_F[p] = (R_PF[p] - R_P[p - n]) / (R_P[p] - R_P[p - n]);//Forward
			if (j != (m - 1)) CHI_B[p] = (R_P[p + n] - R_PB[p]) / (R_P[p + n] - R_P[p]);//Backward

			if (n != 1 && i == (n - 1)) CHI_W[p] = delta_Z[p - (n - 1)] / (delta_Z[p] + delta_Z[p - (n - 1)]);//West(boundary)
			if (n != 1 && i != (n - 1)) CHI_W[p] = delta_Z[p + 1] / (delta_Z[p] + delta_Z[p + 1]);//West
			if (n != 1 && i == 0) CHI_E[p] = delta_Z[p + (n - 1)] / (delta_Z[p] + delta_Z[p + (n - 1)]);//East(boundary)
			if (n != 1 && i != (0)) CHI_E[p] = delta_Z[p - 1] / (delta_Z[p] + delta_Z[p - 1]);//East

			if (k != 0) CHI_N[p] = 0.0;//delta_Z[p - n*m] / (delta_Z[p] + delta_Z[p - n*m]);//North
			if (k != (o - 1)) CHI_S[p] = 1.0l;// delta_Z[p + n*m] / (delta_Z[p] + delta_Z[p + n*m]);//South

			p++;

			}
		}
	}
	
}

void ini_fields(int m, int n, int o, double yo, double rho_ws,
	double *So_coup, double *Sw_coup, double *Sg_coup, 
	double *po_coup, double *pw_coup, double *pg_coup, double *pb_coup, double *X_pb,
	double *So_layer, double *Sw_layer, double *deltaz, double pi, double pbi, double *phi_ref, double phi_ref_0) {
		
	
	int p, k, j, i;
	p = 0;
	for (k = 0; k <= (o - 1); k++) {
		for (j = 0; j <= (m - 1); j++) {
			for (i = 0; i <= (n - 1); i++) {

				Sw_coup[p] = Sw_layer[k];
				So_coup[p] = So_layer[k];		
				Sg_coup[p] = 1.0l - So_coup[p] - Sw_coup[p];
				phi_ref[p] = phi_ref_0;
				p++;

			}
		}
	}

	double grad = 0.0l;
	p = 0;
	for (k = (o-1); k >= 0; k--) {

		if (k < (o - 1)) {
		
			grad = grad + 0.5l*((deltaz[k] * (Sw_layer[k] * rho_ws + So_layer[k] * yo*rho_ws) + deltaz[k + 1] * (Sw_layer[k + 1] * rho_ws + So_layer[k + 1] * yo*rho_ws))) / 144.0l;
		}
		else {
			
			grad = grad + 0.5l* (deltaz[k] * (Sw_layer[k] * rho_ws + So_layer[k] * yo*rho_ws)) / 144.0l;
		}

		for (j = 0; j <= (m - 1); j++) {
			for (i = 0; i <= (n - 1); i++) {

				po_coup[(n*m*o - 1) - p] = pi - grad;
				pb_coup[(n*m*o - 1) - p] = pbi;
				pw_coup[(n*m*o - 1) - p] = po_coup[(n*m*o - 1) - p];// +calc_water_pcw_Sw(Sw[p], Swc, epsilon, pcwmin, pcwmax, Bpcw);
				pg_coup[(n*m*o - 1) - p] = po_coup[(n*m*o - 1) - p];// +calc_gas_pcg(Sw[p], So[p], Swc, Sor, epsilon, pcgmin, pcgmax, Bpcg);

				if (po_coup[p] > pb_coup[p]) {
					X_pb[p] = 1.0l;
				}
				else {
					X_pb[p] = 0.0l;
				}
				

				p++;
			}
		}

	}
}

void reini_fields(int m, int n, int o, double *So_coup, double *Sw_coup, double *Sg_coup,
	double *po_coup, double *pw_coup, double *pg_coup, double *pb_coup, double *X_pb,
	double *phi_ref, double phi_ref_0,
	int res_timesteps, int coup_ts_in, int int_in) {


	int p, k, j, i;
	p = 0;	
	int X1 = 0;
	double  *po_temp = (double *)malloc(2 *m*n*o * sizeof(double));
	double  *pb_temp = (double *)malloc(2 * m*n*o * sizeof(double));
	double  *So_temp = (double *)malloc(2 * m*n*o * sizeof(double));
	double  *Sw_temp = (double *)malloc(2 * m*n*o * sizeof(double));

	
	char arqin1[50];
	sprintf_s(arqin1, "p_%d_%d_%d.m", res_timesteps, coup_ts_in, int_in);
	
	ifstream inFile(arqin1);
	X1 = 0;
	while (!inFile.eof()) {
		inFile >> po_temp[X1++];	
	}
		
	char arqin2[50];
	sprintf_s(arqin2, "pb_%d_%d_%d.m", res_timesteps, coup_ts_in, int_in);

	ifstream inFile2(arqin2);
	X1 = 0;
	while (!inFile2.eof()) {
		inFile2 >> pb_temp[X1++];		
	}
		
	char arqin3[50];
	sprintf_s(arqin3, "So_%d_%d_%d.m", res_timesteps, coup_ts_in, int_in);

	ifstream inFile3(arqin3);
	X1 = 0;
	while (!inFile3.eof()) {
		inFile3 >> So_temp[X1++];
	}
		
	char arqin4[50];
	sprintf_s(arqin4, "Sw_%d_%d_%d.m", res_timesteps, coup_ts_in, int_in);

	ifstream inFile4(arqin4);
	X1 = 0;
	while (!inFile4.eof()) {
		inFile4 >> Sw_temp[X1++];
	}


	for (k = 0; k <= (o - 1); k++) {
		for (j = 0; j <= (m - 1); j++) {
			for (i = 0; i <= (n - 1); i++) {


				po_coup[p] = po_temp[p];
				pb_coup[p] = pb_temp[p];
				pw_coup[p] = po_temp[p];
				pg_coup[p] = po_temp[p];
				Sw_coup[p] = Sw_temp[p];
				So_coup[p] = So_temp[p];
				Sg_coup[p] = 1.l - So_coup[p] - Sw_coup[p];


				if (po_coup[p] > pb_coup[p]) {
					X_pb[p] = 1.0l;
				}
				else {
					X_pb[p] = 0.0l;
				}

				phi_ref[p] = phi_ref_0;

				p++;
			}
		}

	}

	free(po_temp);
	free(pb_temp);
	free(So_temp);
	free(Sw_temp);

}

void grid_R(double re, double rw, int o, int m, int n,
	double *R_P, double *R_PF, double *R_PB) {

	int p = 0;
	double a = pow(re / rw, 1.0l / m);
	double r0 = a*log(a)*rw / (a - 1.0l);
	double RPF = 0.0l;
	double RPB = 0.0l;
	double RPP = 0.0l;

	for (int k = 0; k <= (o - 1); k++) {
		for (int j = 0; j <= (m - 1); j++) {
			
			RPP = pow(a, j)*r0;
			RPF = (RPP - RPP / a) / log(a);
			RPB = (a*RPP-RPP)/log(a);

			for (int i = 0; i <= (n - 1); i++) {

				R_P[p] = RPP;
				R_PF[p] = RPF;
				R_PB[p] = RPB;
				p++;

			}
		}
	}

}

void grid_Z(double L, double *deltaz, int o, int m, int n, 
	double *delta_Z, double  *Z_P) {
	int p = 0;
	for (int k = 0; k <= (o - 1); k++) {
		L = L + 0.5l*deltaz[k];
		
		for (int j = 0; j <= (m - 1); j++) {
			for (int i = 0; i <= (n - 1); i++) {

				Z_P[p] = L;
				delta_Z[p] = deltaz[k];
				p++;

			}
		}
		L = L + 0.5l*deltaz[k];

	}

}

void grid_PHI(int o, int m, int n,
	double *delta_PHI, double  *PHI_P) {
	int p = 0;
	double phip = 0;
	double deltaphi = 2.l*M_PI / n;
	for (int k = 0; k <= (o - 1); k++) {
		for (int j = 0; j <= (m - 1); j++) {
			for (int i = 0; i <= (n - 1); i++) {
				phip = phip + 0.5l*deltaphi;

				PHI_P[p] = phip;
				delta_PHI[p] = deltaphi;
				
				phip = phip + 0.5l*deltaphi;
				p++;

			}
		}
	}

}

double calcula_prop(double x, double *xi, double *yi, int l) {
	int k = 0;
	if (x <= xi[0]) return yi[0];
	while(x>xi[k]) {
		if (k > l) return yi[l];
		k++;
	}	
	
	return (yi[k] - yi[k - 1])*(x - xi[k - 1]) / (xi[k] - xi[k - 1]) + yi[k - 1];
}

double calcula_der_prop(double x, double *xi, double *yi, int l) {
	int k = 0;
	if (x <= xi[0]) return (yi[1] - yi[0]) / (xi[1] - xi[0])/1000.0;
	while (x>xi[k]) {
		if (k > l) return (yi[l] - yi[l - 1]) / (xi[l] - xi[l - 1]) / 1000.0;
		k++;
	}

	return (yi[k] - yi[k - 1]) / (xi[k] - xi[k - 1]);
}

void inicializa_reser(int *res_timesteps,
	double *Re, double *Rw, double *deltaz,
	int *n, int *m, int *o,
	double *pi, double *pbi, 
	double *T, double *g, 
	double *C1, double *C2, double *C3,
	double *rho_air, double *rho_ws,  
	double *rho_os, double *rho_gs, 
	double *yg, double *yo, double *API,
	double *phi_ref, double *Cr,
	double *Swc, double *Sor, double *Sgr, double *epsilon,
	double *pcwmin, double *pcwmax, double *pcgmin, double *pcgmax, 
	double *krwmax, double *nw, double *ng, double *nog, double *now,
	double *parametrosBw, double *parametrosmig, double *ppc, double *Tpc,
	double *Bwi, double *Bpcw, double *Bpcg,
	double *yn2, double *yco2, double *yh2s,
	bool *zonaprodutora,
	double Re_i, double D, double p_i, double pb_i) {

	*Re = Re_i; //raio externo(ft)
	*Rw = 0.5l*in_ft(D); //raio interno(ft)

	//Condição inicial Inicial
	*pi = p_i;//psia
	*pbi = pb_i;//psia	

	*n = 1;//(divisoes do anel)
	*m = 20;//(profundidade)
	*o = 29;//(pra cima)

	*res_timesteps = 15;	

	for (int k = 0; k <= (*o - 1); k++) 
	{	
			if (k < 20) deltaz[k] = 5.0l;
			else deltaz[k] = 40.0l;

			if (k < 12) zonaprodutora[k] = 1;
			else zonaprodutora[k] = 0;

	}
	
	*T = Kel_Far(273.15l+50.0l);//°F
	*g = 32.174049l;//Gravidade ft/s2
	*C1 = 7.324646209006036l*pow(10.0l, -8.0l);
	*C2 = 1.602205l*pow(10.0l, -11.0l);
	*C3 = 1.0l;
	
	//Densidades
	*rho_air = 0.076399332l;//lb/ft³
	*rho_ws = 63.02l; //lb/ft³
	*rho_os = 52.83l;//lb/ft³
	*rho_gs = 0.02804l;//lb/ft³

	*yg = *rho_gs / *rho_air;//
	*yo = *rho_os / *rho_ws;//
	*API = (141.5 / *yo) - 131.5;

	//Parametros phi
	*phi_ref = 25.0l / 100.0l;
	*Cr = 4.700l * pow(10.0l, (-6.0l));//psi-1									  
									   
	//Parametros Saturação
	*Swc = 0.10l;
	*Sor = 0.20l;
	*Sgr = 0.00l;
	*epsilon = 0.00001l;
	*pcwmin = 0.0l;
	*pcwmax = 3.0l;//(psia)
	*pcgmin = 0.0l;
	*pcgmax = 3.9l;//(psia)

	 //Parametros Permeabilidade
	*krwmax = 0.2l;
	*nw = 2.0l;
	*ng = 0.9l;
	*now = 2.0l;
	*nog = 5.5l;
	*yn2 = 0.0l; *yco2 = 0.0l; *yh2s = 0.0l;

	//Variáveis construidas por funções
	calc_const_Bw(parametrosBw);

	calc_const_mig(parametrosmig, *yg, *yn2, *yco2, *yh2s, *T);

	calc_pseudocrit(ppc, Tpc, *yg, *yco2, *yh2s);

	*Bwi = calc_water_Bw_p(*pi, *T, parametrosBw);

	*Bpcw = calc_Bpcw(*Swc, *epsilon, *pcwmin, *pcwmax);

	*Bpcg = calc_Bpcg(*Swc, *Sor, *epsilon, *pcgmin, *pcgmax);
	
}


//FAZ O POSICIONAMENTO DOS VETORES DO PARDISO
void calc_reser_u(int n, int m, int o, int *IA, int *JA,
	int *unW, int *ufW, int *uwbW, int *ueW, int *ucW, int *uwW, int *uebW, int *ubW, int *usW,
	int *unO, int *ufO, int *uwbO, int *ueO, int *ucO, int *uwO, int *uebO, int *ubO, int *usO,
	int *unG, int *ufG, int *uwbG, int *ueG, int *ucG, int *uwG, int *uebG, int *ubG, int *usG, 
	double *Sw, int *neg , bool *NpureW)
{
	int u = 0;
	int p = 0;
	int negtemp = 0;

	for (int k = 0; k <= (o - 1); k++) {
		for (int j = 0; j <= (m - 1); j++) {
			for (int i = 0; i <= (n - 1); i++) {

				neg[p] = negtemp;

				if (Sw[p] >= 0.99l) {
					negtemp = negtemp - 1;
					NpureW[p] = 0;
					
				}
				else {
					NpureW[p] = 1;
				}

				p++;

			}
		}
		neg[p] = negtemp;
	}

	p = 0;
	for (int k = 0; k <= (o - 1); k++) {
		for (int j = 0; j <= (m - 1); j++) {
			for (int i = 0; i <= (n - 1); i++) {

				IA[3 * p + neg[p]] = u + 1;

				//Water		
				{
					//North
					if (k != 0) {
						unW[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3 * n*m + neg[p - n*m];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3 * n*m + neg[p - n*m];
						u = u + 1;
					}

					//Foward
					if (j != 0) {
						ufW[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3 * n + neg[p - n];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3 * n + neg[p - n];
						u = u + 1;
					}

					//West(boundary)
					if (n != 1 && i == (n - 1)) {
						uwbW[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3 * (n - 1) + neg[p - (n - 1)];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3 * (n - 1) + neg[p - (n - 1)];
						u = u + 1;
					}

					//East
					if (n != 1 && i != (0)) {
						ueW[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3 + neg[p - 1];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3 + neg[p - 1];
						u = u + 1;
					}

					//Center
					{
						ucW[p] = u;
						JA[u] = 3 * (p + 1) - 2 + neg[p];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + neg[p];
						u = u + 1;
					}

					//West
					if (n != 1 && i != (n - 1)) {
						uwW[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3 + neg[p + 1];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3 + neg[p + 1];
						u = u + 1;
					}

					//East(boundary)
					if (n != 1 && i == 0) {
						uebW[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3 * (n - 1) + neg[p + (n - 1)];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3 * (n - 1) + neg[p + (n - 1)];
						u = u + 1;
					}

					//Backward
					if (j != (m - 1)) {
						ubW[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3 * n + neg[p + n];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3 * n + neg[p + n];
						u = u + 1;
					}

					//South
					if (k != (o - 1)) {
						usW[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3 * n*m + neg[p + n*m];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3 * n*m + neg[p + n*m];
						u = u + 1;
					}


				}

				IA[3 * p + 1 + neg[p]] = u + 1;
				/*cout << "3 * p + 1 + neg[p]: " << 3 * p + 1 + neg[p] << endl;
				cout << "u: " << u + 1 << endl;*/
				//Oil
				{
					//North
					if (k != 0) {
						unO[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3 * n*m + neg[p - n*m];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3 * n*m + neg[p - n*m];
						u = u + 1;
						if(NpureW[p-n*m]){
							JA[u] = 3 * (p + 1) - 0 - 3 * n*m + neg[p - n*m];
							u = u + 1;
						}
						
					}

					//Foward
					if (j != 0) {
						ufO[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3 * n + neg[p - n];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3 * n + neg[p - n];
						u = u + 1;
						if (NpureW[p - n]){
							JA[u] = 3 * (p + 1) - 0 - 3 * n + neg[p - n];
							u = u + 1;
						}
						
					}

					//West(boundary)
					if (n != 1 && i == (n - 1)) {
						uwbO[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3 * (n - 1) + neg[p - (n - 1)];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3 * (n - 1) + neg[p - (n - 1)];
						u = u + 1;
						if (NpureW[p - (n - 1)]) {
							JA[u] = 3 * (p + 1) - 0 - 3 * (n - 1) + neg[p - (n - 1)];
							u = u + 1;
						}
					}

					//East
					if (n != 1 && i != (0)) {
						ueO[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3 + neg[p - 1];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3 + neg[p - 1];
						u = u + 1;
						if (NpureW[p - 1]) {
							JA[u] = 3 * (p + 1) - 0 - 3 + neg[p - 1];
							u = u + 1;
						}
					}

					//Center
					{
						ucO[p] = u;
						JA[u] = 3 * (p + 1) - 2 + neg[p];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + neg[p];
						u = u + 1;
						if (NpureW[p]) {
							JA[u] = 3 * (p + 1) + neg[p];
							u = u + 1;
						}
					}

					//West
					if (n != 1 && i != (n - 1)) {
						uwO[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3 + neg[p + 1];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3 + neg[p + 1];
						u = u + 1;
						if (NpureW[p + 1]) {
							JA[u] = 3 * (p + 1) - 0 + 3 + neg[p + 1];
							u = u + 1;
						}
					}

					//East(boundary)
					if (n != 1 && i == 0) {
						uebO[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3 * (n - 1) + neg[p + (n - 1)];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3 * (n - 1) + neg[p + (n - 1)];
						u = u + 1;
						if (NpureW[p + (n - 1)]) {
							JA[u] = 3 * (p + 1) - 0 + 3 * (n - 1) + neg[p + (n - 1)];
							u = u + 1;
						}
					}

					//Backward
					if (j != (m - 1)) {
						ubO[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3 * n + neg[p + n];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3 * n + neg[p + n];
						u = u + 1;
						if (NpureW[p + n]) {
							JA[u] = 3 * (p + 1) - 0 + 3 * n + neg[p + n];
							u = u + 1;
						}
					}

					//South
					if (k != (o - 1)) {
						usO[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3 * n*m + neg[p + n*m];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3 * n*m + neg[p + n*m];
						u = u + 1;
						if (NpureW[p + n*m]) {
							JA[u] = 3 * (p + 1) - 0 + 3 * n*m + neg[p + n*m];
							u = u + 1;
						}
					}
				}						

				//Gas	
				if(NpureW[p]){

					IA[3 * p + 2 + neg[p]] = u + 1;					

					//North
					if (k != 0) {
						unG[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3 * n*m + neg[p - n*m];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3 * n*m + neg[p - n*m];
						u = u + 1;
						if (NpureW[p - n*m]) {
							JA[u] = 3 * (p + 1) - 0 - 3 * n*m + neg[p - n*m];
							u = u + 1;
						}
					}

					//Foward
					if (j != 0) {
						ufG[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3 * n + neg[p - n];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3 * n + neg[p - n];
						u = u + 1;
						if (NpureW[p - n]) {
							JA[u] = 3 * (p + 1) - 0 - 3 * n + neg[p - n];
							u = u + 1;
						}
					}

					//West(boundary)
					if (n != 1 && i == (n - 1)) {
						uwbG[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3 * (n - 1) + neg[p - (n - 1)];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3 * (n - 1) + neg[p - (n - 1)];
						u = u + 1;
						if (NpureW[p - (n - 1)]) {
							JA[u] = 3 * (p + 1) - 0 - 3 * (n - 1) + neg[p - (n - 1)];
							u = u + 1;
						}
					}

					//East
					if (n != 1 && i != (0)) {
						ueG[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3 + neg[p - 1];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3 + neg[p - 1];
						u = u + 1;
						if (NpureW[p - 1]) {
							JA[u] = 3 * (p + 1) - 0 - 3 + neg[p - 1];
							u = u + 1;
						}
					}

					//Center
					{
						ucG[p] = u;
						JA[u] = 3 * (p + 1) - 2 + neg[p];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + neg[p];
						u = u + 1;
						if (NpureW[p]) {
							JA[u] = 3 * (p + 1) + neg[p];
							u = u + 1;
						}
					}

					//West
					if (n != 1 && i != (n - 1)) {
						uwG[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3 + neg[p+1];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3 + neg[p+1];
						u = u + 1;
						if (NpureW[p + 1]) {
							JA[u] = 3 * (p + 1) - 0 + 3 + neg[p+1];
							u = u + 1;
						}
					}

					//East(boundary)
					if (n != 1 && i == 0) {
						uebG[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3 * (n - 1) + neg[p + (n - 1)];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3 * (n - 1) + neg[p + (n - 1)];
						u = u + 1;
						if (NpureW[p + 1]) {
							JA[u] = 3 * (p + 1) - 0 + 3 * (n - 1) + neg[p + (n - 1)];
							u = u + 1;
						}
					}

					//Backward
					if (j != (m - 1)) {
						ubG[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3 * n + neg[p + n];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3 * n + neg[p + n];
						u = u + 1;
						if (NpureW[p + n]) {
							JA[u] = 3 * (p + 1) - 0 + 3 * n + neg[p + n];
							u = u + 1;
						}
					}

					//South
					if (k != (o - 1)) {
						usG[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3 * n*m + neg[p + n*m];
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3 * n*m + neg[p + n*m];
						u = u + 1;
						if (NpureW[p + n*m]) {
							JA[u] = 3 * (p + 1) - 0 + 3 * n*m + neg[p + n*m];
							u = u + 1;
						}
					}
				}
				
				p = p + 1;
			}
		}
	}
	IA[3 * p + neg[p]] = u + 1;
	
}


void calc_water_pos_res2(double pw_P, double pw_CHI,
	double *B, double *A,
	double *controladorP, double *controladorSw,
	int y, int u,
	double Ageo, double Adir, double C, double K,
	double CHI,
	double lambda1_P, double lambda1_CHI,
	double d_lambda1_P_po, double d_lambda1_CHI_po,
	double d_lambda1_P_SW, double d_lambda1_CHI_SW, 
	double dpcw_P_SW, double dpcw_CHI_SW)
{
	B[y] = B[y] + C*Ageo*Adir*K*((1.0l-CHI)*lambda1_P + CHI*lambda1_CHI)*(pw_CHI - pw_P);

	A[u] = C*Ageo*Adir*K*((1.0l - CHI)*lambda1_P + CHI*lambda1_CHI + CHI*d_lambda1_CHI_po*(pw_CHI - pw_P));

	*controladorP = *controladorP + C*Ageo*Adir*K*(-(1.0l - CHI)*lambda1_P - CHI*lambda1_CHI + (1.0l - CHI)*d_lambda1_P_po*(pw_CHI - pw_P));

	u = u + 1;

	A[u] = C*Ageo*Adir*K*(-((1.0l - CHI)*lambda1_P + CHI*lambda1_CHI)*dpcw_CHI_SW +CHI*d_lambda1_CHI_SW*(pw_CHI - pw_P));

	*controladorSw = *controladorSw + C*Ageo*Adir*K*(((1.0l - CHI)*lambda1_P + CHI*lambda1_CHI)*dpcw_P_SW + CHI*d_lambda1_P_SW*(pw_CHI - pw_P));
}

void calc_water_grav_res2(double *B, double *A,
	double *controladorP, double *controladorSw,
	int y, int u, double G,
	double Ageo, double C, double K,
	double CHI,
	double lambda2_P, double lambda2_CHI,
	double d_lambda2_P_po, double d_lambda2_CHI_po,
	double d_lambda2_P_SW, double d_lambda2_CHI_SW,
	double dpcw_P_SW, double dpcw_CHI_SW)
{
	B[y] = B[y] + C*Ageo*G*K*((1.0l - CHI)*lambda2_P + CHI*lambda2_CHI);

	A[u] = A[u] + C*Ageo*G*K*(CHI*d_lambda2_CHI_po);

	*controladorP = *controladorP + C*Ageo*G*K*((1.0l - CHI)*d_lambda2_P_po);

	u = u + 1;

	A[u] = A[u] + C*Ageo*G*K*(CHI*d_lambda2_CHI_SW);

	*controladorSw = *controladorSw + C*Ageo*G*K*((1.0l - CHI)*d_lambda2_P_SW);
}


void calc_oil_pos_res2(double po_P, double po_CHI,
	double *B, double *A,
	double *controladorP, double *controladorSw, double *controladorSo,
	int y, int u,
	double Ageo, double Adir, double C, double K,
	double CHI, double X_pb_P, double X_pb_CHI,
	double lambda1_P, double lambda1_CHI,
	double d_lambda1_P_po, double d_lambda1_CHI_po,
	double d_lambda1_P_SW, double d_lambda1_CHI_SW,
	double d_lambda1_P_pb, double d_lambda1_CHI_pb,
	double d_lambda1_P_SO, double d_lambda1_CHI_SO, 
	int neg, bool NpureW)
{
	B[y + neg] = B[y + neg] + C*Ageo*Adir*K*((1.0l - CHI)*lambda1_P + CHI*lambda1_CHI)*(po_CHI - po_P);

	A[u] = C*Ageo*Adir*K*((1.0l - CHI)*lambda1_P + CHI*lambda1_CHI + CHI*d_lambda1_CHI_po*(po_CHI - po_P));

	*controladorP = *controladorP + C*Ageo*Adir*K*(-(1.0l - CHI)*lambda1_P  - CHI*lambda1_CHI + (1.0l - CHI)*d_lambda1_P_po*(po_CHI - po_P));

	u = u + 1;

	A[u] = C*Ageo*Adir*K*(CHI*d_lambda1_CHI_SW*(po_CHI - po_P));

	*controladorSw = *controladorSw + C*Ageo*Adir*K*((1.0l - CHI)*d_lambda1_P_SW*(po_CHI - po_P));
	
	if (NpureW) {
		u = u + 1;

		A[u] = C*Ageo*Adir*K*(CHI*(X_pb_CHI*d_lambda1_CHI_pb + (1.0l - X_pb_CHI)*d_lambda1_CHI_SO)*(po_CHI - po_P));

		*controladorSo = *controladorSo + C*Ageo*Adir*K*((1.0l - CHI)*(X_pb_P*d_lambda1_P_pb + (1.0l - X_pb_P)*d_lambda1_P_SO)*(po_CHI - po_P));
	}
}


void calc_oil_grav_res2(double *B, double *A,
	double *controladorP, double *controladorSw, double *controladorSo,
	int y, int u, double G,
	double Ageo, double C, double K,
	double CHI, double X_pb_P, double X_pb_CHI,
	double lambda2_P, double lambda2_CHI,
	double d_lambda2_P_po, double d_lambda2_CHI_po,
	double d_lambda2_P_SW, double d_lambda2_CHI_SW,
	double d_lambda2_P_pb, double d_lambda2_CHI_pb,
	double d_lambda2_P_SO, double d_lambda2_CHI_SO,
	int neg, bool NpureW)
{
	B[y + neg] = B[y + neg] + C*Ageo*G*K*((1.0l - CHI)*lambda2_P + CHI*lambda2_CHI);

	A[u] = A[u] + C*Ageo*G*K*(CHI*d_lambda2_CHI_po);

	*controladorP = *controladorP + C*Ageo*G*K*((1.0l-CHI)*d_lambda2_P_po);

	u = u + 1;

	A[u] = A[u] + C*Ageo*G*K*(CHI*d_lambda2_CHI_SW);

	*controladorSw = *controladorSw + C*Ageo*G*K*((1.0l - CHI)*d_lambda2_P_SW);

	if (NpureW) {
		u = u + 1;

		A[u] = A[u] + C*Ageo*G*K*((1.0l - CHI)*(X_pb_CHI*d_lambda2_CHI_pb + (1.0l - X_pb_CHI)*d_lambda2_CHI_SO));

		*controladorSo = *controladorSo + C*Ageo*G*K*((1.0l - CHI)*(X_pb_P*d_lambda2_P_pb + (1.0l - X_pb_P)*d_lambda2_P_SO));
	}
}

void calc_gas_pos_res2(double po_P, double po_CHI,
	double pg_P, double pg_CHI,
	double *B, double *A,
	double *controladorP, double *controladorSw, double *controladorSo,
	int y, int u,
	double Ageo, double Adir, double C, double K,
	double CHI, double X_pb_P, double X_pb_CHI,
	double lambda1_P, double lambda1_CHI,
	double d_lambda1_P_po, double d_lambda1_CHI_po,
	double d_lambda1_P_SW, double d_lambda1_CHI_SW,
	double d_lambda1_P_pb, double d_lambda1_CHI_pb,
	double d_lambda1_P_SO, double d_lambda1_CHI_SO,
	double lambda3_P, double lambda3_CHI,
	double d_lambda3_P_po, double d_lambda3_CHI_po,
	double d_lambda3_P_SW, double d_lambda3_CHI_SW,
	double d_lambda3_P_pb, double d_lambda3_CHI_pb,
	double d_lambda3_P_SO, double d_lambda3_CHI_SO,
	double dpcg_P_SW, double dpcg_CHI_SW,
	double dpcg_P_SO, double dpcg_CHI_SO,
	int neg, bool NpureW)
{
	B[y + neg] = B[y + neg] + C*Ageo*Adir*K*(((1.0l - CHI)*lambda1_P + CHI*lambda1_CHI)*(po_CHI - po_P) +
		((1.0l - CHI)*lambda3_P + CHI*lambda3_CHI)*(pg_CHI - pg_P));

	A[u] = C*Ageo*Adir*K*((1.0l - CHI)*lambda1_P + CHI*lambda1_CHI + CHI*d_lambda1_CHI_po*(po_CHI - po_P) +
		(1.0l - CHI)*lambda3_P + CHI*lambda3_CHI + CHI*d_lambda3_CHI_po*(pg_CHI - pg_P));

	*controladorP = *controladorP + C*Ageo*Adir*K*(-(1.0l - CHI)*lambda1_P - CHI*lambda1_CHI + (1.0l - CHI)*d_lambda1_P_po*(po_CHI - po_P) -
		(1.0l - CHI)*lambda3_P - CHI*lambda3_CHI + (1.0l - CHI)*d_lambda3_P_po*(pg_CHI - pg_P));

	u = u + 1;

	A[u] = C*Ageo*Adir*K*(CHI*d_lambda1_CHI_SW*(po_CHI - po_P) +
		((1.0l - CHI)*lambda3_P + CHI*lambda3_CHI)* dpcg_CHI_SW + CHI*d_lambda3_CHI_SW*(pg_CHI - pg_P));

	*controladorSw = *controladorSw + C*Ageo*Adir*K*((1.0l - CHI)*d_lambda1_P_SW*(po_CHI - po_P) -
		((1.0l - CHI)*lambda3_P + CHI*lambda3_CHI)* dpcg_P_SW + (1.0l - CHI)*d_lambda3_P_SW*(pg_CHI - pg_P));

	if (NpureW) {
		u = u + 1;

		A[u] = C*Ageo*Adir*K*(CHI*(X_pb_CHI*d_lambda1_CHI_pb + (1.0l - X_pb_CHI)*d_lambda1_CHI_SO)*(po_CHI - po_P) +
			(1.0l - X_pb_CHI)*((1.0l - CHI)*lambda3_P + CHI*lambda3_CHI)* dpcg_CHI_SO + CHI*(X_pb_CHI*d_lambda3_CHI_pb + (1.0l - X_pb_CHI)*d_lambda3_CHI_SO)*(pg_CHI - pg_P));

		*controladorSo = *controladorSo + C*Ageo*Adir*K*((1.0l - CHI)*(X_pb_P*d_lambda1_P_pb + (1.0l - X_pb_P)*d_lambda1_P_SO)*(po_CHI - po_P) -
			(1.0l - X_pb_P)*((1.0l - CHI)*lambda3_P + CHI*lambda3_CHI)* dpcg_P_SO + (1.0l - CHI)*(X_pb_P*d_lambda3_P_pb + (1.0l - X_pb_P)*d_lambda3_P_SO)*(pg_CHI - pg_P));
	}
}

void calc_gas_grav_res2(double *B, double *A,
	double *controladorP, double *controladorSw, double *controladorSo,
	int y, int u,
	double Ageo, double C, double G, double K,
	double CHI, double X_pb_P, double X_pb_CHI,
	double lambda2_P, double lambda2_CHI,
	double d_lambda2_P_po, double d_lambda2_CHI_po,
	double d_lambda2_P_SW, double d_lambda2_CHI_SW,
	double d_lambda2_P_pb, double d_lambda2_CHI_pb,
	double d_lambda2_P_SO, double d_lambda2_CHI_SO,
	double lambda4_P, double lambda4_CHI,
	double d_lambda4_P_po, double d_lambda4_CHI_po,
	double d_lambda4_P_SW, double d_lambda4_CHI_SW,
	double d_lambda4_P_pb, double d_lambda4_CHI_pb,
	double d_lambda4_P_SO, double d_lambda4_CHI_SO, 
	int neg, bool NpureW)
{
	B[y + neg] = B[y + neg] + C*Ageo*G*K*(((1.0l - CHI)*lambda2_P + CHI*lambda2_CHI) +
		((1.0l - CHI)*lambda4_P + CHI*lambda4_CHI));

	A[u] = A[u] + C*Ageo*G*K*(CHI*d_lambda2_CHI_po +
		CHI*d_lambda4_CHI_po);

	*controladorP = *controladorP + C*Ageo*G*K*((1.0l - CHI)*d_lambda2_P_po +
		(1.0l - CHI)*d_lambda4_P_po);

	u = u + 1;

	A[u] = A[u] + C*Ageo*G*K*(CHI*d_lambda2_CHI_SW +
		CHI*d_lambda4_CHI_SW);

	*controladorSw = *controladorSw + C*Ageo*G*K*((1.0l - CHI)*d_lambda2_P_SW +
		(1.0l - CHI)*d_lambda4_P_SW);

	if (NpureW) {
		u = u + 1;

		A[u] = A[u] + C*Ageo*G*K*(CHI*(X_pb_CHI*d_lambda2_CHI_pb + (1.0l - X_pb_CHI)*d_lambda2_CHI_SO) +
			CHI*(X_pb_CHI*d_lambda4_CHI_pb + (1.0l - X_pb_CHI)*d_lambda4_CHI_SO));

		*controladorSo = *controladorSo + C*Ageo*G*K*((1.0l - CHI)*(X_pb_P*d_lambda2_P_pb + (1.0l - X_pb_P)*d_lambda2_P_SO) +
			(1.0l - CHI)*(X_pb_P*d_lambda4_P_pb + (1.0l - X_pb_P)*d_lambda4_P_SO));
	}
}

void calc_reser_evol(int res_timesteps, int p, int y, int u, int l, int centerP, int centerSw, int centerSo,
	int *IA, int *JA, double *A, double *B, double *X,
	bool dontstopmenow, int i, int k, int j, int o, int m, int n, int dim,
	double controladorP, double controladorSw, double controladorSo,
	double *lambda1_W, double *d_lambda1_po_W, double *d_lambda1_SW_W,
	double *lambda2_W, double *d_lambda2_po_W, double *d_lambda2_SW_W,
	double *lambda1_O, double *d_lambda1_po_O, double *d_lambda1_SW_O, double *d_lambda1_pb_O, double *d_lambda1_SO_O,
	double *lambda2_O, double *d_lambda2_po_O, double *d_lambda2_SW_O, double *d_lambda2_pb_O, double *d_lambda2_SO_O,
	double *lambda1_G, double *d_lambda1_po_G, double *d_lambda1_SW_G, double *d_lambda1_pb_G, double *d_lambda1_SO_G,
	double *lambda2_G, double *d_lambda2_po_G, double *d_lambda2_SW_G, double *d_lambda2_pb_G, double *d_lambda2_SO_G,
	double *lambda3_G, double *d_lambda3_po_G, double *d_lambda3_SW_G, double *d_lambda3_pb_G, double *d_lambda3_SO_G,
	double *lambda4_G, double *d_lambda4_po_G, double *d_lambda4_SW_G, double *d_lambda4_pb_G, double *d_lambda4_SO_G,
	double *Kr, double *Kz, double pi, double T,
	double *PWB, double *PWB_0, double *Pb_WB, int N_wellbore,
	double *po_res, double *pb_res, double *pw_res, double *pg_res, double *Sw_res, double *So_res, double *Sg_res,
	double *po_res_0, double *pb_res_0, double *Sw_res_0, double *So_res_0,	double *Sg_res_0, double *X_pb_0,
	double *phi_res_0, double *Bw_res_0, double *Rs_res_0, double *Bo_res_0, double *Bg_res_0,
	double *po_res_lv0, double *pb_res_lv0, double *So_res_lv0, double *Sw_res_lv0,
	double *varpo, double *varpb, double *varSo, double *varSw,
	double C1, double C2, double C3,
	double *Ageo_R, double *Ageo_PHI, double *Ageo_Z, double *Ageo_t,
	double *Adir_F, double *Adir_B, double *Adir_W, double *Adir_E, double *Adir_N, double *Adir_S,
	double *CHI_F, double *CHI_B, double *CHI_W, double *CHI_E, double *CHI_N, double *CHI_S,
	double *R_P, double *PHI_P, double *Z_P,
	double *R_PB, double *R_PF, double *delta_PHI, double *delta_Z, double deltat_res,
	double g, double Re, double Rw, 
	double rho_air, double rho_ws, double rho_os, double rho_gs,
	double API, double yg, double yo,
	double *phi_ref, double Cr,
	double Qwsinj, double Qosinj, double Qgsinj,
	double *Qw_well, double *Qo_well, double *Qg_well, int coup_ts, int ite_coup,
	double *phi, double *der_phi, 
	double *rho_w, double *der_rho_w, 
	double *mi_w, double *der_mi_w,
	double *Bw, double *der_Bw, 
	double *pcw, double *der_pcw, 
	double *krw, double *der_krw,
	double *Rs, double *der_Rs_p, double *der_Rs_pb, 
	double *kro, double *der_kro_Sw, double *der_kro_So,
	double *Bo, double *der_Bo_p, double *der_Bo_pb,
	double *mi_o, double *der_mi_o_p, double *der_mi_o_pb,
	double *rho_o, double *der_rho_o_p, double *der_rho_o_pb,
	double *pcg, double *der_pcg_Sw, double *der_pcg_So,
	double *krg, double *der_krg, double *mi_g, double *der_mi_g, 
	double *Bg, double *der_Bg, double *rho_g, double *der_rho_g,
	double *K_F, double *K_B, double *K_N, double *K_S, double *K_W, double *K_E,
	double *X_pb, bool *zonaprodutora, int level,
	int *unW, int *ufW, int *uwbW, int *ueW, int *ucW, int *uwW, int *uebW, int *ubW, int *usW,
	int *unO, int *ufO, int *uwbO, int *ueO, int *ucO, int *uwO, int *uebO, int *ubO, int *usO,
	int *unG, int *ufG, int *uwbG, int *ueG, int *ucG, int *uwG, int *uebG, int *ubG, int *usG, 
	int *iparm, double *dparm, int solver, int maxfct, int mnum, int phase, int error,
	int msglvl, int mtype, int nrhs, double ddum, int idum, void *pt,
	double Swc, double Sor, double Sgr, double epsilon,
	double pcwmin, double pcwmax, double pcgmin, double pcgmax,
	double krwmax, double nw, double ng, double nog, double now,
	double *parametrosBw, double *parametrosmig, double ppc, double Tpc,
	double Bwi, double Bpcw, double Bpcg,
	double yn2, double yco2, double yh2s,
	int *neg, bool *NpureW, char directory[])
{
	
	int iy = 0;
	bool nan = 0;	
	double *deltaPWB = (double *)malloc((N_wellbore + 4) * sizeof(double));
	if (level == 0) {
		for (k = 0; k < N_wellbore; k++) {
			deltaPWB[k] = PWB_0[k] - PWB[k];
			if (deltaPWB[k]>2000.l || deltaPWB[k]<(-500.l)) {
				/*cout << "deltaPWB[" << k << "]:" << deltaPWB[k] << endl;
				cout << "PWB_0[" << k << "]:" << PWB_0[k] << endl;
				cout << "PWB[" << k << "]:" << PWB[k] << endl;*/
			
				//system("pause");
			}
		}
	}

	
	for (int t = 0; t < res_timesteps; t++)
	{

		if (level == 0) {		

			for (iy = 0; iy < dim; iy++) {

				po_res_lv0[iy] = po_res_0[iy];
				pb_res_lv0[iy] = pb_res_0[iy];
				So_res_lv0[iy] = So_res_0[iy];
				Sw_res_lv0[iy] = Sw_res_0[iy];

			}

			for (k = 0; k < N_wellbore; k++){
				PWB_0[k] = PWB_0[k] - 0.5l*deltaPWB[k] / res_timesteps;	
				//if(coup_ts==1) cout << "PWB_0[" << k << "]:" << PWB_0[k] << endl;
			}

		}

		calc_coef_grid(Ageo_t, Ageo_R, Ageo_PHI, Ageo_Z,
			Adir_F, Adir_B, Adir_W, Adir_E, Adir_N, Adir_S,
			CHI_F, CHI_B, CHI_W, CHI_E, CHI_N, CHI_S,
			R_P, PHI_P, Z_P,
			R_PB, R_PF, delta_PHI, delta_Z, deltat_res,
			K_F, K_B, K_N, K_S, K_W, K_E,
			Kr, Kz,
			o, m, n);


		l = 0;
		dontstopmenow = 1;
		while (dontstopmenow)
		{
			nan = 0;
			dontstopmenow = 0;
			p = 0;

			calc_props_grid(dim, po_res, pg_res, pw_res, pb_res, Sw_res, So_res, Sg_res, X_pb,
				phi, der_phi,
				rho_w, der_rho_w, mi_w, der_mi_w, Bw, der_Bw, pcw, der_pcw, krw, der_krw,
				Rs, der_Rs_p, der_Rs_pb, kro, der_kro_Sw, der_kro_So, Bo, der_Bo_p, der_Bo_pb,
				mi_o, der_mi_o_p, der_mi_o_pb, rho_o, der_rho_o_p, der_rho_o_pb,
				pcg, der_pcg_Sw, der_pcg_So, krg, der_krg, mi_g, der_mi_g, Bg, der_Bg, rho_g, der_rho_g,
				phi_ref, Cr, Bwi, rho_ws, pi, T, parametrosBw, 
				Bpcw, Bpcg, API, yo, yg, rho_os, rho_gs, 
				parametrosmig, ppc, Tpc, pcwmin, pcwmax,
				lambda1_W, d_lambda1_po_W, d_lambda1_SW_W,
				lambda2_W, d_lambda2_po_W, d_lambda2_SW_W,
				lambda1_O, d_lambda1_po_O, d_lambda1_SW_O, d_lambda1_pb_O, d_lambda1_SO_O,
				lambda2_O, d_lambda2_po_O, d_lambda2_SW_O, d_lambda2_pb_O, d_lambda2_SO_O,
				lambda1_G, d_lambda1_po_G, d_lambda1_SW_G, d_lambda1_pb_G, d_lambda1_SO_G,
				lambda2_G, d_lambda2_po_G, d_lambda2_SW_G, d_lambda2_pb_G, d_lambda2_SO_G,
				lambda3_G, d_lambda3_po_G, d_lambda3_SW_G, d_lambda3_pb_G, d_lambda3_SO_G,
				lambda4_G, d_lambda4_po_G, d_lambda4_SW_G, d_lambda4_pb_G, d_lambda4_SO_G); 

			calc_props_grid_ant(dim, po_res_0, pb_res_0, Sw_res_0, So_res_0, Sg_res_0, X_pb_0,
				phi_res_0, Bw_res_0, Rs_res_0, Bo_res_0, Bg_res_0,
	   			phi_ref, Cr, pi, T, parametrosBw, 
	   			API, yo, yg);

			calc_reser_u(n, m, o, IA, JA,
				unW, ufW, uwbW, ueW, ucW, uwW, uebW, ubW, usW,
				unO, ufO, uwbO, ueO, ucO, uwO, uebO, ubO, usO,
				unG, ufG, uwbG, ueG, ucG, uwG, uebG, ubG, usG,
				Sw_res,	neg, NpureW);

			cout.precision(5);
			p = 0;
			int k_prod = 0;
			

			for (k = 0; k <= (o - 1); k++) {
				for (j = 0; j <= (m - 1); j++) {
					for (i = 0; i <= (n - 1); i++) {

						y = 3 * p;
						controladorP = 0.0;
						controladorSw = 0.0;
						controladorSo = 0.0;
						B[y + neg[p]] = 0.0;

						
						//Water						
						//North
						if (k != 0) {

							calc_water_pos_res2(pw_res[p], pw_res[p - n*m],
								B, A,
								&controladorP, &controladorSw,
								y, unW[p],
								Ageo_Z[p], Adir_N[p], C1, K_N[p],
								CHI_N[p],
								lambda1_W[p], lambda1_W[p - n*m],
								d_lambda1_po_W[p], d_lambda1_po_W[p - n*m],
								d_lambda1_SW_W[p], d_lambda1_SW_W[p - n*m],
								der_pcw[p], der_pcw[p - n*m]);
													
							
							calc_water_grav_res2(B, A,
								&controladorP, &controladorSw,
								y, unW[p], g,
								Ageo_Z[p], C2, K_N[p],
								CHI_N[p],
								lambda2_W[p], lambda2_W[p - n*m],
								d_lambda2_po_W[p], d_lambda2_po_W[p - n*m],
								d_lambda2_SW_W[p], d_lambda2_SW_W[p - n*m],
								der_pcw[p], der_pcw[p - n*m]);

						}
						
						//Foward
						if (j != 0) {

							calc_water_pos_res2(pw_res[p], pw_res[p - n],
								B, A,
								&controladorP, &controladorSw,
								y, ufW[p],
								Ageo_R[p], Adir_F[p], C1, K_F[p],
								CHI_F[p],
								lambda1_W[p], lambda1_W[p - n],
								d_lambda1_po_W[p], d_lambda1_po_W[p - n],
								d_lambda1_SW_W[p], d_lambda1_SW_W[p - n],
								der_pcw[p], der_pcw[p - n]);

						}

						//West(boundary)
						if (n != 1 && i == (n - 1)) {

							calc_water_pos_res2(pw_res[p], pw_res[p - (n - 1)],
								B, A,
								&controladorP, &controladorSw,
								y, uwW[p],
								Ageo_PHI[p], Adir_W[p], C1, K_W[p],
								CHI_W[p],
								lambda1_W[p], lambda1_W[p - (n - 1)],
								d_lambda1_po_W[p], d_lambda1_po_W[p - (n - 1)],
								d_lambda1_SW_W[p], d_lambda1_SW_W[p - (n - 1)],
								der_pcw[p], der_pcw[p - (n - 1)]);

						}

						//East
						if (n != 1 && i != (0)) {

							calc_water_pos_res2(pw_res[p], pw_res[p - 1],
								B, A,
								&controladorP, &controladorSw,
								y, ueW[p],
								Ageo_PHI[p], Adir_E[p], C1, K_E[p],
								CHI_E[p],
								lambda1_W[p], lambda1_W[p - 1],
								d_lambda1_po_W[p], d_lambda1_po_W[p - 1],
								d_lambda1_SW_W[p], d_lambda1_SW_W[p - 1],
								der_pcw[p], der_pcw[p - 1]);

						}

						//Center

						B[y+ neg[p]] = B[y+neg[p]] - C3 * Ageo_t[p] * ((phi[p] * Sw_res[p] / Bw[p]) - (phi_res_0[p] * Sw_res_0[p] / Bw_res_0[p]));


						if (j == (m - 1)) {
							B[y + neg[p]] = B[y+ neg[p]] + Qwsinj / (deltat_res*Bw[p] * n);
						}
						if (j == 0 && zonaprodutora[k]) {
							B[y + neg[p]] = B[y + neg[p]] + C1*Ageo_R[p] * Adir_F[p] * (K_F[p] * lambda1_W[p] * (PWB_0[k_prod] - pw_res[p]));
							Qw_well[k_prod] = -ft3s_bbld(C1*Ageo_R[p] * Adir_F[p] * (K_F[p] * lambda1_W[p] * (PWB_0[k_prod] - pw_res[p])) / (deltat_res));							
						}

						centerP = ucW[p];

						A[ucW[p]] = -C3*Ageo_t[p] * ((Sw_res[p] / Bw[p])*(der_phi[p] - (phi[p] * der_Bw[p] / Bw[p])));

						if (j == (m - 1)) {
							A[ucW[p]] = A[ucW[p]] - Qwsinj*der_Bw[p] / (deltat_res * Bw[p] * Bw[p] * n);
						}
						if (j == 0 && zonaprodutora[k]) {
							A[ucW[p]] = A[ucW[p]] - C1*Ageo_R[p] * Adir_F[p] * (K_F[p] * lambda1_W[p] - K_F[p] * d_lambda1_po_W[p] * (PWB_0[k_prod] - pw_res[p]));
						}

						centerSw = ucW[p] + 1;

						A[ucW[p] + 1] = -Ageo_t[p] * C3* (phi[p] / Bw[p]);

						if (j == 0 && zonaprodutora[k]) {
							A[ucW[p] + 1] = A[ucW[p] + 1] + C1*Ageo_R[p] * Adir_F[p] * (-K_F[p] * lambda1_W[p] * der_pcw[p] + K_F[p] * d_lambda1_SW_W[p] * (PWB_0[k_prod] - pw_res[p]));
						}

						//West
						if (n != 1 && i != (n - 1)) {

							calc_water_pos_res2(pw_res[p], pw_res[p + 1],
								B, A,
								&controladorP, &controladorSw,
								y, uwW[p],
								Ageo_PHI[p], Adir_W[p], C1, K_W[p],
								CHI_W[p],
								lambda1_W[p], lambda1_W[p + 1],
								d_lambda1_po_W[p], d_lambda1_po_W[p + 1],
								d_lambda1_SW_W[p], d_lambda1_SW_W[p + 1],
								der_pcw[p], der_pcw[p + 1]);

						}

						//East(boundary)
						if (n != 1 && i == 0) {

							calc_water_pos_res2(pw_res[p], pw_res[p + (n - 1)],
								B, A,
								&controladorP, &controladorSw,
								y, ueW[p],
								Ageo_PHI[p], Adir_E[p], C1, K_E[p],
								CHI_E[p],
								lambda1_W[p], lambda1_W[p + (n - 1)],
								d_lambda1_po_W[p], d_lambda1_po_W[p + (n - 1)],
								d_lambda1_SW_W[p], d_lambda1_SW_W[p + (n - 1)],
								der_pcw[p], der_pcw[p + (n - 1)]);

						}
												
						//Backward
						if (j != (m - 1)) {

							calc_water_pos_res2(pw_res[p], pw_res[p + n],
								B, A,
								&controladorP, &controladorSw,
								y, ubW[p],
								Ageo_R[p], Adir_B[p], C1, K_B[p],
								CHI_B[p],
								lambda1_W[p], lambda1_W[p + n],
								d_lambda1_po_W[p], d_lambda1_po_W[p + n],
								d_lambda1_SW_W[p], d_lambda1_SW_W[p + n],
								der_pcw[p], der_pcw[p + n]);

						}						
						
						//South
						if (k != (o - 1)) {

							calc_water_pos_res2(pw_res[p], pw_res[p + n*m],
								B, A,
								&controladorP, &controladorSw,
								y, usW[p],
								Ageo_Z[p], Adir_S[p], C1, K_S[p],
								CHI_S[p],
								lambda1_W[p], lambda1_W[p + n*m],
								d_lambda1_po_W[p], d_lambda1_po_W[p + n*m],
								d_lambda1_SW_W[p], d_lambda1_SW_W[p + n*m],
								der_pcw[p], der_pcw[p + n*m]);
														
							calc_water_grav_res2(B, A,
								&controladorP, &controladorSw,
								y, usW[p], -g,
								Ageo_Z[p], C2, K_S[p],
								CHI_S[p],
								lambda2_W[p], lambda2_W[p + n*m],
								d_lambda2_po_W[p], d_lambda2_po_W[p + n*m],
								d_lambda2_SW_W[p], d_lambda2_SW_W[p + n*m],
								der_pcw[p], der_pcw[p + n*m]);

						}						

						A[centerP] = A[centerP] + controladorP;
						A[centerSw] = A[centerSw] + controladorSw;


						controladorP = 0.0;
						controladorSw = 0.0;
						controladorSo = 0.0;
						y = y + 1;
						B[y + neg[p]] = 0.0;


						//Oil
						//North
						if (k != 0) {
							
							calc_oil_pos_res2(po_res[p], po_res[p - n*m],
								B, A,
								&controladorP, &controladorSw, &controladorSo,
								y, unO[p],
								Ageo_Z[p], Adir_N[p], C1, K_N[p],
								CHI_N[p], X_pb[p], X_pb[p - n*m],
								lambda1_O[p], lambda1_O[p - n*m],
								d_lambda1_po_O[p], d_lambda1_po_O[p - n*m],
								d_lambda1_SW_O[p], d_lambda1_SW_O[p - n*m],
								d_lambda1_pb_O[p], d_lambda1_pb_O[p - n*m],
								d_lambda1_SO_O[p], d_lambda1_SO_O[p - n*m],
								neg[p - n*m], NpureW[p - n*m]);


							calc_oil_grav_res2(B, A,
								&controladorP, &controladorSw, &controladorSo,
								y, unO[p], g,
								Ageo_Z[p], C2, K_N[p],
								CHI_N[p], X_pb[p], X_pb[p - n*m],
								lambda2_O[p], lambda2_O[p - n*m],
								d_lambda2_po_O[p], d_lambda2_po_O[p - n*m],
								d_lambda2_SW_O[p], d_lambda2_SW_O[p - n*m],
								d_lambda2_pb_O[p], d_lambda2_pb_O[p - n*m],
								d_lambda2_SO_O[p], d_lambda2_SO_O[p - n*m],
								neg[p - n*m], NpureW[p - n*m]);

						}
						
						//Forward
						if (j != 0) {

							calc_oil_pos_res2(po_res[p], po_res[p - n],
								B, A,
								&controladorP, &controladorSw, &controladorSo,
								y, ufO[p],
								Ageo_R[p], Adir_F[p], C1, K_F[p],
								CHI_F[p], X_pb[p], X_pb[p - n],
								lambda1_O[p], lambda1_O[p - n],
								d_lambda1_po_O[p], d_lambda1_po_O[p - n],
								d_lambda1_SW_O[p], d_lambda1_SW_O[p - n],
								d_lambda1_pb_O[p], d_lambda1_pb_O[p - n],
								d_lambda1_SO_O[p], d_lambda1_SO_O[p - n],
								neg[p - n], NpureW[p - n]);

						}
						
						//West(boundary)
						if (n != 1 && i == (n - 1)) {

							calc_oil_pos_res2(po_res[p], po_res[p - (n - 1)],
								B, A,
								&controladorP, &controladorSw, &controladorSo,
								y, uwO[p],
								Ageo_PHI[p], Adir_W[p], C1, K_W[p],
								CHI_W[p], X_pb[p], X_pb[p - (n - 1)],
								lambda1_O[p], lambda1_O[p - (n - 1)],
								d_lambda1_po_O[p], d_lambda1_po_O[p - (n - 1)],
								d_lambda1_SW_O[p], d_lambda1_SW_O[p - (n - 1)],
								d_lambda1_pb_O[p], d_lambda1_pb_O[p - (n - 1)],
								d_lambda1_SO_O[p], d_lambda1_SO_O[p - (n - 1)], 
								neg[p - (n - 1)], NpureW[p - (n - 1)]);

						}

						//East
						if (n != 1 && i != (0)) {

							calc_oil_pos_res2(po_res[p], po_res[p - 1],
								B, A,
								&controladorP, &controladorSw, &controladorSo,
								y, ueO[p],
								Ageo_PHI[p], Adir_E[p], C1, K_E[p],
								CHI_E[p], X_pb[p], X_pb[p - 1],
								lambda1_O[p], lambda1_O[p - 1],
								d_lambda1_po_O[p], d_lambda1_po_O[p - 1],
								d_lambda1_SW_O[p], d_lambda1_SW_O[p - 1],
								d_lambda1_pb_O[p], d_lambda1_pb_O[p - 1],
								d_lambda1_SO_O[p], d_lambda1_SO_O[p - 1],
								neg[p - 1], NpureW[p - 1]);

						}

						//Center
						B[y + neg[p]] = B[y + neg[p]] - C3 * Ageo_t[p] * (phi[p] * (So_res[p] / Bo[p]) - phi_res_0[p]*(So_res_0[p]/Bo_res_0[p]));
						

						if (j == (m - 1)) {
							B[y + neg[p]] = B[y + neg[p]] + Qosinj / (deltat_res*Bo[p] * n);
						}
						if (j == 0 && zonaprodutora[k]) {
							B[y + neg[p]] = B[y + neg[p]] + C1 * Ageo_R[p] * Adir_F[p] * (K_F[p] * lambda1_O[p] * (PWB_0[k_prod] - po_res[p]));
							Qo_well[k_prod] = -ft3s_bbld(C1 * Ageo_R[p] * Adir_F[p] * (K_F[p] * lambda1_O[p] * (PWB_0[k_prod] - po_res[p])) / (deltat_res));
						}

						centerP = ucO[p];

						A[ucO[p]] = -C3 * Ageo_t[p] * ((So_res[p] / Bo[p])*(der_phi[p] - (phi[p] * der_Bo_p[p] / Bo[p])));
						
						if (j == (m - 1)) {
							A[ucO[p]] = A[ucO[p]] - Qosinj * (der_Bo_p[p] / (deltat_res * n *Bo[p] * Bo[p]));
						}
						if (j == 0 && zonaprodutora[k]) {
							A[ucO[p]] = A[ucO[p]] - C1*Ageo_R[p] * Adir_F[p] * K_F[p] * (lambda1_O[p] - d_lambda1_po_O[p] * (PWB_0[k_prod] - po_res[p]));
						}

						centerSw = ucO[p] + 1;

						A[ucO[p] + 1] = C3 * Ageo_t[p] * X_pb[p] * ((phi[p] / Bo[p]));

						if (j == 0 && zonaprodutora[k]) {
							A[ucO[p] + 1] = A[ucO[p] + 1] + C1*Ageo_R[p] * Adir_F[p] * (K_F[p] * d_lambda1_SW_O[p] * (PWB_0[k_prod] - po_res[p]));
						}

						if (NpureW[p]) {
							centerSo = ucO[p] + 2;

							A[ucO[p] + 2] = -C3 * Ageo_t[p] * ((1.0l - X_pb[p]) * ((phi[p] / Bo[p])) - X_pb[p] * (So_res[p] / Bo[p])*((phi[p] * der_Bo_pb[p] / Bo[p])));
														
							if (j == 0 && zonaprodutora[k]) {
								A[ucO[p] + 2] = A[ucO[p] + 2] + C1*Ageo_R[p] * Adir_F[p] * (K_F[p] * (X_pb[p] * d_lambda1_pb_O[p] + (1.0l - X_pb[p])*d_lambda1_SO_O[p])
									* (PWB_0[k_prod] - po_res[p]));
							}
						}

						//West
						if (n != 1 && i != (n - 1)) {

							calc_oil_pos_res2(po_res[p], po_res[p + 1],
								B, A,
								&controladorP, &controladorSw, &controladorSo,
								y, uwO[p],
								Ageo_PHI[p], Adir_W[p], C1, K_W[p],
								CHI_W[p], X_pb[p], X_pb[p + 1],
								lambda1_O[p], lambda1_O[p + 1],
								d_lambda1_po_O[p], d_lambda1_po_O[p + 1],
								d_lambda1_SW_O[p], d_lambda1_SW_O[p + 1],
								d_lambda1_pb_O[p], d_lambda1_pb_O[p + 1],
								d_lambda1_SO_O[p], d_lambda1_SO_O[p + 1],
								neg[p + 1], NpureW[p + 1]);

						}
						
						//East(boundary)
						if (n != 1 && i == 0) {

							calc_oil_pos_res2(po_res[p], po_res[p + (n - 1)],
								B, A,
								&controladorP, &controladorSw, &controladorSo,
								y, ueO[p],
								Ageo_PHI[p], Adir_E[p], C1, K_E[p],
								CHI_E[p], X_pb[p], X_pb[p + (n - 1)],
								lambda1_O[p], lambda1_O[p + (n - 1)],
								d_lambda1_po_O[p], d_lambda1_po_O[p + (n - 1)],
								d_lambda1_SW_O[p], d_lambda1_SW_O[p + (n - 1)],
								d_lambda1_pb_O[p], d_lambda1_pb_O[p + (n - 1)],
								d_lambda1_SO_O[p], d_lambda1_SO_O[p + (n - 1)],
								neg[p + (n - 1)], NpureW[p + (n - 1)]);

						}
						
						//Backward
						if (j != (m - 1)) {

							calc_oil_pos_res2(po_res[p], po_res[p + n],
								B, A,
								&controladorP, &controladorSw, &controladorSo,
								y, ubO[p],
								Ageo_R[p], Adir_B[p], C1, K_B[p],
								CHI_B[p], X_pb[p], X_pb[p + n],
								lambda1_O[p], lambda1_O[p + n],
								d_lambda1_po_O[p], d_lambda1_po_O[p + n],
								d_lambda1_SW_O[p], d_lambda1_SW_O[p + n],
								d_lambda1_pb_O[p], d_lambda1_pb_O[p + n],
								d_lambda1_SO_O[p], d_lambda1_SO_O[p + n],
								neg[p + n], NpureW[p + n]);

						}
						
						//South						
						if (k != (o - 1)) {

							calc_oil_pos_res2(po_res[p], po_res[p + n*m],
								B, A,
								&controladorP, &controladorSw, &controladorSo,
								y, usO[p],
								Ageo_Z[p], Adir_S[p], C1, K_S[p],
								CHI_S[p], X_pb[p], X_pb[p + n*m],
								lambda1_O[p], lambda1_O[p + n*m],
								d_lambda1_po_O[p], d_lambda1_po_O[p + n*m],
								d_lambda1_SW_O[p], d_lambda1_SW_O[p + n*m],
								d_lambda1_pb_O[p], d_lambda1_pb_O[p + n*m],
								d_lambda1_SO_O[p], d_lambda1_SO_O[p + n*m],
								neg[p + n*m], NpureW[p + n*m]);

							calc_oil_grav_res2(B, A,
								&controladorP, &controladorSw, &controladorSo,
								y, usO[p], -g,
								Ageo_Z[p], C2, K_S[p],
								CHI_S[p], X_pb[p], X_pb[p + n*m],
								lambda2_O[p], lambda2_O[p + n*m],
								d_lambda2_po_O[p], d_lambda2_po_O[p + n*m],
								d_lambda2_SW_O[p], d_lambda2_SW_O[p + n*m],
								d_lambda2_pb_O[p], d_lambda2_pb_O[p + n*m],
								d_lambda2_SO_O[p], d_lambda2_SO_O[p + n*m],
								neg[p + n*m], NpureW[p + n*m]);

						}
						
						A[centerP] = A[centerP] + controladorP;
						A[centerSw] = A[centerSw] + controladorSw;
						if (NpureW[p]) A[centerSo] = A[centerSo] + controladorSo;


						controladorP = 0.0;
						controladorSw = 0.0;
						controladorSo = 0.0;
						y = y + 1;
						//Gas
						if (NpureW[p]) {
							
							B[y + neg[p]] = 0.0;

							//North
							if (k != 0) {

								
								calc_gas_pos_res2(po_res[p], po_res[p - n*m],
									pg_res[p], pg_res[p - n*m],
									B, A,
									&controladorP, &controladorSw, &controladorSo,
									y, unG[p],
									Ageo_Z[p], Adir_N[p], C1, K_N[p],
									CHI_N[p], X_pb[p], X_pb[p - n*m],
									lambda1_G[p], lambda1_G[p - n*m],
									d_lambda1_po_G[p], d_lambda1_po_G[p - n*m],
									d_lambda1_SW_G[p], d_lambda1_SW_G[p - n*m],
									d_lambda1_pb_G[p], d_lambda1_pb_G[p - n*m],
									d_lambda1_SO_G[p], d_lambda1_SO_G[p - n*m],
									lambda3_G[p], lambda3_G[p - n*m],
									d_lambda3_po_G[p], d_lambda3_po_G[p - n*m],
									d_lambda3_SW_G[p], d_lambda3_SW_G[p - n*m],
									d_lambda3_pb_G[p], d_lambda3_pb_G[p - n*m],
									d_lambda3_SO_G[p], d_lambda3_SO_G[p - n*m],
									der_pcg_Sw[p], der_pcg_Sw[p - n*m],
									der_pcg_So[p], der_pcg_So[p - n*m],
									neg[p - n*m], NpureW[p - n*m]);
							
							calc_gas_grav_res2(B, A,
									&controladorP, &controladorSw, &controladorSo,
									y, unG[p],
									Ageo_Z[p], C2, g, K_N[p],
									CHI_N[p], X_pb[p], X_pb[p - n*m],
									lambda2_G[p], lambda2_G[p - n*m],
									d_lambda2_po_G[p], d_lambda2_po_G[p - n*m],
									d_lambda2_SW_G[p], d_lambda2_SW_G[p - n*m],
									d_lambda2_pb_G[p], d_lambda2_pb_G[p - n*m],
									d_lambda2_SO_G[p], d_lambda2_SO_G[p - n*m],
									lambda4_G[p], lambda4_G[p - n*m],
									d_lambda4_po_G[p], d_lambda4_po_G[p - n*m],
									d_lambda4_SW_G[p], d_lambda4_SW_G[p - n*m],
									d_lambda4_pb_G[p], d_lambda4_pb_G[p - n*m],
									d_lambda4_SO_G[p], d_lambda4_SO_G[p - n*m],
									neg[p - n*m], NpureW[p - n*m]);
								
							}
																
							//Forward
							if (j != 0) {

								calc_gas_pos_res2(po_res[p], po_res[p - n],
									pg_res[p], pg_res[p - n],
									B, A,
									&controladorP, &controladorSw, &controladorSo,
									y, ufG[p],
									Ageo_R[p], Adir_F[p], C1, K_F[p],
									CHI_F[p], X_pb[p], X_pb[p - n],
									lambda1_G[p], lambda1_G[p - n],
									d_lambda1_po_G[p], d_lambda1_po_G[p - n],
									d_lambda1_SW_G[p], d_lambda1_SW_G[p - n],
									d_lambda1_pb_G[p], d_lambda1_pb_G[p - n],
									d_lambda1_SO_G[p], d_lambda1_SO_G[p - n],
									lambda3_G[p], lambda3_G[p - n],
									d_lambda3_po_G[p], d_lambda3_po_G[p - n],
									d_lambda3_SW_G[p], d_lambda3_SW_G[p - n],
									d_lambda3_pb_G[p], d_lambda3_pb_G[p - n],
									d_lambda3_SO_G[p], d_lambda3_SO_G[p - n],
									der_pcg_Sw[p], der_pcg_Sw[p - n],
									der_pcg_So[p], der_pcg_So[p - n],
									neg[p - n], NpureW[p - n]);

							}

							//West(boundary)
							if (n != 1 && i == (n - 1)) {

								calc_gas_pos_res2(po_res[p], po_res[p - (n - 1)],
									pg_res[p], pg_res[p - (n - 1)],
									B, A,
									&controladorP, &controladorSw, &controladorSo,
									y, uwG[p],
									Ageo_PHI[p], Adir_W[p], C1, K_W[p],
									CHI_W[p], X_pb[p], X_pb[p - (n - 1)],
									lambda1_G[p], lambda1_G[p - (n - 1)],
									d_lambda1_po_G[p], d_lambda1_po_G[p - (n - 1)],
									d_lambda1_SW_G[p], d_lambda1_SW_G[p - (n - 1)],
									d_lambda1_pb_G[p], d_lambda1_pb_G[p - (n - 1)],
									d_lambda1_SO_G[p], d_lambda1_SO_G[p - (n - 1)],
									lambda3_G[p], lambda3_G[p - (n - 1)],
									d_lambda3_po_G[p], d_lambda3_po_G[p - (n - 1)],
									d_lambda3_SW_G[p], d_lambda3_SW_G[p - (n - 1)],
									d_lambda3_pb_G[p], d_lambda3_pb_G[p - (n - 1)],
									d_lambda3_SO_G[p], d_lambda3_SO_G[p - (n - 1)],
									der_pcg_Sw[p], der_pcg_Sw[p - (n - 1)],
									der_pcg_So[p], der_pcg_So[p - (n - 1)],
									neg[p - (n-1)], NpureW[p - (n - 1)]);

							}

							//East
							if (n != 1 && i != (0)) {

								calc_gas_pos_res2(po_res[p], po_res[p - 1],
									pg_res[p], pg_res[p - 1],
									B, A,
									&controladorP, &controladorSw, &controladorSo,
									y, ueG[p],
									Ageo_PHI[p], Adir_E[p], C1, K_E[p],
									CHI_E[p], X_pb[p], X_pb[p - 1],
									lambda1_G[p], lambda1_G[p - 1],
									d_lambda1_po_G[p], d_lambda1_po_G[p - 1],
									d_lambda1_SW_G[p], d_lambda1_SW_G[p - 1],
									d_lambda1_pb_G[p], d_lambda1_pb_G[p - 1],
									d_lambda1_SO_G[p], d_lambda1_SO_G[p - 1],
									lambda3_G[p], lambda3_G[p - 1],
									d_lambda3_po_G[p], d_lambda3_po_G[p - 1],
									d_lambda3_SW_G[p], d_lambda3_SW_G[p - 1],
									d_lambda3_pb_G[p], d_lambda3_pb_G[p - 1],
									d_lambda3_SO_G[p], d_lambda3_SO_G[p - 1],
									der_pcg_Sw[p], der_pcg_Sw[p - 1],
									der_pcg_So[p], der_pcg_So[p - 1],
									neg[p - 1], NpureW[p - 1]);

							}
							
							//Center
							B[y + neg[p]] = B[y + neg[p]] - C3 * Ageo_t[p] * (phi[p] * ((So_res[p] * Rs[p] / Bo[p]) + (Sg_res[p] / Bg[p]))
								- phi_res_0[p] * ((So_res_0[p] * Rs_res_0[p] / Bo_res_0[p]) + (Sg_res_0[p] / Bg_res_0[p])));
							
							if (j == (m - 1)) {
								B[y + neg[p]] = B[y + neg[p]] + ((Qosinj*Rs[p] / Bo[p]) + (Qgsinj / Bg[p])) / (deltat_res * n);
							}
							if (j == 0 && zonaprodutora[k]) {
								B[y + neg[p]] = B[y + neg[p]] + C1 * Ageo_R[p] * Adir_F[p] * K_F[p] * (lambda1_G[p] * (PWB_0[k_prod] - po_res[p]) + lambda3_G[p] * (PWB_0[k_prod] - pg_res[p]));
								Qg_well[k_prod] = -ft3s_bbld(C1 * Ageo_R[p] * Adir_F[p] * K_F[p] * (lambda1_G[p] * (PWB_0[k_prod] - po_res[p]) + lambda3_G[p] * (PWB_0[k_prod] - pg_res[p])) / (deltat_res));
								Pb_WB[k_prod] = pb_res[p];
							}

							centerP = ucG[p];

							A[ucG[p]] = -C3 * Ageo_t[p] * ((So_res[p] * Rs[p] / Bo[p])*(der_phi[p] - phi[p] * der_Bo_p[p] / Bo[p] + phi[p] * der_Rs_p[p]/ Rs[p])
								+ ((Sg_res[p] / Bg[p])*(der_phi[p] - phi[p] * der_Bg[p] / Bg[p])));


							if (j == (m - 1)) {

								A[ucG[p]] = A[ucG[p]] + (Qosinj*((der_Rs_p[p] / Bo[p]) - Rs[p] * der_Bo_p[p] / (Bo[p] * Bo[p])) - Qgsinj*(der_Bg[p] / (Bg[p] * Bg[p]))) / (deltat_res * n);

							}
							if (j == 0 && zonaprodutora[k]) {

								A[ucG[p]] = A[ucG[p]] - C1*Ageo_R[p] * Adir_F[p] * K_F[p] * (lambda1_G[p] - d_lambda1_po_G[p] * (PWB_0[k_prod] - po_res[p]) + lambda3_G[p] - d_lambda3_po_G[p] * (PWB_0[k_prod] - pg_res[p]));

							}

							centerSw = ucG[p] + 1;

							A[ucG[p] + 1] = C3* Ageo_t[p] * phi[p] * ((1.0l - X_pb[p]) / Bg[p] + (X_pb[p]) * Rs[p] / Bo[p]);
							

							if (j == 0 && zonaprodutora[k]) {
								A[ucG[p] + 1] = A[ucG[p] + 1] + C1*Ageo_R[p] * Adir_F[p] * K_F[p] * (d_lambda1_SW_G[p] * (PWB_0[k_prod] - po_res[p]) + d_lambda3_SW_G[p] * (PWB_0[k_prod] - pg_res[p]) - lambda3_G[p] * der_pcg_Sw[p]);
							}

							if (NpureW[p]) {
								centerSo = ucG[p] + 2;

								A[ucG[p] + 2] = -C3*Ageo_t[p] * (X_pb[p] * (So_res[p] * phi[p] / Bo[p])*(der_Rs_pb[p] - Rs[p] * der_Bo_pb[p] / Bo[p])
									+ (1.0l - X_pb[p])*(Rs[p] * phi[p] / Bo[p])
									- (1.0l - X_pb[p])*(phi[p] / Bg[p]));
							

								if (j == 0 && zonaprodutora[k]) {

									A[ucG[p] + 2] = A[ucG[p] + 2] + C1*Ageo_R[p] * Adir_F[p] * K_F[p] * (X_pb[p] * (d_lambda1_pb_G[p] * (PWB_0[k_prod] - po_res[p]) + d_lambda3_pb_G[p] * (PWB_0[k_prod] - pg_res[p])) +
										(1.0l - X_pb[p])*(d_lambda1_SO_G[p] * (PWB_0[k_prod] - po_res[p]) + d_lambda3_SO_G[p] * (PWB_0[k_prod] - pg_res[p]) - lambda3_G[p] * der_pcg_So[p]));

								}
							}

							//West
							if (n != 1 && i != (n - 1)) {

								calc_gas_pos_res2(po_res[p], po_res[p + 1],
									pg_res[p], pg_res[p + 1],
									B, A,
									&controladorP, &controladorSw, &controladorSo,
									y, uwG[p],
									Ageo_PHI[p], Adir_W[p], C1, K_W[p],
									CHI_W[p], X_pb[p], X_pb[p + 1],
									lambda1_G[p], lambda1_G[p + 1],
									d_lambda1_po_G[p], d_lambda1_po_G[p + 1],
									d_lambda1_SW_G[p], d_lambda1_SW_G[p + 1],
									d_lambda1_pb_G[p], d_lambda1_pb_G[p + 1],
									d_lambda1_SO_G[p], d_lambda1_SO_G[p + 1],
									lambda3_G[p], lambda3_G[p + 1],
									d_lambda3_po_G[p], d_lambda3_po_G[p + 1],
									d_lambda3_SW_G[p], d_lambda3_SW_G[p + 1],
									d_lambda3_pb_G[p], d_lambda3_pb_G[p + 1],
									d_lambda3_SO_G[p], d_lambda3_SO_G[p + 1],
									der_pcg_Sw[p], der_pcg_Sw[p + 1],
									der_pcg_So[p], der_pcg_So[p + 1],
									neg[p + 1], NpureW[p + 1]);

							}

							//East(boundary)
							if (n != 1 && i == 0) {

								calc_gas_pos_res2(po_res[p], po_res[p + (n - 1)],
									pg_res[p], pg_res[p + (n - 1)],
									B, A,
									&controladorP, &controladorSw, &controladorSo,
									y, ueG[p],
									Ageo_PHI[p], Adir_E[p], C1, K_E[p],
									CHI_E[p], X_pb[p], X_pb[p + (n - 1)],
									lambda1_G[p], lambda1_G[p + (n - 1)],
									d_lambda1_po_G[p], d_lambda1_po_G[p + (n - 1)],
									d_lambda1_SW_G[p], d_lambda1_SW_G[p + (n - 1)],
									d_lambda1_pb_G[p], d_lambda1_pb_G[p + (n - 1)],
									d_lambda1_SO_G[p], d_lambda1_SO_G[p + (n - 1)],
									lambda3_G[p], lambda3_G[p + (n - 1)],
									d_lambda3_po_G[p], d_lambda3_po_G[p + (n - 1)],
									d_lambda3_SW_G[p], d_lambda3_SW_G[p + (n - 1)],
									d_lambda3_pb_G[p], d_lambda3_pb_G[p + (n - 1)],
									d_lambda3_SO_G[p], d_lambda3_SO_G[p + (n - 1)],
									der_pcg_Sw[p], der_pcg_Sw[p + (n - 1)],
									der_pcg_So[p], der_pcg_So[p + (n - 1)],
									neg[p + (n - 1)], NpureW[p + (n - 1)]);

							}

							//Backward
							if (j != (m - 1)) {

								calc_gas_pos_res2(po_res[p], po_res[p + n],
									pg_res[p], pg_res[p + n],
									B, A,
									&controladorP, &controladorSw, &controladorSo,
									y, ubG[p],
									Ageo_R[p], Adir_B[p], C1, K_B[p],
									CHI_B[p], X_pb[p], X_pb[p + n],
									lambda1_G[p], lambda1_G[p + n],
									d_lambda1_po_G[p], d_lambda1_po_G[p + n],
									d_lambda1_SW_G[p], d_lambda1_SW_G[p + n],
									d_lambda1_pb_G[p], d_lambda1_pb_G[p + n],
									d_lambda1_SO_G[p], d_lambda1_SO_G[p + n],
									lambda3_G[p], lambda3_G[p + n],
									d_lambda3_po_G[p], d_lambda3_po_G[p + n],
									d_lambda3_SW_G[p], d_lambda3_SW_G[p + n],
									d_lambda3_pb_G[p], d_lambda3_pb_G[p + n],
									d_lambda3_SO_G[p], d_lambda3_SO_G[p + n],
									der_pcg_Sw[p], der_pcg_Sw[p + n],
									der_pcg_So[p], der_pcg_So[p + n],
									neg[p + n], NpureW[p + n]);

							}
							
							//South
							if (k != (o - 1)) {

								calc_gas_pos_res2(po_res[p], po_res[p + n*m],
									pg_res[p], pg_res[p + n*m],
									B, A,
									&controladorP, &controladorSw, &controladorSo,
									y, usG[p],
									Ageo_Z[p], Adir_S[p], C1, K_S[p],
									CHI_S[p], X_pb[p], X_pb[p + n*m],
									lambda1_G[p], lambda1_G[p + n*m],
									d_lambda1_po_G[p], d_lambda1_po_G[p + n*m],
									d_lambda1_SW_G[p], d_lambda1_SW_G[p + n*m],
									d_lambda1_pb_G[p], d_lambda1_pb_G[p + n*m],
									d_lambda1_SO_G[p], d_lambda1_SO_G[p + n*m],
									lambda3_G[p], lambda3_G[p + n*m],
									d_lambda3_po_G[p], d_lambda3_po_G[p + n*m],
									d_lambda3_SW_G[p], d_lambda3_SW_G[p + n*m],
									d_lambda3_pb_G[p], d_lambda3_pb_G[p + n*m],
									d_lambda3_SO_G[p], d_lambda3_SO_G[p + n*m],
									der_pcg_Sw[p], der_pcg_Sw[p + n*m],
									der_pcg_So[p], der_pcg_So[p + n*m],
									neg[p + n*m], NpureW[p + n*m]);


								calc_gas_grav_res2(B, A,
									&controladorP, &controladorSw, &controladorSo,
									y, usG[p],
									Ageo_Z[p], C2, -g, K_S[p],
									CHI_S[p], X_pb[p], X_pb[p + n*m],
									lambda2_G[p], lambda2_G[p + n*m],
									d_lambda2_po_G[p], d_lambda2_po_G[p + n*m],
									d_lambda2_SW_G[p], d_lambda2_SW_G[p + n*m],
									d_lambda2_pb_G[p], d_lambda2_pb_G[p + n*m],
									d_lambda2_SO_G[p], d_lambda2_SO_G[p + n*m],
									lambda4_G[p], lambda4_G[p + n*m],
									d_lambda4_po_G[p], d_lambda4_po_G[p + n*m],
									d_lambda4_SW_G[p], d_lambda4_SW_G[p + n*m],
									d_lambda4_pb_G[p], d_lambda4_pb_G[p + n*m],
									d_lambda4_SO_G[p], d_lambda4_SO_G[p + n*m],
									neg[p + n*m], NpureW[p + n*m]);


							}							
							
							A[centerP] = A[centerP] + controladorP;
							A[centerSw] = A[centerSw] + controladorSw;
							if (NpureW[p]) A[centerSo] = A[centerSo] + controladorSo;
						}				

						p++;
						if (j == 0 && zonaprodutora[k]) k_prod++;
						
					}
				}

			}

			y = y + 1 + neg[p];
			for (int wer = 0; wer < y; wer++)
			{
				B[wer] = -B[wer];
			}			

			grava_ia(IA, y + 1, 0, 1);
			grava_b(B, y, 0, 1);
			grava_ja(JA, centerSw, 0, 1);
			grava_a(A, centerSw, 0, 1);

			system("pause");

			//SOLVER
			phase = 11;
			pardiso(pt, &maxfct, &mnum, &mtype, &phase,
				&y, A, IA, JA, &idum, &nrhs,
				iparm, &msglvl, &ddum, &ddum, &error, dparm);
			if (error != 0) {
				cout << "\nERROR during symbolic factorization: %d" << endl;
				dontstopmenow = 1;
			}

			phase = 33;
			pardiso(pt, &maxfct, &mnum, &mtype, &phase,
				&y, A, IA, JA, &idum, &nrhs,
				iparm, &msglvl, B, X, &error, dparm);
			if (error != 0) {
				printf("\nERROR during solution: %d", error);
				dontstopmenow = 1;
			}
			
			cout.precision(19);
			p = 0;
			for (k = 0; k <= (o - 1); k++) {
				for (j = 0; j <= (m - 1); j++) {
					for (i = 0; i <= (n - 1); i++) {

												
						po_res[p] = po_res[p] + X[3 * p + neg[p]];
						Sw_res[p] = min(0.999999999l, max(0.000000001l, Sw_res[p] + X[3 * p + 1 + neg[p]]));
						So_res[p] = X_pb[p] * (1.0l - Sw_res[p]) + (1.0l - X_pb[p])*max(So_res[p] + X[3 * p + 2 + neg[p]], 0.0000000001l);
						if(NpureW[p]) pb_res[p] = max(X_pb[p] * (pb_res[p] + X[3 * p + 2 + neg[p]]) + (1.0l - X_pb[p])*po_res[p], 15.0l);
						else pb_res[p]= pb_res[p];
																	

						pw_res[p] = po_res[p];//- calcula_prop(Sw_res[p], Sw_pcow, pcow_tab, length_pcow);
						Sg_res[p] = (1.0l - So_res[p] - Sw_res[p]);
						pg_res[p] = po_res[p];// +calcula_prop(Sg_res[p], Sg_pcgo, pcgo_tab, length_pcgo);

					
						if (fabs(X[3 * p + neg[p]] / po_res[p]) > 0.0002 || fabs(X[3 * p + 1 + neg[p]]) > 0.0002 || (1.0l - X_pb[p])*fabs(X[3 * p + 2 + neg[p]]) > 0.0002
								|| X_pb[p]*fabs(X[3 * p + 2 + neg[p]] / pb_res[p]) > 0.0002) {
								dontstopmenow = 1;						
						}									
						
						if ((po_res[p]<=pb_res[p]) && X_pb[p] > 0.5l) {
							X_pb[p] = 0.0l;							
							pb_res[p] = po_res[p];
							So_res[p] = So_res[p] - 0.0001l;//0.01l
							Sg_res[p] = 0.0001l;
							//dontstopmenow = 1;
						}
						if (Sg_res[p] <= -0.00000000000001l && X_pb[p] < 0.5l) {
							X_pb[p] = 1.0l;
							pb_res[p] = po_res[p] - 0.01l;//0.1l
							So_res[p] = 1 - Sw_res[p];
						//	dontstopmenow = 1;
						}


						if (isnan(po_res[p]) || isnan(So_res[p]) || isnan(Sw_res[p]) || isnan(pb_res[p]) ||
							po_res[p] < 0.000011l  || Sw_res[p] > 1.0l || So_res[p] > 1.00000000000000001l || fabs(X[3 * p + neg[p]])>100.l || X_pb[p] * fabs(X[3 * p + 2 + neg[p]])>50.l) {
							//|| fabs(X[3 * p + neg[p]])>100.l|| X_pb[p] * fabs(X[3 * p + 2 + neg[p]])>50.l
							//|| po_res[p]>(po_res_0[p] + 10.l) || fabs(pb_res[p] - pb_res_0[p])>50.0l
							//||po_res[p]>(po_res_0[p] + 10.l)|||| pb_res[p] < 500l


						/*	if (isnan(po_res[p])) cout << "erro: " << 1 << endl;
							if (isnan(So_res[p])) cout << "erro: " << 2 << endl;
							if (isnan(Sw_res[p])) cout << "erro: " << 3 << endl;
							if (isnan(pb_res[p])) cout << "erro: " << 4 << endl;
							if (po_res[p] < 0.000011l) cout << "erro: " << 5 << endl;
							if (pb_res[p] < 500l) cout << "erro: " << 6 << endl;
							if (Sw_res[p] > 1.0l) cout << "erro: " << 7 << endl;
							if (So_res[p] > 1.00000000000000001l) cout << "erro: " << 8 << endl;*/
							//if (fabs(X[3 * p + neg[p]]) > 100.l) {
							//	cout << "erro: " << 9 << endl;
							//	cout << "deltapb: " << X[3 * p + neg[p]] << endl;
							//if (level>5) {
							//	
							//		grava_p(po_res, dim, 0, (t + 1), ite_coup, coup_ts);
							//		grava_pb(pb_res, dim, 0, (t + 1), ite_coup, coup_ts);
							//		grava_Sw(Sw_res, dim, 0, (t + 1), ite_coup, coup_ts);
							//		grava_So(So_res, dim, 0, (t + 1), ite_coup, coup_ts);
							//		cout << "po: " << po_res[p] << endl;
							//		cout << "SO: " << So_res[p] << endl;
							//		cout << "Sg: " << Sg_res[p] << endl;
							//		cout << "SW: " << Sw_res[p] << endl;

							//		cout << "po0: " << po_res_0[p] << endl;
							//		cout << "So0: " << So_res_0[p] << endl;
							//		cout << "Sg0: " << Sg_res_0[p] << endl;
							//		cout << "Sw0: " << Sw_res_0[p] << endl;									

							//		cout << "PWB: " << PWB_0[0] << endl;
							//		cout << "K_F: " << K_F[p]<<endl;
							//		cout << "kro: " << kro[p] << endl;
							//		cout << "krg: " << krg[p] << endl;
							//		cout << "krw: " << krw[p] << endl;
							//		cout << "Rs: " << Rs[p] << endl;
							//		cout << "Bo: " << Bo[p] << endl;
							//		cout << "Bg: " << Bg[p] << endl;
							//		cout << "Bw: " << Bw[p] << endl;
							//		cout << "muo: " << mi_o[p] << endl;
							//		cout << "mug: " << mi_g[p] << endl;
							//		cout << "muw: " << mi_w[p] << endl;
							//		cout << "K_F: " << K_F[p] << endl;

							//		cout << "Rs0: " << Rs_res_0[p] << endl;
							//		cout << "Bo0: " << Bo_res_0[p] << endl;
							//		cout << "Bg0: " << Bg_res_0[p] << endl;
							//		cout << "Bw0: " << Bw_res_0[p] << endl;
							//		
							//		system("pause");
							//}
							////}
							//if (X_pb[p] * fabs(X[3 * p + 2 + neg[p]])>50.l) {
							//	cout << "erro: " << 10 << endl;
							//	cout << "deltapb: " << X[3 * p + 2] << endl;								
							//}


							nan = 1;

							i = n*o*m;
							j = n*o*m;
							k = n*o*m;
							dontstopmenow = 0;

						}
						p++;

					}
				}
			}
			
			l++;			

			if (l > 5) {
				dontstopmenow = 0;
			}

		}

		
		p = 0;
		for (k = 0; k <= (o - 1); k++) {
			for (j = 0; j <= (m - 1); j++) {
				for (i = 0; i <= (n - 1); i++) {
					if (abs(So_res[p] - So_res_0[p]) >= 0.05 || abs(Sw_res[p] - Sw_res_0[p]) >= 0.05) {
						nan = 1;
						p++;
					}
				}
			}
		}

		
		p = 0;
		if (nan == 1 || l > 5) {
						
			double factor = 0;
			if (level == 0) factor = pow(0.5l, 2.0l);
			if (level > 0) factor = 0.5l;


			for (k = 0; k <= (o - 1); k++) {
				for (j = 0; j <= (m - 1); j++) {
					for (i = 0; i <= (n - 1); i++) {
						varpo[p] = factor*varpo[p];
						varpb[p] = factor*varpb[p];
						varSo[p] = factor*varSo[p];
						varSw[p] = factor*varSw[p];
						po_res[p] = po_res_0[p] - varpo[p];
						Sw_res[p] = Sw_res_0[p] - varSw[p];
						So_res[p] = So_res_0[p] - varSo[p];
						pb_res[p] = pb_res_0[p] - varpb[p];
						Sg_res[p] = 1.0l - So_res[p] - Sw_res[p];
						pw_res[p] = po_res[p]; //- calcula_prop(Sw_res[p], Sw_pcow, pcow_tab, length_pcow);
						pg_res[p] = po_res[p]; //+ calcula_prop(Sg_res[p], Sg_pcgo, pcgo_tab, length_pcgo);
						if (po_res[p] <= pb_res[p]) X_pb[p] = 0.0l;
						else X_pb[p] = 1.0l;

						p++;

					}
				}
			}

			deltat_res = factor*deltat_res;

			if (isnan(factor)) {
				cout << "factor: " << factor << endl;
				system("pause");
			}

			calc_coef_grid(Ageo_t, Ageo_R, Ageo_PHI, Ageo_Z,
				Adir_F, Adir_B, Adir_W, Adir_E, Adir_N, Adir_S,
				CHI_F, CHI_B, CHI_W, CHI_E, CHI_N, CHI_S,
				R_P, PHI_P, Z_P,
				R_PB, R_PF, delta_PHI, delta_Z, deltat_res,
				K_F, K_B, K_N, K_S, K_W, K_E,
				Kr, Kz,
				o, m, n);

			calc_reser_evol(int(1 / factor), p, y, u, l, centerP, centerSw, centerSo,
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
				Qw_well, Qo_well, Qg_well, coup_ts, ite_coup,
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
				mi_g, der_mi_g,
				Bg, der_Bg, 
				rho_g, der_rho_g,
				K_F, K_B, K_N, K_S, K_W, K_E,
				X_pb, zonaprodutora, (level + 1),
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


			nan = 0;
			deltat_res = deltat_res / factor;

			calc_coef_grid(Ageo_t, Ageo_R, Ageo_PHI, Ageo_Z,
				Adir_F, Adir_B, Adir_W, Adir_E, Adir_N, Adir_S,
				CHI_F, CHI_B, CHI_W, CHI_E, CHI_N, CHI_S,
				R_P, PHI_P, Z_P,
				R_PB, R_PF, delta_PHI, delta_Z, deltat_res,
				K_F, K_B, K_N, K_S, K_W, K_E,
				Kr, Kz,
				o, m, n);

			for (int iy = 0; iy<dim; iy++) {
				
				varpo[p] = varpo[p] / factor;
				varpb[p] = varpb[p] / factor;
				varSo[p] = varSo[p] / factor;
				varSw[p] = varSw[p] / factor;

			}
		}
		
		if (level >= 2) {

			Sleep(3 * 100);
		}
			
		p = 0;
		for (k = 0; k <= (o - 1); k++) {
			for (j = 0; j <= (m - 1); j++) {
				for (i = 0; i <= (n - 1); i++) {
					
					po_res_0[p] = po_res[p];
					Sw_res_0[p] = Sw_res[p];
					So_res_0[p] = So_res[p];
					Sg_res_0[p] = Sg_res[p];
					pb_res_0[p] = pb_res[p];

					p++;
					
				}
			}
		}
				
		if (level == 0) {

			if ((t == 0 || t == (res_timesteps - 1)) && coup_ts % 50 == 0) {
				//
				grava_var(po_res, dim, 0, (t + 1), ite_coup, directory, 1);
				grava_var(po_res, dim, 0, (t + 1), ite_coup, directory, 2);
				grava_var(po_res, dim, 0, (t + 1), ite_coup, directory, 3);
				grava_var(po_res, dim, 0, (t + 1), ite_coup, directory, 4);
				
			}

			for (k = 0; k < N_wellbore; k++) {
				PWB_0[k] = PWB_0[k] - 0.5l*deltaPWB[k] / res_timesteps;
			}

			for (iy = 0; iy < dim; iy++) {
				varpo[iy] = po_res_lv0[iy] - po_res[iy];
				varpb[iy] = pb_res_lv0[iy] - pb_res[iy];
				varSo[iy] = So_res_lv0[iy] - So_res[iy];
				varSw[iy] = Sw_res_lv0[iy] - Sw_res[iy];
				po_res[iy] = po_res[iy] - (varpo[iy]);
				pb_res[iy] = pb_res[iy] - (varpb[iy]);
				So_res[iy] = So_res[iy] - (varSo[iy]);
				Sw_res[iy] = Sw_res[iy] - (varSw[iy]);
				Sg_res[iy] = 1.l - So_res[iy] - Sw_res[iy];

			}

			for (iy = 0; iy < dim; iy++) {

				po_res_lv0[iy] = po_res_0[iy];
				pb_res_lv0[iy] = pb_res_0[iy];
				So_res_lv0[iy] = So_res_0[iy];
				Sw_res_lv0[iy] = Sw_res_0[iy];
			}

		}

		dontstopmenow = 1;

	}

	free(deltaPWB);

}