#include "stdafx.h"
#include "reser_func.h"

double calcula_THP(double tempo, double *tempo_THP, double *THP) {
	int k = 0;
	while(tempo>tempo_THP[k]) {
		k++;
	}
	return (THP[k] - THP[k - 1])*(tempo - tempo_THP[k - 1]) / (tempo_THP[k] - tempo_THP[k - 1]) + THP[k - 1];
}

void ini_fields(int m, int n, int o, double yo, double rho_ws,
	double *So_coup, double *Sw_coup, double *Sg_coup, 
	double *po_coup, double *pw_coup, double *pg_coup, double *pb_coup, bool *acimapb,
	double *So_layer, double *Sw_layer, double *deltaz, double pi, double pbi) {
		
	
	int p, k, j, i;
	p = 0;
	for (k = 0; k <= (o - 1); k++) {
		for (j = 0; j <= (m - 1); j++) {
			for (i = 0; i <= (n - 1); i++) {

				Sw_coup[p] = Sw_layer[k];
				So_coup[p] = So_layer[k];		
				Sg_coup[p] = 1.0l - So_coup[p] - Sw_coup[p];				
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
					acimapb[p] = 1;
				}
				else {
					acimapb[p] = 0;
				}
				

				p++;
			}
		}

	}
}

void inicializa_reser(int *res_timesteps, double *H, double *Re, double *Rw, int *n, int *m, int *o, double *deltaphi, double *deltaz, double *a, double *R0, double *pi, double *pbi, double *Swi,
	double *Soi, double *T, double *g, double *C1, double *C2, double *C3, double *C4, double *C5, double *C6, double *Qwsinj, double *Qosinj, double *Qgsinj, double *rho_air, double *rho_ws, double *API, double *yg,
	double *yo, double *rho_os, double *rho_gs, double *yn2, double *yco2, double *yh2s, double *phi_ref, double *Cr, double *pgr, double *Swc, double *Sor, double *Sgr, double *epsilon, double *pcwmin, double *pcwmax, double *pcgmin,
	double *pcgmax, double *krwmax, double *nw, double *ng, double *now, double *nog, double *parametrosBw, double *parametrosmig, double *ppc, double *Tpc, double *parametrosZ, double *Bwi, double *Bpcw, double *Bpcg, bool *zonaprodutora) {

	*H = 40l;//altura do reser(ft)
	*Re = 1500;//raio externo(ft)
	*Rw = 0.145833333333333l;//raio interno(ft)

	*n = 1;//(divisoes do anel)
	*m = 20;//(profundidade)
	*o = 29;//(pra cima)

	*res_timesteps = 15;	
		
	*deltaphi = 2.0l * 3.14159265358979323846l / *n;

	for (int k = 0; k <= (*o - 1); k++) 
	{	
			if (k < 20) deltaz[k] = 5.0l;
			else deltaz[k] = 40.0l;

			if (k < 4) zonaprodutora[k] = 1;
			else zonaprodutora[k] = 0;

	}

	*a = pow((*Re / *Rw), (1.0 / *m)) - 1;//(fator de expansão)
	*R0 = log(*a + 1)**Rw / *a;

	//Condição inicial Inicial
	*pi = 3600;
	*pbi = 2000;
	*Swi = 0.2000000l;
	*Soi = 0.8000000l;
	*T = Kel_Far(273.15l+75.0l);//°F
	*g = 32.174049l;//ft/s2


	*C1 = 7.324646209006036l*pow(10.0l, -8.0l);
	*C2 = 1.602205l*pow(10.0l, -11.0l);
	*C3 = 1.304575164438286l*pow(10.0l, -8.0l);
	*C4 = 2.85365l*pow(10.0l, -12.0l);
	*C5 = 1.0l / 5.61458l;
	*C6 = 1.l;


	//Injeção
	*Qwsinj = 0.0l;//ft³/s
	*Qosinj = 0.0l*0.0000000005l*(2.0l * 3.14159265358979323846l**H**Re);//ft³/s
	*Qgsinj = 0.0l*0.0000000000000000000001l*(2.0l * 3.14159265358979323846l**H**Re);//ft³/s
	

	//Densidades
	*rho_air = 0.076399332l;//lb/ft3
	*rho_ws = 62.30959092l; //lb/ft³
	*API = 35.3884l;//
	*yg = 0.73670l;//
	*yo = 141.5l / (131.5l + *API);//
	*rho_os = *yo**rho_ws;//lb/ft³
	*rho_gs = *yg**rho_air / 5.61458l;//lb*bbl/ft^6 (sim, exatamente)
	*yn2 = 0.0l;
	*yco2 = 0.0l;
	*yh2s = 0.0l;

	//Parametros phi
	*phi_ref =25.0l / 100.0l;
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

	//Variável auxiliar do Z
	*pgr = 666.l;

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
	int *unG, int *ufG, int *uwbG, int *ueG, int *ucG, int *uwG, int *uebG, int *ubG, int *usG)
{
	int u = 0;
	int p = 0;

	for (int k = 0; k <= (o - 1); k++) {
		for (int j = 0; j <= (m - 1); j++) {
			for (int i = 0; i <= (n - 1); i++) {

				IA[3 * p] = u + 1;

				//Water		
				{
					//North
					if (k != 0) {
						unW[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3 * n*m;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3 * n*m;
						u = u + 1;
					}

					//Foward
					if (j != 0) {
						ufW[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3 * n;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3 * n;
						u = u + 1;
					}

					//West(boundary)
					if (n != 1 && i == (n - 1)) {
						uwbW[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3 * (n - 1);
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3 * (n - 1);
						u = u + 1;
					}

					//East
					if (n != 1 && i != (0)) {
						ueW[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3;
						u = u + 1;
					}

					//Center
					{
						ucW[p] = u;
						JA[u] = 3 * (p + 1) - 2;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1;
						u = u + 1;
					}

					//West
					if (n != 1 && i != (n - 1)) {
						uwW[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3;
						u = u + 1;
					}

					//East(boundary)
					if (n != 1 && i == 0) {
						uebW[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3 * (n - 1);
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3 * (n - 1);
						u = u + 1;
					}

					//Backward
					if (j != (m - 1)) {
						ubW[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3 * n;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3 * n;
						u = u + 1;
					}

					//South
					if (k != (o - 1)) {
						usW[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3 * n*m;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3 * n*m;
						u = u + 1;
					}


				}

				IA[3 * p + 1] = u + 1;


				//Oil
				{
					//North
					if (k != 0) {
						unO[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3 * n*m;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3 * n*m;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 0 - 3 * n*m;
						u = u + 1;
					}

					//Foward
					if (j != 0) {
						ufO[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3 * n;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3 * n;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 0 - 3 * n;
						u = u + 1;
					}

					//West(boundary)
					if (n != 1 && i == (n - 1)) {
						uwbO[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3 * (n - 1);
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3 * (n - 1);
						u = u + 1;
						JA[u] = 3 * (p + 1) - 0 - 3 * (n - 1);
						u = u + 1;
					}

					//East
					if (n != 1 && i != (0)) {
						ueO[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 0 - 3;
						u = u + 1;
					}

					//Center
					{
						ucO[p] = u;
						JA[u] = 3 * (p + 1) - 2;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1;
						u = u + 1;
						JA[u] = 3 * (p + 1);
						u = u + 1;
					}

					//West
					if (n != 1 && i != (n - 1)) {
						uwO[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 0 + 3;
						u = u + 1;
					}

					//East(boundary)
					if (n != 1 && i == 0) {
						uebO[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3 * (n - 1);
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3 * (n - 1);
						u = u + 1;
						JA[u] = 3 * (p + 1) - 0 + 3 * (n - 1);
						u = u + 1;
					}

					//Backward
					if (j != (m - 1)) {
						ubO[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3 * n;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3 * n;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 0 + 3 * n;
						u = u + 1;
					}

					//South
					if (k != (o - 1)) {
						usO[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3 * n*m;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3 * n*m;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 0 + 3 * n*m;
						u = u + 1;
					}
				}

				IA[3 * p + 2] = u + 1;

				//Gas	
				{
					//North
					if (k != 0) {
						unG[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3 * n*m;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3 * n*m;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 0 - 3 * n*m;
						u = u + 1;
					}

					//Foward
					if (j != 0) {
						ufG[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3 * n;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3 * n;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 0 - 3 * n;
						u = u + 1;
					}

					//West(boundary)
					if (n != 1 && i == (n - 1)) {
						uwbG[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3 * (n - 1);
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3 * (n - 1);
						u = u + 1;
						JA[u] = 3 * (p + 1) - 0 - 3 * (n - 1);
						u = u + 1;
					}

					//East
					if (n != 1 && i != (0)) {
						ueG[p] = u;
						JA[u] = 3 * (p + 1) - 2 - 3;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 - 3;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 0 - 3;
						u = u + 1;
					}

					//Center
					{
						ucG[p] = u;
						JA[u] = 3 * (p + 1) - 2;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1;
						u = u + 1;
						JA[u] = 3 * (p + 1);
						u = u + 1;
					}

					//West
					if (n != 1 && i != (n - 1)) {
						uwG[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 0 + 3;
						u = u + 1;
					}

					//East(boundary)
					if (n != 1 && i == 0) {
						uebG[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3 * (n - 1);
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3 * (n - 1);
						u = u + 1;
						JA[u] = 3 * (p + 1) - 0 + 3 * (n - 1);
						u = u + 1;
					}

					//Backward
					if (j != (m - 1)) {
						ubG[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3 * n;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3 * n;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 0 + 3 * n;
						u = u + 1;
					}

					//South
					if (k != (o - 1)) {
						usG[p] = u;
						JA[u] = 3 * (p + 1) - 2 + 3 * n*m;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 1 + 3 * n*m;
						u = u + 1;
						JA[u] = 3 * (p + 1) - 0 + 3 * n*m;
						u = u + 1;
					}
				}

				p = p + 1;
			}
		}

	}
	IA[3 * p] = u + 1;
}

//BLZ
void calc_water_pos_res(double pwp,double pwnm,
	double *B, double *A, 
	double *controladorP, double *controladorSw, 
	int y, int u,
	double pesop, double pesocon,double Ageo, double C,
	double lambda1, double lambda2, double lambda3, double lambda4, double lambda5,
	double lambda1int, double lambda2int, double lambda3int, double lambda4int, double lambda5int)
{
	B[y] = B[y] - Ageo*C*(pesop*lambda1 + pesocon*lambda1int)*(pwp - pwnm);
		
	A[u] = Ageo*C*(pesop*lambda1+ pesocon*lambda1int - pesocon*lambda2int*(pwp - pwnm));
	
	*controladorP = *controladorP + Ageo*C*(-pesocon*lambda1int - pesop*lambda1 - pesop*lambda2*(pwp - pwnm));

	u = u + 1;
	
	A[u] = Ageo*C*(-pesocon*lambda5int*(pwp - pwnm) - (pesocon*lambda4int+ pesop*lambda1*lambda3int));
	
	*controladorSw = *controladorSw + Ageo*C*(-pesop*lambda5*(pwp - pwnm)+ pesop*lambda4 + pesocon*lambda1int*lambda3);
}

//BLZ
void calc_water_grav_res(double rhonm, double rhop, double drhonm, double drhop,
	double *B, double *A, 
	double *controladorP, double *controladorSw, 
	int y, int u,
	double pesop, double pesocon, double Ageo, double C,
	double lambda1, double lambda2, double lambda5,
	double lambda1int, double lambda2int, double lambda5int,
	double deltaz, double g)
{	
	B[y] = B[y] + Ageo*C*deltaz*g*(pesocon*lambda1int*rhonm + pesop*lambda1*rhop);
	
	A[u] = A[u] + Ageo*C*deltaz*g*pesocon*(lambda1int*drhonm + lambda2int*rhonm);

	*controladorP = *controladorP + Ageo*C*deltaz*g*pesop*(lambda1*drhop + lambda2*rhop);

	u = u + 1;
	
	A[u] = A[u]+Ageo*C*deltaz*g*pesocon*(rhonm*lambda5int);

	*controladorSw = *controladorSw+ Ageo*C*deltaz*g*pesop*(rhop*lambda5);
}


//BLZ
void calc_oil_pos_res(double pop, double ponm, 
	double *B, double *A,
	double *controladorP, double *controladorSw, double *controladorSo,
	int y, int u, 
	double pesop, double pesocon, double Ageo, double C,
	double lambda1, double lambda2, double lambda3, double lambda4,
	double lambda1int, double lambda2int, double lambda3int, double lambda4int)
{
	B[y] = B[y] - Ageo*C*(pesocon*lambda1int + pesop*lambda1)*(pop - ponm);
		
	A[u] = Ageo*C*(pesop*lambda1 + pesocon*lambda1int - pesocon*lambda2int*(pop - ponm));

	*controladorP = *controladorP + Ageo*C*(-pesop*lambda1 - pesocon*lambda1int - pesop*lambda2*(pop - ponm));

	u = u + 1;
	
	A[u] = -Ageo*C*pesocon*lambda3int*(pop - ponm);

	*controladorSw = *controladorSw - Ageo*C*pesop*lambda3*(pop - ponm);

	u = u + 1;
	
	A[u] = -Ageo*C*pesocon*lambda4int*(pop - ponm);

	*controladorSo = *controladorSo - Ageo*C*pesop*lambda4*(pop - ponm);
}

//BLZ
void calc_oil_grav_res(double rhop, double rhonm, double drhop_p, double drhonm_p, double drhop_pb, double drhonm_pb,
	double *B, double *A, 
	double *controladorP, double *controladorSw, double *controladorSo, 
	int y, int u,
	double pesop, double pesocon, double Ageo, double C,
	double lambda1, double lambda2, double lambda3, double lambda4,
	double lambda1int, double lambda2int, double lambda3int, double lambda4int,
	double deltaz, double g)
{	
	B[y] = B[y] + Ageo*C*deltaz*g*(pesocon*lambda1int*rhonm + pesop*rhop*lambda1);

	
	A[u] = A[u] + Ageo*C*deltaz*g*pesocon*(lambda1int*drhonm_p + lambda2int*rhonm);

	*controladorP = *controladorP + Ageo*C*deltaz*g*pesop*(lambda1*drhop_p + lambda2*rhop);
	
	u = u + 1;
	
	A[u] = A[u] + Ageo*C*deltaz*g*pesocon*lambda3int*rhonm;

	*controladorSw = *controladorSw + Ageo*C*deltaz*g*pesop*lambda3*rhop;

	u = u + 1;
	
	A[u] = A[u] + Ageo*C*deltaz*g*pesocon*(lambda1int*drhonm_pb + lambda4int*rhonm);
		
	*controladorSo = *controladorSo + Ageo*C*deltaz*g*pesop*(lambda1*drhop_pb + lambda4*rhop);
}


//BLZ
void calc_gas_pos_res(double pop, double ponm, double pgp, double pgnm,
	double *B, double *A, 
	double *controladorP, double *controladorSw, double *controladorSo, 
	int y, int u,
	double pesop, double pesocon, double Ageo, double C,
	double lambda1, double lambda2, double lambda3, double lambda4, double lambda5, double lambda6, double lambda7, double lambda8, double dpcgp,
	double lambda1int, double lambda2int, double lambda3int, double lambda4int, double lambda5int, double lambda6int, double lambda7int, double lambda8int, double dpcgnm,
	bool acimapbp, bool acimapbnm)
{	
	B[y] = B[y] - Ageo*C*((pesocon*lambda1int + pesop*lambda1)*(pop - ponm));
	
	//if (acimapbp == 0 && acimapbnm == 0) {
		B[y] = B[y] - Ageo*C*((pesocon*lambda3int + pesop*lambda3)*(pgp - pgnm));
	//}

	A[u] = Ageo*C*((pesocon*lambda1int+ pesop*lambda1) - pesocon*lambda2int*(pop - ponm));

	//if (acimapbp == 0 && acimapbnm == 0) {
		A[u] = A[u]+ Ageo*C*((pesocon*lambda3int + pesop*lambda3) - pesocon*lambda4int*(pgp - pgnm));
	//}

	*controladorP = *controladorP + Ageo*C*(-(pesocon*lambda1int + pesop*lambda1) - pesop*lambda2*(pop - ponm));

	//if (acimapbp == 0 && acimapbnm == 0) {
		*controladorP = *controladorP + Ageo*C*( - (pesocon*lambda3int + pesop*lambda3) - pesop*lambda4*(pgp - pgnm));
	//}

	u = u + 1;

	A[u] = -Ageo*C*pesocon*lambda5int*(pop - ponm);

	//if (acimapbp == 0 && acimapbnm == 0) {
		A[u] = A[u] - Ageo*C*pesocon*lambda7int*(pgp - pgnm) + Ageo*C*(pesocon*lambda3int + pesop*lambda3)*dpcgnm;
	//}
	

	*controladorSw = *controladorSw - Ageo*C*(pesop*lambda5*(pop - ponm));

	//if (acimapbp == 0 && acimapbnm == 0) {
		*controladorSw = *controladorSw + Ageo*C*(- pesop*lambda7*(pgp - pgnm) - (pesop*lambda6 + pesocon*lambda3int*dpcgp));
	//}

	u = u + 1;

	A[u] = Ageo*C*(-pesocon*lambda8int*(pop - ponm));

	//if (acimapbp == 0 && acimapbnm == 0) {
		A[u] = A[u]+ Ageo*C*( - pesocon*lambda7int*(pgp - pgnm) + (pesocon*lambda3int + pesop*lambda3)*dpcgnm);
	//}

	*controladorSo = *controladorSo + Ageo*C*(-pesop*lambda8*(pop - ponm));
	
	//if (acimapbp == 0 && acimapbnm == 0) {
		*controladorSo = *controladorSo + Ageo*C*(-(pesop*lambda6 + pesocon*lambda3int*dpcgp) - pesop*lambda7*(pgp - pgnm));
	//}

}

void calc_gas_grav_res(double rhoop, double drhoop, double rhoonm, double drhoonm, double rhogp, double drhogp, double rhognm, double drhognm, 
	double drhoopbp, double drhoopbnm,
	double *B, double *A,
	double *controladorP, double *controladorSw, double *controladorSo, 
	int y, int u,
	double pesop, double pesocon, double Ageo, double C,
	double lambda1, double lambda2, double lambda3, double lambda4, double lambda5, double lambda7, double lambda8, 
	double lambda1int, double lambda2int, double lambda3int, double lambda4int, double lambda5int, double lambda7int, double lambda8int,
	double deltaz, double g, bool acimapbp, bool acimapbnm)
{
	B[y] = B[y] + Ageo*C*deltaz*g*(pesocon*lambda1int*rhoonm + pesop*lambda1*rhoop);

	//if (acimapbp == 0 && acimapbnm == 0) {
		B[y] = B[y] + Ageo*C*deltaz*g*(pesocon*lambda3int*rhognm + pesop*lambda3*rhognm);
	//}

	A[u] = A[u] + Ageo*C*deltaz*g*pesocon*(lambda1int*drhoonm + rhoonm*lambda2int);
	
	//if (acimapbp == 0 && acimapbnm == 0) {
		A[u] = A[u] + Ageo*C*deltaz*g*pesocon*(lambda3int*drhognm + lambda4int*rhognm);
	//}

	*controladorP = *controladorP + Ageo*C*deltaz*g*pesop*(lambda1*drhoop + lambda2*rhoop);
	
	//if (acimapbp == 0 && acimapbnm == 0) {
		*controladorP = *controladorP + Ageo*C*deltaz*g*pesocon*(lambda3*drhogp + lambda4*rhogp);
	//}

	u = u + 1;

	A[u] = A[u] + Ageo*C*deltaz*g*pesocon*lambda5int*rhoonm;

	//if (acimapbp == 0 && acimapbnm == 0) {
		A[u] = A[u] + Ageo*C*deltaz*g*pesocon*lambda7int*rhognm;
	//}

	*controladorSw = *controladorSw + Ageo*C*deltaz*g*pesop*lambda5*rhoop;

	//if (acimapbp == 0 && acimapbnm == 0) {
		*controladorSw = *controladorSw + Ageo*C*deltaz*g*pesop*lambda7*rhogp;
	//}

	u = u + 1;

	A[u] = A[u] + Ageo*C*deltaz*g*pesocon*lambda8int*rhoonm+ Ageo*C*deltaz*g*pesocon*lambda1int*drhoopbnm;

	//if (acimapbp == 0 && acimapbnm == 0) {
		A[u] = A[u] + Ageo*C*deltaz*g*pesocon*lambda7int*rhognm;
	//}
	
	*controladorSo = *controladorSo + Ageo*C*deltaz*g*pesop*lambda8*rhoop + Ageo*C*deltaz*g*pesop*lambda1*drhoopbp;

	//if (acimapbp == 0 && acimapbnm == 0) {
		*controladorSo = *controladorSo +  Ageo*C*deltaz*g*pesop*lambda7*rhogp;
	//}
}


void calc_reser_evol(int res_timesteps, int p, int y, int u, int l, int centerP, int centerSw, int centerSo,
	int *IA, int *JA, double *A, double *B, double *X,
	bool dontstopmenow, int i, int k, int j, int o, int m, int n, int dim, 
	double controladorP, double controladorSw, double controladorSo,
	double *lambda1_w, double *lambda2_w, double *lambda4_w, double *lambda5_w,
	double *lambda1_o, double *lambda2_o, double *lambda3_o, double *lambda4_o,
	double *lambda1_g, double *lambda2_g, double *lambda3_g, double *lambda4_g, 
	double *lambda5_g, double *lambda6_g, double *lambda7_g, double *lambda8_g,
	double *Kx, double *Kz,double pi, double T, double *PWB ,double *PWB_0, double *Pb_WB,
	double *po_res, double *pb_res,double *pw_res, double *pg_res, double *Sw_res, double *So_res, double *Sg_res,
	double *po_res_0, double *pb_res_0, double *Sw_res_0, double *So_res_0, 
	double *po_res_lv0, double *pb_res_lv0, double *So_res_lv0, double *Sw_res_lv0,
	double *varpo, double *varpb, double *varSo, double *varSw,
	double C1, double C2, double C3, double C4, double C5, double C6,
	double deltaphi, double *deltaz, double deltat, double a, 
	double *Ax, double *Aphi, double *Axint, double *Az, double *At, 
	double g, double H, double Re, double Rw, double R0, double pesoF, double pesoB,
	double rho_air, double rho_ws, double API, double yg, double yo, double rho_os, double rho_gs, double yn2, double yco2, double yh2s,
	double pgr, double krwmax, double nw, double ng, double now, double nog, double phi_ref, double Cr,
	double Swc, double Sor, double Sgr, double epsilon, double pcwmin, double pcwmax, double pcgmin, double pcgmax,
	double *parametrosBw, double *parametrosmig, double *parametrosZ, double ppc, double Tpc, double Bwi, double Bpcw, double Bpcg, 
	double Qwsinj, double Qosinj, double Qgsinj,
	int *iparm, double *dparm, int solver, int maxfct, int mnum, int phase, int error,
	int msglvl, int mtype, int nrhs, double ddum, int idum, void *pt,
	double *Qw_well, double *Qo_well, double *Qg_well, int coup_ts, int ite_coup,
	double *phi, double *der_phi, double *rho_w, double *der_rho_w, double *mi_w, double *der_mi_w,
	double *Bw, double *der_Bw, double *pcw, double *der_pcw, double *krw, double *der_krw,
	double *Rs, double *der_Rs_p, double *der_Rs_pb, double *kro, double *der_kro_Sw, double *der_kro_So, 
	double *Bo, double *der_Bo_p, double *der_Bo_pb,
	double *mi_o, double *der_mi_o_p, double *der_mi_o_pb, double *rho_o, double *der_rho_o_p, double *der_rho_o_pb,
	double *pcg, double *der_pcg_Swo, double *krg,double *der_krg, double *mi_g, double *der_mi_g, double *Z,
	double *Bg, double *der_Bg, double *rho_g, double *der_rho_g,
	bool *acimapb, bool *zonaprodutora, int level,
	int *unW, int *ufW, int *uwbW, int *ueW, int *ucW, int *uwW, int *uebW, int *ubW, int *usW,
	int *unO, int *ufO, int *uwbO, int *ueO, int *ucO, int *uwO, int *uebO, int *ubO, int *usO,
	int *unG, int *ufG, int *uwbG, int *ueG, int *ucG, int *uwG, int *uebG, int *ubG, int *usG)
{	
	
	
	bool teste = 0;
	int errou = 0;
	int iy = 0;
	bool nan = 0;
	double *Azinterno = (double *)malloc(dim * sizeof(double));
	
	
	for (int iy = 0; iy<dim; iy++) {
		Azinterno[iy] = Az[iy];		
	}

	if (level == 0) {
		for (k = 0; k <= (o - 1); k++) PWB_0[k] = PWB_0[k] - 0.5l*(PWB_0[k] - PWB[k]) / res_timesteps;
	}
	
	for (int t = 0; t < res_timesteps; t++) 
	{
		//cout << "1" << endl;
		if (level == 0) {

			for (iy = 0; iy < dim; iy++) {

				po_res_lv0[iy] = po_res_0[iy];
				pb_res_lv0[iy] = pb_res_0[iy];
				So_res_lv0[iy] = So_res_0[iy];
				Sw_res_lv0[iy] = Sw_res_0[iy];
			}

		}

		//cout << "2" << endl;
		l = 0;
		dontstopmenow = 1;
		while (dontstopmenow)
		{
			nan = 0;
			dontstopmenow = 0;
			p = 0;					
			
		
			calc_props_grid(dim, po_res, pb_res, Sw_res, So_res, Sg_res, phi, der_phi, rho_w, der_rho_w, mi_w, der_mi_w, Bw, der_Bw, pcw, der_pcw, krw, der_krw,
				Rs, der_Rs_p, der_Rs_pb, kro, der_kro_Sw, der_kro_So, Bo, der_Bo_p, der_Bo_pb,
				mi_o, der_mi_o_p, der_mi_o_pb, rho_o, der_rho_o_p, der_rho_o_pb,
				pcg, der_pcg_Swo, krg, der_krg, mi_g, der_mi_g, Z, Bg, der_Bg, rho_g, der_rho_g,
				phi_ref, Cr, Bwi, rho_ws, pi, T, parametrosBw,
				Bpcw, Bpcg, API, yo, yg, rho_os, rho_gs,
				parametrosmig, T / Tpc, ppc, Tpc, pcwmin, pcwmax, acimapb);
					
			//WATER
			{
				for (p = 0; p <= (dim - 1); p++)
				{
					lambda1_w[p] = krw[p] / (mi_w[p] * Bw[p]);
					if (isnan(lambda1_w[p])) {
						cout << "lambda1_w: " << lambda1_w[p] << endl;
						cout << "po_res: " << po_res[p] << endl;
						cout << "Sw_res: " << Sw_res[p] << endl;
						cout << "So_res: " << So_res[p] << endl;
						cout << "pb_res: " << pb_res[p] << endl;
						//system("pause");;
					}

				}

				for (p = 0; p <= (dim - 1); p++)
				{
					lambda2_w[p] = -((der_mi_w[p] / mi_w[p]) + (der_Bw[p] / Bw[p]))*lambda1_w[p];
					if (isnan(lambda2_w[p])) {
						cout << "lambda2_w: " << lambda2_w[p] << endl;
						cout << "po_res: " << po_res[p] << endl;
						cout << "Sw: " << Sw_res[p] << endl;
						cout << "So_res: " << So_res[p] << endl;
						cout << "pb_res: " << pb_res[p] << endl;						
					}
					
				}

				//lambda3 é a maldita derivada da pressão capilar
				for (p = 0; p <= (dim - 1); p++)
				{
					lambda4_w[p] = der_pcw[p] * lambda1_w[p];
					if (isnan(lambda4_w[p])) {
						cout << "lambda4_w: " << lambda4_w[p] << endl;
						cout << "po_res: " << po_res[p] << endl;
						cout << "Sw_res: " << Sw_res[p] << endl;
						cout << "So_res: " << So_res[p] << endl;
						cout << "pb_res: " << pb_res[p] << endl;						
					}
				}

				for (p = 0; p <= (dim - 1); p++)
				{
					lambda5_w[p] = der_krw[p] / (mi_w[p] * Bw[p]);
					if (isnan(lambda5_w[p])) {
						cout << "lambda5_w: " << lambda5_w[p] << endl;
						cout << "po_res: " << po_res[p] << endl;
						cout << "Sw_res: " << Sw_res[p] << endl;
						cout << "So_res: " << So_res[p] << endl;
						cout << "pb_res: " << pb_res[p] << endl;						
					}
				}
			}

			//OIL
			{
				for (p = 0; p <= (dim - 1); p++)
				{
					lambda1_o[p] = kro[p] / (mi_o[p] * Bo[p]);
					if (isnan(lambda1_o[p])) {
						cout << "lambda1_o: " << lambda1_o[p] << endl;
						cout << "po_res: " << po_res[p] << endl;
						cout << "Sw_res: " << Sw_res[p] << endl;
						cout << "So_res: " << So_res[p] << endl;
						cout << "pb_res: " << pb_res[p] << endl;
					}
				}

				for (p = 0; p <= (dim - 1); p++)
				{
					lambda2_o[p] = -((der_mi_o_p[p] / mi_o[p]) + (der_Bo_p[p] / Bo[p]))*lambda1_o[p];
					if (isnan(lambda2_o[p])) {
						cout << "lambda2_o: " << lambda2_o[p] << endl;
						cout << "po_res: " << po_res[p] << endl;
						cout << "Sw_res: " << Sw_res[p] << endl;
						cout << "So_res: " << So_res[p] << endl;
						cout << "pb_res: " << pb_res[p] << endl;
					}
				}

				for (p = 0; p <= (dim - 1); p++)
				{
					lambda3_o[p] = der_kro_Sw[p] / (mi_o[p] * Bo[p]);
					if (isnan(lambda3_o[p])) {
						cout << "lambda3_o: " << lambda3_o[p] << endl;
						cout << "po_res: " << po_res[p] << endl;
						cout << "Sw_res: " << Sw_res[p] << endl;
						cout << "So_res: " << So_res[p] << endl;
						cout << "pb_res: " << pb_res[p] << endl;
					}
				}

				for (p = 0; p <= (dim - 1); p++)
				{
					lambda4_o[p] = der_kro_So[p] / (mi_o[p] * Bo[p]);
					if (isnan(lambda4_o[p])) {
						cout << "lambda4_o: " << lambda4_o[p] << endl;
						cout << "po_res: " << po_res[p] << endl;
						cout << "Sw_res: " << Sw_res[p] << endl;
						cout << "So_res: " << So_res[p] << endl;
						cout << "pb_res: " << pb_res[p] << endl;
					}
				}
			}

			//GAS
			{
				for (p = 0; p <= (dim - 1); p++)
				{
					lambda1_g[p] = Rs[p] * lambda1_o[p];						
					if (isnan(lambda1_g[p])) {
						cout << "lambda1_g: " << lambda1_g[p] << endl;
						cout << "po_res: " << po_res[p] << endl;
						cout << "Sw_res: " << Sw_res[p] << endl;
						cout << "So_res: " << So_res[p] << endl;
						cout << "pb_res: " << pb_res[p] << endl;
					}
				}

				for (p = 0; p <= (dim - 1); p++)
				{
					lambda2_g[p] = -(Rs[p] * ((der_mi_o_p[p] / mi_o[p]) + (der_Bo_p[p] / Bo[p])) - der_Rs_p[p])* lambda1_o[p];
					if (isnan(lambda2_g[p])) {
						cout << "lambda2_g: " << lambda2_g[p] << endl;
						cout << "po_res: " << po_res[p] << endl;
						cout << "Sw_res: " << Sw_res[p] << endl;
						cout << "So_res: " << So_res[p] << endl;
						cout << "pb_res: " << pb_res[p] << endl;
					}
				}

				for (p = 0; p <= (dim - 1); p++)
				{
					lambda3_g[p] = krg[p] / (mi_g[p] * Bg[p]);					
					if (isnan(lambda3_g[p])) {
						cout << "lambda3_g: " << lambda3_g[p] << endl;
						cout << "po_res: " << po_res[p] << endl;
						cout << "Sw_res: " << Sw_res[p] << endl;
						cout << "So_res: " << So_res[p] << endl;
						cout << "pb_res: " << pb_res[p] << endl;
					}
				}

				for (p = 0; p <= (dim - 1); p++)
				{
					lambda4_g[p] = -((der_mi_g[p] / mi_g[p]) + (der_Bg[p] / Bg[p]))*lambda3_g[p];
					if (isnan(lambda4_g[p])) {
						cout << "lambda4_g: " << lambda4_g[p] << endl;
						cout << "po_res: " << po_res[p] << endl;
						cout << "Sw_res: " << Sw_res[p] << endl;
						cout << "So_res: " << So_res[p] << endl;
						cout << "pb_res: " << pb_res[p] << endl;
					}
				}

				for (p = 0; p <= (dim - 1); p++)
				{
					lambda5_g[p] = Rs[p] * lambda3_o[p];
					if (isnan(lambda5_g[p])) {
						cout << "lambda5_g: " << lambda5_g[p] << endl;
						cout << "po_res: " << po_res[p] << endl;
						cout << "Sw_res: " << Sw_res[p] << endl;
						cout << "So_res: " << So_res[p] << endl;
						cout << "pb_res: " << pb_res[p] << endl;
					}
				}

				for (p = 0; p <= (dim - 1); p++)
				{
					lambda6_g[p] = der_pcg_Swo[p] * lambda3_g[p];
					if (isnan(lambda6_g[p])) {
						cout << "lambda6_g: " << lambda6_g[p] << endl;
						cout << "po_res: " << po_res[p] << endl;
						cout << "Sw_res: " << Sw_res[p] << endl;
						cout << "So_res: " << So_res[p] << endl;
						cout << "pb_res: " << pb_res[p] << endl;
					}
				}

				for (p = 0; p <= (dim - 1); p++)
				{
					lambda7_g[p] = der_krg[p] / (mi_g[p] * Bg[p]);
					if (isnan(lambda7_g[p])) {
						cout << "lambda7_g: " << lambda7_g[p] << endl;
						cout << "po_res: " << po_res[p] << endl;
						cout << "Sw_res: " << Sw_res[p] << endl;
						cout << "So_res: " << So_res[p] << endl;
						cout << "pb_res: " << pb_res[p] << endl;
					}
				}

				for (p = 0; p <= (dim - 1); p++)
				{
					lambda8_g[p] = Rs[p] * lambda4_o[p];
					if (isnan(lambda8_g[p])) {
						cout << "lambda8_g: " << lambda8_g[p] << endl;
						cout << "po_res: " << po_res[p] << endl;
						cout << "Sw_res: " << Sw_res[p] << endl;
						cout << "So_res: " << So_res[p] << endl;
						cout << "pb_res: " << pb_res[p] << endl;
					}
				}
			}					

			//acima pb (óleo)
			for (p = 0; p <= (dim - 1); p++)
			{
				if (acimapb[p])	lambda2_o[p] = -((der_mi_o_p[p] / mi_o[p]))*lambda1_o[p];
				if (isnan(lambda2_o[p])) {
					cout << "lambda2_o2: " << lambda2_o[p] << endl;
					cout << "po_res: " << po_res[p] << endl;
					cout << "Sw_res: " << Sw_res[p] << endl;
					cout << "So_res: " << So_res[p] << endl;
					cout << "pb_res: " << pb_res[p] << endl;
				}
			}

			for (p = 0; p <= (dim - 1); p++)
			{
				if (acimapb[p]) lambda4_o[p] = -((der_mi_o_pb[p] / mi_o[p]) + (der_Bo_pb[p] / Bo[p]))*lambda1_o[p];
				if (isnan(lambda4_o[p])) {
					cout << "lambda4_o2: " << lambda4_o[p] << endl;
					cout << "po_res: " << po_res[p] << endl;
					cout << "Sw_res: " << Sw_res[p] << endl;
					cout << "So_res: " << So_res[p] << endl;
					cout << "pb_res: " << pb_res[p] << endl;
				}
			}
			//
						
			//acima pb (gás)
			for (p = 0; p <= (dim - 1); p++)
			{
				if (acimapb[p]) lambda2_g[p] = -((der_mi_o_p[p] / mi_o[p]))*lambda1_g[p];
				if (isnan(lambda2_g[p])) {
					cout << "lambda2_g2: " << lambda2_g[p] << endl;
					cout << "po_res: " << po_res[p] << endl;
					cout << "Sw_res: " << Sw_res[p] << endl;
					cout << "So_res: " << So_res[p] << endl;
					cout << "pb_res: " << pb_res[p] << endl;
				}
			}
					
			for (p = 0; p <= (dim - 1); p++)
			{
				if (acimapb[p]) lambda8_g[p] = -(Rs[p] * ((der_mi_o_pb[p] / mi_o[p]) + (der_Bo_pb[p] / Bo[p])) - der_Rs_pb[p])* lambda1_o[p];
				if (isnan(lambda8_g[p])) {
					cout << "lambda8_g2: " << lambda8_g[p] << endl;
					cout << "po_res: " << po_res[p] << endl;
					cout << "Sw_res: " << Sw_res[p] << endl;
					cout << "So_res: " << So_res[p] << endl;
					cout << "pb_res: " << pb_res[p] << endl;
				}
			}
			//

			//CALC COEFS
			p = 0;
			for (k = 0; k <= (o - 1); k++) {
				for (j = 0; j <= (m - 1); j++) {
					for (i = 0; i <= (n - 1); i++) {
					
						y = 3 * p;
						controladorP = 0.0;
						controladorSw = 0.0;
						controladorSo = 0.0;
						B[y] = 0.0;					

						//WATER

								//North
								if (k != 0) {

									calc_water_pos_res(pw_res[p], pw_res[p - n*m],
										B, A,
										&controladorP, &controladorSw,
										y, unW[p],
										1.0l, 0.0l, Azinterno[p] / (0.5l*(deltaz[k] + deltaz[k - 1])), C1,
										Kz[p]*lambda1_w[p], Kz[p]*lambda2_w[p], der_pcw[p], Kz[p] * lambda4_w[p], Kz[p] * lambda5_w[p],
										Kz[p - n*m] * lambda1_w[p - n*m], Kz[p - n*m] * lambda2_w[p - n*m], der_pcw[p - n*m], Kz[p - n*m] * lambda4_w[p - n*m], Kz[p-n*m] * lambda5_w[p - n*m]);
									

									calc_water_grav_res(rho_w[p - n*m], rho_w[p], der_rho_w[p - n*m], der_rho_w[p],
										B, A,
										&controladorP, &controladorSw,
										y, unW[p],
										1.0l, 0.0l, Azinterno[p] / (0.5l*(deltaz[k] + deltaz[k - 1])), C2,
										Kz[p] * lambda1_w[p], Kz[p] * lambda2_w[p], Kz[p] * lambda5_w[p],
										Kz[p - n*m] * lambda1_w[p - n*m], Kz[p - n*m] *lambda2_w[p - n*m], Kz[p - n*m] * lambda5_w[p - n*m],
										(deltaz[k] + deltaz[k - 1]) / 2.0l, g);

								}

								//Foward
								if (j != 0) {

									calc_water_pos_res(pw_res[p], pw_res[p - n],
										B, A,
										&controladorP, &controladorSw,
										y, ufW[p],
										1.0l, 0.0l, Ax[p], C1,
										Kx[p] * lambda1_w[p], Kx[p] * lambda2_w[p], der_pcw[p], Kx[p] * lambda4_w[p], Kx[p] * lambda5_w[p],
										Kx[p - n] * lambda1_w[p - n], Kx[p - n] * lambda2_w[p - n], der_pcw[p - n], Kx[p - n] * lambda4_w[p - n], Kx[p-n] * lambda5_w[p - n]);

								}

								//West(boundary)
								if (n != 1 && i == (n - 1)) {

									calc_water_pos_res(pw_res[p], pw_res[p - (n - 1)],
										B, A,
										&controladorP, &controladorSw,
										y, uwbW[p],
										0.5, 0.5, Aphi[p], C1,
										lambda1_w[p], lambda2_w[p], der_pcw[p], lambda4_w[p], lambda5_w[p],
										lambda1_w[p - (n - 1)], lambda2_w[p - (n - 1)], der_pcw[p - (n - 1)], lambda4_w[p - (n - 1)], lambda5_w[p - (n - 1)]);
								
								}

								//East
								if (n != 1 && i != (0)) {

									calc_water_pos_res(pw_res[p], pw_res[p - 1],
										B, A,
										&controladorP, &controladorSw,
										y, ueW[p],
										0.5, 0.5, Aphi[p], C1,
										lambda1_w[p], lambda2_w[p], der_pcw[p], lambda4_w[p], lambda5_w[p],
										lambda1_w[p - 1], lambda2_w[p - 1], der_pcw[p - 1], lambda4_w[p - 1], lambda5_w[p - 1]);
								
								}

								//Center
								
									B[y] = B[y] - C6*At[p] * ((phi[p] * Sw_res[p] / Bw[p]) - (calc_phi_p(phi_ref, Cr, po_res_0[p], pi)*Sw_res_0[p] / calc_water_Bw_p(po_res_0[p], T, parametrosBw)));
									if (j == (m - 1)) {
										B[y] = B[y] + C6*At[p] * deltat*Qwsinj / Bw[p];
									}
									if (j == 0 && zonaprodutora[k]) {
										B[y] = B[y] + Axint[p] * C1*(Kx[p]*lambda1_w[p] * (PWB_0[k] - pw_res[p]));
										Qw_well[k] = -ft3s_bbld(Axint[p] * C1*(Kx[p] * lambda1_w[p] * (PWB_0[k] - pw_res[p])) / (deltat));										
									}
																		
									centerP = ucW[p];
									
									A[ucW[p]] = -C6*At[p] * ((Sw_res[p] / Bw[p])*(der_phi[p] - (phi[p] * der_Bw[p] / Bw[p])));
									
									if (j == (m - 1)) {
										A[ucW[p]] = A[ucW[p]] - C6*At[p] * deltat*Qwsinj*(der_Bw[p] / (Bw[p] * Bw[p]));
									}
									if (j == 0 && zonaprodutora[k]) {
										A[ucW[p]] = A[ucW[p]] - Axint[p] * C1*(Kx[p] * lambda1_w[p] - Kx[p] * lambda2_w[p] * (PWB_0[k] - pw_res[p]));
									}
									
									centerSw = ucW[p] + 1;

									A[ucW[p] + 1] = -At[p] * C6* (phi[p] / Bw[p]);
									if (j == 0 && zonaprodutora[k]) {
										A[ucW[p] + 1] = A[ucW[p] + 1] + Axint[p] * C1*(Kx[p] * lambda5_w[p] * (PWB_0[k] - pw_res[p]));
									}
																									

								//West
								if (n != 1 && i != (n - 1)) {

									calc_water_pos_res(pw_res[p], pw_res[p + 1],
										B, A,
										&controladorP, &controladorSw,
										y, uwW[p],
										0.5, 0.5, Aphi[p], C1,
										lambda1_w[p], lambda2_w[p], der_pcw[p], lambda4_w[p], lambda5_w[p],
										lambda1_w[p + 1], lambda2_w[p + 1], der_pcw[p + 1], lambda4_w[p + 1], lambda5_w[p + 1]);
									
								}

								//East(boundary)
								if (n != 1 && i == 0) {

									calc_water_pos_res(pw_res[p], pw_res[p + (n - 1)],
										B, A,
										&controladorP, &controladorSw,
										y, uebW[p],
										0.5, 0.5, Aphi[p], C1,
										lambda1_w[p], lambda2_w[p], der_pcw[p], lambda4_w[p], lambda5_w[p],
										lambda1_w[p + (n - 1)], lambda2_w[p + (n - 1)], der_pcw[p + (n - 1)], lambda4_w[p + (n - 1)], lambda5_w[p + (n - 1)]);
								
								}

								//Backward
								if (j != (m - 1)) {

									calc_water_pos_res(pw_res[p], pw_res[p + n],
										B, A,
										&controladorP, &controladorSw,
										y, ubW[p],
										0.0l, 1.0l, Ax[p], C1,
										Kx[p] * lambda1_w[p], Kx[p] * lambda2_w[p], der_pcw[p], Kx[p] * lambda4_w[p], Kx[p] * lambda5_w[p],
										Kx[p + n] * lambda1_w[p + n], Kx[p + n] * lambda2_w[p + n], der_pcw[p + n], Kx[p + n] * lambda4_w[p + n], Kx[p+n] * lambda5_w[p + n]);
								
								}

								//South
								if (k != (o - 1)) {

									calc_water_pos_res(pw_res[p], pw_res[p + n*m],
										B, A,
										&controladorP, &controladorSw,
										y, usW[p],
										0.0l, 1.0l, Azinterno[p] / (0.5l*(deltaz[k] + deltaz[k + 1])), C1,
										Kz[p] * lambda1_w[p], Kz[p]*lambda2_w[p], der_pcw[p], Kz[p]*lambda4_w[p], Kz[p] * lambda5_w[p],
										Kz[p + n*m] * lambda1_w[p + n*m], Kz[p + n*m] * lambda2_w[p + n*m], der_pcw[p + n*m], Kz[p + n*m] * lambda4_w[p + n*m], Kz[p+n*m] * lambda5_w[p + n*m]);

									calc_water_grav_res(rho_w[p + n*m], rho_w[p], der_rho_w[p + n*m], der_rho_w[p],
										B, A,
										&controladorP, &controladorSw,
										y, usW[p],
										0.0l, 1.0l, Azinterno[p] / (0.5l*(deltaz[k] + deltaz[k + 1])), C2,
										Kz[p] * lambda1_w[p], Kz[p] * lambda2_w[p], Kz[p] * lambda5_w[p],
										Kz[p + n*m] * lambda1_w[p + n*m], Kz[p + n*m] * lambda2_w[p + n*m], Kz[p + n*m] * lambda5_w[p + n*m],
										(deltaz[k] + deltaz[k + 1]) / 2.0l, -g);

								}

								A[centerP] = A[centerP] + controladorP;
								A[centerSw] = A[centerSw] + controladorSw;

						
						controladorP = 0.0;
						controladorSw = 0.0;
						controladorSo = 0.0;
						y = y + 1;
						B[y] = 0.0;
						

						//OIL
													
								//North
								if (k != 0) {

									calc_oil_pos_res(po_res[p], po_res[p - n*m],
										B, A,
										&controladorP, &controladorSw, &controladorSo,
										y, unO[p],
										1.0l, 0.0l, Azinterno[p] / (0.5l*(deltaz[k] + deltaz[k - 1])), C1,
										Kz[p] * lambda1_o[p], Kz[p] * lambda2_o[p], Kz[p] * lambda3_o[p], Kz[p] * lambda4_o[p],
										Kz[p - n*m] * lambda1_o[p - n*m], Kz[p - n*m] * lambda2_o[p - n*m], Kz[p - n*m] * lambda3_o[p - n*m], Kz[p - n*m] * lambda4_o[p - n*m]);

									calc_oil_grav_res(rho_o[p], rho_o[p - n*m], der_rho_o_p[p], der_rho_o_p[p - n*m], der_rho_o_pb[p], der_rho_o_pb[p - n*m],
										B, A,
										&controladorP, &controladorSw, &controladorSo,
										y, unO[p],
										1.0l, 0.0l, Azinterno[p] / (0.5l*(deltaz[k] + deltaz[k - 1])), C2,
										Kz[p] * lambda1_o[p], Kz[p] * lambda2_o[p], Kz[p] * lambda3_o[p], Kz[p] * lambda4_o[p],
										Kz[p - n*m] * lambda1_o[p - n*m], Kz[p - n*m] * lambda2_o[p - n*m], Kz[p - n*m] * lambda3_o[p - n*m], Kz[p - n*m] * lambda4_o[p - n*m],
										(deltaz[k]+ deltaz[k-1])/2.0l, g);

								}

								//Forward
								if (j != 0) {

									calc_oil_pos_res(po_res[p], po_res[p - n],
										B, A,
										&controladorP, &controladorSw, &controladorSo,
										y, ufO[p],
										1.0l, 0.0l, Ax[p], C1,
										Kx[p] * lambda1_o[p], Kx[p] * lambda2_o[p], Kx[p] * lambda3_o[p], Kx[p] * lambda4_o[p],
										Kx[p - n] * lambda1_o[p - n], Kx[p - n] * lambda2_o[p - n], Kx[p - n] * lambda3_o[p - n], Kx[p-n] * lambda4_o[p - n]);
								
								}

								//West(boundary)
								if (n != 1 && i == (n - 1)) {

									calc_oil_pos_res(po_res[p], po_res[p - (n - 1)],
										B, A,
										&controladorP, &controladorSw, &controladorSo,
										y, uwbO[p],
										0.5l, 0.5l, Aphi[p], C1,
										lambda1_o[p], lambda2_o[p], lambda3_o[p], lambda4_o[p],
										lambda1_o[p - (n - 1)], lambda2_o[p - (n - 1)], lambda3_o[p - (n - 1)], lambda4_o[p - (n - 1)]);
								
								}

								//East
								if (n != 1 && i != (0)) {

									calc_oil_pos_res(po_res[p], po_res[p - 1],
										B, A,
										&controladorP, &controladorSw, &controladorSo,
										y, ueO[p],
										0.5l, 0.5l, Aphi[p], C1,
										lambda1_o[p], lambda2_o[p], lambda3_o[p], lambda4_o[p],
										lambda1_o[p - 1], lambda2_o[p - 1], lambda3_o[p - 1], lambda4_o[p - 1]);
								
								}
								
								//Center
								
									B[y] = B[y] - C6*At[p] * phi[p] * (So_res[p] / Bo[p]) + C6*At[p] * calc_phi_p(phi_ref, Cr, po_res_0[p], pi)*(So_res_0[p] / calc_oil_Bo_p(yo, yg, po_res_0[p], pb_res_0[p], API, T));

									if (j == (m - 1)) {
										B[y] = B[y] + C6*At[p] * deltat*Qosinj / Bo[p];
									}
									if (j == 0 && zonaprodutora[k]) {
										B[y] = B[y] + Axint[p] * C1*(Kx[p] * lambda1_o[p] * (PWB_0[k] - po_res[p]));
										Qo_well[k] = -ft3s_bbld(Axint[p] * C1*(Kx[p] * lambda1_o[p] * (PWB_0[k] - po_res[p])) / (deltat));
									}
																		
									centerP = ucO[p];

									A[ucO[p]] = -C6*At[p] * ((So_res[p] / Bo[p])*(der_phi[p] - (phi[p] * der_Bo_p[p] / Bo[p])));
									
									if (j == (m - 1)) {
										A[ucO[p]] = A[ucO[p]] - C6*At[p] * deltat*Qosinj*(der_Bo_p[p] / (Bo[p] * Bo[p]));
									}
									if (j == 0 && zonaprodutora[k]) {
										A[ucO[p]] = A[ucO[p]] - Axint[p] * C1*Kx[p] * (lambda1_o[p] - lambda2_o[p] * (PWB_0[k] - po_res[p]));
									}

									centerSw = ucO[p] + 1;

									if (acimapb[p]==0) {
										A[ucO[p] + 1] = 0.0l;
									}
									if (acimapb[p]) {
										A[ucO[p] + 1] = At[p] * C6* ((phi[p] / Bo[p]));
									}

									if (j == 0 && zonaprodutora[k]) {
										A[ucO[p] + 1] = A[ucO[p] + 1] + Axint[p] * C1*(Kx[p] * lambda3_o[p] * (PWB_0[k] - po_res[p]));
									}

									centerSo = ucO[p] + 2;

									if (acimapb[p] == 0) {
										A[ucO[p] + 2] = -At[p] * C6* ((phi[p] / Bo[p]));
									}
									if (acimapb[p]) {
										A[ucO[p] + 2] = At[p] * C6* (So_res[p] / Bo[p])*((phi[p] * der_Bo_pb[p] / Bo[p]));
									}

									if (j == 0 && zonaprodutora[k]) {
										A[ucO[p] + 2] = A[ucO[p] + 2] + Axint[p] * C1*(Kx[p] * lambda4_o[p] * (PWB_0[k] - po_res[p]));
									}
								

								//West
								if (n != 1 && i != (n - 1)) {

									calc_oil_pos_res(po_res[p], po_res[p + 1],
										B, A,
										&controladorP, &controladorSw, &controladorSo,
										y, uwO[p],
										0.5l, 0.5l, Aphi[p], C1,
										lambda1_o[p], lambda2_o[p], lambda3_o[p], lambda4_o[p],
										lambda1_o[p + 1], lambda2_o[p + 1], lambda3_o[p + 1], lambda4_o[p + 1]);

								}

								//East(boundary)
								if (n != 1 && i == 0) {

									calc_oil_pos_res(po_res[p], po_res[p + (n - 1)],
										B, A,
										&controladorP, &controladorSw, &controladorSo,
										y, uebO[p],
										0.5l, 0.5l, Aphi[p], C1,
										lambda1_o[p], lambda2_o[p], lambda3_o[p], lambda4_o[p],
										lambda1_o[p + (n - 1)], lambda2_o[p + (n - 1)], lambda3_o[p + (n - 1)], lambda4_o[p + (n - 1)]);
							
								}

								//Backward
								if (j != (m - 1)) {

									calc_oil_pos_res(po_res[p], po_res[p + n],
										B, A,
										&controladorP, &controladorSw, &controladorSo,
										y, ubO[p],
										0.0l, 1.0l, Ax[p], C1,
										Kx[p] * lambda1_o[p], Kx[p] * lambda2_o[p], Kx[p] * lambda3_o[p], Kx[p] * lambda4_o[p],
										Kx[p + n] * lambda1_o[p + n], Kx[p + n] * lambda2_o[p + n], Kx[p + n] * lambda3_o[p + n], Kx[p+n] * lambda4_o[p + n]);

								}

								//South
								if (k != (o - 1)) {
									
									calc_oil_pos_res(po_res[p], po_res[p + n*m],
										B, A,
										&controladorP, &controladorSw, &controladorSo,
										y, usO[p],
										0.0l, 1.0l, Azinterno[p] / (0.5l*(deltaz[k] + deltaz[k + 1])), C1,
										Kz[p] * lambda1_o[p], Kz[p] * lambda2_o[p], Kz[p] * lambda3_o[p], Kz[p] * lambda4_o[p],
										Kz[p+n*m] * lambda1_o[p + n*m], Kz[p + n*m]*lambda2_o[p + n*m], Kz[p + n*m] * lambda3_o[p + n*m], Kz[p + n*m] * lambda4_o[p + n*m]);
									

									calc_oil_grav_res(rho_o[p], rho_o[p + n*m], der_rho_o_p[p], der_rho_o_p[p + n*m], der_rho_o_pb[p], der_rho_o_pb[p + n*m],
										B, A,
										&controladorP, &controladorSw, &controladorSo,
										y, usO[p],
										0.0l, 1.0l, Azinterno[p] / (0.5l*(deltaz[k] + deltaz[k + 1])), C2,
										Kz[p] * lambda1_o[p], Kz[p] * lambda2_o[p], Kz[p] * lambda3_o[p], Kz[p] * lambda4_o[p],
										Kz[p + n*m] * lambda1_o[p + n*m], Kz[p + n*m] * lambda2_o[p + n*m], Kz[p + n*m] * lambda3_o[p + n*m], Kz[p + n*m] *lambda4_o[p + n*m],
										(deltaz[k] + deltaz[k + 1]) / 2.0l, -g);
										
								}

								A[centerP] = A[centerP] + controladorP;
								A[centerSw] = A[centerSw] + controladorSw;
								A[centerSo] = A[centerSo] + controladorSo;

					

						controladorP = 0.0;
						controladorSw = 0.0;
						controladorSo = 0.0;
						y = y + 1;
						B[y] = 0.0;

						
						//GAS
						
							
								//North
								if (k != 0) {
									
									calc_gas_pos_res(po_res[p], po_res[p - n*m], pg_res[p], pg_res[p - n*m],
										B, A,
										&controladorP, &controladorSw, &controladorSo,
										y, unG[p],
										1.0l, 0.0l, Azinterno[p] / (0.5l*(deltaz[k] + deltaz[k - 1])), C3,
										Kz[p]*lambda1_g[p], Kz[p] * lambda2_g[p], Kz[p] * lambda3_g[p], Kz[p] * lambda4_g[p], Kz[p] * lambda5_g[p], Kz[p] * lambda6_g[p], Kz[p] * lambda7_g[p], Kz[p] * lambda8_g[p], der_pcg_Swo[p],
										Kz[p - n*m] * lambda1_g[p - n*m], Kz[p - n*m] * lambda2_g[p - n*m], Kz[p - n*m] * lambda3_g[p - n*m], Kz[p - n*m] * lambda4_g[p - n*m], Kz[p - n*m] * lambda5_g[p - n*m], Kz[p - n*m] * lambda6_g[p - n*m], Kz[p - n*m] * lambda7_g[p - n*m], Kz[p-n*m] * lambda8_g[p - n*m], der_pcg_Swo[p - n*m],
										acimapb[p], acimapb[p-n*m]);

									calc_gas_grav_res(rho_o[p], der_rho_o_p[p], rho_o[p - n*m], der_rho_o_p[p - n*m], rho_g[p], der_rho_g[p], rho_g[p - n*m], der_rho_g[p - n*m],
										der_rho_o_pb[p], der_rho_o_pb[p - n*m],
										B, A,
										&controladorP, &controladorSw, &controladorSo,
										y, unG[p],
										1.0l, 0.0l, Azinterno[p] / (0.5l*(deltaz[k] + deltaz[k - 1])), C4,
										Kz[p] * lambda1_g[p], Kz[p] * lambda2_g[p], Kz[p] * lambda3_g[p], Kz[p] * lambda4_g[p], Kz[p] * lambda5_g[p], Kz[p] * lambda7_g[p], Kz[p] * lambda8_g[p],
										Kz[p - n*m] * lambda1_g[p - n*m], Kz[p - n*m] * lambda2_g[p - n*m], Kz[p - n*m] * lambda3_g[p - n*m], Kz[p - n*m] * lambda4_g[p - n*m], Kz[p - n*m] * lambda5_g[p - n*m], Kz[p - n*m] * lambda7_g[p - n*m], Kz[p - n*m] * lambda8_g[p - n*m],
										(deltaz[k] + deltaz[k - 1]) / 2.0l, g, acimapb[p], acimapb[p - n*m]);
																	
								}

								//Forward
								if (j != 0) {

									calc_gas_pos_res(po_res[p], po_res[p - n], pg_res[p], pg_res[p - n],
										B, A,
										&controladorP, &controladorSw, &controladorSo,
										y, ufG[p],
										1.0l, 0.0l, Ax[p], C3,
										Kx[p]*lambda1_g[p], Kx[p] * lambda2_g[p], Kx[p] * lambda3_g[p], Kx[p] * lambda4_g[p], Kx[p] * lambda5_g[p], Kx[p] * lambda6_g[p], Kx[p] * lambda7_g[p], Kx[p] * lambda8_g[p], der_pcg_Swo[p],
										Kx[p-n] * lambda1_g[p - n], Kx[p - n] * lambda2_g[p - n], Kx[p - n] * lambda3_g[p - n], Kx[p - n] * lambda4_g[p - n], Kx[p - n] * lambda5_g[p - n], Kx[p - n] * lambda6_g[p - n], Kx[p - n] * lambda7_g[p - n], Kx[p - n] * lambda8_g[p - n], der_pcg_Swo[p - n],
										acimapb[p], acimapb[p - n]);

								}

								//West(boundary)
								if (n != 1 && i == (n - 1)) {

									calc_gas_pos_res(po_res[p], po_res[p - (n - 1)], pg_res[p], pg_res[p - (n - 1)],
										B, A,
										&controladorP, &controladorSw, &controladorSo,
										y, uwbG[p],
										0.5l, 0.5l, Aphi[p], C3,
										lambda1_g[p], lambda2_g[p], lambda3_g[p], lambda4_g[p], lambda5_g[p], lambda6_g[p], lambda7_g[p], lambda8_g[p], der_pcg_Swo[p],
										lambda1_g[p - (n - 1)], lambda2_g[p - (n - 1)], lambda3_g[p - (n - 1)], lambda4_g[p - (n - 1)], lambda5_g[p - (n - 1)], lambda6_g[p - (n - 1)], lambda7_g[p - (n - 1)], lambda8_g[p - (n - 1)], der_pcg_Swo[p - (n - 1)],
										acimapb[p], acimapb[p - (n - 1)]);

								}

								//East
								if (n != 1 && i != (0)) {

									calc_gas_pos_res(po_res[p], po_res[p - 1], pg_res[p], pg_res[p - 1],
										B, A,
										&controladorP, &controladorSw, &controladorSo,
										y, ueG[p],
										0.5l, 0.5l, Aphi[p], C3,
										lambda1_g[p], lambda2_g[p], lambda3_g[p], lambda4_g[p], lambda5_g[p], lambda6_g[p], lambda7_g[p], lambda8_g[p], der_pcg_Swo[p],
										lambda1_g[p - 1], lambda2_g[p - 1], lambda3_g[p - 1], lambda4_g[p - 1], lambda5_g[p - 1], lambda6_g[p - 1], lambda7_g[p - 1], lambda8_g[p - 1], der_pcg_Swo[p - 1],
										acimapb[p], acimapb[p - 1]);

								}

								//Center
								
									B[y] = B[y] - At[p] * C5*phi[p] * ((So_res[p] * Rs[p] / Bo[p]) + (Sg_res[p] / Bg[p]))
										+ At[p] * C5*calc_phi_p(phi_ref, Cr, po_res_0[p], pi)*
										((So_res_0[p] * calc_gas_Rs_p(yg, po_res_0[p], pb_res_0[p], API, T) / calc_oil_Bo_p(yo, yg, po_res_0[p], pb_res_0[p], API, T)) + ((1 - Sw_res_0[p] - So_res_0[p]) / calc_gas_Bg_p(po_res_0[p])));
									
									if (j == (m - 1)){
										B[y] = B[y] + At[p] * C5*deltat*((Qosinj*Rs[p] / Bo[p]) + (Qgsinj / Bg[p]));
									}									
									if (j == 0 && zonaprodutora[k]) {
										B[y] = B[y] + Axint[p] * C3*Kx[p] * (lambda1_g[p] * (PWB_0[k] - po_res[p]) + lambda3_g[p] * (PWB_0[k] - pg_res[p]));
										Qg_well[k] = -bbl_ft3(ft3s_bbld(Axint[p] * C3*Kx[p] * (lambda1_g[p] * (PWB_0[k] - po_res[p]) + lambda3_g[p] * (PWB_0[k] - pg_res[p])) / (deltat)));
									}

									centerP = ucG[p];									
									
									if (acimapb[p] == 0) {

										A[ucG[p]] = -At[p] * C5*((Sg_res[p] / Bg[p])*(der_phi[p] - phi[p] * der_Bg[p] / Bg[p]))
											- At[p] * C5*(So_res[p] * Rs[p] / Bo[p])*(der_phi[p] - phi[p] * der_Bo_p[p] / Bo[p])
											- At[p] * C5*So_res[p] * phi[p] * der_Rs_p[p] / Bo[p];

									}

									if (acimapb[p]) {
									
										A[ucG[p]] = -At[p] * C5*(So_res[p] * Rs[p] / Bo[p])*(der_phi[p]);
									
									}

									if (j == (m - 1)){
										
										A[ucG[p]] = A[ucG[p]] + At[p] * C5*deltat*(Qosinj*((der_Rs_p[p] / Bo[p]) - Rs[p] * der_Bo_p[p] / (Bo[p] * Bo[p])) - Qgsinj*(der_Bg[p] / (Bg[p] * Bg[p])));
									
									}
									
									if (j == 0 && zonaprodutora[k]) {
									
										A[ucG[p]] = A[ucG[p]] - Axint[p] * C3*Kx[p] * (lambda1_g[p] - lambda2_g[p] * (PWB_0[k] - po_res[p]) + lambda3_g[p] - lambda4_g[p] * (PWB_0[k] - pg_res[p]));
									
									}

									centerSw = ucG[p] + 1;
									
									if (acimapb[p] == 0) {

										A[ucG[p] + 1] = At[p] * C5*(phi[p] / Bg[p]);

									}
									
									if (acimapb[p]) {

										A[ucG[p] + 1] = At[p] * C5*phi[p] * Rs[p] / Bo[p];

									}
									
									if (j == 0 && zonaprodutora[k]) {

										A[ucG[p] + 1] = A[ucG[p] + 1] + Axint[p] * C3*Kx[p] * (lambda5_g[p] * (PWB_0[k] - po_res[p]) + lambda7_g[p] * (PWB_0[k] - pg_res[p]) - lambda6_g[p]);
									
									}

									
									u = u + 1;
									centerSo = ucG[p] + 2;

									if (acimapb[p]==0) {
									
										A[ucG[p] + 2] = A[ucG[p] + 1] - At[p] * C5*phi[p] * Rs[p] / Bo[p] + At[p] * C5*(phi[p] / Bg[p]);
									
									}

									if (acimapb[p]){
									
										A[ucG[p] + 2] = At[p] * C5* (So_res[p] * phi[p] / Bo[p])*((Rs[p]  * der_Bo_pb[p] / Bo[p]) - der_Rs_pb[p]);
									
									}

									if (j == 0 && zonaprodutora[k]){
 
										A[ucG[p] + 2] = A[ucG[p] + 2] + Axint[p] * C3*Kx[p] * (lambda8_g[p] * (PWB_0[k] - po_res[p]) + lambda7_g[p] * (PWB_0[k] - pg_res[p]) - lambda6_g[p]);
									
									}
									
									u = u + 1;

								

								//West
								if (n != 1 && i != (n - 1)) {

									calc_gas_pos_res(po_res[p], po_res[p + 1], pg_res[p], pg_res[p + 1],
										B, A,
										&controladorP, &controladorSw, &controladorSo,
										y, uwG[p],
										0.5l, 0.5l, Aphi[p], C3,
										lambda1_g[p], lambda2_g[p], lambda3_g[p], lambda4_g[p], lambda5_g[p], lambda6_g[p], lambda7_g[p], lambda8_g[p], der_pcg_Swo[p],
										lambda1_g[p + 1], lambda2_g[p + 1], lambda3_g[p + 1], lambda4_g[p + 1], lambda5_g[p + 1], lambda6_g[p + 1], lambda7_g[p + 1], lambda8_g[p + 1], der_pcg_Swo[p + 1],
										acimapb[p], acimapb[p + 1]);

								}

								//East(boundary)
								if (n != 1 && i == 0) {

									calc_gas_pos_res(po_res[p], po_res[p + (n - 1)], pg_res[p], pg_res[p + (n - 1)],
										B, A,
										&controladorP, &controladorSw, &controladorSo,
										y, uebG[p],
										0.5l, 0.5l, Aphi[p], C3,
										lambda1_g[p], lambda2_g[p], lambda3_g[p], lambda4_g[p], lambda5_g[p], lambda6_g[p], lambda7_g[p], lambda8_g[p], der_pcg_Swo[p],
										lambda1_g[p + (n - 1)], lambda2_g[p + (n - 1)], lambda3_g[p + (n - 1)], lambda4_g[p + (n - 1)], lambda5_g[p + (n - 1)], lambda6_g[p + (n - 1)], lambda7_g[p + (n - 1)], lambda8_g[p + (n - 1)], der_pcg_Swo[p + (n - 1)],
										acimapb[p], acimapb[p + (n - 1)]);

								}

								//Backward
								if (j != (m - 1)) {

									calc_gas_pos_res(po_res[p], po_res[p + n], pg_res[p], pg_res[p + n],
										B, A,
										&controladorP, &controladorSw, &controladorSo,
										y, ubG[p],
										0.0l, 1.0l, Ax[p], C3,
										Kx[p] * lambda1_g[p], Kx[p]*lambda2_g[p], Kx[p] * lambda3_g[p], Kx[p] * lambda4_g[p], Kx[p] * lambda5_g[p], Kx[p] * lambda6_g[p], Kx[p] * lambda7_g[p], Kx[p] * lambda8_g[p], der_pcg_Swo[p],
										Kx[p + n] * lambda1_g[p + n], Kx[p + n] * lambda2_g[p + n], Kx[p + n] * lambda3_g[p + n], Kx[p + n] * lambda4_g[p + n], Kx[p + n] * lambda5_g[p + n], Kx[p + n] * lambda6_g[p + n], Kx[p + n] * lambda7_g[p + n], Kx[p+n] * lambda8_g[p + n], der_pcg_Swo[p + n],
										acimapb[p], acimapb[p + n]);

								}
								
								//South
								if (k != (o - 1)) {
									
									calc_gas_pos_res(po_res[p], po_res[p + n*m], pg_res[p], pg_res[p + n*m],
										B, A,
										&controladorP, &controladorSw, &controladorSo,
										y, usG[p],
										0.0l, 1.0l, Azinterno[p]/ (0.5l*(deltaz[k] + deltaz[k + 1])), C3,
										Kz[p] * lambda1_g[p], Kz[p] * lambda2_g[p], Kz[p] * lambda3_g[p], Kz[p] * lambda4_g[p], Kz[p] * lambda5_g[p], Kz[p] * lambda6_g[p], Kz[p] * lambda7_g[p], Kz[p] * lambda8_g[p], der_pcg_Swo[p],
										Kz[p + n*m] * lambda1_g[p + n*m], Kz[p + n*m] * lambda2_g[p + n*m], Kz[p + n*m] * lambda3_g[p + n*m], Kz[p + n*m] * lambda4_g[p + n*m], Kz[p + n*m] * lambda5_g[p + n*m], Kz[p + n*m] * lambda6_g[p + n*m], Kz[p + n*m] * lambda7_g[p + n*m], Kz[p + n*m] * lambda8_g[p + n*m], der_pcg_Swo[p + n*m],
										acimapb[p], acimapb[p + n*m]);

									calc_gas_grav_res(rho_o[p], der_rho_o_p[p], rho_o[p + n*m], der_rho_o_p[p + n*m], rho_g[p], der_rho_g[p], rho_g[p + n*m], der_rho_g[p + n*m],
										der_rho_o_pb[p], der_rho_o_pb[p + n*m],
										B, A,
										&controladorP, &controladorSw, &controladorSo,
										y, usG[p],
										0.0l, 1.0l, Azinterno[p] / (0.5l*(deltaz[k] + deltaz[k + 1])), C4,
										Kz[p] * lambda1_g[p], Kz[p] * lambda2_g[p], Kz[p] * lambda3_g[p], Kz[p] * lambda4_g[p], Kz[p] * lambda5_g[p], Kz[p] * lambda7_g[p], Kz[p] * lambda8_g[p],
										Kz[p + n*m] * lambda1_g[p + n*m], Kz[p + n*m] * lambda2_g[p + n*m], Kz[p + n*m] * lambda3_g[p + n*m], Kz[p + n*m] * lambda4_g[p + n*m], Kz[p + n*m] * lambda5_g[p + n*m], Kz[p + n*m] * lambda7_g[p + n*m], Kz[p + n*m] * lambda8_g[p + n*m],
										(deltaz[k] + deltaz[k + 1]) / 2.0l, -g, acimapb[p], acimapb[p + n*m]);
									
								}


								A[centerP] = A[centerP] + controladorP;
								A[centerSw] = A[centerSw] + controladorSw;
								A[centerSo] = A[centerSo] + controladorSo;						
								


						p++;
					}
				}

			}
						
			y = y + 1;		
			for (int wer = 0; wer < y; wer++)
			{
				B[wer] = -B[wer];
			}
			//
			
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
			//

									
			cout.precision(19);
			p = 0;
			for (k = 0; k <= (o - 1); k++) {
				for (j = 0; j <= (m - 1); j++) {
					for (i = 0; i <= (n - 1); i++) {

						teste = 0;

						if (acimapb[p] == 0) {

							po_res[p] = po_res[p] + X[3 * p];
							Sw_res[p] = Sw_res[p] + X[3 * p + 1];
							So_res[p] = So_res[p] + X[3 * p + 2];
							if (So_res[p] < 0.0000000001l){
								So_res[p] = 0.0000000001l;								
							}
							pb_res[p] = po_res[p];

						}

						if (acimapb[p]) {

							po_res[p] = po_res[p] + X[3 * p];
							Sw_res[p] = Sw_res[p] + X[3 * p + 1];
							pb_res[p] = pb_res[p] + X[3 * p + 2];
							So_res[p] = 1.0l - Sw_res[p];

						}

						pw_res[p] = po_res[p];// +calc_water_pcw_Sw(Sw_res[p], Swc, epsilon, pcwmin, pcwmax, Bpcw);
						pg_res[p] = po_res[p];// +calc_gas_pcg(Sw_res[p], So_res[p], Swc, Sor, epsilon, pcgmin, pcgmax, Bpcg);
						Sg_res[p] = (1.0l - So_res[p] - Sw_res[p]);

												
						if (acimapb[p]==0) {
							if (fabs(X[3 * p] / po_res[p]) > 0.0002 || fabs(X[3 * p + 1] / Sw_res[p]) > 0.0002 || fabs(X[3 * p + 2] / So_res[p]) > 0.0002 || teste) {
								dontstopmenow = 1;								
							}
						}

						if (acimapb[p]) {
							if (fabs(X[3 * p]/po_res[p]) > 0.0002 || fabs(X[3 * p + 1] / Sw_res[p]) > 0.0002 || fabs(X[3 * p + 2] / pb_res[p]) > 0.0002 || teste) {
								dontstopmenow = 1;						
							}
						}					
						

						if ((po_res[p]<pb_res[p]) && acimapb[p]) {
							acimapb[p] = 0;
							pb_res[p] = po_res[p];
							So_res[p] = So_res[p] - 0.001l;//0.01l
							Sg_res[p] = 0.001l;							
						}
						if (Sg_res[p] <= -0.00000000000001l && acimapb[p] == 0) {
							acimapb[p] = 1;
							pb_res[p] = po_res[p] - 2.0l;//0.1
							So_res[p] = 1 - Sw_res[p];							
						}

						
						if (isnan(po_res[p]) || isnan(So_res[p]) || isnan(Sw_res[p]) || isnan(pb_res[p]) || 
							po_res[p] < 0.000011l || pb_res[p] < 500l || Sw_res[p] > 1.0l || So_res[p] > 1.00000000000000001l) {
							//|| po_res[p]>(po_res_0[p] + 10.l) || fabs(pb_res[p] - pb_res_0[p])>50.0l
							//||po_res[p]>(po_res_0[p] + 10.l)||
							
							if (level > 250) {
							
								grava_p(po_res, dim, 0, (t + 1), ite_coup, coup_ts);
								grava_pb(pb_res, dim, 0, (t + 1), ite_coup, coup_ts);
								grava_Sw(Sw_res, dim, 0, (t + 1), ite_coup, coup_ts);
								grava_So(So_res, dim, 0, (t + 1), ite_coup, coup_ts);
							
							}

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

		//cout << "3" << endl;
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

		//cout << "4" << endl;
		p = 0;
		if (nan == 1 || l>5) {

			double factor = 0;
			if (level==0) factor = pow(0.5l, 6.0 - 2.0*level);
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
						pw_res[p] = po_res[p];//+ calc_water_pcw_Sw(Sw_res[p], Swc, epsilon, pcwmin, pcwmax, Bpcw);
						pg_res[p] = po_res[p];//+ calc_gas_pcg(Sw_res[p], So_res[p], Swc, Sor, epsilon, pcgmin, pcgmax, Bpcg);
						if (po_res[p] <= pb_res[p]) acimapb[p] = 0;
						else acimapb[p] = 1;
						if(j==0) Pb_WB[k] = pb_res[p];
						p++;

					}
				}
			}			

			deltat = factor*deltat;
			
			for (iy = 0; iy<dim; iy++) {
				Azinterno[iy] = factor*Azinterno[iy];
				Ax[iy] = factor*Ax[iy];
				Aphi[iy] = factor*Aphi[iy];
				Axint[iy] = factor*Axint[iy];
		
			}			
			
			calc_reser_evol(int(1 / factor), p, y, u, l, centerP, centerSw, centerSo,
				IA, JA, A, B, X,
				dontstopmenow, i, k, j, o, m, n, dim,
				controladorP, controladorSw, controladorSo,
				lambda1_w, lambda2_w, lambda4_w, lambda5_w,
				lambda1_o, lambda2_o, lambda3_o, lambda4_o,
				lambda1_g, lambda2_g, lambda3_g, lambda4_g,
				lambda5_g, lambda6_g, lambda7_g, lambda8_g,
				Kx, Kz, pi, T, PWB, PWB_0, Pb_WB,
				po_res, pb_res, pw_res, pg_res, Sw_res, So_res, Sg_res,
				po_res_0, pb_res_0, Sw_res_0, So_res_0,
				po_res_lv0, pb_res_lv0, So_res_lv0, Sw_res_lv0,
				varpo, varpb, varSo, varSw,
				C1, C2, C3, C4, C5, C6,
				deltaphi, deltaz, deltat, a,
				Ax, Aphi, Axint, Azinterno, At,
				g, H, Re, Rw, R0, pesoF, pesoB,
				rho_air, rho_ws, API, yg, yo, rho_os, rho_gs, yn2, yco2, yh2s,
				pgr, krwmax, nw, ng, now, nog, phi_ref, Cr,
				Swc, Sor, Sgr, epsilon, pcwmin, pcwmax, pcgmin, pcgmax,
				parametrosBw, parametrosmig, parametrosZ, ppc, Tpc, Bwi, Bpcw, Bpcg,
				Qwsinj, Qosinj, Qgsinj,
				iparm, dparm, solver, maxfct, mnum, phase, error,
				msglvl, mtype, nrhs, ddum, idum, pt,
				Qw_well, Qo_well, Qg_well, coup_ts, ite_coup,
				phi, der_phi, rho_w, der_rho_w, mi_w, der_mi_w,
				Bw, der_Bw, pcw, der_pcw, krw, der_krw,
				Rs, der_Rs_p, der_Rs_pb, kro, der_kro_Sw, der_kro_So,
				Bo, der_Bo_p, der_Bo_pb,
				mi_o, der_mi_o_p, der_mi_o_pb, rho_o, der_rho_o_p, der_rho_o_pb,
				pcg, der_pcg_Swo, krg, der_krg, mi_g, der_mi_g, Z,
				Bg, der_Bg, rho_g, der_rho_g,
				acimapb, zonaprodutora, (level + 1),
				unW, ufW, uwbW, ueW, ucW, uwW, uebW, ubW, usW,
				unO, ufO, uwbO, ueO, ucO, uwO, uebO, ubO, usO,
				unG, ufG, uwbG, ueG, ucG, uwG, uebG, ubG, usG);
			
			nan = 0;
			deltat = deltat/ factor;			

			for (int iy = 0; iy<dim; iy++) {
				Azinterno[iy] = Azinterno[iy]/ factor;
				Ax[iy] = Ax[iy] / factor;
				Aphi[iy] = Aphi[iy] / factor;
				Axint[iy] = Axint[iy] / factor;
				varpo[p] = varpo[p] / factor;
				varpb[p] = varpb[p] / factor;
				varSo[p] = varSo[p] / factor;
				varSw[p] = varSw[p] / factor;
			}		
		}

		//cout << "5" << endl;
		if (level >= 2) {
				
				Sleep(3 * 100);
		}

		//cout << "6" << endl;
		p = 0;
		for (k = 0; k <= (o - 1); k++) {
			for (j = 0; j <= (m - 1); j++) {
				for (i = 0; i <= (n - 1); i++) {
					po_res_0[p] = po_res[p];
					Sw_res_0[p] = Sw_res[p];
					So_res_0[p] = So_res[p];
					pb_res_0[p] = pb_res[p];
					p++;
				}
			}
		}

		//cout << "7" << endl;
		if (level == 0) {

			for (k = 0; k <= (o - 1); k++) PWB_0[k] = PWB_0[k] - (PWB_0[k] - PWB[k]) / res_timesteps;
		
			if ((t == 0 || t == (res_timesteps-1)) && coup_ts % 20 == 0) {
				
				grava_p(po_res, dim, 0, (t + 1), ite_coup, coup_ts);
				grava_pb(pb_res, dim, 0, (t + 1), ite_coup, coup_ts);
				grava_Sw(Sw_res, dim, 0, (t + 1), ite_coup, coup_ts);
				grava_So(So_res, dim, 0, (t + 1), ite_coup, coup_ts);
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
			}	

		}

		dontstopmenow = 1;

	}

}


