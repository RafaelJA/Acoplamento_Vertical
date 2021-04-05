#include "stdafx.h"
#include "props_res.h"
#include <omp.h>

//#include <math.h>

//Press�o(psia)(na verdade foda-se)
double calc_phi_p(double phi_ref, double Cr, double p, double pi) {
	return (phi_ref*exp(Cr*(p - pi)));
}

//Press�o(psia)(na verdade foda-se)
double calc_der_phi_p(double phi_ref, double Cr, double p, double pi) {
	return (Cr*calc_phi_p(phi_ref, Cr, p, pi));
}

//Press�o(psia), Temperatura(~F}, Bw(Rb/STB), rho_ws e sa�da (g/cm3)
double calc_water_rho_p(double rho_ws, double Bw) {
	return (rho_ws / Bw);
}

//Press�o(psia), Temperatura(~F}, Bw(Rb/STB), rho_ws e sa�da (g/cm3)
double calc_der_water_rho_p(double rho_ws, double Bw, double der_Bw) {
	return rho_ws * (-der_Bw/(Bw*Bw));
}

//
double calc_oil_rho_p(double X_pb, double po, double pb, double co, double rho_os, double rho_gs,
	double *p_Bo, double *Bo_tab, int length_Bo,
	double *p_Rso, double *Rso_tab, int length_Rso){
	return X_pb*((rho_os / calcula_prop(pb, p_Bo, Bo_tab, length_Bo)) + (rho_gs*ft3_bbl(calcula_prop(pb, p_Rso, Rso_tab, length_Rso)) / calcula_prop(pb, p_Bo, Bo_tab, length_Bo)))*(1.l+co*(po-pb))
		+(1.l-X_pb)*((rho_os / calcula_prop(po, p_Bo, Bo_tab, length_Bo)) + (rho_gs*ft3_bbl(calcula_prop(po, p_Rso, Rso_tab, length_Rso)) / calcula_prop(po, p_Bo, Bo_tab, length_Bo)));
}

//
double calc_der_oil_rho_p(double X_pb, double po, double pb, double co, double rho_os, double rho_gs,
	double *p_Bo, double *Bo_tab, int length_Bo,
	double *p_Rso, double *Rso_tab, int length_Rso) {

	return X_pb*((rho_os / calcula_prop(pb, p_Bo, Bo_tab, length_Bo)) + (rho_gs*ft3_bbl(calcula_prop(pb, p_Rso, Rso_tab, length_Rso)) / calcula_prop(pb, p_Bo, Bo_tab, length_Bo)))*co
	+ (1.l - X_pb)*(((rho_os + (rho_gs * ft3_bbl(calcula_prop(po, p_Rso, Rso_tab, length_Rso))))*(-calcula_der_prop(po, p_Bo, Bo_tab, length_Bo) / (calcula_prop(po, p_Bo, Bo_tab, length_Bo)*calcula_prop(po, p_Bo, Bo_tab, length_Bo)))) + (rho_gs*(ft3_bbl(calcula_der_prop(po, p_Rso, Rso_tab, length_Rso)) / calcula_prop(po, p_Bo, Bo_tab, length_Bo))));
}

double calc_der_oil_rho_pb(double X_pb, double po, double pb, double co, double rho_os, double rho_gs,
	double *p_Bo, double *Bo_tab, int length_Bo,
	double *p_Rso, double *Rso_tab, int length_Rso) {

	return -X_pb*((rho_os / calcula_prop(pb, p_Bo, Bo_tab, length_Bo)) + (rho_gs*ft3_bbl(calcula_prop(pb, p_Rso, Rso_tab, length_Rso)) / calcula_prop(pb, p_Bo, Bo_tab, length_Bo)))*co
		+ X_pb*(((rho_os + (rho_gs * ft3_bbl(calcula_prop(pb, p_Rso, Rso_tab, length_Rso))))*(-calcula_der_prop(pb, p_Bo, Bo_tab, length_Bo) / (calcula_prop(pb, p_Bo, Bo_tab, length_Bo)*calcula_prop(pb, p_Bo, Bo_tab, length_Bo)))) + (rho_gs*(ft3_bbl(calcula_der_prop(pb, p_Rso, Rso_tab, length_Rso)) / calcula_prop(pb, p_Bo, Bo_tab, length_Bo))))*(1.l + co*(po - pb));
}

//
double calc_gas_rho_p(double rho_gs, double Bg) {
	return rho_gs / Bg;
}

//
double calc_der_gas_rho_p(double rho_gs, double Bg, double der_Bg) {
	return -rho_gs*der_Bg / (Bg*Bg);
}

double calc_oil_Bo(double X_pb, double po, double pb, double co, double *p_Bo, double *Bo_tab, int length_Bo) {

	return (1.l - X_pb)*calcula_prop(po, p_Bo, Bo_tab, length_Bo)
		+ X_pb * calcula_prop(pb, p_Bo, Bo_tab, length_Bo)*(1.l - co * (po - pb));

}

double calc_der_oil_Bo_p(double X_pb, double po,double pb, double co, double *p_Bo, double *Bo_tab, int length_Bo) {

	return (1.l - X_pb)*calcula_der_prop(po, p_Bo, Bo_tab, length_Bo) - X_pb*co*calcula_prop(pb, p_Bo, Bo_tab, length_Bo);

}

double calc_der_oil_Bo_pb(double X_pb, double po, double pb, double co, double *p_Bo, double *Bo_tab, int length_Bo) {
	
	return X_pb * calcula_der_prop(pb, p_Bo, Bo_tab, length_Bo)*(1.l - co*X_pb*(po - pb)) + X_pb*co*calcula_prop(pb, p_Bo, Bo_tab, length_Bo);

}

double calc_oil_Rs(double X_pb, double po, double pb, double *p_Rso, double *Rso_tab, int length_Rso) {

	return ft3_bbl((1.l - X_pb)*calcula_prop(po, p_Rso, Rso_tab, length_Rso)
		+ X_pb * calcula_prop(pb, p_Rso, Rso_tab, length_Rso));

}

double calc_der_oil_Rs_p(double X_pb, double po, double pb, double *p_Rso, double *Rso_tab, int length_Rso) {

	return ft3_bbl((1.l - X_pb)*calcula_der_prop(po, p_Rso, Rso_tab, length_Rso));

}

double calc_der_oil_Rs_pb(double X_pb, double po, double pb, double *p_Rso, double *Rso_tab, int length_Rso) {

	return ft3_bbl((X_pb)*calcula_der_prop(pb, p_Rso, Rso_tab, length_Rso));

}

double calc_oil_muo(double X_pb, double po, double pb, double c_mu, double *p_muo, double *muo_tab, int length_muo) {

	return (1.l - X_pb)*calcula_prop(po, p_muo, muo_tab, length_muo)
		+ X_pb * calcula_prop(pb, p_muo, muo_tab, length_muo)*(1.l + c_mu * (po - pb));

}

double calc_der_oil_muo_p(double X_pb, double po, double pb, double c_mu, double *p_muo, double *muo_tab, int length_muo) {

	return (1.l - X_pb)*calcula_der_prop(po, p_muo, muo_tab, length_muo) + X_pb*c_mu*calcula_prop(pb, p_muo, muo_tab, length_muo);

}

double calc_der_oil_muo_pb(double X_pb, double po, double pb, double c_mu, double *p_muo, double *muo_tab, int length_muo) {

	return X_pb * calcula_der_prop(pb, p_muo, muo_tab, length_muo)*(1.l + c_mu*X_pb*(po - pb)) - X_pb*c_mu*calcula_prop(pb, p_muo, muo_tab, length_muo);

}


void calc_props_grid(int N, double *po, double *pg, double *pw, double *pb, double *Sw, double *So, double *Sg, double *X_pb,
	double *phi, double *der_phi,
	double *rho_w, double *der_rho_w, double *mi_w, double *der_mi_w, double *Bw, double *der_Bw, double *pcw, double *der_pcw, double *krw, double *der_krw,
	double *Rs, double *der_Rs_p, double *der_Rs_pb, double *kro, double *der_kro_Sw, double *der_kro_So, double *Bo, double *der_Bo_p, double *der_Bo_pb,
	double *mi_o, double *der_mi_o_p, double *der_mi_o_pb, double *rho_o, double *der_rho_o_p, double *der_rho_o_pb,
	double *pcg, double *der_pcg_Sw, double *der_pcg_So, double *krg, double *der_krg, double *mi_g, double *der_mi_g, double *Bg, double *der_Bg, double *rho_g, double *der_rho_g,
	double *phi_ref, double Cr, double pi,
	double rho_ws, double rho_os, double rho_gs,
	double *p_Bg, double *Bg_tab, int length_Bg,
	double *p_Bo, double *Bo_tab, int length_Bo,
	double *p_Bw, double *Bw_tab, int length_Bw,
	double *Sg_krg, double *krg_tab, int length_krg,
	double *Sg_krog, double *krog_tab, int length_krog,
	double *Sw_krw, double *krw_tab, int length_krw,
	double *Sw_krow, double *krow_tab, int length_krow,
	double *p_mug, double *mug_tab, int length_mug,
	double *p_muo, double *muo_tab, int length_muo,
	double *Sg_pcgo, double *pcgo_tab, int length_pcgo,
	double *Sw_pcow, double *pcow_tab, int length_pcow,
	double *p_Rso, double *Rso_tab, int length_Rso,
	double *lambda1_W, double *d_lambda1_po_W, double *d_lambda1_SW_W,
	double *lambda2_W, double *d_lambda2_po_W, double *d_lambda2_SW_W,
	double *lambda1_O, double *d_lambda1_po_O, double *d_lambda1_SW_O, double *d_lambda1_pb_O, double *d_lambda1_SO_O,
	double *lambda2_O, double *d_lambda2_po_O, double *d_lambda2_SW_O, double *d_lambda2_pb_O, double *d_lambda2_SO_O,
	double *lambda1_G, double *d_lambda1_po_G, double *d_lambda1_SW_G, double *d_lambda1_pb_G, double *d_lambda1_SO_G,
	double *lambda2_G, double *d_lambda2_po_G, double *d_lambda2_SW_G, double *d_lambda2_pb_G, double *d_lambda2_SO_G,
	double *lambda3_G, double *d_lambda3_po_G, double *d_lambda3_SW_G, double *d_lambda3_pb_G, double *d_lambda3_SO_G,
	double *lambda4_G, double *d_lambda4_po_G, double *d_lambda4_SW_G, double *d_lambda4_pb_G, double *d_lambda4_SO_G) {

	int p;
	double krog, krow, der_krog, der_krow;
	double co = 3.e-6;
	double cw = 4.e-6;
	double c_mu = 0.l;

#pragma omp parallel private(p) 
	{
#pragma omp for schedule(dynamic) nowait 
		for (p = 0; p <= (N - 1); p++)
		{

			if (po[p] > pb[p]) X_pb[p] = 1.0l;
			else X_pb[p] = 0.0l;

			phi[p] = calc_phi_p(phi_ref[p], Cr, po[p], pi);
			der_phi[p] = calc_der_phi_p(phi_ref[p], Cr, po[p], pi);

			mi_w[p] = 0.96l;
			der_mi_w[p] = 0.0l;

			Bw[p] = calcula_prop(po[p], p_Bw, Bw_tab, length_Bw);
			der_Bw[p] = calcula_der_prop(po[p], p_Bw, Bw_tab, length_Bw);

			rho_w[p] = calc_water_rho_p(rho_ws, Bw[p]);
			der_rho_w[p] = calc_der_water_rho_p(rho_ws, Bw[p], der_Bw[p]);

			pcw[p] = calcula_prop(Sw[p], Sw_pcow, pcow_tab, length_pcow);
			der_pcw[p] = calcula_der_prop(Sw[p], Sw_pcow, pcow_tab, length_pcow);

			pcg[p] = calcula_prop(Sg[p], Sg_pcgo, pcgo_tab, length_pcgo);
			der_pcg_Sw[p] = -(1.0l - X_pb[p])*calcula_der_prop(Sg[p], Sg_pcgo, pcgo_tab, length_pcgo);
			der_pcg_So[p] = -(1.0l- X_pb[p])*calcula_der_prop(Sg[p], Sg_pcgo, pcgo_tab, length_pcgo);

			pw[p] = po[p] - pcw[p];
			pg[p] = po[p] + pcg[p];

			krw[p] = calcula_prop(Sw[p], Sw_krw, krw_tab, length_krw);
			der_krw[p] = calcula_der_prop(Sw[p], Sw_krw, krw_tab, length_krw);

			krow = calcula_prop(Sw[p], Sw_krow, krow_tab, length_krow);
			der_krow = calcula_der_prop(Sw[p], Sw_krow, krow_tab, length_krow);

			krog = calcula_prop(Sg[p], Sg_krog, krog_tab, length_krog);
			der_krog = -calcula_der_prop(Sg[p], Sg_krog, krog_tab, length_krog);

			krg[p] = calcula_prop(Sg[p], Sg_krg, krg_tab, length_krg);
			der_krg[p] = -calcula_der_prop(Sg[p], Sg_krg, krg_tab, length_krg);

			kro[p] = ((krow + krw[p])*(krog + krg[p]) - (krw[p] + krg[p]));
			der_kro_Sw[p] = ((der_krow + der_krw[p])*(krog + krg[p]) + (krow + krw[p])*(der_krog + der_krg[p]) - (der_krw[p] + der_krg[p]));
			der_kro_So[p] = ((krow + krw[p])*(der_krog + der_krg[p]) - (der_krg[p]));
			if (kro[p] <= 0.0l) {
				kro[p] = 0.0l;
				der_kro_Sw[p] = -0.0l;
				der_kro_So[p] = 0.0l;
			}

			Rs[p] = calc_oil_Rs(X_pb[p], po[p], pb[p], p_Rso, Rso_tab, length_Rso);
			der_Rs_p[p] = calc_der_oil_Rs_p(X_pb[p], po[p], pb[p], p_Rso, Rso_tab, length_Rso);
			der_Rs_pb[p] = calc_der_oil_Rs_pb(X_pb[p], po[p], pb[p], p_Rso, Rso_tab, length_Rso);

			Bo[p] = calc_oil_Bo(X_pb[p], po[p], pb[p], co, p_Bo, Bo_tab, length_Bo);
			der_Bo_p[p] = calc_der_oil_Bo_p(X_pb[p], po[p], pb[p], co, p_Bo, Bo_tab, length_Bo);
			der_Bo_pb[p] = calc_der_oil_Bo_pb(X_pb[p], po[p], pb[p], co, p_Bo, Bo_tab, length_Bo);

			mi_o[p] = calc_oil_muo(X_pb[p], po[p], pb[p], c_mu, p_muo, muo_tab, length_muo);
			der_mi_o_p[p] = calc_der_oil_muo_p(X_pb[p], po[p], pb[p], c_mu, p_muo, muo_tab, length_muo);
			der_mi_o_pb[p] = calc_der_oil_muo_pb(X_pb[p], po[p], pb[p], c_mu, p_muo, muo_tab, length_muo);

			rho_o[p] = calc_oil_rho_p(X_pb[p], po[p], pb[p], co, rho_os, rho_gs,
				p_Bo, Bo_tab, length_Bo,
				p_Rso, Rso_tab, length_Rso);
			der_rho_o_p[p] = calc_der_oil_rho_p(X_pb[p], po[p], pb[p], co, rho_os, rho_gs,
				p_Bo, Bo_tab, length_Bo,
				p_Rso, Rso_tab, length_Rso);
			der_rho_o_pb[p] = calc_der_oil_rho_pb(X_pb[p], po[p], pb[p], co, rho_os, rho_gs,
				p_Bo, Bo_tab, length_Bo,
				p_Rso, Rso_tab, length_Rso);/////

			mi_g[p] = calcula_prop(po[p], p_mug, mug_tab, length_mug);
			der_mi_g[p] = calcula_der_prop(po[p], p_mug, mug_tab, length_mug);

			Bg[p] = calcula_prop(po[p], p_Bg, Bg_tab, length_Bg);
			der_Bg[p] = calcula_der_prop(po[p], p_Bg, Bg_tab, length_Bg);

			rho_g[p] = calc_gas_rho_p(rho_gs, Bg[p]);
			der_rho_g[p] = calc_der_gas_rho_p(rho_gs, Bg[p], der_Bg[p]);

			lambda1_W[p] = krw[p] / (mi_w[p] * Bw[p]);
			d_lambda1_po_W[p] = -(der_mi_w[p] / mi_w[p] + der_Bw[p] / Bw[p])*lambda1_W[p];
			d_lambda1_SW_W[p] = der_krw[p] / (mi_w[p] * Bw[p]);

			lambda2_W[p] = rho_w[p] * krw[p] / (mi_w[p] * Bw[p]);
			d_lambda2_po_W[p] = der_rho_w[p] * lambda1_W[p] + rho_w[p] * d_lambda1_po_W[p];
			d_lambda2_SW_W[p] = rho_w[p] * d_lambda1_SW_W[p];

			lambda1_O[p] = kro[p] / (mi_o[p] * Bo[p]);
			d_lambda1_po_O[p] = -(der_mi_o_p[p] / mi_o[p] + der_Bo_p[p] / Bo[p])*lambda1_O[p];
			d_lambda1_SW_O[p] = der_kro_Sw[p] / (mi_o[p] * Bo[p]);
			d_lambda1_pb_O[p] = -(der_mi_o_pb[p] / mi_o[p] + der_Bo_pb[p] / Bo[p])*lambda1_O[p];
			d_lambda1_SO_O[p] = der_kro_So[p] / (mi_o[p] * Bo[p]);

			lambda2_O[p] = rho_o[p] * kro[p] / (mi_o[p] * Bo[p]);
			d_lambda2_po_O[p] = der_rho_o_p[p] * lambda1_O[p] + rho_o[p] * d_lambda1_po_O[p];
			d_lambda2_SW_O[p] = rho_o[p] * d_lambda1_SW_O[p];
			d_lambda2_pb_O[p] = der_rho_o_pb[p] * lambda1_O[p] + rho_o[p] * d_lambda1_pb_O[p];
			d_lambda2_SO_O[p] = rho_o[p] * d_lambda1_SO_O[p];

			lambda1_G[p] = Rs[p]*kro[p] / (mi_o[p] * Bo[p]);
			d_lambda1_po_G[p] = (der_Rs_p[p]/Rs[p] - der_mi_o_p[p] / mi_o[p] - der_Bo_p[p] / Bo[p])*lambda1_G[p];
			d_lambda1_SW_G[p] = Rs[p] * der_kro_Sw[p] / (mi_o[p] * Bo[p]);
			d_lambda1_pb_G[p] = (der_Rs_pb[p] / Rs[p] - der_mi_o_pb[p] / mi_o[p] - der_Bo_pb[p] / Bo[p])*lambda1_G[p];
			d_lambda1_SO_G[p] = Rs[p] * der_kro_So[p] / (mi_o[p] * Bo[p]);

			lambda2_G[p] = Rs[p] * rho_o[p] * kro[p] / (mi_o[p] * Bo[p]);
			d_lambda2_po_G[p] = der_rho_o_p[p] * lambda1_G[p] + rho_o[p] * d_lambda1_po_G[p];
			d_lambda2_SW_G[p] = rho_o[p] * d_lambda1_SW_G[p];
			d_lambda2_pb_G[p] = der_rho_o_pb[p] * lambda1_G[p] + rho_o[p] * d_lambda1_pb_G[p];
			d_lambda2_SO_G[p] = rho_o[p] * d_lambda1_SO_G[p];

			lambda3_G[p] = krg[p] / (mi_g[p] * Bg[p]);
			d_lambda3_po_G[p] = -(der_mi_g[p] / mi_g[p] + der_Bg[p] / Bg[p])*lambda3_G[p];
			d_lambda3_SW_G[p] = der_krg[p] / (mi_g[p] * Bg[p]);
			d_lambda3_pb_G[p] = 0.0l;
			d_lambda3_SO_G[p] = der_krg[p] / (mi_g[p] * Bg[p]);

			lambda4_G[p] = rho_g[p]*krg[p] / (mi_g[p] * Bg[p]);
			d_lambda4_po_G[p] = der_rho_g[p] * lambda3_G[p] + rho_g[p] * d_lambda3_po_G[p];
			d_lambda4_SW_G[p] = rho_g[p] * d_lambda3_SW_G[p];
			d_lambda4_pb_G[p] = 0.0l;
			d_lambda4_SO_G[p] = rho_g[p] * d_lambda3_SO_G[p];
			
		}
	}

	//grava_p(X_pb, N, 0, 0, 0, 0);
	//system("pause");


}

double perm_harm(double Kp, double Knm) {
	
	if (Kp < 0.000000001 || Knm < 0.0000000001) return 0.0;
	
		return 2.0l / ((1.0l / Kp) + (1.0l / Knm));
}

//https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
double generateGaussianNoise(double mu, double sigma)
{
	const double eps = 0.00000000000005;
	const double two_pi = 2.0*3.14159265358979323846;

	static double z0, z1;
	static bool generate;
	generate = !generate;

	if (!generate)
		return z1 * sigma + mu;

	double u1, u2;
	do
	{
		u1 = rand() * (1.0 / RAND_MAX);
		u2 = rand() * (1.0 / RAND_MAX);
	} while (u1 <= eps);

	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
	return z0 * sigma + mu;
}

void get_perm(double *k, int dimx, int dimy, int dimz)
{
	char nome_in[] = "permeabilidade.txt";
	int pos = 0;
	ifstream fpin(nome_in);
	for (int z = 0; z<dimz; z++) for (int y = 0; y<dimy; y++) for (int x = 0; x<dimx; x++) fpin >> k[pos], pos++;
}

void gera_permeabilidades(double *k, int dimx, int dimy, int dimz, double perm, double varperm)
{
	int pos = 0;
	char nome_out[20], nome_out2[20];
	sprintf_s(nome_out, "permeabilidade.m");
	sprintf_s(nome_out2, "permeabilidade.txt");
	ofstream fper(nome_out);
	ofstream fper2(nome_out2);
	fper2.precision(17);
	fper.precision(17);
	fper << "k=[" << endl;
	for (int z = 0; z<dimz; z++) for (int y = 0; y<dimy; y++) for (int x = 0; x<dimx; x++) k[pos] = generateGaussianNoise(perm, varperm), fper << k[pos] << endl, fper2 << k[pos] << endl, pos++;
	fper << "];" << endl;
	fper << "k=reshape(k," << dimx << "," << dimy << "," << dimz << ");";
	fper.close();
}

void grava_p(double *p, int fim, int inicio, int j, int controlador, int zy)
{
	char nomep[50];


	sprintf_s(nomep, "p_%d_%d_%d.m", j, zy, controlador);

	ofstream fp(nomep);
	fp.precision(14);

	fp << "p = [ " << endl;
	for (int i = inicio; i<fim; i++)
	{
		fp << p[i];
		fp << endl;
	}

	cout << "\n\nGravando pressao." << endl;

	fp << "];" << endl;

	fp.close();
}

void grava_pb(double *p, int fim, int inicio, int j, int controlador, int zy)
{
	char nomep[50];


	sprintf_s(nomep, "pb_%d_%d_%d.m", j, zy, controlador);

	ofstream fp(nomep);
	fp.precision(14);

	fp << "p = [ " << endl;
	for (int i = inicio; i<fim; i++)
	{
		fp << p[i];
		fp << endl;
	}

	cout << "\n\nGravando pressao de bolha." << endl;

	fp << "];" << endl;

	fp.close();
}

void grava_Sw(double *p, int fim, int inicio, int j, int controlador, int zy)
{
	char nomep[50];


	sprintf_s(nomep, "Sw_%d_%d_%d.m", j, zy, controlador);

	ofstream fp(nomep);
	fp.precision(14);

	fp << "Sw = [ " << endl;
	for (int i = inicio; i<fim; i++)
	{
		fp << p[i];
		fp << endl;
	}


	cout << "\n\nGravando Saturacao agua." << endl;

	fp << "];" << endl;

	fp.close();
}

void grava_So(double *p, int fim, int inicio, int j, int controlador, int zy)
{
	char nomep[50];


	sprintf_s(nomep, "So_%d_%d_%d.m", j, zy, controlador);

	ofstream fp(nomep);
	fp.precision(14);

	fp << "So = [ " << endl;
	for (int i = inicio; i<fim; i++)
	{
		fp << p[i];
		fp << endl;
	}


	cout << "\n\nGravando Saturacao oleo." << endl;

	fp << "];" << endl;

	fp.close();
}

void grava_ia(int *p, int fim, int inicio, int j)
{
	char nomep[50];


	sprintf_s(nomep, "ia_%d.m", j);

	ofstream fp(nomep);
	fp.precision(14);

	fp << "ia = [ " << endl;
	for (int i = inicio; i<fim; i++)
	{
		fp << p[i];
		fp << endl;
	}

	cout << "\n\nGravando IA." << endl;

	fp << "];" << endl;

	fp.close();
}

void grava_a(double *p, int fim, int inicio, int j)
{
	char nomep[50];


	sprintf_s(nomep, "a_%d.m", j);

	ofstream fp(nomep);
	fp.precision(14);

	fp << "a = [ " << endl;
	for (int i = inicio; i<fim; i++)
	{
		fp << p[i];
		fp << endl;
	}


	cout << "\n\nGravando A." << endl;

	fp << "];" << endl;

	fp.close();
}

void grava_b(double *p, int fim, int inicio, int j)
{
	char nomep[50];


	sprintf_s(nomep, "b_%d.m", j);

	ofstream fp(nomep);
	fp.precision(14);

	fp << "b = [ " << endl;
	for (int i = inicio; i<fim; i++)
	{
		fp << p[i];
		fp << endl;
	}


	cout << "\n\nGravando B." << endl;

	fp << "];" << endl;

	fp.close();
}

void grava_e(double *p, int fim, int inicio, int j)
{
	char nomep[50];


	sprintf_s(nomep, "e_%d.m", j);

	ofstream fp(nomep);
	fp.precision(14);

	fp << "e = [ " << endl;
	for (int i = inicio; i<fim; i++)
	{
		fp << p[i];
		fp << endl;
	}


	cout << "\n\nGravando E." << endl;

	fp << "];" << endl;

	fp.close();
}

void grava_ja(int *p, int fim, int inicio, int j)
{
	char nomep[50];


	sprintf_s(nomep, "ja_%d.m", j);

	ofstream fp(nomep);
	fp.precision(14);

	fp << "ja = [ " << endl;
	for (int i = inicio; i<fim; i++)
	{
		fp << p[i];
		fp << endl;
	}


	cout << "\n\nGravando JA." << endl;

	fp << "];" << endl;

	fp.close();
}

