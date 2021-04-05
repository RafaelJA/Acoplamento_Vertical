#include "stdafx.h"
#include "Acoplamento_poco.h"

double converge_p_sup(double L_in, int nz_in, double dT_in, double P_in_in, double T_in_in, double OFR_in, double GOR_in, double API_in,
	double SG_G_in, double WC_in, double P_surf_in, double P_initial, double rug_in, int coup_ts, int iter_coup, __int8 correlation, char directory[])
{


	if (coup_ts % 50 == 0) {

		char arc[200];
		char arc2[30];
		sprintf_s(arc, directory);

		sprintf_s(arc2, "\\well_data_inlet_%d_%d.txt", coup_ts, iter_coup);

		strcat_s(arc, arc2);

		fstream data(arc, ios::out);
		data.precision(8);

		data << "---parâmetros iniciais---" << endl;
		data << "L=\t" << L_in << "\t|m" << endl;
		data << "dT=\t" << dT_in << "\t|m" << endl;
		data << "nz=\t" << nz_in << "\t|" << endl;
		data << "P_in=\t" << P_in_in << "\t|Pa" << endl;
		data << "T_in=\t" << T_in_in << "\t|K" << endl;
		data << "OFR=\t" << OFR_in << "\t|bbl/d" << endl;
		data << "GOR=\t" << GOR_in << "\t|SCF/stb" << endl;
		data << "API=\t" << API_in << "\t|" << endl;
		data << "SG_G=\t" << SG_G_in << "\t|" << endl;
		data << "WC=\t" << WC_in << "\t|" << endl;

		data.close();
	}


	double deltap = 1000.l * 6894.l;
	cout.precision(4);
	double P_sup, P_inf;
	P_sup = P_in_in;

	double P_surf;
	double P_chut;
	double P_chut_n = 0.;
	bool error_SL_CH = 1;
	long double error;
	int count = 0;
	Well chutes;


	//cout << "Começo das iterações"<< endl;
	for (double P = P_sup - deltap; P > 0; P = P_sup - deltap)
	{
		chutes.clear_mem_well();

		P_inf = P;
		//cout << "Começo dentro do intervalo" << endl;
		do
		{
			count++;

			P_chut = P_chut_n;
			P_chut_n = (P_sup + P_inf) * .5l;


			chutes.inicializa_poco(L_in, nz_in, dT_in, P_chut_n, T_in_in, OFR_in, GOR_in, API_in, SG_G_in, WC_in, &error_SL_CH, P_initial, rug_in);
			if (correlation == 1) chutes.integration_HB_m(1, nz_in);
			else if (correlation == 2) chutes.integration_GRAMP(1, nz_in, &error_SL_CH);
			else if (correlation == 3) chutes.integration_ChFr(1, nz_in);
			else if (correlation == 4) chutes.integration_BeB(1, nz_in,1);

			P_surf = chutes.get_surf_P();
			error = P_surf - P_surf_in;

			if (error > 0.) P_sup = P_chut_n;
			else if (error < 0. || isnan(P_surf)) P_inf = P_chut_n;
			else break;

			cout.precision(4);

		} while ((fabs(error) >= 6894.l || isnan(P_surf)) && count < 20);

		count = 0;

		if (error >= 6894.l)
		{
			P_sup = P_inf;
			error_SL_CH = 1;

		}
		else if (fabs(error) >= 6894.l || isnan(P_surf))
		{
			P_sup += deltap;
			deltap = deltap / 2;
			error_SL_CH = 1;

		}
		else if (error_SL_CH == 0)
		{
			P_sup += deltap;
			deltap = deltap / 2;
			error_SL_CH = 1;

		}
		else break;

		if (abs(deltap) < 10e-10) break;

	}

	if (coup_ts % 50 == 0) chutes.grav_tudo(coup_ts, iter_coup, directory);

	double retorno = chutes.get_bot_P();

	chutes.~Well();

	return retorno;

}
