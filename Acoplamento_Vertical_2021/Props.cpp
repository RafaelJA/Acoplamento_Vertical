#include "stdafx.h"
#include "Props.h"


double Bubble_Point(double T, double API, double SG_G, double GOR)
{
	double T_far, GOR_SCF, PB_psia, Yg;
	T_far = Kel_Far(T);
	GOR_SCF = bbl_ft3(GOR);

	
	Yg = 0.00091l*T_far - 0.0125l*API;
	PB_psia = 18.2l*(pow(GOR_SCF / SG_G, pow(1.204819l, -1.l))*pow(10.l, Yg)-1.4l);
	
	return Psia_Pa(PB_psia);
}

double Tensuper(double rhol, double rhog)
{
	return pow(560. / 200.*(rhol - rhog)*1.e-3, 4.f)*1.e-3;
}

void H_CP_Black(double T, double API, double * Hl, double * Cp)
{
	double T_Far = Kel_Far(T);
	double Tref_Far = Kel_Far(273.15l);
	double SG_O = 141.5l / (131.5l + API);
	double Watsonk = 11.9l;
	double A1, A2, A3, Cp_Btu, Hl_Btu;
	A1 = 0.055l*Watsonk + 0.35l;//aqui
	A2 = 0.6811l - 0.308l*SG_O;
	A3 = (8.15l - 3.06l*SG_O)*pow(10.l, -4.l);
	//Calcula a entalpia e o Cp no sistema britanico de unidades térmicas
	Cp_Btu = A1*(A2 + A3*T_Far);
	Hl_Btu = A1*A2*(T_Far - Tref_Far) + A1*A3 / 2.*(pow(T_Far, 2.) - pow(Tref_Far, 2.));
	//convertendo do sistema britânico para SI
	*Cp = Cp_Btu*4.1868*1e+3;
	*Hl = Hl_Btu*2.326*1e+3;
}

double Rs_Quality(double GOR, double x, double rhol_st, double rhog_st)
{
	return GOR - x*(rhol_st / rhog_st + GOR);
}

double Quality_Rs(double GOR, double Rs, double rhol_st, double rhog_st)
{
	return max(0.0l, (GOR - Rs) / (rhol_st / rhog_st + GOR));
}

void water_H_CP(double P, double T, double * Hlw, double * Hlgw, double * Cplw, double * Cpgw)
{
	double T_C = T - 273.2l;

	*Cplw = 4179. + 1.282e-2l*(T_C - 40.)*(T_C - 40.);

	*Hlw = (*Cpgw)*T_C;

	*Cpgw = 1864. + .0802l*pow(T_C, 1.5l)*exp(T_C / 133.);

	*Hlgw = 2.5018e+6l - 2.386e+3*T_C - 5.254e-5l*pow(T_C, 4.);
}

double tens_water(double P, double T)
{
	double t = T - 273.2l;
	double ta, tc, tx;
	ta = t + 273.15;
	tc = 374.15 + 273.15;
	tx = 1. - ta / tc;
	return 0.23580*pow(tx, 1.256)*(1. - 0.625*tx);
}

double Z_factor(double P, double T, double SG_G)
{
	double P_psia, T_Far, PPC, PPR, TPC, TPR, A, B, C, D, E, F, G, err, Dr, Fdr, Dfdr, Drn;
	int count = 0;

	P_psia = Pa_Psia(P);
	T_Far = Kel_Far(T);

	PPC = 667. + 15.*SG_G - 37.5f*SG_G*SG_G;
	TPC = 168. + 325.* SG_G - 12.5f*SG_G*SG_G;

	PPR = P_psia / PPC;
	TPR = (T_Far + 460.) / TPC;

	A = 0.06423f;
	B = 0.5353f*TPR - 0.6123f;
	C = 0.3151f*TPR - 1.0467f - 0.5783f / (TPR*TPR);
	D = TPR;
	E = 0.6816f / (TPR*TPR);
	F = 0.6845f;
	G = 0.27f*PPR;

	Dr = G / TPR;

	do
	{
		Fdr = A*pow(Dr, 6.f) + B*pow(Dr, 3.f) + C*pow(Dr, 2.f) + D*Dr + E*pow(Dr, 3.f)*(1 + pow(F, 2.f))*exp(-F*pow(Dr, 2.f)) - G;
		Dfdr = 6.*A*pow(Dr, 5.f) + 3.f*B*pow(Dr, 2.f) + 2.f*C*Dr + D + E*pow(Dr, 2.f)*(3.f + F*pow(Dr, 2.f)*(3.f - 2.f*F*pow(Dr, 2.f)))*exp(-F*pow(Dr, 2.f));
		Drn = Dr - Fdr / Dfdr;
		err = fabs((Drn - Dr) / Drn);
		Dr = Drn;
		count++;
	} while (err > (1e-6) && count < 50);

	return 2.7e-1f*PPR / Dr / TPR;
}
