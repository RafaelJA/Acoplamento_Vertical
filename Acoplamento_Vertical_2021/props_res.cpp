#include "stdafx.h"
#include "props_res.h"
#include <omp.h>

//#include <math.h>

//Press�o(psia)(na verdade foda-se)
double calc_phi_p(double phi_ref, double Cr, double p, double pi) {
	return (phi_ref * exp(Cr * (p - pi)));
}

//Press�o(psia)(na verdade foda-se)
double calc_der_phi_p(double phi_ref, double Cr, double p, double pi) {
	return (Cr * calc_phi_p(phi_ref, Cr, p, pi));
}

//Press�o(psia), Temperatura(~F}
double calc_water_Rsw_p(double T, double p) {
	double Arsw = 2.12l + 3.45l * pow(10.0l, (-3.0l)) * T - 3.59l * pow(10.0l, (-5.0l)) * T * T;
	double Brsw = 0.0107l - 5.26l * pow(10.0l, (-5.0l)) * T + 1.48l * pow(10, (-7.0l)) * T * T;
	double Crsw = -8.75 * pow(10, (-7)) + 3.9 * pow(10, (-9)) * T - 1.02 * pow(10, (-11)) * T * T;
	return (Arsw + Brsw * p + Crsw * p * p);
}

//Press�o(psia), Temperatura(~F}
double calc_der_water_Rsw_p(double T, double p) {
	double Brsw = 0.0107 - 5.26 * pow(10, (-5)) * T + 1.48 * pow(10, (-7)) * T * T;
	double Crsw = -8.75 * pow(10, (-7)) + 3.9 * pow(10, (-9)) * T - 1.02 * pow(10, (-11)) * T * T;
	return (Brsw + 2. * Crsw * p);
}

//Press�o(psia), Temperatura(~F}
double calc_water_hat_p(double T, double p) {
	double Ahat = 3.8546 - 1.34 * pow(10, (-4)) * p;
	double Bhat = -0.01052 + 4.77 * pow(10, (-7)) * p;
	double Chat = 3.9267 * pow(10, (-5)) - 8.8 * pow(10, (-10)) * p;
	return (Ahat + Bhat * T + Chat * T * T);
}

//Press�o(psia), Temperatura(~F}
double calc_der_water_hat_p(double T) {
	return (-1.34 * pow(10, (-4)) + 4.77 * pow(10, (-7)) * T - 8.8 * pow(10, (-10)) * T * T);
}

//Press�o(psia), Temperatura(~F}
double calc_water_compr_p(double Tf, double p) {
	return ((1.l + 0.0089l * calc_water_Rsw_p(Tf, p)) * calc_water_hat_p(Tf, p) * pow(10, (-6)));
}

//Press�o(psia), Temperatura(~F}
double calc_der_water_compr_p(double Tf, double p) {
	return (((1.l + 0.0089l * calc_water_Rsw_p(Tf, p)) * calc_der_water_hat_p(Tf) + calc_water_hat_p(Tf, p) * calc_der_water_Rsw_p(Tf, p)) * pow(10, (-6)));
}

//Press�o(psia), Temperatura(~F}, Bw(Rb/STB), rho_ws e sa�da (g/cm3)
double calc_water_rho_p(double Bwi, double rho_ws, double T, double p, double pi) {
	return ((rho_ws / Bwi) * (1.0l + calc_water_compr_p(T, p) * (p - pi)));
}

//Press�o(psia), Temperatura(~F}, Bw(Rb/STB), rho_ws e sa�da (g/cm3)
double calc_der_water_rho_p(double Bwi, double rho_ws, double T, double p, double pi) {
	return ((rho_ws / Bwi) * (calc_water_compr_p(T, p) + calc_der_water_compr_p(T, p) * (p - pi)));
}

//Press�o(psia), Temperatura(~F} Sa�da cp
double calc_water_mi_p(double T, double p) {
	double Tk = 273.15 + ((T - 32.0) / 1.8);
	return(0.02414 * pow(10, (247.8 / (Tk - 140.))) * (1 + 3.5 * pow(10, (-12)) * (T - 40.0l) * p * p));
}

//Press�o(psia), Temperatura(~F}, Sa�da cp
double calc_der_water_mi_p(double T, double p) {
	double Tk = 273.15 + ((T - 32.0) / 1.8);
	return(2 * 0.02414 * pow(10, (247.8 / (Tk - 140.))) * 3.5 * (1.e-12l) * (T - 40.0l) * p);
}

//Press�o(psia), Temperatura(~F}, Bw(Rb/STB)
double calc_water_Bw_p(double p, double T, double* parametrosBw) {
	double A = parametrosBw[0] + parametrosBw[1] * T + parametrosBw[2] * T * T;
	double B = parametrosBw[3] + parametrosBw[4] * T + parametrosBw[5] * T * T;
	double C = parametrosBw[6] + parametrosBw[7] * T + parametrosBw[8] * T * T;
	return (A + B * p + C * p * p);
}

//Press�o(psia), Temperatura(~F}, Bw(Rb/STB)
double calc_der_water_Bw_p(double p, double T, double* parametrosBw) {
	double B = parametrosBw[3] + parametrosBw[4] * T + parametrosBw[5] * T * T;
	double C = parametrosBw[6] + parametrosBw[7] * T + parametrosBw[8] * T * T;
	return (B + 2 * C * p);
}

//
void calc_const_Bw(double* parametrosBw) {
	parametrosBw[0] = 0.9911;
	parametrosBw[1] = 6.35 * pow(10, (-5));
	parametrosBw[2] = 8.5 * pow(10, (-7));
	parametrosBw[3] = -1.093 * pow(10, (-6));
	parametrosBw[4] = -3.497 * pow(10, (-9));
	parametrosBw[5] = 4.57 * pow(10, (-12));
	parametrosBw[6] = -5.0 * pow(10, (-11));
	parametrosBw[7] = -6.429 * pow(10, (-13));
	parametrosBw[8] = -1.43 * pow(10, (-15));
}

//Press�o(psia), Resto adimensionaal
double calc_water_pcw_Sw(double Sw, double Swc, double epsilon, double pcwmin, double pcwmax, double B) {
	if (Swc >= Sw) return pcwmax;
	else return (pcwmin + B * log((Sw - Swc + epsilon) / (1 - Swc)));
}

//Press�o(psia), Resto adimensionaal
double calc_der_water_pcw_Sw(double Sw, double Swc, double epsilon, double B) {
	if (Swc >= Sw)  0.0;
	return (B / (Sw - Swc + epsilon));
}

//Press�o(psia), Resto adimensionaal
double calc_Bpcw(double Swc, double epsilon, double pcwmin, double pcwmax) {
	return((pcwmax - pcwmin) / (log(epsilon / (1. - Swc))));
}

//Press�o(psia), Resto adimensionaal
double calc_gas_pcg(double Sw, double So, double Swc, double Sor, double epsilon, double pcgmin, double pcgmax, double B) {
	if ((Sw + So) <= (Swc + Sor)) return (pcgmin + B * log((epsilon) / (1. - Swc - Sor)));
	else return (pcgmin + B * log((Sw + So - Sor - Swc + epsilon) / (1. - Swc - Sor)));
}

//Press�o(psia), Resto adimensionaal
double calc_der_gas_pcg_Swo(double Sw, double So, double Swc, double Sor, double epsilon, double B) {
	if ((Sw + So) <= (Swc + Sor))  return 0.0;
	else if ((1 - So - Sw) <= 0) return 0.0;
	else return (B / (So + Sw - Swc - Sor + epsilon));
}

//Press�o(psia), Resto adimensionaal
double calc_Bpcg(double Swc, double Sor, double epsilon, double pcgmin, double pcgmax) {
	return((pcgmax - pcgmin) / (log(epsilon / (1. - Swc - Sor))));
}

//Press�o(psia), Resto adimensionaal (t� correto)
double calc_water_krw_Sw(double krwmax, double Sw, double Swc, double Sor, double nw) {
	if (Swc >= Sw) return 0.0;
	else if (Sor >= (1 - Sw)) return krwmax;
	else return krwmax * pow(((Sw - Swc) / (1. - Sor - Swc)), nw);
}

double calc_water_krw_Sw(double Sw) {
	if (0.2116228 >= Sw) return 0.0l;
	else if (Sw >= 0.98387) return 1.0l;
	else return 7.18578687 * 0.00001 * pow(1.626403520 / 0.0001, Sw);
}

//Press�o(psia), Resto adimensionaal (t� certo)
double calc_der_water_krw_Sw(double krwmax, double Sw, double Swc, double Sor, double nw) {
	if (Swc >= Sw) return 0.00000000001;
	else if (Sor >= (1 - Sw)) return -0.000000000001;
	else return krwmax * nw * pow(((Sw - Swc) / (1. - Sor - Swc)), (nw - 1)) * (1. / (1. - Swc - Sor));
}

double calc_der_water_krw_Sw(double Sw) {
	if (0.2116228 >= Sw) return 0.00000000001l;
	else if (Sw >= 0.98387) return -0.0000000001l;
	else return 7.18578687 * 0.00001 * pow(1.626403520 / 0.0001, Sw) * log(1.626403520);
}

//Press�o(psia), Resto adimensionaal (t� certo)
double calc_gas_krg_Swo(double So, double Sw, double Sgr, double Sor, double Swc, double ng) {
	if (Sgr >= (1 - So - Sw)) return 0.00000;
	else if ((1 - So - Sw) >= (1 - Sor - Swc)) return pow(((1 - Sor - Swc - Sgr) / (1 - Sor - Swc - Sgr)), ng);
	else return pow(((1 - So - Sw - Sgr) / (1 - Sor - Swc - Sgr)), ng);
}

double calc_gas_krg_Swo(double Sg) {
	if (0.08675 >= Sg) return 0.0l;
	else if (Sg >= 0.78843) return 1.0l;
	else return 158.33844154 * pow(Sg, 6.0l) - 387.05577263 * pow(Sg, 5.0l) + 345.8214412 * pow(Sg, 4.0l) - 139.25334951 * pow(Sg, 3.0l) + 27.45435678 * pow(Sg, 2.0l) - 2.03393602 * Sg + 0.04299577;
}

double calc_der_gas_krg_Swo(double Sg) {
	if (0.08675 >= Sg) return -0.0000000001l;
	else if (Sg >= 0.78843) return 0.00000000001l;
	else return -6.0l * 158.33844154 * pow(Sg, 5.0l) + 5.0l * 387.05577263 * pow(Sg, 4.0l) - 4.0l * 345.8214412 * pow(Sg, 3.0l) + 3.0l * 139.25334951 * pow(Sg, 2.0l) - 2.0l * 27.45435678 * Sg + 2.03393602;
}

//Press�o(psia), Resto adimensionaal  (t� certo)
double calc_der_gas_krg_Swo(double So, double Sw, double Sgr, double Sor, double Swc, double ng) {
	if (Sgr >= (1 - So - Sw)) return -0.0000000001l;
	else if ((1 - So - Sw) >= (1 - Sor - Swc)) return 0.00000000001l;
	else	return (-ng * pow(((1 - So - Sw - Sgr) / (1 - Sor - Swc - Sgr)), (ng - 1))) * (1 / (1 - Sor - Swc - Sgr));
}

//Press�o(psia), Resto adimensionaal (t� certo)
double calc_water_krow_Sw(double Sw, double Swc, double Sor, double now) {
	if (Swc >= Sw) return 1.0;
	else if (Sor >= (1 - Sw)) return 0.0;
	else return pow(((1. - Sw - Sor) / (1. - Sor - Swc)), (now));
}

double calc_water_krow_Sw(double Sw) {
	if (0.2467 >= Sw) return 1.0;
	else if (Sw >= 0.843) return 0.0;
	else return 89.25490222 * pow(Sw, 5.0l) - 287.31217267 * pow(Sw, 4.0l) + 351.80268542 * pow(Sw, 3.0l) - 199.17941299 * pow(Sw, 2.0l) + 48.74763796 * Sw - 3.20427908;
}

//Press�o(psia), Resto adimensionaal (t� certo)
double calc_der_water_krow_Sw(double Sw, double Swc, double Sor, double now) {
	if (Swc >= Sw) return 0.000000001;
	else if (Sor >= (1 - Sw)) return -0.0000000001;
	else return ((-now) * (1. / (1. - Swc - Sor)) * pow(((1. - Sw - Sor) / (1. - Sor - Swc)), (now - 1)));
}

double calc_der_water_krow_Sw(double Sw) {
	if (0.2467 >= Sw) return -0.0000000001;
	else if (Sw >= 0.843) return 0.000000001;
	else return 5.0l * 89.25490222 * pow(Sw, 4.0l) - 4.0l * 287.31217267 * pow(Sw, 3.0l) + 3.0l * 351.80268542 * pow(Sw, 2.0l) - 2.0l * 199.17941299 * Sw + 48.74763796;
}

//Press�o(psia), Resto adimensionaal (t� certo)
double calc_gas_krog_Swo(double Sw, double So, double Swc, double Sor, double Sgr, double nog) {
	if (Sgr >= (1 - So - Sw)) return pow(((1.0 - Sgr - Sor - Swc) / (1. - Sor - Swc)), (nog));
	else if ((1 - So - Sw) >= (1 - Sor - Swc)) return 0.0;
	else return pow(((So + Sw - Sor - Swc) / (1. - Sor - Swc)), (nog));
}

//Press�o(psia), Resto adimensionaal (t� certo)
double calc_gas_krog_Swo(double Sg) {
	if (0.0142 >= Sg) return 1.0l;
	else if (Sg >= 0.6995) return 0.0l;
	return 9.96741762 * pow(Sg, 5.0l) - 16.55870922 * pow(Sg, 4.0l) + 4.26589291 * pow(Sg, 3.0l) + 6.86968292 * pow(Sg, 2.0l) - 5.04975874 * Sg + 1.00679655;
}

//Press�o(psia), Resto adimensionaal (t� certo)
double calc_der_gas_krog_Swo(double Sw, double So, double Swc, double Sor, double Sgr, double nog) {
	if (Sgr >= (1 - So - Sw)) return 0.0000000001;
	else if ((1 - So - Sw) >= (1 - Sor - Swc)) return -0.000000001;
	return ((nog) * (1. / (1. - Swc - Sor)) * pow(((So + Sw - Swc - Sor) / (1. - Sor - Swc)), (nog - 1)));
}

double calc_der_gas_krog_Swo(double Sg) {
	if (0.0142 >= Sg) return -0.000000001l;
	else if (Sg >= 0.6995) return 0.0000000001l;
	return -5.0l * 9.96741762 * pow(Sg, 4.0l) + 4.0l * 16.55870922 * pow(Sg, 3.0l) - 3.0l * 4.26589291 * pow(Sg, 2.0l) - 2.0l * 6.86968292 * Sg + 5.04975874;
}

double calc_oil_kro_Swo(double Sw, double Sg) {
	return ((calc_water_krow_Sw(Sw) + calc_water_krw_Sw(Sw)) * (calc_gas_krog_Swo(Sg) + calc_gas_krg_Swo(Sg))
		- (calc_water_krw_Sw(Sw) + calc_gas_krg_Swo(Sg)));
}

double calc_der_oil_kro_Sw(double Sw, double Sg) {
	return ((calc_der_water_krow_Sw(Sw) + calc_der_water_krw_Sw(Sw)) * (calc_gas_krog_Swo(Sg) + calc_gas_krg_Swo(Sg))
		+ (calc_water_krow_Sw(Sw) + calc_water_krw_Sw(Sw)) * (calc_der_gas_krog_Swo(Sg) + calc_der_gas_krg_Swo(Sg))
		- (calc_der_water_krw_Sw(Sw) + calc_der_gas_krg_Swo(Sg)));
}

double calc_der_oil_kro_So(double Sw, double Sg) {
	return ((calc_water_krow_Sw(Sw) + calc_water_krw_Sw(Sw)) * (calc_der_gas_krog_Swo(Sg) + calc_der_gas_krg_Swo(Sg))
		- (calc_der_gas_krg_Swo(Sg)));
}

double calc_oil_kro_Swo(double krwmax, double Sw, double So, double Sgr, double Swc, double Sor, double nw, double ng, double now, double nog) {
	if (So <= Sor) return 0.00;
	else {
		double Sno = (So - Sor) / (1.0 - Swc - Sor);
		double Snw = (Sw - Swc) / (1.0 - Swc - Sor);
		double Sng = (1.0 - So - Sw) / (1.0 - Swc - Sor);
		double betaw = calc_water_krow_Sw(Sw, Swc, Sor, now) / (1.0 - Snw);
		double betag = calc_gas_krog_Swo(Sw, So, Swc, Sor, Sgr, nog) / (1.0 - Sng);
		return Sno * betag * betaw;
	}
}

double calc_der_oil_kro_Sw(double krwmax, double Sw, double So, double Sgr, double Swc, double Sor, double nw, double ng, double now, double nog) {
	if (So <= Sor) return 0.0000001;
	else {
		double Sno = (So - Sor) / (1.0 - Swc - Sor);
		double Snw = (Sw - Swc) / (1.0 - Swc - Sor);
		double Sng = (1.0 - So - Sw) / (1.0 - Swc - Sor);
		double betaw = calc_water_krow_Sw(Sw, Swc, Sor, now) / (1.0 - Snw);
		double betag = calc_gas_krog_Swo(Sw, So, Swc, Sor, Sgr, nog) / (1.0 - Sng);
		double dbetawdSw = (calc_der_water_krow_Sw(Sw, Swc, Sor, now) / (1.0 - Snw)) + (1.0 / ((1.0 - Swc - Sor) * (1.0 - Snw) * (1.0 - Snw)));
		double dbetagdSw = (calc_der_gas_krog_Swo(Sw, So, Swc, Sor, Sgr, nog) / (1.0 - Sng)) - (1.0 / ((1.0 - Swc - Sor) * (1.0 - Sng) * (1.0 - Sng)));
		return Sno * (betag * dbetawdSw + betaw * dbetagdSw);
	}
}

double calc_der_oil_kro_So(double krwmax, double Sw, double So, double Sgr, double Swc, double Sor, double nw, double ng, double now, double nog) {
	if (So <= Sor) return 0.0000001;
	else {
		double Sno = (So - Sor) / (1.0 - Swc - Sor);
		double Snw = (Sw - Swc) / (1.0 - Swc - Sor);
		double Sng = (1.0 - So - Sw) / (1.0 - Swc - Sor);
		double betaw = calc_water_krow_Sw(Sw, Swc, Sor, now) / (1.0 - Snw);
		double betag = calc_gas_krog_Swo(Sw, So, Swc, Sor, Sgr, nog) / (1.0 - Sng);
		double dbetagdSo = (calc_der_gas_krog_Swo(Sw, So, Swc, Sor, Sgr, nog) / (1.0 - Sng)) - (1.0 / ((1.0 - Swc - Sor) * (1.0 - Sng) * (1.0 - Sng)));
		return betaw * (betag / (1.0 - Swc - Sor) + Sno * dbetagdSo);
	}
}

//Press�o(psia), Temperatura(~F}, Rs(SCF/STB)
double calc_gas_Rs_p(double yg, double p, double pb, double API, double T) {
	if (p <= pb) {
		return (yg * pow(pow(10.l, 0.0125l * API - 0.00091l * T) / 18.2l, 1.204819l) * pow(p + 1.4l * 18.2l, 1.204819l));
	}

	else {
		return (yg * pow(pow(10.l, 0.0125l * API - 0.00091l * T) / 18.2l, 1.204819l) * pow(pb + 1.4l * 18.2l, 1.204819l));
	}
}

double calc_der_gas_Rs_p(double yg, double p, double API, double T) {

	return (1.2048 * yg * pow(pow(10., 0.0125 * API - 0.00091 * T) / 18.4, 1.2048) * pow(p + 1.4 * 18.4, 0.2048));

}


//Press�o(psia), Temperatura(~F}, Rs(SCF/STB)
double calc_der_gas_Rs_p(double yg, double p, double pb, double API, double T) {
	if (p <= pb) {
		return (1.2048 * yg * pow(pow(10., 0.0125 * API - 0.00091 * T) / 18.4, 1.2048) * pow(p + 1.4 * 18.4, 0.2048));
	}

	else {
		return 3e-6 * (yg * pow(pow(10.l, 0.0125l * API - 0.00091l * T) / 18.2l, 1.204819l) * pow(pb + 1.4l * 18.2l, 1.204819l));
	}
}

double calc_der_gas_Rs_pb(double yg, double p, double pb, double API, double T) {
	return (1.2048 * yg * pow(pow(10., 0.0125 * API - 0.00091 * T) / 18.4, 1.2048) * pow(p + 1.4 * 18.4, 0.2048)) * (1.l + 3e-6 * (p - pb)) - 3e-6 * (yg * pow(pow(10.l, 0.0125l * API - 0.00091l * T) / 18.2l, 1.204819l) * pow(pb + 1.4l * 18.2l, 1.204819l));
}

//Press�o(psia), Temperatura(~F}, Bo(Rb/STB)
double calc_oil_Bo_p(double yo, double yg, double p, double pb, double API, double T) {
	return (0.9759 + 0.000147 * pow((calc_gas_Rs_p(yg, p, pb, API, T) * pow((yg / yo), 0.5) + 1.25 * T), 1.175));
}

//Press�o(psia), Temperatura(~F}, Bo(Rb/STB)
double calc_der_oil_Bo_p(double yo, double yg, double p, double API, double T) {
	return (0.000147 * 1.175 * pow((yg / yo), 0.5) * pow((calc_gas_Rs_p(yg, p, p, API, T) * pow((yg / yo), 0.5) + 1.25 * T), 0.175) * calc_der_gas_Rs_p(yg, p, API, T));
}

//Press�o(psia), Temperatura(~F}, mi(cp)
double calc_oil_mi_p(double Rs, double T, double API, double p, double pb) {
	double parametrosmi[6];
	//A
	parametrosmi[0] = pow(10.0l, (0.43l + (8.33l / API)));//A
														  //a
	parametrosmi[1] = Rs * (2.2l * pow(10.l, (-7.0l)) * Rs - 7.4 * pow(10.l, (-4.0l)));//a
																			   //c
	parametrosmi[2] = 8.62l * pow(10.l, (-5.0l)) * Rs;//c
												  //d
	parametrosmi[3] = 1.1l * pow(10.l, (-3.0l)) * Rs;//d
												 //e
	parametrosmi[4] = 3.74l * pow(10.l, (-3.0l)) * Rs;//e
												  //b
	parametrosmi[5] = (0.68l / pow(10.0l, parametrosmi[2])) + (0.25l / pow(10.0l, parametrosmi[3])) + (0.062l / pow(10.0l, parametrosmi[4]));//b
	double miod = (0.32l + (1.8l * pow(10.0l, 7.0l) / pow(API, 4.53l))) * pow((360 / (T + 200)), parametrosmi[0]);
	if (p < pb) {
		return pow(10.0, parametrosmi[1]) * pow(miod, parametrosmi[5]);
	}
	else {
		double miob = pow(10.0, parametrosmi[1]) * pow(miod, parametrosmi[5]);
		return miob + 0.001l * (p - pb) * (0.024l * pow(miob, 1.6l) + 0.38l * pow(miob, 0.56l));
	}
}

//Press�o(psia), Temperatura(~F}, mi(cp)
double calc_der_oil_mi_p(double Rs, double dRs, double T, double API, double p, double pb) {
	double parametrosmi[6];
	parametrosmi[0] = pow(10.0l, (0.43l + (8.33l / API)));//A
														  //a
	parametrosmi[1] = Rs * (2.2l * pow(10.l, (-7.0l)) * Rs - 7.4 * pow(10.l, (-4.0l)));//a
																			   //c
	parametrosmi[2] = 8.62l * pow(10.l, (-5.0l)) * Rs;//c
												  //d
	parametrosmi[3] = 1.1l * pow(10.l, (-3.0l)) * Rs;//d
												 //e
	parametrosmi[4] = 3.74l * pow(10.l, (-3.0l)) * Rs;//e
												  //b
	parametrosmi[5] = (0.68l / pow(10.0l, parametrosmi[2])) + (0.25l / pow(10.0l, parametrosmi[3])) + (0.062l / pow(10.0l, parametrosmi[4]));//b
	double miod = (0.32l + (1.8l * pow(10.0, 7.0) / pow(API, 4.53l))) * pow((360 / (T + 200)), parametrosmi[0]);
	if (p < pb) {
		double db = (-0.0001349683278109390 / pow(10.0, parametrosmi[2]) - 0.0006332109005733627 / pow(10.0, parametrosmi[3]) - 0.0005339234313634593 / pow(10.0, parametrosmi[4])) * dRs;
		double da = dRs * (0.00000044 * Rs - 0.00074);
		return 2.30258509299 * pow(miod, parametrosmi[5]) * pow(10.0, parametrosmi[1]) * da + pow(10.0, parametrosmi[1]) * pow(miod, parametrosmi[5]) * log(miod) * db;
	}
	else {
		double miob = pow(10.0, parametrosmi[1]) * pow(miod, parametrosmi[5]);
		return 0.001l * (0.024l * pow(miob, 1.6l) + 0.38l * pow(miob, 0.56l));
	}
}

double calc_der_oil_mi_pb(double Rs, double dRs, double T, double API, double p, double pb) {
	double parametrosmi[6];
	parametrosmi[0] = pow(10.0l, (0.43l + (8.33l / API)));//A
														  //a
	parametrosmi[1] = Rs * (2.2l * pow(10.l, (-7.0l)) * Rs - 7.4 * pow(10.l, (-4.0l)));//a
																			   //c
	parametrosmi[2] = 8.62l * pow(10.l, (-5.0l)) * Rs;//c
												  //d
	parametrosmi[3] = 1.1l * pow(10.l, (-3.0l)) * Rs;//d
												 //e
	parametrosmi[4] = 3.74l * pow(10.l, (-3.0l)) * Rs;//e
												  //b
	parametrosmi[5] = (0.68l / pow(10.0l, parametrosmi[2])) + (0.25l / pow(10.0l, parametrosmi[3])) + (0.062l / pow(10.0l, parametrosmi[4]));//b
	double miod = (0.32l + (1.8l * pow(10.0, 7.0) / pow(API, 4.53l))) * pow((360 / (T + 200)), parametrosmi[0]);
	if (p > pb) {
		double db = (-0.0001349683278109390 / pow(10.0, parametrosmi[2]) - 0.0006332109005733627 / pow(10.0, parametrosmi[3]) - 0.0005339234313634593 / pow(10.0, parametrosmi[4])) * dRs;
		double da = dRs * (0.00000044 * Rs - 0.00074);
		double der_miob = 2.30258509299 * pow(miod, parametrosmi[5]) * pow(10.0, parametrosmi[1]) * da + pow(10.0, parametrosmi[1]) * pow(miod, parametrosmi[5]) * log(miod) * db;
		double miob = pow(10.0, parametrosmi[1]) * pow(miod, parametrosmi[5]);
		return der_miob - 0.001l * (0.024l * pow(miob, 1.6l) + 0.38l * pow(miob, 0.56l)) + 0.001l * (p - pb) * der_miob * (0.024l * 1.6l * pow(miob, 0.6l) + 0.38l * 0.56l * pow(miob, -0.44l));
	}
	else {
		return 0.0l;
	}
}

//
double calc_oil_rho_p(double Rs, double Bo, double rho_os, double rho_gs) {
	return ((rho_os / Bo) + (rho_gs * Rs / Bo));
}


double calc_der_oil_rho_p(double Rs, double dRs, double Bo, double dBo, double rho_os, double rho_gs) {
	return (((rho_os + (rho_gs * Rs)) * (-dBo / (Bo * Bo))) + (rho_gs * (dRs / Bo)));
}

//
void calc_const_mig(double* parametrosmig, double yg, double yn2, double yco2, double yh2s, double Tf) {
	//A0 to A3
	parametrosmig[0] = -2.4621182l;
	parametrosmig[1] = 2.97054714l;
	parametrosmig[2] = -0.286264054l;
	parametrosmig[3] = 8.05420522l * pow(10.l, (-3.l));
	//B0 to B3
	parametrosmig[4] = 2.80860949l;
	parametrosmig[5] = -3.49803305l;
	parametrosmig[6] = 0.36037302l;
	parametrosmig[7] = -1.04432413l * pow(10.l, (-2.l));
	//C0 to C3
	parametrosmig[8] = -0.793385684l;
	parametrosmig[9] = 1.39643306l;
	parametrosmig[10] = -0.149144925l;
	parametrosmig[11] = 4.41015512l * pow(10.l, (-3.l));
	//D0 to D3
	parametrosmig[12] = 0.0839387178l;
	parametrosmig[13] = -0.186408848l;
	parametrosmig[14] = 0.0203367881l;
	parametrosmig[15] = -6.095799263l * pow(10.l, (-4.l));
	//mi_c
	parametrosmig[16] = (1.709l * pow(10.l, (-5.l)) - 2.062l * pow(10.l, (-6.l)) * yg) * Tf
		+ pow(10.l, (-3.l)) * (8.188l - 6.15l * log(yg)
			+ yn2 * (9.59l + 8.48l * log(yg)) + yco2 * (6.24l + 9.08l * log(yg))
			+ yh2s * (3.73l + 8.49l * log(yg)));
}

//Tpc(R), p(psia)
void calc_pseudocrit(double* ppc, double* Tpc, double yg, double yco2, double yh2s) {
	*ppc = 706.l - 51.7l * yg - 11.1l * yg * yg;
	*Tpc = 187.l + 330.l * yg - 71.5l * yg * yg;
	double WA = 120.0l * (pow((yco2 + yh2s), 0.9l) - pow((yco2 + yh2s), 1.6l)) - 15.0 * (pow(yh2s, 0.5l) - pow(yh2s, 4.0l));
	*ppc = *ppc * ((*Tpc) - WA) / ((*Tpc) + yh2s * (1.0l - yh2s) * WA);
	*Tpc = *Tpc - WA;
}

//Sa�da cp
double calc_gas_mig_p(double* parametrosmig, double Tred, double pred) {

	//A,B,C e D
	parametrosmig[17] = parametrosmig[0] + parametrosmig[1] * pred + parametrosmig[2] * pred * pred + parametrosmig[3] * pred * pred * pred;
	parametrosmig[18] = parametrosmig[4] + parametrosmig[5] * pred + parametrosmig[6] * pred * pred + parametrosmig[7] * pred * pred * pred;
	parametrosmig[19] = parametrosmig[8] + parametrosmig[9] * pred + parametrosmig[10] * pred * pred + parametrosmig[11] * pred * pred * pred;
	parametrosmig[20] = parametrosmig[12] + parametrosmig[13] * pred + parametrosmig[14] * pred * pred + parametrosmig[15] * pred * pred * pred;

	return exp(parametrosmig[17] + parametrosmig[18] * Tred + parametrosmig[19] * Tred * Tred + parametrosmig[20] * Tred * Tred * Tred) * parametrosmig[16] / Tred;
}

//Sa�da cp
double calc_der_gas_mig_p(double* parametrosmig, double Tred, double pred, double ppc, double Tpc) {
	//A,B,C e D
	parametrosmig[17] = parametrosmig[0] + parametrosmig[1] * pred + parametrosmig[2] * pred * pred + parametrosmig[3] * pred * pred * pred;
	parametrosmig[18] = parametrosmig[4] + parametrosmig[5] * pred + parametrosmig[6] * pred * pred + parametrosmig[7] * pred * pred * pred;
	parametrosmig[19] = parametrosmig[8] + parametrosmig[9] * pred + parametrosmig[10] * pred * pred + parametrosmig[11] * pred * pred * pred;
	parametrosmig[20] = parametrosmig[12] + parametrosmig[13] * pred + parametrosmig[14] * pred * pred + parametrosmig[15] * pred * pred * pred;

	//dA,dB,dC e dD /dp
	parametrosmig[21] = (parametrosmig[1] / ppc) + 2 * (parametrosmig[2] / ppc) * pred + 3 * (parametrosmig[3] / ppc) * pred * pred;
	parametrosmig[22] = (parametrosmig[5] / ppc) + 2 * (parametrosmig[6] / ppc) * pred + 3 * (parametrosmig[7] / ppc) * pred * pred;
	parametrosmig[23] = (parametrosmig[9] / ppc) + 2 * (parametrosmig[10] / ppc) * pred + 3 * (parametrosmig[11] / ppc) * pred * pred;
	parametrosmig[24] = (parametrosmig[13] / ppc) + 2 * (parametrosmig[14] / ppc) * pred + 3 * (parametrosmig[15] / ppc) * pred * pred;

	return exp(parametrosmig[17] + parametrosmig[18] * Tred + parametrosmig[19] * Tred * Tred + parametrosmig[20] * Tred * Tred * Tred) * (parametrosmig[21] + parametrosmig[22] * Tred + parametrosmig[23] * Tred * Tred + parametrosmig[24] * Tred * Tred * Tred) * parametrosmig[16] / Tred;
}

//Sa�da sem unidade
double calc_gas_Z_p(double pred, double Tred, double* parametrosZ, double pgr) {
	parametrosZ[0] = 0.06423l;
	parametrosZ[1] = 0.5353l * Tred - 0.6123l;
	parametrosZ[2] = 0.3151l * Tred - 1.0467l - (0.5783l / (Tred * Tred));
	parametrosZ[3] = Tred;
	parametrosZ[4] = 0.6816l / (Tred * Tred);
	parametrosZ[5] = 0.6845l;
	parametrosZ[6] = 0.27l * pred;
	pgr = parametrosZ[6] / Tred;
	double F;
	double dF;
	double erro = 1;
	while (erro > 0.00000005) {
		F = parametrosZ[0] * (pgr * pgr * pgr * pgr * pgr * pgr) + parametrosZ[1] * (pgr * pgr * pgr)
			+ parametrosZ[2] * (pgr * pgr) + parametrosZ[3] * pgr + parametrosZ[4] * (pgr * pgr * pgr) * (1 + parametrosZ[5] * pgr * pgr) * exp(-parametrosZ[5] * pgr * pgr) - parametrosZ[6];
		dF = 6 * parametrosZ[0] * (pgr * pgr * pgr * pgr * pgr) + 3 * parametrosZ[1] * (pgr * pgr)
			+ 2 * parametrosZ[2] * (pgr)+parametrosZ[3] + parametrosZ[4] * (pgr * pgr) * (3 + parametrosZ[5] * pgr * pgr * (3 - 2 * parametrosZ[5] * pgr * pgr)) * exp(-parametrosZ[5] * pgr * pgr);
		pgr = pgr - F / dF;
		erro = abs(F / dF);
	}
	return (parametrosZ[6] / (Tred * pgr));
}

//Sa�da bbl/SCF
double calc_gas_Bg_p(double p) {

	return 2.314078032892617l * pow(p, -0.95528532l);

}

//Sa�da bbl/SCF
double calc_der_gas_Bg_p(double p) {

	return -0.95528532 * 2.314078032892617 * pow(p, -1.95528532);
}

double calc_gas_rho_p(double rho_gs, double Bg) {
	return rho_gs / Bg;
}

double calc_der_gas_rho_p(double rho_gs, double Bg, double der_Bg) {
	return -rho_gs * der_Bg / (Bg * Bg);
}

void calc_props_grid_ant(int N, double* po_0, double* pb_0, double* Sw_0, double* So_0, double* Sg_0, double* X_pb_0,
	double* phi_0, double* Bw_0, double* Rs_0, double* Bo_0, double* Bg_0,
	double* phi_ref, double Cr, double pi, double T, double* parametrosBw,
	double API, double yo, double yg) {

	int p;
#pragma omp parallel private(p) 
	{
#pragma omp for schedule(dynamic) nowait 
		for (p = 0; p <= (N - 1); p++)
		{

			if (po_0[p] > pb_0[p]) X_pb_0[p] = 1.0l;
			else X_pb_0[p] = 0.0l;

			phi_0[p] = calc_phi_p(phi_ref[p], Cr, po_0[p], pi);
			Bw_0[p] = calc_water_Bw_p(po_0[p], T, parametrosBw);
			Bo_0[p] = calc_oil_Bo_p(yo, yg, po_0[p], pb_0[p], API, T);
			Bg_0[p] = bbl_ft3(calc_gas_Bg_p(po_0[p]));
			Rs_0[p] = ft3_bbl(calc_gas_Rs_p(yg, po_0[p], pb_0[p], API, T));

		}
	}
}


void calc_props_grid(int N, double* po, double* pg, double* pw, double* pb, double* Sw, double* So, double* Sg, double* X_pb,
	double* phi, double* der_phi,
	double* rho_w, double* der_rho_w, double* mi_w, double* der_mi_w, double* Bw, double* der_Bw, double* pcw, double* der_pcw, double* krw, double* der_krw,
	double* Rs, double* der_Rs_p, double* der_Rs_pb, double* kro, double* der_kro_Sw, double* der_kro_So, double* Bo, double* der_Bo_p, double* der_Bo_pb,
	double* mi_o, double* der_mi_o_p, double* der_mi_o_pb, double* rho_o, double* der_rho_o_p, double* der_rho_o_pb,
	double* pcg, double* der_pcg_Sw, double* der_pcg_So, double* krg, double* der_krg, double* mi_g, double* der_mi_g, double* Bg, double* der_Bg, double* rho_g, double* der_rho_g,
	double* phi_ref, double Cr, double Bwi, double rho_ws, double pi, double T, double* parametrosBw,
	double Bpcw, double Bpcg, double API, double yo, double yg, double rho_os, double rho_gs,
	double* parametrosmig, double ppc, double Tpc, double pcwmin, double pcwmax,
	double* lambda1_W, double* d_lambda1_po_W, double* d_lambda1_SW_W,
	double* lambda2_W, double* d_lambda2_po_W, double* d_lambda2_SW_W,
	double* lambda1_O, double* d_lambda1_po_O, double* d_lambda1_SW_O, double* d_lambda1_pb_O, double* d_lambda1_SO_O,
	double* lambda2_O, double* d_lambda2_po_O, double* d_lambda2_SW_O, double* d_lambda2_pb_O, double* d_lambda2_SO_O,
	double* lambda1_G, double* d_lambda1_po_G, double* d_lambda1_SW_G, double* d_lambda1_pb_G, double* d_lambda1_SO_G,
	double* lambda2_G, double* d_lambda2_po_G, double* d_lambda2_SW_G, double* d_lambda2_pb_G, double* d_lambda2_SO_G,
	double* lambda3_G, double* d_lambda3_po_G, double* d_lambda3_SW_G, double* d_lambda3_pb_G, double* d_lambda3_SO_G,
	double* lambda4_G, double* d_lambda4_po_G, double* d_lambda4_SW_G, double* d_lambda4_pb_G, double* d_lambda4_SO_G) {

	int p;
	double Tred = T / Tpc;
#pragma omp parallel private(p) 
	{
#pragma omp for schedule(dynamic) nowait 
		for (p = 0; p <= (N - 1); p++)
		{

			if (po[p] > pb[p]) X_pb[p] = 1.0l;
			else X_pb[p] = 0.0l;

			phi[p] = calc_phi_p(phi_ref[p], Cr, po[p], pi);
			der_phi[p] = calc_der_phi_p(phi_ref[p], Cr, po[p], pi);

			rho_w[p] = calc_water_rho_p(Bwi, rho_ws, T, po[p], pi);
			der_rho_w[p] = calc_der_water_rho_p(Bwi, rho_ws, T, po[p], pi);

			mi_w[p] = calc_water_mi_p(T, po[p]);
			der_mi_w[p] = calc_der_water_mi_p(T, po[p]);

			Bw[p] = calc_water_Bw_p(po[p], T, parametrosBw);
			der_Bw[p] = calc_der_water_Bw_p(po[p], T, parametrosBw);

			pcw[p] = 0.0l;// calc_water_pcw_Sw(Sw[p], Swc, epsilon, pcwmin, pcwmax, Bpcw);
			der_pcw[p] = 0.0l;//calc_der_water_pcw_Sw(Sw[p], Swc, epsilon, Bpcw);

			pcg[p] = 0.0l;//calc_gas_pcg(Sw[p], So[p], Swc, Sor, epsilon, pcgmin, pcgmax, Bpcg);
			der_pcg_Sw[p] = 0.0l;//calc_der_gas_pcg_Swo(Sw[p], So[p], Swc, Sor, epsilon, Bpcg);
			der_pcg_So[p] = 0.0l;//calc_der_gas_pcg_Swo(Sw[p], So[p], Swc, Sor, epsilon, Bpcg);

			pw[p] = po[p];
			pg[p] = po[p];

			krw[p] = calc_water_krw_Sw(Sw[p]);
			der_krw[p] = calc_der_water_krw_Sw(Sw[p]);

			krg[p] = calc_gas_krg_Swo(Sg[p]);
			der_krg[p] = calc_der_gas_krg_Swo(Sg[p]);

			kro[p] = calc_oil_kro_Swo(Sg[p], Sw[p]);
			der_kro_Sw[p] = calc_der_oil_kro_Sw(Sg[p], Sw[p]);
			der_kro_So[p] = calc_der_oil_kro_So(Sg[p], Sw[p]);

			Rs[p] = ft3_bbl(calc_gas_Rs_p(yg, po[p], pb[p], API, T));
			der_Rs_p[p] = (1 - X_pb[p]) * ft3_bbl(calc_der_gas_Rs_p(yg, po[p], API, T));
			der_Rs_pb[p] = X_pb[p] * ft3_bbl(calc_der_gas_Rs_p(yg, pb[p], API, T));

			Bo[p] = calc_oil_Bo_p(yo, yg, po[p], pb[p], API, T);
			der_Bo_p[p] = (1.0l - X_pb[p]) * calc_der_oil_Bo_p(yo, yg, po[p], API, T);
			der_Bo_pb[p] = X_pb[p] * calc_der_oil_Bo_p(yo, yg, pb[p], API, T);

			mi_o[p] = calc_oil_mi_p(bbl_ft3(Rs[p]), T, API, po[p], pb[p]);
			der_mi_o_p[p] = calc_der_oil_mi_p(bbl_ft3(Rs[p]), bbl_ft3(der_Rs_p[p]), T, API, po[p], pb[p]);
			der_mi_o_pb[p] = calc_der_oil_mi_pb(bbl_ft3(Rs[p]), bbl_ft3(der_Rs_pb[p]), T, API, po[p], pb[p]);

			rho_o[p] = calc_oil_rho_p(Rs[p], Bo[p], rho_os, rho_gs);
			der_rho_o_p[p] = calc_der_oil_rho_p(Rs[p], der_Rs_p[p], Bo[p], der_Bo_p[p], rho_os, rho_gs);
			der_rho_o_pb[p] = calc_der_oil_rho_p(Rs[p], der_Rs_pb[p], Bo[p], der_Bo_pb[p], rho_os, rho_gs);

			mi_g[p] = calc_gas_mig_p(parametrosmig, Tred, po[p] / ppc);
			der_mi_g[p] = calc_der_gas_mig_p(parametrosmig, Tred, po[p] / ppc, ppc, Tpc);

			Bg[p] = bbl_ft3(calc_gas_Bg_p(po[p]));
			der_Bg[p] = bbl_ft3(calc_der_gas_Bg_p(po[p]));

			rho_g[p] = calc_gas_rho_p(rho_gs, Bg[p]);
			der_rho_g[p] = calc_der_gas_rho_p(rho_gs, Bg[p], der_Bg[p]);

			lambda1_W[p] = krw[p] / (mi_w[p] * Bw[p]);
			d_lambda1_po_W[p] = -(der_mi_w[p] / mi_w[p] + der_Bw[p] / Bw[p]) * lambda1_W[p];
			d_lambda1_SW_W[p] = der_krw[p] / (mi_w[p] * Bw[p]);

			lambda2_W[p] = rho_w[p] * krw[p] / (mi_w[p] * Bw[p]);
			d_lambda2_po_W[p] = der_rho_w[p] * lambda1_W[p] + rho_w[p] * d_lambda1_po_W[p];
			d_lambda2_SW_W[p] = rho_w[p] * d_lambda1_SW_W[p];

			lambda1_O[p] = kro[p] / (mi_o[p] * Bo[p]);
			d_lambda1_po_O[p] = -(der_mi_o_p[p] / mi_o[p] + der_Bo_p[p] / Bo[p]) * lambda1_O[p];
			d_lambda1_SW_O[p] = der_kro_Sw[p] / (mi_o[p] * Bo[p]);
			d_lambda1_pb_O[p] = -(der_mi_o_pb[p] / mi_o[p] + der_Bo_pb[p] / Bo[p]) * lambda1_O[p];
			d_lambda1_SO_O[p] = der_kro_So[p] / (mi_o[p] * Bo[p]);

			lambda2_O[p] = rho_o[p] * kro[p] / (mi_o[p] * Bo[p]);
			d_lambda2_po_O[p] = der_rho_o_p[p] * lambda1_O[p] + rho_o[p] * d_lambda1_po_O[p];
			d_lambda2_SW_O[p] = rho_o[p] * d_lambda1_SW_O[p];
			d_lambda2_pb_O[p] = der_rho_o_pb[p] * lambda1_O[p] + rho_o[p] * d_lambda1_pb_O[p];
			d_lambda2_SO_O[p] = rho_o[p] * d_lambda1_SO_O[p];

			lambda1_G[p] = Rs[p] * kro[p] / (mi_o[p] * Bo[p]);
			d_lambda1_po_G[p] = (der_Rs_p[p] / Rs[p] - der_mi_o_p[p] / mi_o[p] - der_Bo_p[p] / Bo[p]) * lambda1_G[p];
			d_lambda1_SW_G[p] = Rs[p] * der_kro_Sw[p] / (mi_o[p] * Bo[p]);
			d_lambda1_pb_G[p] = (der_Rs_pb[p] / Rs[p] - der_mi_o_pb[p] / mi_o[p] - der_Bo_pb[p] / Bo[p]) * lambda1_G[p];
			d_lambda1_SO_G[p] = Rs[p] * der_kro_So[p] / (mi_o[p] * Bo[p]);

			lambda2_G[p] = Rs[p] * rho_o[p] * kro[p] / (mi_o[p] * Bo[p]);
			d_lambda2_po_G[p] = der_rho_o_p[p] * lambda1_G[p] + rho_o[p] * d_lambda1_po_G[p];
			d_lambda2_SW_G[p] = rho_o[p] * d_lambda1_SW_G[p];
			d_lambda2_pb_G[p] = der_rho_o_pb[p] * lambda1_G[p] + rho_o[p] * d_lambda1_pb_G[p];
			d_lambda2_SO_G[p] = rho_o[p] * d_lambda1_SO_G[p];

			lambda3_G[p] = krg[p] / (mi_g[p] * Bg[p]);
			d_lambda3_po_G[p] = -(der_mi_g[p] / mi_g[p] + der_Bg[p] / Bg[p]) * lambda3_G[p];
			d_lambda3_SW_G[p] = der_krg[p] / (mi_g[p] * Bg[p]);
			d_lambda3_pb_G[p] = 0.0l;
			d_lambda3_SO_G[p] = der_krg[p] / (mi_g[p] * Bg[p]);

			lambda4_G[p] = rho_g[p] * krg[p] / (mi_g[p] * Bg[p]);
			d_lambda4_po_G[p] = der_rho_g[p] * lambda3_G[p] + rho_g[p] * d_lambda3_po_G[p];
			d_lambda4_SW_G[p] = rho_g[p] * d_lambda3_SW_G[p];
			d_lambda4_pb_G[p] = 0.0l;
			d_lambda4_SO_G[p] = rho_g[p] * d_lambda3_SO_G[p];

		}
	}


}

double perm_harm(double Kp, double Knm) {

	if (Kp < 0.000000001 || Knm < 0.0000000001) return 0.0;

	return 2.0l / ((1.0l / Kp) + (1.0l / Knm));
}

//https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
double generateGaussianNoise(double mu, double sigma)
{
	const double eps = 0.00000000000005;
	const double two_pi = 2.0 * 3.14159265358979323846;

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

void get_perm(double* k, int dimx, int dimy, int dimz)
{
	char nome_in[] = "permeabilidade.txt";
	int pos = 0;
	ifstream fpin(nome_in);
	for (int z = 0; z < dimz; z++) for (int y = 0; y < dimy; y++) for (int x = 0; x < dimx; x++) fpin >> k[pos], pos++;
}

void gera_permeabilidades(double* k, int dimx, int dimy, int dimz, double perm, double varperm)
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
	for (int z = 0; z < dimz; z++) for (int y = 0; y < dimy; y++) for (int x = 0; x < dimx; x++) k[pos] = generateGaussianNoise(perm, varperm), fper << k[pos] << endl, fper2 << k[pos] << endl, pos++;
	fper << "];" << endl;
	fper << "k=reshape(k," << dimx << "," << dimy << "," << dimz << ");";
	fper.close();
}

void grava_var(double* p, int fim, int inicio, int res_step, int ite_coup, char directory[], __int8 vartoprint)
{
	char nomep[200];
	char nomep2[20];
	sprintf_s(nomep, directory);

	if (vartoprint == 1) sprintf_s(nomep2, "\\p_%d_%d.csv", res_step, ite_coup);
	else if (vartoprint == 2) sprintf_s(nomep2, "\\pb_%d_%d.csv", res_step, ite_coup);
	else if (vartoprint == 3) sprintf_s(nomep2, "\\Sw_%d_%d.csv", res_step, ite_coup);
	else if (vartoprint == 4) sprintf_s(nomep2, "\\So_%d_%d.csv", res_step, ite_coup);
	else if (vartoprint == 5) sprintf_s(nomep2, "\\PWB_%d_%d.csv", res_step, ite_coup);
	else if (vartoprint == 6) sprintf_s(nomep2, "\\Qo_%d_%d.csv", res_step, ite_coup);
	else if (vartoprint == 7) sprintf_s(nomep2, "\\Qg_%d_%d.csv", res_step, ite_coup);

	strcat_s(nomep, nomep2);

	ofstream fp(nomep);
	fp.precision(14);

	if (vartoprint == 1) fp << "p" << endl;
	else if (vartoprint == 2) fp << "pb" << endl;
	else if (vartoprint == 3) fp << "Sw" << endl;
	else if (vartoprint == 4) fp << "So" << endl;
	else if (vartoprint == 5) fp << "PWB" << endl;
	else if (vartoprint == 6) fp << "Qo" << endl;
	else if (vartoprint == 7) fp << "Qg" << endl;

	for (int i = inicio; i < fim; i++)
	{
		fp << p[i] << "\t";
		fp << endl;
	}

	fp.close();
}

//
//void grava_Qo(double* p, int fim, int inicio, int j, int controlador, int zy)
//{
//	char nomep[100];
//
//
//	sprintf_s(nomep, "Qo_%d_%d_%d.m", j, zy, controlador);
//
//	ofstream fp(nomep);
//	fp.precision(14);
//
//	fp << "Qo = [ " << endl;
//	for (int i = inicio; i < fim; i++)
//	{
//		fp << p[i];
//		fp << endl;
//	}
//
//	//cout << "\n\nGravando pressao." << endl;
//
//	fp << "];" << endl;
//
//	fp.close();
//}

//functions for debug
void grava_ia(int* p, int fim, int inicio, int j)
{
	char nomep[50];


	sprintf_s(nomep, "ia_%d.m", j);

	ofstream fp(nomep);
	fp.precision(14);

	fp << "ia = [ " << endl;
	for (int i = inicio; i < fim; i++)
	{
		fp << p[i];
		fp << endl;
	}

	cout << "\n\nGravando IA." << endl;

	fp << "];" << endl;

	fp.close();
}

void grava_a(double* p, int fim, int inicio, int j)
{
	char nomep[50];


	sprintf_s(nomep, "a_%d.m", j);

	ofstream fp(nomep);
	fp.precision(14);

	fp << "a = [ " << endl;
	for (int i = inicio; i < fim; i++)
	{
		fp << p[i];
		fp << endl;
	}


	cout << "\n\nGravando A." << endl;

	fp << "];" << endl;

	fp.close();
}

void grava_b(double* p, int fim, int inicio, int j)
{
	char nomep[50];


	sprintf_s(nomep, "b_%d.m", j);

	ofstream fp(nomep);
	fp.precision(14);

	fp << "b = [ " << endl;
	for (int i = inicio; i < fim; i++)
	{
		fp << p[i];
		fp << endl;
	}


	cout << "\n\nGravando B." << endl;

	fp << "];" << endl;

	fp.close();
}

void grava_e(double* p, int fim, int inicio, int j)
{
	char nomep[50];


	sprintf_s(nomep, "e_%d.m", j);

	ofstream fp(nomep);
	fp.precision(14);

	fp << "e = [ " << endl;
	for (int i = inicio; i < fim; i++)
	{
		fp << p[i];
		fp << endl;
	}


	cout << "\n\nGravando E." << endl;

	fp << "];" << endl;

	fp.close();
}

void grava_ja(int* p, int fim, int inicio, int j)
{
	char nomep[50];


	sprintf_s(nomep, "ja_%d.m", j);

	ofstream fp(nomep);
	fp.precision(14);

	fp << "ja = [ " << endl;
	for (int i = inicio; i < fim; i++)
	{
		fp << p[i];
		fp << endl;
	}


	cout << "\n\nGravando JA." << endl;

	fp << "];" << endl;

	fp.close();
}

