#include "stdafx.h"
#include "Well.h"


Well::Well()
{
}

Well::Well(double L_in, int nz_in, double dT_in, double P_in_in, double T_in_in, double OFR_in,
	double GOR_in, double API_in, double SG_G_in, double WC_in, bool *error_SL_CH, double P_initial_in, double rug_in)
{
	L = L_in;
	nz = nz_in;
	dz = L / (nz - 1);
	dT = dT_in;
	P_in = P_in_in;
	T_in = T_in_in;
	OFR = OFR_in;
	GOR = GOR_in;
	API = API_in;
	SG_G = SG_G_in;
	rug = rug_in;
	WC = WC_in;
	P_initial = P_initial_in;
	calc_pseudocrit(&ppc, &Tpc, SG_G, y_co2, y_h2s);	

	double parametersBw[9];
	calc_const_Bw(parametersBw);
	Bwi = calc_water_Bw_p(Pa_Psia(P_initial), Kel_Far(T_in), parametersBw);
	error_SL_CH_in = 1;
	sprintf_s(escreve_pasta, "mkdir ");
	sprintf_s(pasta, "\L-%f_OFR-%f_GOR-%f", L, OFR, bbl_ft3(GOR));
	strcat_s(escreve_pasta, pasta);
	//cout << "teste=" << escreve_pasta << endl;
	system(escreve_pasta);


	cout << "teste=" << endl;
	dp = (double*)calloc(4, sizeof(double));
	dpout = (double*)calloc(4, sizeof(double));
	

	P_duct[0] = P_in;
	T_duct[0] = T_in;


	physical_prop(P_in, T_in);
	Pb0 = min(Bubble_Point(T_in, API, SG_G, GOR), P_in);
	Pb_duct[0] = Pb0;
	Rso0 = ft3_bbl(calc_gas_Rs_p(SG_G, Pa_Psia(P_in), Pa_Psia(Pb0), API, Kel_Far(T_in)));//Rs_GOR(Pb0, T_in, API, SG_G);
	Rsw0 = 0.0l;// ft3_bbl(calc_water_Rsw_p(Kel_Far(T_in), Pa_Psia(P_in)));
	Rs0 = Rso0*(1 - WC) + Rsw0*WC;
	Rs_duct[0] = Rs0;
	dp_duct[0] = 0.l;
	rhol_duct[0] = rhol;
	rhog_duct[0] = rhog;
	mul_duct[0] = mul;
	mug_duct[0] = mug;
	*error_SL_CH = error_SL_CH_in;
	
	inlet(0);

	//inlet_data_write();
	//Titulo_saida_HBm();
}

void Well::inicializa_poco(double L_in, int nz_in, double dT_in, double P_in_in, double T_in_in, double OFR_in,
	double GOR_in, double API_in, double SG_G_in, double WC_in, bool *error_SL_CH, double P_initial_in, double rug_in)
{
	L = L_in;
	nz = nz_in;
	dz = L / (nz - 1);
	dT = dT_in;
	P_in = P_in_in;
	T_in = T_in_in;
	OFR = OFR_in;
	GOR = GOR_in;
	API = API_in;
	SG_G = SG_G_in;
	rug = rug_in;
	WC = WC_in;
	P_initial = P_initial_in;
	calc_pseudocrit(&ppc, &Tpc, SG_G, y_co2, y_h2s);
	double parametersBw[9];
	calc_const_Bw(parametersBw);
	Bwi = calc_water_Bw_p(Pa_Psia(P_initial), Kel_Far(T_in), parametersBw);
	error_SL_CH_in = 1;	

	x_duct = new double[nz];
	hold_up_duct = new double[nz];
	dpdz_fric_duct = new double[nz];
	dpdz_grav_duct = new double[nz];
	dpdz_acc_duct = new double[nz];
	dp_duct = new double[nz];
	rhol_duct = new double[nz];
	rhog_duct = new double[nz];
	mul_duct = new double[nz];
	mug_duct = new double[nz];
	Rs_duct = new double[nz];
	reg_duct = new int[nz];
	vel_gas_sup = new double[nz];
	vel_liq_sup = new double[nz];
	P_duct = new double[nz]; 
	T_duct = new double[nz]; 
	Pb_duct = new double[nz];
	
	P_duct[0] = P_in;
	T_duct[0] = T_in;
		
	physical_prop(P_in, T_in);
	Pb0 = min(Bubble_Point(T_in, API, SG_G, GOR), P_in);
	Pb_duct[0] = Pb0;
	Rso0 = ft3_bbl(calc_gas_Rs_p(SG_G, Pa_Psia(P_in), Pa_Psia(Pb0), API, Kel_Far(T_in)));
	Rsw0 = 0.0l;// ft3_bbl(calc_water_Rsw_p(Kel_Far(T_in), Pa_Psia(P_in)));
	Rs0 = Rso0*(1 - WC) + Rsw0*WC;
	Rs_duct[0] = Rs0;
	dp_duct[0] = 0.;
	rhol_duct[0] = rhol;
	rhog_duct[0] = rhog;
	mul_duct[0] = mul;
	mug_duct[0] = mug;
	*error_SL_CH = error_SL_CH_in;

	inlet(0);

	//inlet_data_write();
}

double Well::fric_pure()
{
	Re = G*dT / mul;
	//cout << "Re=" << Re << endl;
	root_f = 4.e-2;
	do
	{
		root_fn = 1 / (3.48 - 4.*log10(2 * rug + 9.35 / Re / root_f));
		root_f = root_fn;
	} while (fabs(root_f - root_fn)>1e-6);

	f = pow(root_f, 2.);

	return -2 *f*G*G / rhol / dT;
}

double Well::grav_pure(double theta)
{
	//cout << "rhol----grav=" << rhol << "\tsin(theta)=" << sin(theta) << endl;;
	return -rhol*g*sin(theta);
}

double Well::fric_HBm(void)
{
	return -ff*Mt*Mt / (7.413*pow(10, 10) *pow((dT * 100 / 2.54) / 12., 5)*(rhom / 16.0185f)) / 144. / 14.7 / 10 * 1e+6 / 0.3048;;
}

double Well::fric_Friedel(void)
{
	double relo, dpdz_lo, fflo, fgo, flo, denstp, Fr, We, A1, A2, A3, PHI_lo;

	relo = G*dT / mul;

	if (relo < 2.3e+3l) fflo = 16. / relo;
	else if (relo <= 8.e+3)fflo = (1.4338e+2l /
		(pow(log(relo), 4.5l))) + 1.069e-3l;
	else if (relo <= 2.4e+4l) fflo = 7.9e-2l*pow(relo, -.25);
	else fflo = 4.6e-2l*pow(relo, -2e-1l);

	dpdz_lo = 2 * fflo*G*G / (rhol*dT);

	fgo = .079l*pow(G*dT / mug, -.25);
	flo = .079l*pow(G*dT / mul, -.25);
	denstp = 1. / (x / rhog + (1.l - x) / rhol);
	Fr = G*G / (g*dT*denstp*denstp);
	We = G*G*dT / (denstp*tens);
	A3 = pow(rhol / rhog, .91l)*pow(mug / mul, .19l)*pow(1. - mug / mul, .7);
	A2 = pow(x, .78l)*pow(1. - x, .224l);
	A1 = (1.l - x)*(1.l - x) + x*x*rhol*fgo / rhog / flo;

	PHI_lo = A1 + 3.24l*A2*A3*pow(Fr, -.045l)*pow(We, -.035l);


	return -PHI_lo*dpdz_lo;
}

double Well::fric_BeB()
{
	double f_ns, f_s, Re_ns, S, y;

	//Reynolds no slip
	Re_ns = G*dT / (mul*lambda_l + mug*(1.l - lambda_l));
	//rela��o entre os hold_up slip e no slip
	y = lambda_l / (Yl*Yl);
	//cout << "y=" << y << endl;
	//calculo do fator S
	if (y <= 1.2l && y >= 1.l) S = log(2.2l*y - 1.2l);
	else S = log(y) / (-.0523l + 3.182l*log(y) -
		.8725l*log(y)*log(y) + .01853l*pow(log(y), 4.l));
	//cout << "S=" << S << endl;
	f_ns = pow(2 * log10(Re_ns / (4.5223l*log10(Re_ns) - 3.8215l)), -2.l);

	f_s = f_ns*exp(S);
	//if (abs(Pa_Psia(-f_s*G*vm / (2.l * dT))) > 10.l || isnan(-f_s*G*vm / (2.l * dT))) {
	//	cout << "y=" << y << endl;
	//	cout << "S=" << S << endl;
	//	cout << "f_s=" << f_s << endl;
	//	cout << "G=" << G << endl;
	//	cout << "v_m=" << vm << endl;
	//	system("pause");
	//}
	//cout << "f_s=" << f_s << endl;
	return -f_s*G*vm / (2.l * dT);
}

double Well::grav_HBm(double theta)
{
	return -rhom / 16.0185f / 144. / 14.7 / 10 * 1e+6 / 0.3048*sin(theta);
}

double Well::grav_Friedel(double theta)
{
	return -(Yl*rhol + (1.l - Yl)*rhog)*g*sin(theta);
}

void Well::CHURN_model(int churn_model, double * tauw)
{
	churn_new(tauw);
}

void Well::churn_new(double * Tauw)
{
	unsigned int nmax_it = 500;
	unsigned int itera = 0;
	double tol = 1.e-12L;
	double Void = 1.L - Yl;
	double a, b, c, d, e, f, h, i;
	double dp = -1000.L;
	double d_void = .87L;
	double errf;
	double errx;
	double GG, GL;
	double dI, fIW, fIB, fI, fW, dT_star, DfIDX2;
	GG = G*x;
	GL = G*(1.L - x);
	double TAUI, TAUW;
	//la�o newton-rapshon
	dT_star = dT*sqrt((rhol - rhog)*g / tens);
	do
	{
		errf = 0.L;
		errx = 0.L;

		//processo Jayanti & Brauner
		{
			//tens�o interfacial
			dI = dT*sqrt(Void);

			fIB = .005L + pow(10.L, -.56L + 9.07L / dT_star)*pow(dT_star*(1.L - Void)*.25L, 1.63L + 4.74L / dT_star);

			fIW = .005L + .375L*(1 - Void);

			fI = .5L*(fIB + fIW);

			TAUI = fI*GG*GG*.5L / rhog / (Void*Void);

			//tens�o da parede
			Re = GL*dT / mul;
			if (Re <= 2100.L)fW = 16.L / Re;
			else fW = .079L*pow(Re, -.25L);

			TAUW = fW*GL*GL*.5L / (rhol*(1.L - Void)*(1.L - Void));

			h = 4.L / dI*TAUI + rhog*g - dp;
			i = 4.L / dT*TAUW + (rhol*(1.L - Void) + rhog*Void)*g - dp;

			a = 1.L;
			c = 1.L;
			DfIDX2 = -.5L*(pow(10.L, -.56L + 9.07L / dT_star)*pow(dT_star*(1.L - Void)*.25L, .63L + 4.74L / dT_star)*(1.63L + 4.74L / dT_star)*dT_star*.25L - .375L);

			b = 2.L*GG*GG / (rhog*dT*Void)*(DfIDX2 - fI / Void);
			d = 4.L / dT*fW / rhol*GL*GL / ((1 - Void)*(1 - Void)*(1 - Void)) - (rhol - rhog);
		}
		//termino do m�todo jayanti...

		errf = abs(h) + abs(i);
		if (errf <= tol)break;


		/*
		__  __    _ _
		|a	b| h  |h|
		|c  d| i =|i|
		--  --    - -
		*/



		//resolver sistema linear
		{
			h = -(b*i - d*h) / (a*d - b*c);
			i = (a*i - c*h) / (a*d - b*c);
		}
		//termino do sistema
		errx = abs(h) + abs(i);
		dp += h;
		Void += i;
		itera++;
	} while (errx>tol&&itera<nmax_it);
	Yl = 1.L - Void;
	*Tauw = (dp - (rhol*Yl + rhog*Void)*g)*dT*.25L;
}

void Well::ANNULAR_model()
{
	double GL, GG, GLE, GLF, GLFi, USTAR, FT, FTs,
		FTlo, FThi, dpdz_loc, TAUW, UC, rhoc, fSC, fI, TAUI,
		XCORE, R, Ri, TAUn;

	//EF = entrained_fraction();
	GL = G*(1.l - x);
	GG = G*x;
	GLE = GL*EF;
	GLF = GL*(1.l - EF);
	//cout << "GLE=" << GLE << "\tGLF=" << GLF << "\tEF=" << EF << endl;
	//system("pause");
	//chute inicial para a tens�o na parede
	dpdz_loc = fric_Friedel();
	TAUn = -dT*.25*dpdz_loc;
	//cout << "Taun=" << TAUn << endl;

	//la�o iterativo para melhor termo de fric��o com a parede
	int count1 = 0;
	R = dT*.5l;
	do
	{
		int count2 = 0;
		TAUW = TAUn;

		USTAR = sqrt(TAUW / rhol);

		//cout << "Ustar=" << USTAR << endl;
		FT = 1.e-2l;
		FTs = FT*rhol / mul*USTAR;

		do
		{
			GLFi = film_mass(FTs, dT, USTAR);
			FT = 2.l*FT / (1.l + GLFi / GLF);
			FTs = FT*rhol / mul*USTAR;
			count2++;
		} while (fabs((GLFi - GLF) / GLF) > 1.e-5l&&count2<30);
		FT = mul / rhol*FTs / USTAR;

		UC = (GLE / rhol + GG / rhog);
		rhoc = rhog + GLE / UC;
		fSC = .079l*pow((GG + GLE)*dT / mug, -.25l);
		fI = fSC*(1.l + 360.l*FT / dT);
		TAUI = .5l*rhoc*UC*UC*fI;


		Ri = R - FT;
		XCORE = x / (x*(1.l - EF) + EF);
		Yl = 1.l - (1.l - 4.l*FT / dT - (1.l - XCORE)*rhog / (rhol*XCORE + (1.l - XCORE)*rhog));
		//cout << "Yl=" << Yl << endl;
		TAUn = ((rhol - rhoc)*g + dpdz_acc + 2 * R*R*TAUI / Ri / (R*R - Ri*Ri))*(R*R - Ri*Ri)*.5l / R;
		//cout << "Taun=" << TAUn << endl;
		//system("pause");
		count1++;
	} while (fabs(TAUn - TAUW) > 1.e-6l&&count1<50);

	dpdz_fric = -4.l*TAUW / dT;
	dpdz_grav = -g*(Yl*rhol + (1.l - Yl)*rhog);
}

double Well::film_mass(double fts, double diam, double us)
{
	double R, Rstar, MLFstar;

	R = dT*.5l;
	Rstar = R*rhol / mul*us;

	if (fts <= 5.l) MLFstar = .5l*fts*fts;
	else if (fts <= 30.l) MLFstar = 12.45l - 8.05l*fts + 5.l*fts*log(fts);
	else MLFstar = 2.5l*fts*log(fts) + 2.956l*fts - 62.5l;

	return MLFstar*4.l*mul / dT;
}

void Well::SLUG_model()
{
	double BO, NF, M, GAMA, V0, VGB, VLB, VOIDB,
		DC, FM, VOIDS, BETA, FT, VLBN, REF, XGS,
		RHOSS, MUS, REUM, fUM, dpdz_fs, dpdz_fb;
	//Numero de BOND e viscosidade admensional (WALLIS, 1969 p. 288)
	BO = (rhol - rhog)*g*dT*dT / tens;
	NF = sqrt(rhol*(rhol - rhog)*g*dT*dT*dT / (mul*mul));
	//velocidade da bolha em liquido estagnado
	if (NF >= 250.l) M = 10.l;
	else if (NF > 18.l) M = 69.*pow(NF, -.35L);
	else M = 25.l;

	GAMA = .345l*(1.l - exp(-.01l*NF / .345l))*
		(1.l - exp((3.37l - BO) / M));
	V0 = GAMA*sqrt(g*dT);

	VGB = 1.2l*vm + V0;

	VLBN = -.5l;
	int iter = 0;
	do
	{
		VLB = VLBN;
		//fra��o da bolha de taylor
		VOIDB = (vm - VLB) / (VGB - VLB);
		//fra��o do pist�o liquido
		DC = 2.l * sqrt(.4l*tens / g / (rhol - rhog));
		FM = .046l*pow(rhol*vm*dT / mul, -.2l);
		VOIDS = .058l*pow(DC*pow(2.l * FM*vm*vm*vm / dT, .4l)*
			pow(rhol / tens, .6l) - .725l, 2.l);

		if (BO <= 140.) BETA = vsg / VOIDB / VGB;
		else BETA = (vsg - VOIDS*vm) / (VOIDB*VGB - VOIDS*vm);

		FT = dT*.5l*(1.l - sqrt(VOIDB));
		REF = -VLB*FT*rhol / mul;
		if (REF < 750.l) VLBN = -.333l*g*dT*dT*rhol / mul*
			pow(1.L - sqrt(VOIDB), 2.l);
		else VLBN = -11.2L*sqrt(g*dT*(1.l - sqrt(VOIDB)));
		iter++;
	} while (fabs(VLBN - VLB) > 1.e-6l&&iter < 2500);
	//termo friccional
	XGS = rhog*VOIDS / (rhog*VOIDS + rhol*(1.l - VOIDS));
	RHOSS = pow(XGS / rhog + (1.l - XGS) / rhol, -1.l);
	MUS = pow(XGS / mug + (1.l - XGS) / mul, -1.l);

	REUM = RHOSS*vm*dT / MUS;
	fUM = .079l*pow(REUM, -.25l);
	dpdz_fs = fUM*RHOSS*vm*vm*.5l / dT;

	dpdz_fb = -rhol*g*(1.l - VOIDB);

	dpdz_fric = -((1.l - BETA)*dpdz_fs + BETA*dpdz_fb);
	Yl = 1.l - ((1.l - BETA)*VOIDS + BETA*VOIDB);
	dpdz_grav = -(Yl*rhol + (1.l - Yl)*rhog)*g;
}

void Well::dpdz_GRAMP(double * dp_in, double * dpdz_out)
{
	EF = dp_in[3];

	padrao();
	if (x)
	{
		vsl = G*(1. - x) / rhol;
		vsg = G*x / rhog;
		vm = vsg + vsl;

		if (!strcmp(regime, "BUBBLE"))
		{
			Yl = Yl_Chexaullelouche(P_integral + dp_in[0], vsl, vsg);
			dpdz_fric = fric_Friedel();
			dpdz_grav = grav_Friedel(M_PI_2);
		}
		else if (!strcmp(regime, "SLUG"))
		{
			SLUG_model();
		}
		else if (!strcmp(regime, "ANNULAR"))ANNULAR_model();
		else if (!strcmp(regime, "CHURN"))
		{
			double tauw;
			CHURN_model(3, &tauw);
			dpdz_fric = -4.l / dT*tauw;
			dpdz_grav = -g*(Yl*rhol + (1.l - Yl)*rhog);
		}
		else
		{
			strcpy_s(regime, "SLUG");
			SLUG_model();
		}
	}
	else
	{
		Yl = 1.;
		dpdz_fric = fric_pure();
		dpdz_grav = grav_pure(M_PI_2);
	}
	

	hlocal = -g + dhdz_kine;
	cpmix = x*CPg + (1. - x)*CPl;
	cpj = (1.f - x) / rhol - T_integral / Z_gas*dzdtp(P_integral, T_integral)*x / rhog;

	if (!strcmp(regime, "ANNULAR"))
	{
		ENT = entrainment();
		DEP = deposition();
		dpdz_out[3] = (4.l*(ENT - DEP) / dT / G) / (1.l - x);
	}
	else if (!strcmp(regime, "CHURN"))
	{
		EF = (.95L + 342.55L*sqrt(rhol / rhog*(1.l - x) / x)*dT*dT)*1e-2L;
	}

	dpdz_out[0] = dpdz_fric + dpdz_grav + dpdz_acc;
	dpdz_out[2] = ((hlocal + cpj*dpdz_out[0]) / cpmix);
}

double Well::Yl_HBm(double P, double Uls, double Ugs)
{
	double N_lv, N_gv, N_d, N_l, psi;
	double aux_CNl, CNL, X1, X2, X3;
	double H_l;
	double auxiliar;
	double index, modified;

	auxiliar = g / (rhol*pow(tens, 3.f));
	N_l = mul*(pow(auxiliar, .25f));

	

	auxiliar = rhol / (g*tens);
	N_lv = Uls*(pow(auxiliar, .25f));
	N_gv = Ugs*(pow(auxiliar, .25f));


	auxiliar = rhol*g / tens;
	N_d = dT*sqrt(auxiliar);
	

	X1 = log10(N_l + 3.f);
	aux_CNl = -2.69851f + 1.584095e-1f*X1 - 5.509976e-1f*X1*X1 +
		5.478492e-1f*pow(X1, 3.) - 1.219458e-1f*pow(X1, 4.);
	CNL = pow(10.f, aux_CNl);
	

	X3 = N_gv * pow(N_l, 0.38f) / (pow(N_d, 2.14f));
	
	X2 = (N_lv / pow(N_gv, 0.575f)) * pow(Pa_Psia(P) / 14.7f, .1f) * (CNL / N_d);
		
	index = (X3 - .012) / fabs(X3 - .012);
	modified = (1. - index)*.006 + (1 + index)*.5*X3;
	psi = 9.116257e-1f - 4.821756e+0f*modified + 1.232250e+3f*modified*modified -
		2.225358e+4f*pow(modified, 3.) + 1.161743e+5f*pow(modified, 4.);

	auxiliar = log10(X2) + 6.;
	H_l = psi*(-1.030658e-1 + 6.177740e-1*auxiliar - 6.329460e-1*auxiliar*auxiliar +
		2.959800e-1*pow(auxiliar, 3.) - 4.010000e-2*pow(auxiliar, 4.));

	return H_l;
}

double Well::Yl_BeB(double P, double Uls, double Ugs, double theta, double Cg, double Cl)
{
	double Fr, L1, L2, L3, L4, Ylo, A, B, B0, Beta, E, F, I, J;//I e J pq o alfabeto ta foda...
	int flow_pattern;
	//Cg = Ugs / (Ugs + Uls);
	//Cl = Uls / (Ugs + Uls);

	L1 = 316 * pow(Cl, .302);
	L2 = .0009252*pow(Cl, -2.4684);
	L3 = 1.e-1*pow(Cl, -1.4516);
	L4 = 5.e-1*pow(Cl, -6.738);

	Fr = (Uls + Ugs)*(Uls + Ugs) / (g*dT);

	/*
	1-Segregado
	2-Intermitente
	3-Distribuido
	4-Transi��o
	*/

	/*cout << "Cl=" << Cl << endl;
	cout << "Fr=" << Fr << endl;*/

	if ((Cl < 1.e-2&&Fr < L1) || (Cl >= 1e-2&&Fr < L2)) flow_pattern = 1;
	else if ((Cl >= 1e-2&&Cl<4e-1 && Fr>L3&&Fr <= L1) || (Cl >= 4e-1&&Fr>L3&&Fr <= L4)) flow_pattern = 2;
	else if ((Cl < 4e-1&&Fr >= L1) || (Cl >= 4e-1&&Fr > L4)) flow_pattern = 3;
	else if ((Cl > 1e-2&&Fr > L2&&Fr < L3))flow_pattern = 4;
	else flow_pattern = 3;
	//cout << "Flow=" << flow_pattern << endl;

	//Segregado
	if (flow_pattern == 1) Ylo = .98*pow(Cl, .4846) / pow(Fr, .0868);

	//Intermitente
	else if (flow_pattern == 2) Ylo = .845*pow(Cl, .5351) / pow(Fr, .0173);

	//Distribuido
	else if (flow_pattern == 3) Ylo = 1.065*pow(Cl, .5824) / pow(Fr, .0609);

	double Nlv, YB;

	if (flow_pattern == 1)
	{
		E = 11.e-3;
		F = -3.7608;
		I = 3.5390;
		J = -1.6140;
		Nlv = Uls*(pow(rhol / (g*tens), .25f));
		Beta = (1 - Cl)*log(E*pow(Cl, F)*pow(Nlv, I)*pow(Fr, J));
		if (Beta < 0)Beta = 0.;
		B0 = 1 + Beta*(sin(1.8*theta) - pow(sin(1.8*theta), 3) / 3);
		YB = Ylo*B0;
	}
	else if (flow_pattern == 2)
	{
		E = 2.960;
		F = .3050;
		I = -.4473;
		J = .0978;
		Nlv = Uls*(pow(rhol / (g*tens), .25f));
		Beta = (1 - Cl)*log(E*pow(Cl, F)*pow(Nlv, I)*pow(Fr, J));
		if (Beta < 0)Beta = 0.;
		B0 = 1 + Beta*(sin(1.8*theta) - pow(sin(1.8*theta), 3) / 3);
		YB = Ylo*B0;
	}
	else if (flow_pattern == 3)
	{
		Beta = 0.;
		YB = Ylo;
	}
	else if (flow_pattern == 4)
	{
		A = (L3 - Fr) / (L3 - L2);
		B = 1. - A;

		E = 11.e-3;
		F = -3.7608;
		I = 3.5390;
		J = -1.6140;
		Nlv = Uls*(pow(rhol / (g*tens), .25f));
		Beta = (1 - Cl)*log(E*pow(Cl, F)*pow(Nlv, I)*pow(Fr, J));
		if (Beta < 0)Beta = 0.;
		B0 = 1 + Beta*(sin(1.8*theta) - pow(sin(1.8*theta), 3) / 3);
		Ylo = 9.8e-1*pow(Cl, .4846) / pow(Fr, .0868);
		YB = A*Ylo*B0;

		E = 2.960;
		F = .3050;
		I = -.4473;
		J = .0978;
		if (Beta < 0)Beta = 0.;
		Beta = (1 - Cl)*log(E*pow(Cl, F)*pow(Nlv, I)*pow(Fr, J));
		B0 = 1 + Beta*(sin(1.8*theta) - pow(sin(1.8*theta), 3) / 3);
		Ylo = 8.45e-1*pow(Cl, .5351) / pow(Fr, .0173);
		YB += B*Ylo*B0;
	}
	/*cout << "Ylo=" << Ylo << endl;
	cout << "YB=" << YB << endl;*/
	return YB;
}

double Well::Yl_Chexaullelouche(double P, double Uls, double Ugs)
{
	double reg, rel, pcrit, C1, C2, C3, C4, C5, C6, C7, C8, C9,
		A1, B1, L, KO, CO, r, _voidn, _void, re, D2, ugjo, ugj;

	reg = G*x*dT / mug;
	rel = G*(1. - x) / mul;

	_voidn = 5.e-1l;

	if (reg > rel) re = reg;
	else re = rel;

	pcrit = 2.212e+7l;
	C1 = 4.*pcrit*pcrit / (P*(pcrit - P));
	A1 = 1. / (1. + exp(-re / 6e+4l));
	//if (8.e-1l < A1)B1 = 8.e-1l;
	//else B1 = A1;
	B1 = min(8.e-1l, A1);
	D2 = 9.144e-2l;
	C7 = pow(D2 / dT, 6.e-1l);
	C8 = C7 / (1. - C7);
	if (C7 >= 1.l) C4 = 1.l;
	else C4 = 1.l / (1.l - exp(-C8));
	C5 = sqrt(150 * rhog / rhol);
	if ((rhog / rhol) <= 1.8e+1l) C2 = 4.757e-1l*pow(log(rhol / rhog), 7.e-1l);
	else
	{
		if (C5 >= 1.l) C2 = 1.l;
		else C2 = 1.l / (1.l - exp(1.l - C5));
	}
	C3 = max(5.e-1l, 2.*exp(-fabs(rel / 6.e+4l)));
	r = (1.l + 1.57l*rhog / rhol) / (1.l - B1);
	KO = B1 + (1.l - B1)*pow(rhog / rhol, .25);
	int iter = 0;
	do
	{
		_void = _voidn;

		C9 = pow(1.l - _void, B1);

		ugjo = 1.41l*C2*C3*C4*pow((rhol - rhog)*tens*g / rhol / rhol, .25);
		ugj = ugjo*C9;

		L = pow(_void, 2.5e-2l*(1. + 1.e+1*_void))*exp(5.e-1*(1.l - _void));

		CO = L / (KO + (1.l - KO)*pow(_void, .25l));

		_voidn = Ugs / (CO*(Uls + Ugs) + ugj);
		iter++;
	} while (fabs(_voidn-_void)>1e-10l&&iter<20);
	return 1.l-_void;
}

void Well::derivs_HBm(double * dp_in, double * dpdz_out)
{
	padrao();
	
	if (x>10e-14)
	{
		vsl = G*(1. - x) / rhol;
		vsg = G*x / rhog;
		vm  = vsg + vsl;
		lambda_g = vsg / vm;
		lambda_l = vsl / vm;
		
		
		
		Yl = Yl_HBm(P_integral + dp_in[0], vsl, vsg);

		
			
		rhom = rhol*Yl + (1. - Yl)*rhog;
		mum = pow(mul, Yl)*pow(mug, 1. - Yl);

		Rem = G*dT / mum;
		ff = pow(-4.*log10(rug / dT) / 3.7065f - 5.0452f / Rem*log10(pow(rug / dT, 1.1098f) / 2.8257f + pow(7.149f / Rem, 0.8981f)), -2.f);

		dpdz_fric = fric_HBm();
		dpdz_grav = grav_HBm(M_PI_2);		

		
	}	
	else
	{
		Yl = 1.;
		dpdz_fric = fric_pure();
		dpdz_grav = grav_pure(M_PI_2);
			
	}
	

	hlocal = -g + dhdz_kine;
	cpmix = x*CPg + (1. - x)*CPl;
	cpj= (1.f - x) / rhol - T_integral / Z_gas*dzdtp(P_integral,T_integral)*x / rhog;


	dpdz_out[0] = dpdz_fric + dpdz_grav + dpdz_acc;
	dpdz_out[2] = ((hlocal + cpj*dpdz_out[0]) / cpmix);
	dpdz_out[1] = x;
}

void Well::derivs_ChFr(double * dp_in, double * dpdz_out)
{
	padrao();
	if (x)
	{
		vsl = G*(1. - x) / rhol;
		vsg = G*x / rhog;

		Yl = Yl_Chexaullelouche(P_integral + dp_in[0], vsl, vsg);
		

		dpdz_fric = fric_Friedel();
		dpdz_grav = grav_Friedel(M_PI_2);
		

	}
	else
	{
		Yl = 1.;
		dpdz_fric = fric_pure();
		dpdz_grav = grav_pure(M_PI_2);
		
	}
	hlocal = -g + dhdz_kine;
	cpmix = x*CPg + (1. - x)*CPl;
	cpj = (1.f - x) / rhol - T_integral / Z_gas*dzdtp(P_integral, T_integral)*x / rhog;


	dpdz_out[0] = dpdz_fric + dpdz_grav + dpdz_acc;
	dpdz_out[2] = ((hlocal + cpj*dpdz_out[0]) / cpmix);
	dpdz_out[1] = x;
}

void Well::derivs_BeB(double * dp_in, double * dpdz_out)
{
	
	Rs = ft3_bbl(calc_gas_Rs_p(SG_G, Pa_Psia(P_integral + dp_in[0]), Pa_Psia(Pb), API, Kel_Far(dp_in[2])));
	x = Quality_Rs(GOR, Rs, rhol_st, rhog_st);
	
	if (x>10e-14)
	{
		vsl = G*(1. - x) / rhol;
		vsg = G*x / rhog;
		vm = vsl + vsg;
		lambda_g = vsg / vm;
		lambda_l = vsl / vm;
		Yl = Yl_BeB(P_integral + dp_in[0], vsl, vsg, M_PI_2, lambda_g, lambda_l);
		//cout << "Yl: " << Yl << endl;
		dpdz_fric = fric_BeB();
		rhom = Yl*rhol + (1.l - Yl)*rhog;
		dpdz_grav = grav_HBm(M_PI_2);
		//cout << "fric=" << dpdz_fric << endl;
		//cout << "grav=" << dpdz_grav << endl;
	}
	else
	{
		Yl = 1.;
		dpdz_fric = fric_pure();
		dpdz_grav = grav_pure(M_PI_2);
	}
	hlocal = -g + dhdz_kine;
	cpmix = x*CPg + (1. - x)*CPl;
	cpj = (1.f - x) / rhol - T_integral / Z_gas*dzdtp(P_integral, T_integral)*x / rhog;


	dpdz_out[0] = dpdz_fric + dpdz_grav + dpdz_acc;
	dpdz_out[2] = ((hlocal + cpj*dpdz_out[0]) / cpmix);
	dpdz_out[1] = x;
}

void Well::RK4_HBm(double * dp_in, double * dpz, double * dp_out)
{
	double *dy = (double*)calloc(4, sizeof(double));
	double *dym = (double*)calloc(4, sizeof(double));
	double *dyo = (double*)calloc(4, sizeof(double));
	double *yt = (double*)calloc(4, sizeof(double));
	double dx2 = .5*dz;
	double dx6 = dz / 6.;
	for (int i = 0; i < 3; i++) yt[i] = dp_in[i] + dx2*dpz[i];

	derivs_HBm(yt, dy);
	for (int i = 0; i < 3; i++) yt[i] = dp_in[i] + dx2*dy[i];
	derivs_HBm(yt, dym);
	for (int i = 0; i < 3; i++)
	{
		yt[i] = dp_in[i] + dz*dym[i];
		dym[i] += dy[i];
	}
	derivs_HBm(yt, dyo);
	//update para a saida
	for (int i = 0; i < 3; i++) dpout[i] = dp_in[i] + dx6*(dpz[i] + dyo[i] + 2.*dym[i]);

	free(dyo);
	free(dym);
	free(dy);
	free(yt);
}

void Well::RK4_ChFr(double * dp_in, double * dpz, double * dp_out)
{
	double *dy = (double*)calloc(3, sizeof(double));
	double *dym = (double*)calloc(3, sizeof(double));
	double *dyo = (double*)calloc(3, sizeof(double));
	double *yt = (double*)calloc(3, sizeof(double));
	double dx2 = .5*dz;
	double dx6 = dz / 6.;
	for (int i = 0; i < 3; i++) yt[i] = dp_in[i] + dx2*dpz[i];

	derivs_ChFr(yt, dy);
	for (int i = 0; i < 3; i++) yt[i] = dp_in[i] + dx2*dy[i];
	derivs_ChFr(yt, dym);
	for (int i = 0; i < 3; i++)
	{
		yt[i] = dp_in[i] + dz*dym[i];
		dym[i] += dy[i];
	}
	derivs_ChFr(yt, dyo);
	//update para a saida
	for (int i = 0; i < 3; i++) dpout[i] = dp_in[i] + dx6*(dpz[i] + dyo[i] + 2.*dym[i]);

	free(dyo);
	free(dym);
	free(dy);
	free(yt);
}

void Well::RK4_GRAMP(double * dp_in, double * dpz, double * dp_out)
{
	double *dy = (double*)calloc(4, sizeof(double));
	double *dym = (double*)calloc(4, sizeof(double));
	double *dyo = (double*)calloc(4, sizeof(double));
	double *yt = (double*)calloc(4, sizeof(double));
	double dx2 = .5*dz;
	double dx6 = dz / 6.;
	for (int i = 0; i < 4; i++) yt[i] = dp_in[i] + dx2*dpz[i];

	dpdz_GRAMP(yt, dy);
	for (int i = 0; i < 4; i++) yt[i] = dp_in[i] + dx2*dy[i];
	dpdz_GRAMP(yt, dym);
	for (int i = 0; i < 4; i++)
	{
		yt[i] = dp_in[i] + dz*dym[i];
		dym[i] += dy[i];
	}
	dpdz_GRAMP(yt, dyo);
	//update para a saida
	for (int i = 0; i < 4; i++) dpout[i] = dp_in[i] + dx6*(dpz[i] + dyo[i] + 2.*dym[i]);

	free(dyo);
	free(dym);
	free(dy);
	free(yt);
}

void Well::RK4_BeB(double * dp_in, double * dpz, double * dp_out)
{
	double *dy = (double*)calloc(5, sizeof(double));
	double *dym = (double*)calloc(5, sizeof(double));
	double *dyo = (double*)calloc(5, sizeof(double));
	double *yt = (double*)calloc(5, sizeof(double));
	double dx2 = .5*dz;
	double dx6 = dz / 6.;
	for (int i = 0; i < 3; i++) yt[i] = dp_in[i] + dx2*dpz[i];

	derivs_BeB(yt, dy);
	for (int i = 0; i < 3; i++) yt[i] = dp_in[i] + dx2*dy[i];
	derivs_BeB(yt, dym);
	for (int i = 0; i < 4; i++)
	{
		yt[i] = dp_in[i] + dz*dym[i];
		dym[i] += dy[i];
	}
	derivs_BeB(yt, dyo);
	//update para a saida
	for (int i = 0; i < 3; i++) dpout[i] = dp_in[i] + dx6*(dpz[i] + dyo[i] + 2.*dym[i]);

	free(dyo);
	free(dym);
	free(dy);
	free(yt);
}

void Well::integration_HB_m(int init, int end)
{
	cout.precision(15);
	dp[0] = 0.l;
	dp[1] =  x_duct[init - 1];
	dp[2] =  T_duct[init - 1];
	
	int count;
	double *dpdz = (double*)calloc(3, sizeof(double));

	P_integral = P_duct[init - 1];
	T_integral = T_duct[init - 1];


	//cout << "P_integral: " << P_integral << endl;
	//cout << "T_integral: " << T_integral << endl;


	physical_prop(P_integral, T_integral);

	dpdz_acc_n = 1.e-4;
	dhdz_kine_n = 1.e-4;

	for (int i = init; i < end; i++)
	{
		rhoho = pow(x / rhog + (1. - x) / rhol, -1.0);

		//inicio da itera��o para calculo dos termos acelerativos e perdas de carga no
		//geral.	
		//cout << "i: " << i << endl;
		
		count = 0;
		do
		{
			//cout << "Count: " << count << endl;
			dpdz_acc = dpdz_acc_n;
			dhdz_kine = dhdz_kine_n;
			
			derivs_HBm(dp, dpdz);
			RK4_HBm(dp, dpdz, dpout);


			//avan�o de press�o
			P_integral += dpout[0];			

			//avan�o de temperatura
			T_integral += dpout[2] - dp[2];		
			
			//update da press�o de bolha
			Pb = min(P_integral, Bubble_Point(T_integral, API, SG_G, GOR));
			calc_quality(P_integral, T_integral);
			//titulo da saida
			dpout[1] = x;

			//recalculo das propriedades na nova condi��o
			physical_prop(P_integral, T_integral);
			//atualiza��o da densidade homogenea
			rhoh = pow(x / rhog + (1. - x) / rhol, -1.);
			//calculo dos termos acelerativos
			//cout << "rhoh=" << rhoh << "\trhoho=" << rhoho << endl;
			dpdz_acc_n = -pow(G, 2.)*(1. / rhoh - 1. / rhoho) / dz;
			dhdz_kine_n = dpdz_acc_n / rhoh;

			//update itera��o
			count++;
			//desatualiza as condi��es
			P_integral += -dpout[0];
			T_integral += -dpout[2] + dp[2];
			physical_prop(P_integral, T_integral);

			//cout << "--------------count=" << count << endl;
		} while ((fabs((dpdz_acc_n - dpdz_acc) / dpdz_acc_n)) > 1.e-10 &&
			fabs((dhdz_kine_n - dhdz_kine) / dhdz_kine_n) > 1.e-10 &&
			count < 200);		

		
		P_integral += dpout[0];
		//avan�o de temperatura
		T_integral += dpout[2] - dp[2];
		
		//recalculo das propriedades na nova condi��o
		physical_prop(P_integral, T_integral);
		Rs = ft3_bbl(calc_gas_Rs_p(SG_G, Pa_Psia(P_integral), Pa_Psia(Pb), API, Kel_Far(T_integral)));
		x = Quality_Rs(GOR, Rs, rhol_st, rhog_st);

		P_duct[i] = P_integral;
		T_duct[i] = T_integral;
		Pb_duct[i] = Pb;
		x_duct[i] = x;
		Rs_duct[i] = Rs;
		dp_duct[i] = dpout[0];
		rhol_duct[i] = rhol;
		rhog_duct[i] = rhog;
		mul_duct[i] = mul;
		mug_duct[i] = mug;
		dpdz_fric_duct[i] = dpdz_fric;
		dpdz_grav_duct[i] = dpdz_grav;
		dpdz_acc_duct[i]  = dpdz_acc;		
		hold_up_duct[i] = Yl;

		for (int j = 1; j < 3; j++) dp[j] = dpout[j];
		//linha_saida_HBm(i);
	}
	
	free(dpdz);
}

void Well::integration_ChFr(int init, int end)
{
	cout.precision(15);
	dp[0] = 0.l;
	dp[1] = x_duct[init - 1];
	dp[2] = T_duct[init - 1];
	int count;
	double *dpdz = (double*)calloc(3, sizeof(double));

	P_integral = P_duct[init - 1];
	T_integral = T_duct[init - 1];

	physical_prop(P_integral, T_integral);

	dpdz_acc_n = 1.e-4;
	dhdz_kine_n = 1.e-4;

	for (int i = init; i < end; i++)
	{
		rhoho = pow(x / rhog + (1.l - x) / rhol, -1.0l);

		//inicio da itera��o para calculo dos termos acelerativos e perdas de carga no
		//geral.
		count = 0;
		do
		{
			
			dpdz_acc = dpdz_acc_n;
			dhdz_kine = dhdz_kine_n;

			derivs_ChFr(dp, dpdz);
			RK4_ChFr(dp, dpdz, dpout);

			//avan�o de press�o
			P_integral += dpout[0];
			//avan�o de temperatura
			T_integral += dpout[2] - dp[2];			
			//update da press�o de bolha
			Pb = min(P_integral,Bubble_Point(T_integral, API, SG_G, GOR));
			calc_quality(P_integral, T_integral);
			//titulo da saida
			dpout[1] = x;

			//recalculo das propriedades na nova condi��o
			physical_prop(P_integral, T_integral);
			//atualiza��o da densidade homogenea
			rhoh = pow(x / rhog + (1. - x) / rhol, -1.);

			//calculo dos termos acelerativos			
			dpdz_acc_n = -pow(G, 2.)*(1. / rhoh - 1. / rhoho) / dz;
			dhdz_kine_n = dpdz_acc_n / rhoh;

			//update itera��o
			count++;

			//desatualiza as condi��es
			P_integral += -dpout[0];
			T_integral += -dpout[2] + dp[2];
			physical_prop(P_integral, T_integral);

			//cout << "--------------count=" << count << endl;
		} while ((fabs((dpdz_acc_n - dpdz_acc) / dpdz_acc_n)) > 1.e-14 &&
			fabs((dhdz_kine_n - dhdz_kine) / dhdz_kine_n) > 1.e-14 &&
			count < 50);

		
		//atualiza para pr�ximo ponto ap�s sair da itera��o
		//avan�o de press�o
		P_integral += dpout[0];
		//avan�o de temperatura
		T_integral += dpout[2] - dp[2];

		//recalculo das propriedades na nova condi��o
		physical_prop(P_integral, T_integral);
		Rs = ft3_bbl(calc_gas_Rs_p(SG_G, Pa_Psia(P_integral), Pa_Psia(Pb), API, Kel_Far(T_integral)));
		x = Quality_Rs(GOR, Rs, rhol_st, rhog_st);

		P_duct[i] = P_integral;
		T_duct[i] = T_integral;
		Pb_duct[i] = Pb;
		hold_up_duct[i] = Yl;
		dpdz_fric_duct[i] = dpdz_fric;
		dpdz_grav_duct[i] = dpdz_grav;
		dpdz_acc_duct[i] = dpdz_acc;
		dp_duct[i] = P_duct[i] - P_duct[i - 1];
		rhol_duct[i] = rhol;		
		rhog_duct[i] = rhog;
		mul_duct[i] = mul;
		mug_duct[i] = mug;
		Rs_duct[i] = Rs;
		x_duct[i] = x;
		//regime
		if (strcmp(regime, "BUBBLE") == 0) reg_duct[i] = 1;
		else if (strcmp(regime, "SLUG") == 0) reg_duct[i] = 2;
		else if (strcmp(regime, "CHURN")==0) reg_duct[i] = 3;
		else if (strcmp(regime, "ANNULAR")==0) reg_duct[i] = 4;
		else reg_duct[i] = NAN;
		vel_gas_sup[i] = vsg;
		vel_liq_sup[i] = vsl;	
		

		for (int j = 1; j < 3; j++) dp[j] = dpout[j];
		//linha_saida_HBm(i);

	}

	free(dpdz);

}

void Well::integration_GRAMP(int init, int end, bool *error_SL_CH)
{
	cout.precision(4);
	dp[0] = 0;
	dp[1] = x_duct[init - 1];
	dp[2] = T_duct[init - 1];
	int count;
	double *dpdz = (double*)calloc(4, sizeof(double));

	P_integral = P_duct[init - 1];
	T_integral = T_duct[init - 1];

	physical_prop(P_integral, T_integral);

	dpdz_acc_n = 1.e-4;
	dhdz_kine_n = 1.e-4;

	for (int i = init; i < end; i++)
	{
		rhoho = pow(x / rhog + (1. - x) / rhol, -1.0);
		//inicio da itera��o para calculo dos termos acelerativos e perdas de carga no
		//geral.
		count = 0;
		do
		{
			dpdz_acc = dpdz_acc_n;
			dhdz_kine = dhdz_kine_n;
						
			dpdz_GRAMP(dp, dpdz);
			RK4_GRAMP(dp, dpdz, dpout);

			//avan�o de press�o
			P_integral += dpout[0];
			//avan�o de temperatura
			T_integral += dpout[2] - dp[2];

			//update da press�o de bolha
			Pb = min(P_integral, Bubble_Point(T_integral, API, SG_G, GOR));
			calc_quality(P_integral, T_integral);
			//titulo da saida
			dpout[1] = x;

			//recalculo das propriedades na nova condi��o
			physical_prop(P_integral, T_integral);
			//atualiza��o da densidade homogenea
			rhoh = pow(x / rhog + (1. - x) / rhol, -1.);
			//calculo dos termos acelerativos			
			dpdz_acc_n = -pow(G, 2.)*(1. / rhoh - 1. / rhoho) / dz;
			dhdz_kine_n = dpdz_acc_n / rhoh;

			//update itera��o
			count++;
			//desatualiza as condi��es
			P_integral += -dpout[0];
			T_integral += -dpout[2] + dp[2];
			physical_prop(P_integral, T_integral);

			//cout << "--------------count=" << count << endl;
		} while ((fabs((dpdz_acc_n - dpdz_acc) / dpdz_acc_n)) > 1.e-14 &&
			fabs((dhdz_kine_n - dhdz_kine) / dhdz_kine_n) > 1.e-14 &&
			count < 250);

		//atualiza para pr�ximo ponto ap�s sair da itera��o
		//avan�o de press�o
		P_integral += dpout[0];
		//avan�o de temperatura
		T_integral += dpout[2] - dp[2];
		if (isnan(P_integral))break;

		//recalculo das propriedades na nova condi��o
		physical_prop(P_integral, T_integral);
		Rs = ft3_bbl(calc_gas_Rs_p(SG_G, Pa_Psia(P_integral), Pa_Psia(Pb), API, Kel_Far(T_integral)));
		x = Quality_Rs(GOR, Rs, rhol_st, rhog_st);

		P_duct[i] = P_integral;
		T_duct[i] = T_integral;
		Pb_duct[i] = Pb;
		hold_up_duct[i] = Yl;
		dpdz_fric_duct[i] = dpdz_fric;
		dpdz_grav_duct[i] = dpdz_grav;
		dpdz_acc_duct[i] = dpdz_acc;
		dp_duct[i] = P_duct[i] - P_duct[i - 1];
		rhol_duct[i] = rhol;
		rhog_duct[i] = rhog;
		mul_duct[i] = mul;
		mug_duct[i] = mug;
		Rs_duct[i] = Rs;
		x_duct[i] = x;

		//regime
		if (!strcmp(regime, "BUBBLE")) reg_duct[i] = 1;
		else if (!strcmp(regime, "SLUG")) reg_duct[i] = 2;
		else if (!strcmp(regime, "CHURN")) reg_duct[i] = 3;
		else if (!strcmp(regime, "ANNULAR")) reg_duct[i] = 4;
		else reg_duct[i] = NAN;
		vel_gas_sup[i] = vsg;
		vel_liq_sup[i] = vsl;		

		for (int j = 1; j < 4; j++) dp[j] = dpout[j];
		
	}

	*error_SL_CH = error_SL_CH_in;

	free(dpdz);
	
}

void Well::integration_BeB(int init, int end, int cor)
{
	cout.precision(4);
	dp[0] = 0;
	dp[1] = x_duct[init - 1];
	dp[2] = T_duct[init - 1];
	int count;
	double *dpdz = (double*)calloc(3, sizeof(double));

	P_integral = P_duct[init - 1];
	T_integral = T_duct[init - 1];
	//cout << "P_integral: " << Pa_Psia(P_integral) << endl;
	physical_prop(P_integral, T_integral);

	dpdz_acc_n = 1.e-4;
	dhdz_kine_n = 1.e-4;

	for (int i = init; i < end; i++)
	{
		rhoho = pow(x / rhog + (1. - x) / rhol, -1.0);
		//inicio da itera��o para calculo dos termos acelerativos e perdas de carga no
		//geral.
		count = 0;
		do
		{			
			dpdz_acc = dpdz_acc_n;
			dhdz_kine = dhdz_kine_n;
			
			derivs_BeB(dp, dpdz);
			RK4_BeB(dp, dpdz, dpout);

			//avan�o de press�o
			P_integral += dpout[0];
			//avan�o de temperatura
			T_integral += dpout[2] - dp[2];

			//update da press�o de bolha
			Pb = min(P_integral, Bubble_Point(T_integral, API, SG_G, GOR));
			calc_quality(P_integral, T_integral);
			//titulo da saida
			dpout[1] = x;

			//recalculo das propriedades na nova condi��o
			physical_prop(P_integral, T_integral);
			//atualiza��o da densidade homogenea
			rhoh = pow(x / rhog + (1. - x) / rhol, -1.);

			//calculo dos termos acelerativos
			//cout << "rhoh=" << rhoh << "\trhoho=" << rhoho << endl;
			dpdz_acc_n = -pow(G, 2.)*(1. / rhoh - 1. / rhoho) / dz;
			dhdz_kine_n = dpdz_acc_n / rhoh;

			//update itera��o
			count++;
			//desatualiza as condi��es
			P_integral += -dpout[0];
			T_integral += -dpout[2] + dp[2];
			physical_prop(P_integral, T_integral);

			//cout << "--------------count=" << count << endl;
		} while ((fabs((dpdz_acc_n - dpdz_acc) / dpdz_acc_n)) > 1.e-14 &&
			fabs((dhdz_kine_n - dhdz_kine) / dhdz_kine_n) > 1.e-14 &&
			count < 250);

		//atualiza para pr�ximo ponto ap�s sair da itera��o
		//avan�o de press�o
		P_integral += dpout[0];
		//avan�o de temperatura
		T_integral += dpout[2] - dp[2];
		if (isnan(P_integral))break;

		//if (i == init) system("pause");
		//recalculo das propriedades na nova condi��o
		physical_prop(P_integral, T_integral);
		Rs = ft3_bbl(calc_gas_Rs_p(SG_G, Pa_Psia(P_integral), Pa_Psia(Pb), API, Kel_Far(T_integral)));
		x = Quality_Rs(GOR, Rs, rhol_st, rhog_st);
		P_duct[i] = P_integral;
		T_duct[i] = T_integral;
		Pb_duct[i] = Pb;
		hold_up_duct[i] = Yl;
		dpdz_fric_duct[i] = dpdz_fric;
		dpdz_grav_duct[i] = dpdz_grav;
		dpdz_acc_duct[i] = dpdz_acc;
		dp_duct[i] = P_duct[i] - P_duct[i - 1];
		rhol_duct[i] = rhol;
		rhog_duct[i] = rhog;
		mul_duct[i] = mul;
		mug_duct[i] = mug;
		Rs_duct[i] = Rs;
		x_duct[i] = x;
		//regime
		if (!strcmp(regime, "BUBBLE")) reg_duct[i] = 1;
		else if (!strcmp(regime, "SLUG")) reg_duct[i] = 2;
		else if (!strcmp(regime, "CHURN")) reg_duct[i] = 3;
		else if (!strcmp(regime, "ANNULAR")) reg_duct[i] = 4;
		else reg_duct[i] = NAN;
		vel_gas_sup[i] = vsg;
		vel_liq_sup[i] = vsl;

		for (int j = 1; j < 3; j++) dp[j] = dpout[j];

		}


	free(dpdz);

}

void Well::inlet(int pos)
{
	cout.precision(4);
	double OFR_SI, flow_oil, flow_solution_gas, flow_free_gas;
	Area = M_PI*dT*dT*.25;

	//densidades na condi��o std
	rhol_st = 141.5f / (API + 131.5f)*(1 - WC)*rhow_std + WC*rhow_std;
	Z_gas = Z_factor(P_std, T_std, SG_G);
	rhog_st = SG_G*rho_air;

	//conver��o de bbl/d para m�/s
	OFR_SI = bbld_m3s(OFR);

	//vaz�o de �leo morto-com agua
	flow_oil = rhol_st*OFR_SI;

	xini = Quality_Rs(GOR, Rs0, rhol_st, rhog_st);
	x = xini;
		
	x_duct[pos] = x;

	//vaz�o de g�s dissolvido
	flow_solution_gas = (rhog_st / rhol_st)*Rs0*flow_oil;

	//vaz�o de g�s livre
	flow_free_gas = (rhog_st / rhol_st)*(GOR - Rs0)*flow_oil;
	
	G = (flow_oil + flow_solution_gas + flow_free_gas) / Area;

	//if (model == "HB")
	Mt = OFR*((GOR*5.619)*SG_G*.0765f +
		(141.5f / (131.5f + API))*62.4f*5.615);
	
}

void Well::physical_prop(double P, double T)
{
	double parametrosmug[25];
	Pb = min(P, Bubble_Point(T_integral, API, SG_G, GOR));
	Rsb = ft3_bbl(calc_gas_Rs_p(SG_G, Pa_Psia(P_in), Pa_Psia(Pb), API, Kel_Far(T_in)));
	double SG_O = 141.5f / (API + 131.5f);
	double Bo=calc_oil_Bo_p(SG_O, SG_G, Pa_Psia(P), Pa_Psia(Pb), API, Kel_Far(T));
	rhool = calc_oil_rho_p(Rsb, Bo, SG_O*rhow_std, SG_G*rho_air);
	muol= cP_Pas(calc_oil_mi_p(bbl_ft3(Rsb), Kel_Far(T), API, Pa_Psia(P), Pa_Psia(Pb)));

	double Bg = bbl_ft3(calc_gas_Bg_p(Pa_Psia(P)));
	rhog = calc_gas_rho_p(SG_G*rho_air, Bg);
	calc_const_mig(parametrosmug, SG_G, y_n2, y_co2, y_h2s, Kel_Far(T));
	mug = cP_Pas(calc_gas_mig_p(parametrosmug, Kel_Rank(T)/Tpc, Pa_Psia(P)/ppc));

	//CHECAR AQUI PRA VER SE BATE COM O WELLBORE
	tens = Tensuper(rhol, rhog);
	H_CP_Black(T, API, &Hl, &CPol);
	Hlog = 0.;
	//

	rhow = calc_water_rho_p(Bwi, rhow_std, Kel_Far(T), Pa_Psia(P), Pa_Psia(P_initial));
	muw = cP_Pas(calc_water_mi_p(Kel_Far(T), Pa_Psia(P)));

	water_H_CP(P, T, &Hlw, &Hlgw, &CPwl, &CPwg);

	rhol = rhow*WC + rhool*(1 - WC);
	mul = muw*WC + muol*(1 - WC);
	CPl = CPwl*WC + CPol*(1 - WC);
	CPg = 3.e+3;
	
}

double Well::dzdtp(double P, double T)
{
	double deltat = 1.e-6;
	double Z1, Z2;
	Z1 = Z_factor(P, T, SG_G);
	Z2 = Z_factor(P, T + deltat, SG_G);
	return (Z2 - Z1) / deltat;
}

void Well::calc_quality(double P, double T)
{	
	
	Rso = ft3_bbl(calc_gas_Rs_p(SG_G, Pa_Psia(P), Pa_Psia(Pb), API, Kel_Far(T)));
	Rsw = 0.0l;// ft3_bbl(calc_water_Rsw_p(Kel_Far(T_in), Pa_Psia(P_in)));
	Rs = Rso*(1.l - WC) + Rsw*WC;
	x = Quality_Rs(GOR, Rs, rhol_st, rhog_st);

	if (x < 0.l) x = 0.l;

}

void Well::padrao(void)
{
	double UGSstar, Void, flood;

	DRIFT_VOID_BUBBLE(&Void);
	//cout << "void: " << Void << endl;

	if (x == 0) UGSstar = 0.;
	else UGSstar = G*x / sqrt(g*dT*(rhol - rhog)*rhog);

	if (UGSstar < 1.&&UGSstar>0) SLUG_TO_CHURN(&flood);
	else flood = 10.;

	if (x <= 0.&&UGSstar <= 0) strcpy_s(regime, "SPHASE");
	else if (Void < .25) strcpy_s(regime, "BUBBLE");
	else
	{
		if (UGSstar >= 1.) strcpy_s(regime, "ANNULAR");
		else
		{
			if (flood >= 1.) strcpy_s(regime, "CHURN");
			else strcpy_s(regime, "SLUG");
		}
	}
}

void Well::DRIFT_VOID_BUBBLE(double * Void)
{
	double UGJ, voidn;
	voidn = .1;
	Co = 1.;
	vsl = G*(1. - x) / rhol;
	vsg = G*x / rhog;
	vm = vsl + vsg;

	do
	{
		*Void = voidn;
		UGJ = 1.53f*pow(tens*g*(rhol - rhog) / (rhol*rhol), .25);
		voidn = vsg / (Co*(vm)+UGJ);
	} while (fabs(voidn - *Void) >= 1e-6);
}

void Well::SLUG_TO_CHURN(double * flooding)
{
	double FTSLU, UBS, UFS, REFS, FTSLUn, UFSstar, UBSstar, LTAYLOR, LD, MCOEF;

	FTSLUn = 2e-4;
	FTSLU = 2e-4;
	int count = 0;

	do
	{
		if (count>0) FTSLU = 0.6l*FTSLUn+0.4l*FTSLU;

		UBS = (1.l - 4.l*FTSLU / dT)*(1.2F*(vm)+.35l*sqrt(g*dT*(rhol - rhog) / rhol));
		UFS = UBS - vm;
		REFS = rhol*UFS*FTSLU / mul;

		FTSLUn = .1719l*pow(REFS, 2.l / 3.l) /
			pow(g*(rhol - rhog) / pow(mul / rhol, 2.l) / rhol, 1.l / 3.l);

		if (count % 50==0 && count!=0 ) {
			FTSLUn = FTSLUn + 0.00015;
		}
				
		count = count + 1;

	} while (fabs(FTSLU-FTSLUn)>1e-6 && count<500);

	UBS = UFS + vm;

	UFSstar = UFS*sqrt(rhol / (g*dT*(rhol - rhog)));
	UBSstar = UBS*sqrt(rhog / (g*dT*(rhol - rhog)));

	LTAYLOR = 6.l*dT*(vsg - vsl) / (UFS + vsl);

	LD = LTAYLOR / dT;
	if (LD <= 1.2e+2) MCOEF = 1.92e-1l + 1.089e-2l*LD - 3.754e-5l*LD*LD;
	else MCOEF = 9.6e-1l;

	if(count>=500) error_SL_CH_in=0;


	*flooding = sqrt(UBSstar) + MCOEF*sqrt(UFSstar);
}

void Well::clear_mem_well() {

	delete[] x_duct;
	delete[] hold_up_duct;
	delete[] dpdz_fric_duct;
	delete[] dpdz_grav_duct;
	delete[] dpdz_acc_duct;
	delete[] dp_duct;
	delete[] rhol_duct;
	delete[] rhog_duct;
	delete[] mul_duct;
	delete[] mug_duct;
	delete[] Rs_duct;
	delete[] reg_duct;
	delete[] vel_gas_sup;
	delete[] vel_liq_sup;
	delete[] P_duct;
	delete[] T_duct;
	delete[] Pb_duct;

}

Well::~Well()
{

	delete[] x_duct;
	delete[] hold_up_duct;
	delete[] dpdz_fric_duct;
	delete[] dpdz_grav_duct;
	delete[] dpdz_acc_duct;
	delete[] dp_duct;
	delete[] rhol_duct;
	delete[] rhog_duct;
	delete[] mul_duct;
	delete[] mug_duct;
	delete[] Rs_duct;
	delete[] reg_duct;
	delete[] vel_gas_sup;
	delete[] vel_liq_sup;
	delete[] P_duct;
	delete[] T_duct;
	delete[] Pb_duct;
	free(dp);
	free(dpout);

}

void Well::inlet_data_write()
{
	setlocale(LC_ALL, "");
	char arq[600];
	char write[600];
	sprintf_s(write, pasta);
	sprintf_s(arq, "\\data_inlet.txt");
	strcat_s(write, arq);

	fstream data(write, ios::out);
	data.precision(5);

	data << "---par�metros iniciais---" << endl;
	data << "L=\t" << L << "\t|m" << endl;
	data << "dT=\t" << dT << "\t|m" << endl;
	data << "nz=\t" << nz << "\t|" << endl;
	data << "dz=\t" << dz << "\t|m" << endl;
	data << "P_in=\t" << P_in << "\t|Pa" << endl;
	data << "T_in=\t" << T_in << "\t|K" << endl;
	data << "OFR=\t" << OFR << "\t|bbl/d" << endl;
	data << "GOR=\t" << GOR*5.619 << "\t|SCF/stb" << endl;
	data << "API=\t" << API << "\t|" << endl;
	data << "SG_G=\t" << SG_G << "\t|" << endl;
	data << "-------calculado---------" << endl;
	data << "G=\t" << G << "\t\t|kg/m2s" << endl;
	data << "Mt=\t" << Mt << "\t|lbm/d" << endl;
}

void Well::pres_depth_mat_write()
{
	setlocale(LC_ALL, "");
	char arq[50];
	char write[600];
	sprintf_s(write, pasta);
	sprintf_s(arq, "\\pres.m");
	strcat_s(write, arq);

	fstream data(write, ios::out);
	data.precision(14);

	data << "p=[";

	for (int i = 0; i < nz; i++) data << P_duct[i] << "\t" << i*dz << endl;

	data << "];";



}

void Well::buble_pres_mat_write()
{
	setlocale(LC_ALL, "");
	char arq[50];
	char write[500];
	sprintf_s(write, pasta);
	sprintf_s(arq, "\\bub_pres.m");
	strcat_s(write, arq);

	fstream data(write, ios::out);
	data.precision(14);

	data << "pb=[";

	for (int i = 0; i < nz; i++) data << Pb_duct[i] << "\t" << i*dz << endl;

	data << "];";
}

void Well::Titulo_saida_HBm()
{
	char arq[50];
	char write[500];
	sprintf_s(write, pasta);
	sprintf_s(arq, "\\data_out.csv");
	strcat_s(write, arq);

	fstream tit(write, ios::out);

	tit << "ponto,altura,profundidade,pressao,temperatura,titulo,hold_up,pb,dpdz_fric,dpdz_grav,dpdz_acc,dp,rhol,rhog,viscl,viscg,regime,Rs" << endl;
}

void Well::linha_saida_HBm(int pos)
{
	char arq[50];
	char write[500];
	sprintf_s(write, pasta);
	sprintf_s(arq, "\\data_out.csv");
	strcat_s(write, arq);

	fstream lin(write, ios::out | ios::app);
	lin.precision(6);

	//tit << "ponto,altura,profundidade,pressao,temperatura,titulo,hold_up,pb,dpdz_fric,dpdz_grav,dpdz_acc,rhol,rhog,viscl,viscg" << endl;
	lin << pos << "," << pos*dz << "," << L - pos*dz << "," << Pa_Psia(P_duct[pos]) << "," << T_duct[pos] << "," << x << "," << Yl << "," << Pa_Psia(Pb) << "," <<
		Pa_Psia(dpdz_fric) << "," << Pa_Psia(dpdz_grav) << "," << Pa_Psia(dpdz_acc) << "," << dpout[0] << "," << rhol << "," << rhog << "," << mul << "," << mug << "," << regime
		<< "," << bbl_ft3(Rs) << endl;
}

double Well::get_bot_P()
{
	return P_duct[0];
}

double Well::get_bot_P1()
{
	return P_duct[1];
}

double Well::get_surf_P()
{
	return P_duct[nz - 1];
}

double Well::get_bot_T()
{
	return T_duct[0];
}

double Well::get_surf_T()
{
	return T_duct[nz - 1];
}

double Well::entrainment()
{
	double GLF, GLE, GG;
	GLF = G*(1.l - x)*(1.l - EF);
	GLE = G*(1.l - x)*EF;
	GG = G*x;

	double group,
		glfc = mul / dT*exp(5.8504l + .4249l*(mug / mul)*pow(rhol / rhog, .5l));
	if (GLF < glfc) return 0.l;
	else
	{
		group = 5.75e-5l*pow(dT*rhol / tens / (rhog*rhog), .316l);
		return group*GG*pow(GLF - glfc, .632l);
	}
}

double Well::deposition()
{
	double GLF, GLE, GG;
	GLF = G*(1.l - x)*(1.l - EF);
	GLE = G*(1.l - x)*EF;
	GG = G*x;
	double C, group, K;
	C = GLE / (GG / rhog + GLE / rhol);
	group = pow(rhog*dT / tens, -.5l);
	if ((C / rhog) <= .3l)K = .18l*group;
	else K = .083l*group*pow(rhog / C, .65l);

	return K*C;
}

void Well::grav_tudo(int coup_ts,int iter_coup, char directory[])
{

	char arc[200];
	char arc2[30];
	sprintf_s(arc, directory);

	sprintf_s(arc2, "\\well_data_out_duct_%d_%d.csv", coup_ts, iter_coup);
	
	strcat_s(arc, arc2);

	int i = 0;

	fstream	data(arc, ios::out);
	data.precision(12);

	data << "ponto;altura [m];profundidade [m];pressao [psi];temperatura [K];titulo ;hold_up ;pb [psi];dpdz_fric [psi];dpdz_grav [psi];dpdz_acc [psi];dp [psi];rhol;rhog;viscl;viscg;Rs [scf/stb];vsl;vsg;regime" << endl;	
	for (int i = 1; i < nz; i++)
		data << i << ";" << i*dz << ";" << (L - i*dz) << ";" << Pa_Psia(P_duct[i]) << ";" << T_duct[i] << ";" << x_duct[i] <<
		";" << hold_up_duct[i] << ";" << Pa_Psia(Pb_duct[i]) << ";" << Pa_Psia(dpdz_fric_duct[i]) << ";" << Pa_Psia(dpdz_grav_duct[i]) <<
		";" << dpdz_acc_duct[i] << ";" << dp_duct[i] << ";" << rhol_duct[i] << ";" << mul_duct[i] <<
		";" << rhog_duct[i] << ";" << mug_duct[i] << ";" << bbl_ft3(Rs_duct[i]) << ";" << vel_liq_sup[i] << ";" << vel_gas_sup[i] <<
		";" << reg_duct[i] << endl;

	data.close();
}
