#include "stdafx.h"
#include "pressure_drop_wellbore.h"


double cal_fric_fact_wellbore(double Re, double Rew, double D, double eps) {
	if (Re <= 3000) {
		if (Rew >= 0) {
			return cal_fric_fact_lam_prod_mono(Re, Rew);
		}
		else {

			return cal_fric_fact_lam_inj_mono(Re, Rew);
		}
	}

	else {
		double f0 = cal_fric_fact_turb_colle(Re, D, eps);
		if (Rew >= 0) {
			//cout << "f0: " << f0 << endl;
			return cal_fric_fact_turb_prod_mono(Re, Rew, f0);
		}
		else {
			return cal_fric_fact_turb_inj_mono(Re, Rew, f0);
		}
	}
}

double cal_fric_fact_lam_prod_mono(double Re, double Rew) {
	return 16.l / Re * (1.0l + 0.04304l * pow(Rew, 0.6142l));
}

double cal_fric_fact_lam_inj_mono(double Re, double Rew) {

	if (abs(Rew) < 4.626l) {
		return (16.l / Re) * (1.0l - 0.0625l * pow(-Rew, 1.3056l) / pow((4.626l + Rew), -0.2724l));
	}

	return (16.l / Re) * (1.0l - 0.0625l * pow(-Rew, 1.3056l) / pow((0.000000000001l), -0.2724l));
}

double cal_fric_fact_turb_prod_mono(double Re, double Rew, double f0) {
	return f0 * (1.0l - 0.0153l * pow(Rew, 0.3978l));
}

double cal_fric_fact_turb_inj_mono(double Re, double Rew, double f0) {
	return f0 * (1.0l - 17.5l * Rew / pow(Re, 0.75l));
}

double cal_fric_fact_turb_colle(double Re, double D, double eps) {
	double A = -2 * log10(eps / (3.7l * D) + 12.0l / Re);
	double B = -2 * log10(eps / (3.7l * D) + 2.51l * A / Re);
	double C = -2 * log10(eps / (3.7l * D) + 2.51l * B / Re);
	return 0.25l / pow(A - (B - A) * (B - A) / (C - 2 * B + A), 2);
}

double cal_Reynolds(double rho, double Q, double mu, double D) {
	return 4.0l * rho * Q / (M_PI * mu * D);//1.4776*rho*Q/ (mu*D);
}

double cal_Reynolds_w(double rho, double Q_I, double mu) {
	return 4.0l * rho * Q_I / (M_PI * mu);//1.4776*rho*J*drawdown/(mu*D_perf*perf_dens);
}

double calc_pres_drop_ouyang_mono(double rho, double mu, double D, double D_perf, double perf_dens, double teta,
	double g, double e, double gama, double deltaL, double Q_med, double U_med, double Q_I) {

	double Re = cal_Reynolds(rho, Q_med, mu, D);


	double Re_w = cal_Reynolds_w(rho, Q_I, mu);


	double f = cal_fric_fact_wellbore(abs(Re), Re_w, D, D * e);


	double R_af = Q_I * D / (f * Q_med);
	double R_gf = g * D * sin(teta) / (2.0l * f * U_med * U_med);
	double R_da = 0.25l * (4.0l * Q_I / (perf_dens * M_PI * D_perf * D_perf)) * sin(2.l * gama) / U_med;


	return 2.0l * f * (deltaL / D) * rho * U_med * U_med * (1.l + R_af * (1.l - R_da) + R_gf);
}

double calc_f_0_Chen(double eps_0, double Re) {
	if (Re < 2500) {
		return 64.0l / Re;
	}

	double A = -1.81l * log10(pow(eps_0 / 3.7l, 1.11l) + 6.9l / Re);
	return 1 / (A * A);
}

double calc_pres_drop_zhang(double q_w, double q, double deltax, double S, double phi, double A, double D,
	double rho, double mu, double eps_0) {
	//q_w wall influx m3/s
	//rho kg/m³
	//deltax m
	//S perimeter m
	//phi wall opening ratio
	//eps_0 relative roughness of a regular pipe

	double v_w = q_w / (S * deltax);
	double v = q / A;
	double xi_w = v_w / v;
	double Re = rho * v * D / mu;

	double f_0 = calc_f_0_Chen(eps_0, Re);
	/*cout << "f_0: " << f_0 << endl;
	system("pause");*/

	double f_w_i = 0.98l * f_0;//initial guess
	double f_w = 8001;//mais de 8000
	double A1, A2;

	while (abs(f_w_i - f_w) > (1e-7)) {
		A1 = sqrt(8.0l / f_0) - log(Re) * (phi - xi_w) + 1.25 * log(f_w_i / f_0);
		A2 = A1 * A1;
		f_w = f_w_i;
		f_w_i = 8.0l / A2;
	}

	return -rho * v * v * (0.5l * f_w_i * (1 - phi) + 8.0l * xi_w + 32.0l * xi_w * xi_w) * deltax / D;
}

double calc_pres_drop_siwon(double q_w, double q, double deltax, double S, double phi, double A, double D,
	double rho, double mu, double eps_0) {

	double v_p = q_w / (phi * S * deltax);
	double q_s = q_w / deltax;
	double v = q / A;
	double xi_p = v_p / v;
	double Re = rho * v * D / mu;

	double f_p = 0.0106l * pow(phi, 0.413l);
	double eps = eps_0 + 0.282l * pow(phi, 2.4l);

	double  b = 10.l * pow(1000l * phi, -4.2l) + 4e-7;
	double beta = 1.05l * (1.l + 1.175 * pow(b / (xi_p * xi_p) + 1.235l, -2.0l));

	double f_0 = 0.0;

	if (Re < 2500) {
		f_0 = 64.0l / Re;
	}
	else {
		f_0 = 0.11l * pow(68.0l / Re + eps, 0.25l);
	}

	double f_w = f_p + f_0;

	return -(0.5l * f_w * rho * v * v / D + beta * rho * v * q_s / A) * deltax;

}

double calc_pres_drop_asheim(double q_w, double q, double deltax, double S, double phi, double A, double D,
	double rho, double mu, double eps_0, double n) {

	//n perforation density

	double v_p = q_w / (phi * S * deltax);
	double q_s = q_w / deltax;
	double v = q / A;
	double Re = rho * v * D / mu;

	double f_w = 0.16l * pow(Re, -0.19l);//-0.19 meixmo
	double ksi_s = q_s / q;

	double f_acc = 2.l * D * ksi_s * (2.l + ksi_s / n);

	double f_t = f_acc + f_w;

	return -(0.5l * f_t * rho * v * v / D) * deltax;

}

double calc_pres_drop_yalniz(double q_m, double q_i, double D, double A, double d_i, double deltax, double rho, double mu) {

	//q_m main flow rate
	//q_i injection flow rate
	//d_i perforation diameter 

	double q_t = q_m + q_i;
	double v = q_t / A;
	double v2 = 0.5l * (q_t + q_m) / A;
	double Re = rho * v * D / mu;

	double d_i_mm = d_i * 1000.l;

	double f_w = (8.0l - d_i_mm) * (98.52l / Re) / (8.0l - 4.0l)
		+ (d_i_mm - 4.0l) * (57.22l / Re) / (8.0l - 4.0l);

	double f_p;
	double R = q_i / q_t;
	double H = sqrt(D / d_i);

	if (R > 0.0933l * H && R <= 0.184l * H) {
		// cout << "1" << endl;
		f_p = (H * H * H * H * D / deltax) * (0.045444l - 0.00011457l * H * H / (R * R));
	}

	if (R > 0.0129l * H && R <= 0.0933l * H) {
		// cout << "2" << endl;
		f_p = (H * H * H * H * D / deltax) * (0.5373l * pow(H * H, -0.52321l) * pow(R, 1.184642l));
	}

	if (R >= 0.00732l * H && R <= 0.0129l * H) {
		//cout << "3" << endl;
		f_p = (H * H * H * H * D / deltax) * (0.019889l - 0.000001927l * log(H * H) / (R * R));
		/*cout << "X: " << R*R / (H*H) << endl;
		cout << "Y: " << f_p*deltax / (D*H*H*H*H) << endl;*/
		cout << "H: " << H << endl;
		cout << "R: " << R << endl;
	}

	double f_t = f_w + f_p;

	/*cout << "f_w: " << f_w << endl;
	cout << "f_p: " << f_p << endl;
	cout << "f_t: " << f_t << endl;
	system("pause");*/

	return -(0.5l * f_t * rho * v2 * v2 * deltax / D);

}

double calc_B_yue(double Re, double teta, double phi) {
	//teta Completion shot phasing °
	//phi completion shot density shots/ft

	if (phi >= 5.l && phi <= 10.l) {

		double phi1 = 5.l;
		double phi2 = 10.l;

		if (teta >= 90.l && teta <= 180.l) {
			double teta1 = 90.l;
			double teta2 = 180.l;
			double B11 = 0.873l * pow(Re, -0.335l);
			double B12 = 0.622l * pow(Re, -0.295l);
			double B21 = 0.160l * pow(Re, -0.170l);
			double B22 = 0.360l * pow(Re, -0.265l);

			return 1.l / ((teta2 - teta1) * (phi2 - phi1)) * (B11 * (teta2 - teta) * (phi2 - phi) + B21 * (teta - teta1) * (phi2 - phi) +
				B12 * (teta2 - teta) * (phi - phi1) + B22 * (teta - teta1) * (phi - phi1));
		}

		if (teta >= 180.l && teta <= 360.l) {
			double teta1 = 180.l;
			double teta2 = 360.l;
			double B11 = 0.622l * pow(Re, -0.295l);
			double B12 = 0.651l * pow(Re, -0.314l);
			double B21 = 0.360l * pow(Re, -0.265l);
			double B22 = 0.755l * pow(Re, -0.290l);

			return 1.l / ((teta2 - teta1) * (phi2 - phi1)) * (B11 * (teta2 - teta) * (phi2 - phi) + B21 * (teta - teta1) * (phi2 - phi) +
				B12 * (teta2 - teta) * (phi - phi1) + B22 * (teta - teta1) * (phi - phi1));
		}

	}

	if (phi >= 10.l && phi <= 20.l) {

		double phi1 = 10.l;
		double phi2 = 20.l;

		if (teta >= 90.l && teta <= 180.l) {
			double teta1 = 90.l;
			double teta2 = 180.l;
			double B11 = 0.160l * pow(Re, -0.170l);
			double B12 = 0.360l * pow(Re, -0.265l);
			double B21 = 1.471l * pow(Re, -0.435l);
			double B22 = 0.452l * pow(Re, -0.270l);

			return 1.l / ((teta2 - teta1) * (phi2 - phi1)) * (B11 * (teta2 - teta) * (phi2 - phi) + B21 * (teta - teta1) * (phi2 - phi) +
				B12 * (teta2 - teta) * (phi - phi1) + B22 * (teta - teta1) * (phi - phi1));
		}

		if (teta >= 180.l && teta <= 360.l) {
			double teta1 = 180.l;
			double teta2 = 360.l;
			double B11 = 0.360l * pow(Re, -0.265l);
			double B12 = 0.755l * pow(Re, -0.290l);
			double B21 = 0.452l * pow(Re, -0.270l);
			double B22 = 1.078l * pow(Re, -0.350l);

			return 1.l / ((teta2 - teta1) * (phi2 - phi1)) * (B11 * (teta2 - teta) * (phi2 - phi) + B21 * (teta - teta1) * (phi2 - phi) +
				B12 * (teta2 - teta) * (phi - phi1) + B22 * (teta - teta1) * (phi - phi1));
		}

	}

}

double calc_pres_drop_yue(double Q, double q, double teta, double phi, double D, double deltax, double rho, double mu) {
	//Q main flow rate m³/s
	//q Volumetric influx flow rate from each injection opening, m³/s/shoot
	//teta Completion shot phasing °
	//phi completion shot density shots/m

	double A = (0.0042l * teta + 1.866l) * phi + 1298.98l / teta - 19.7185l;

	double deltaQ = q * (phi * deltax);
	double u_avg = (Q + 0.5 * deltaQ) / (0.25l * M_PI * D * D);

	double B = calc_B_yue(rho * u_avg * D / mu, teta, 0.3048l * phi);
	double f_t = B + 2 * D * A * q / (Q + 0.5 * deltaQ);

	/*cout << "A: " << 2 * D*A*q / (Q + 0.5*deltaQ) << endl;
	cout << "B: " << B << endl;
	system("pause");*/

	return -0.5 * f_t * rho * u_avg * u_avg * deltax / D;

}

