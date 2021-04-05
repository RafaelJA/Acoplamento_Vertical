#include "stdafx.h"
#include "calc_wellbore.h"

void wellbore_mono(int correlation, int N, double *PWB, double *Q_I,
	double BHP, double PWB_0, double *deltaL,
	double rho, double mu,
	double D, double D_perf, double perf_dens, 
	double teta, double g, double e, double gama, double teta_yue) {
		
	/*  1- Ouyang
		2 - Zhang
		3 - Siwon
		4 - Asheim
		5 - Yalniz-Not working
		6 - Yue*/

	double DELTAP,error, BHP_CALC;

	if(correlation == 1) {

		double Q_med, U_med, Q;
		double f1, f0, PWB_1,a;
		error = 10;
		double cont = 0;

		while (error>1e-7) {					

			PWB[N-1] = PWB_0;
			Q = 0.0l;
			for (int i = N-2; i>=-1; i--) {	

				Q_med = Q + 0.5l*Q_I[i + 1] * deltaL[i + 1];
				U_med = 4.0l*Q_med / (M_PI*D*D);
				Q = Q + Q_I[i + 1] *deltaL[i + 1];
				
				DELTAP = calc_pres_drop_ouyang_mono(rho, mu, D, D_perf, perf_dens, teta, g, e, gama, deltaL[i+1], Q_med, U_med, Q_I[i+1]);
		
				if (i == -1) {					
					BHP_CALC= PWB[i + 1] - Pa_Psia(DELTAP);
					
				}
				else{
					PWB[i] = PWB[i + 1] - Pa_Psia(DELTAP);					
				}
				
								
			}	

			error = abs(BHP_CALC - BHP);
			

			if (cont == 0) {
				f1 = BHP_CALC - BHP;
				PWB_1 = PWB_0;
				if (BHP_CALC > BHP) {
					PWB_0 = PWB_0 - 0.5l;
				}
				else {
					PWB_0 = PWB_0 + 0.5l;
				}
			}
			else {
				f0 = BHP_CALC - BHP;
				a = (PWB_1 - PWB_0) / (f1 - f0);
				PWB_1 = PWB_0;
				f1 = f0;
				PWB_0 = PWB_0-f0/a;
			}

			cont++;

		}

	}

	else if (correlation == 2) {

		//cout << "ZHANG" << endl;
		//q_w wall influx m3/s
		//rho kg/m�
		//deltax m
		//S perimeter m
		//phi wall opening ratio
		//eps_0 relative roughness of a regular pipe

		double phi = 0.l;
		double A = 0.25l*M_PI*D*D;
		double S = M_PI*D;
		
		double Q;
		double f1, f0, PWB_1, a;
		error = 10;
		double cont = 0;

		while (error>1e-7) {

			PWB[N - 1] = PWB_0;
			Q = 0.0l;

			for (int i = N - 2; i >= -1; i--) {
				
				phi = perf_dens*D_perf*D_perf / (D*deltaL[i + 1]);

				Q = Q + Q_I[i + 1] * deltaL[i + 1];

				DELTAP = calc_pres_drop_zhang(Q_I[i + 1] * deltaL[i + 1], Q, deltaL[i + 1], S, phi, A, D, rho, mu, e);
				
				if (i == -1) {
					BHP_CALC = PWB[i + 1] - Pa_Psia(DELTAP);
				}
				else {
					PWB[i] = PWB[i + 1] - Pa_Psia(DELTAP);
				}							

			}

			error = abs(BHP_CALC - BHP);

			if (cont == 0) {
				f1 = BHP_CALC - BHP;
				PWB_1 = PWB_0;
				if (BHP_CALC > BHP) {
					PWB_0 = PWB_0 - 0.5l;
				}
				else {
					PWB_0 = PWB_0 + 0.5l;
				}
			}
			else {
				f0 = BHP_CALC - BHP;
				a = (PWB_1 - PWB_0) / (f1 - f0);
				PWB_1 = PWB_0;
				f1 = f0;
				PWB_0 = PWB_0 - f0 / a;
			}

			cont++;

		}
		
	}

	else if (correlation == 3) {

		//cout << "SIWON" << endl;
		//q_w wall influx m3/s
		//rho kg/m�
		//deltax m
		//S perimeter m
		//phi wall opening ratio
		//eps_0 relative roughness of a regular pipe
		
		double phi = 0.0l;
		double A = 0.25l*M_PI*D*D;
		double S = M_PI*D;

		double Q;
		double f1, f0, PWB_1, a;
		error = 10;
		double cont = 0;

		while (error>1e-7) {
		

			PWB[N - 1] = PWB_0;
			Q = 0.0l;

			for (int i = N - 2; i >= -1; i--) {

				phi = perf_dens*D_perf*D_perf / (D*deltaL[i + 1]);
				Q = Q + Q_I[i + 1] * deltaL[i + 1];

				DELTAP = calc_pres_drop_siwon(Q_I[i + 1] * deltaL[i + 1], Q, deltaL[i + 1], S, phi, A, D, rho, mu, e);
				
				if (i == -1) {
					BHP_CALC = PWB[i + 1] - Pa_Psia(DELTAP);
				}
				else {
					PWB[i] = PWB[i + 1] - Pa_Psia(DELTAP);
				}

			}

			error = abs(BHP_CALC - BHP);

			if (cont == 0) {
				f1 = BHP_CALC - BHP;
				PWB_1 = PWB_0;
				if (BHP_CALC > BHP) {
					PWB_0 = PWB_0 - 0.5l;
				}
				else {
					PWB_0 = PWB_0 + 0.5l;
				}
			}
			else {
				f0 = BHP_CALC - BHP;
				a = (PWB_1 - PWB_0) / (f1 - f0);
				PWB_1 = PWB_0;
				f1 = f0;
				PWB_0 = PWB_0 - f0 / a;
			}

			cont++;


		}

	}

	else if (correlation == 4) {
		
		//cout << "ASHEIM" << endl;
		//q_w wall influx m3/s
		//rho kg/m�
		//deltax m
		//S perimeter m
		//phi wall opening ratio
		//eps_0 relative roughness of a regular pipe
		//n perf dens

		double phi = 0.0l;
		double A = 0.25l*M_PI*D*D;
		double S = M_PI*D;

		double Q;
		double f1, f0, PWB_1, a;
		error = 10;
		double cont = 0;

		while (error>1e-7) {

			PWB[N - 1] = PWB_0;
			Q = 0.0l;

			for (int i = N - 2; i >= -1; i--) {

				phi = perf_dens*D_perf*D_perf / (D*deltaL[i + 1]);
				Q = Q + Q_I[i + 1] * deltaL[i + 1];

				DELTAP = calc_pres_drop_asheim(Q_I[i + 1] * deltaL[i + 1], Q, deltaL[i + 1], S, phi, A, D, rho, mu, e, perf_dens);
				
				if (i == -1) {
					BHP_CALC = PWB[i + 1] - Pa_Psia(DELTAP);
				}
				else {
					PWB[i] = PWB[i + 1] - Pa_Psia(DELTAP);
				}

			}

			error = abs(BHP_CALC - BHP);

			if (cont == 0) {
				f1 = BHP_CALC - BHP;
				PWB_1 = PWB_0;
				if (BHP_CALC > BHP) {
					PWB_0 = PWB_0 - 0.5l;
				}
				else {
					PWB_0 = PWB_0 + 0.5l;
				}
			}
			else {
				f0 = BHP_CALC - BHP;
				a = (PWB_1 - PWB_0) / (f1 - f0);
				PWB_1 = PWB_0;
				f1 = f0;
				PWB_0 = PWB_0 - f0 / a;
			}

			cont++;

		}

	}

	else if (correlation == 5) {
		cout << "ESSA CORRELA��O N�O FAZ NENHUM SENTIDO.PAU NO CU DO YALNIZ(TALVEZ DA YALNIZ) E DO OSKAN" << endl;
		//cout << "YALNIZ" << endl;
		//q_w wall influx m3/s
		//rho kg/m�
		//deltax m
		//S perimeter m
		//phi wall opening ratio
		//eps_0 relative roughness of a regular pipe
		//n perf dens

		/*
		double phi = perf_dens*D_perf*D_perf / (D*deltaL);
		double A = 0.25l*M_PI*D*D;
		double S = M_PI*D;

		double Q_med, Q;
		error = 10;

		while (error>1e-7) {

			
			PWB[0] = PWB_0;
			Q = 0.0l;

			for (int i = N - 2; i >= -1; i--) {
				
				
				DELTAP = calc_pres_drop_yalniz(Q, Q_I[i + 1] * deltaL, D, A, D_perf, deltaL, rho, mu);

				Q = Q + Q_I[i + 1] * deltaL;

				if (i == -1) {
					BHP_CALC= PWB[i + 1] - Pa_Psia(DELTAP);
				}
				else{
					PWB[i] = PWB[i + 1] - Pa_Psia(DELTAP);
				}

			}

			error = abs(BHP_CALC - BHP);

			if (cont == 0) {
				f1 = BHP_CALC - BHP;
				PWB_1 = PWB_0;
				if (BHP_CALC > BHP) {
					PWB_0 = PWB_0 - 0.5l;
				}
				else {
					PWB_0 = PWB_0 + 0.5l;
				}
			}
			else {
				f0 = BHP_CALC - BHP;
				a = (PWB_1 - PWB_0) / (f1 - f0);
				PWB_1 = PWB_0;
				f1 = f0;
				PWB_0 = PWB_0-f0/a;
			}

		}
		*/




	}

	else if (correlation == 6) {
		
		//Q main flow rate m�/s
		//q Volumetric influx flow rate from each injection opening, m�/s/shoot
		//teta Completion shot phasing �
		//phi completion shot density shots/m
		
		double phi = perf_dens;		

		double Q;
		double f1, f0, PWB_1, a;
		error = 10;
		double cont = 0;

		while (error>1e-7) {

			PWB[N - 1] = PWB_0;
			Q = 0.0l;

			for (int i = N - 2; i >= -1; i--) {
						

				DELTAP = calc_pres_drop_yue(Q, Q_I[i + 1]/phi, teta_yue, phi, D, deltaL[i + 1], rho, mu);

				Q = Q + Q_I[i + 1] * deltaL[i + 1];

				if (i == -1) {
					BHP_CALC = PWB[i + 1] - Pa_Psia(DELTAP);
				}
				else {
					PWB[i] = PWB[i + 1] - Pa_Psia(DELTAP);
				}

			}

			error = abs(BHP_CALC - BHP);

			if (cont == 0) {
				f1 = BHP_CALC - BHP;
				PWB_1 = PWB_0;
				if (BHP_CALC > BHP) {
					PWB_0 = PWB_0 - 0.5l;
				}
				else {
					PWB_0 = PWB_0 + 0.5l;
				}
			}
			else {
				f0 = BHP_CALC - BHP;
				a = (PWB_1 - PWB_0) / (f1 - f0);
				PWB_1 = PWB_0;
				f1 = f0;
				PWB_0 = PWB_0 - f0 / a;
			}

			cont++;


		}

	}
	
	else {
		cout << "ERROR 001- THE WELLBORE CORRELATION CANNOT BE FOUND" << endl;
	}

}