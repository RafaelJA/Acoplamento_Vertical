#pragma once
class Well
{
private:
	//Tubbing
	double L, rug, dT, Area;

	//Well integration dimensions
	double dz;
	int nz;

	//Inlet and Outlet conditions
	double T_in, P_in, P_out;
	double G, OFR;

	//Base well and flow condition
	double GOR, API, SG_G;
	double WC, Rso0, Rsw0;

	//Oil base condition
	double Rso, Rsb, xini, Pb0, Rs0;
	double rhol_st, rhog_st;

	//Variaveis do poço
	/*double x_duct[1000];
	double hold_up_duct[1000];
	double dpdz_fric_duct[1000];
	double dpdz_grav_duct[1000];
	double dpdz_acc_duct[1000];
	double dp_duct[1000];
	double rhol_duct[1000];
	double rhog_duct[1000];
	double mul_duct[1000];
	double mug_duct[1000];
	double Rs_duct[1000];
	int reg_duct[1000];
	double vel_gas_sup[1000];
	double vel_liq_sup[1000];*/
	double* x_duct;
	double* hold_up_duct;
	double* dpdz_fric_duct;
	double* dpdz_grav_duct;
	double* dpdz_acc_duct;
	double* dp_duct;
	double* rhol_duct;
	double* rhog_duct;
	double* mul_duct;
	double* mug_duct;
	double* Rs_duct;
	int* reg_duct;
	double* vel_gas_sup;
	double* vel_liq_sup;

	//propriedades Black Oil
	double rhol, mul, CPl, Hl, Co;
	double rhog, mug, Z_gas, CPg, Hlg;
	double rhool, muol, CPol, Hlog;
	double Pb, Rs, tens, x;
	double cpmix, dzdt_p, cpj, hlocal;
	double ppc, Tpc;
	double y_co2 = 0.0l, y_h2s = 0.0l, y_n2 = 0.0l;

	//água
	double rhow, muw, Hlw, Hlgw, CPwl, CPwg;
	double Rsw;
	double P_initial, Bwi;

	//vetor de pressão e temperatura
	//double P_duct[1000], T_duct[1000], Pb_duct[1000];
	double* P_duct, * T_duct, * Pb_duct;

	//padrão de escoamento
	char regime[30];

	//condição std
	const long double P_std = 101325.l, T_std = 288.71l, rhow_std = 998.97l, rho_air = 1.2928l;

	//gravidade
	const long double g = 9.81l;

	//termos de perda de carga
	double dpdz_fric, dpdz_grav, dpdz_acc;

	//termo de perda de cinética
	double dhdz_kine;

	//termos iterativos
	double dpdz_acc_n, dhdz_kine_n;

	//Temperatura e pressão dentro da integral
	double P_integral, T_integral;

	//termos da mistura com slip
	double rhom, mum, Rem, fn, ff;

	//termos da mistura no slip
	double lambda_l, lambda_g, rhos, mus;

	//termos da mistura homogenea
	double rhoh, rhoho;

	//vetores de propragação
	double* dp = (double*)calloc(5, sizeof(double));
	double* dpout = (double*)calloc(5, sizeof(double));

	//correlações adotadas/
	//velocidades superficiais necessárias
	double vsl = 0.0, vsg = 0.0, vm, Yl;

	//HBm
	double Mt, Lb;

	//Beggs&Brill

	//void drift

	//modelo fenomenologico GRAMP
	double EF, DEP, ENT;

	//Petalas&Aziz

	//termo de fricção para monofásico
	double f, Re, root_f, root_fn;

	//variavel de arquivos e pasta
	char escreve_pasta[300];
	char pasta[300];

	//remendo
	bool error_SL_CH_in;

public:
	//construtor padrão
	Well();
	//construtor para parametros de entrada
	Well(double L_in, int nz_in, double dT_in, double P_in_in, double T_in_in, double OFR_in, double GOR_in, double API_in, double SG_G_in, double WC_in, bool* error_SL_CH, double P_initial_in, double rug_in);//terminar-ok

	//inicializador caso nao seja usado o construtor
	void inicializa_poco(double L_in, int nz_in, double dT_in, double P_in_in, double T_in_in, double OFR_in,
		double GOR_in, double API_in, double SG_G_in, double WC_in, bool* error_SL_CH, double P_initial_in, double rug_in);

	//perdas de carga S_PHASE
	double fric_pure();
	double grav_pure(double theta);

	//fricção cada modelo
	double fric_HBm(void);
	double fric_Friedel(void);//continuar_linha 2525 gramp -> n DPDZ_FRIELEL
	double fric_BeB();

	//gravidade cada modelo
	double grav_HBm(double theta);
	double grav_Friedel(double theta);

	//GRAMP
	void CHURN_model(int churn_model, double* tauw);
	void churn_new(double* Tauw);
	void ANNULAR_model();
	double film_mass(double fts, double diam, double us);
	void SLUG_model();
	void dpdz_GRAMP(double* dp_in, double* dpdz_out);

	//hold up de cada modelo
	double Yl_HBm(double P, double Uls, double Ugs);
	double Yl_BeB(double P, double Uls, double Ugs, double theta, double Cg, double Cl);
	double Yl_Chexaullelouche(double P, double Uls, double Ugs);

	//derivadas
	void derivs_HBm(double* dp_in, double* dpdz_out);
	void derivs_ChFr(double* dp_in, double* dpdz_out);
	void derivs_BeB(double* dp_in, double* dpdz_out);

	//RK4
	void RK4_HBm(double* dp_in, double* dpz, double* dp_out);
	void RK4_ChFr(double* dp_in, double* dpz, double* dp_out);
	void RK4_GRAMP(double* dp_in, double* dpz, double* dp_out);
	void RK4_BeB(double* dp_in, double* dpz, double* dp_out);

	//integrais
	void integration_HB_m(int init, int end);
	void integration_ChFr(int init, int end);
	void integration_GRAMP(int init, int end, bool* error_SL_CH);
	void integration_BeB(int init, int end, int cor);

	//inlet
	void inlet(int pos);

	//propriedades
	void physical_prop(double P, double T);
	double dzdtp(double P, double T);
	void calc_quality(double P, double T);

	//determinação do padrão de escoamento
	void padrao(void);
	void DRIFT_VOID_BUBBLE(double* Void);
	void SLUG_TO_CHURN(double* flooding);

	void clear_mem_well();
	~Well();

	void inlet_data_write();

	void pres_depth_mat_write();
	void buble_pres_mat_write();

	void Titulo_saida_HBm();
	void linha_saida_HBm(int pos);

	double get_bot_P();
	double get_bot_P1();
	double get_surf_P();
	double get_bot_T();
	double get_surf_T();

	//funções do entranhamento
	double entrainment();
	double deposition();

	void grav_tudo(int t, int it, char directory[]);



};

