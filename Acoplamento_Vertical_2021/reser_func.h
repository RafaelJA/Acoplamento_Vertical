#pragma once

//void Sleep(float s);

double calcula_THP(double tempo, double *tempo_THP, double *THP);

void ini_fields(int m, int n, int o, double yo, double rho_ws,
	double *So_coup, double *Sw_coup, double *Sg_coup,
	double *po_coup, double *pw_coup, double *pg_coup, double *pb_coup, bool *acimapb,
	double *So_layer, double *Sw_layer, double *deltaz, double pi, double pbi);

void inicializa_reser(int *res_timesteps, double *H, double *Re, double *Rw, int *n, int *m, int *o, double *deltaphi, double *deltaz, double *a, double *R0, double *pi, double *pbi, double *Swi,
	double *Soi, double *T, double *g, double *C1, double *C2, double *C3, double *C4, double *C5, double *C6, double *Qwsinj, double *Qosinj, double *Qgsinj, double *rho_air, double *rho_ws, double *API, double *yg,
	double *yo, double *rho_os, double *rho_gs, double *yn2, double *yco2, double *yh2s, double *phi_ref, double *Cr, double *pgr, double *Swc, double *Sor, double *Sgr, double *epsilon, double *pcwmin, double *pcwmax, double *pcgmin,
	double *pcgmax, double *krwmax, double *nw, double *ng, double *now, double *nog, double *parametrosBw, double *parametrosmig, double *ppc, double *Tpc, double *parametrosZ, double *Bwi, double *Bpcw, double *Bpcg, bool *zonaprodutora);

void calc_reser_u(int n, int m, int o, int *IA, int *JA,
	int *unW, int *ufW, int *uwbW, int *ueW, int *ucW, int *uwW, int *uebW, int *ubW, int *usW,
	int *unO, int *ufO, int *uwbO, int *ueO, int *ucO, int *uwO, int *uebO, int *ubO, int *usO,
	int *unG, int *ufG, int *uwbG, int *ueG, int *ucG, int *uwG, int *uebG, int *ubG, int *usG);


void calc_water_pos_res(double pwp, double pwnm,
	double *B, double *A,
	double *controladorP, double *controladorSw,
	int y, int u,
	double pesop, double pesocon, double Ageo, double C,
	double lambda1, double lambda2, double lambda3, double lambda4, double lambda5,
	double lambda1int, double lambda2int, double lambda3int, double lambda4int, double lambda5int);


void calc_water_grav_res(double rhonm, double rhop, double drhonm, double drhop,
	double *B, double *A,
	double *controladorP, double *controladorSw,
	int y, int u,
	double pesop, double pesocon, double Ageo, double C,
	double lambda1, double lambda2, double lambda5,
	double lambda1int, double lambda2int, double lambda5int,
	double deltaz, double g);

void calc_oil_pos_res(double pop, double ponm,
	double *B, double *A,
	double *controladorP, double *controladorSw, double *controladorSo,
	int y, int u,
	double pesop, double pesocon, double Ageo, double C,
	double lambda1, double lambda2, double lambda3, double lambda4,
	double lambda1int, double lambda2int, double lambda3int, double lambda4int);

void calc_oil_grav_res(double rhop, double rhonm, double drhop_p, double drhonm_p, double drhop_pb, double drhonm_pb,
	double *B, double *A,
	double *controladorP, double *controladorSw, double *controladorSo,
	int y, int u,
	double pesop, double pesocon, double Ageo, double C,
	double lambda1, double lambda2, double lambda3, double lambda4,
	double lambda1int, double lambda2int, double lambda3int, double lambda4int,
	double deltaz, double g);

void calc_gas_pos_res(double pop, double ponm, double pgp, double pgnm,
	double *B, double *A,
	double *controladorP, double *controladorSw, double *controladorSo,
	int y, int u,
	double pesop, double pesocon, double Ageo, double C,
	double lambda1, double lambda2, double lambda3, double lambda4, double lambda5, double lambda6, double lambda7, double lambda8, double dpcgp,
	double lambda1int, double lambda2int, double lambda3int, double lambda4int, double lambda5int, double lambda6int, double lambda7int, double lambda8int, double dpcgnm,
	bool acimapbp, bool acimapbnm);

void calc_gas_grav_res(double rhoop, double drhoop, double rhoonm, double drhoonm, double rhogp, double drhogp, double rhognm, double drhognm,
	double drhoopbp, double drhoopbnm,
	double *B, double *A,
	double *controladorP, double *controladorSw, double *controladorSo,
	int y, int u,
	double pesop, double pesocon, double Ageo, double C,
	double lambda1, double lambda2, double lambda3, double lambda4, double lambda5, double lambda7, double lambda8,
	double lambda1int, double lambda2int, double lambda3int, double lambda4int, double lambda5int, double lambda7int, double lambda8int,
	double deltaz, double g, bool acimapbp, bool acimapbnm);

void calc_reser_evol(int passos, int p, int y, int u, int l, int centerP, int centerSw, int centerSo,
	int *IA, int *JA, double *A, double *B, double *X,
	bool dontstopmenow, int i, int k, int j, int o, int m, int n, int dim, 
	double controladorP, double controladorSw, double controladorSo,
	double *lambda1_w, double *lambda2_w, double *lambda4_w, double *lambda5_w,
	double *lambda1_o, double *lambda2_o, double *lambda3_o, double *lambda4_o,
	double *lambda1_g, double *lambda2_g, double *lambda3_g, double *lambda4_g,
	double *lambda5_g, double *lambda6_g, double *lambda7_g, double *lambda8_g,
	double *Kx, double *Kz, double pi, double T, double *PWB, double *PWB_0, double *Pb_WB,
	double *po_res, double *pb_res, double *pw_res, double *pg_res, double *Sw_res, double *So_res, double *Sg_res,
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
	double *Qw_well, double *Qo_well, double *Qg_well, int zy, int cont,
	double *phi, double *der_phi, double *rho_w, double *der_rho_w, double *mi_w, double *der_mi_w,
	double *Bw, double *der_Bw, double *pcw, double *der_pcw, double *krw, double *der_krw,
	double *Rs, double *der_Rs_p, double *der_Rs_pb, double *kro, double *der_kro_Sw, double *der_kro_So,
	double *Bo, double *der_Bo_p, double *der_Bo_pb,
	double *mi_o, double *der_mi_o_p, double *der_mi_o_pb, double *rho_o, double *der_rho_o_p, double *der_rho_o_pb,
	double *pcg, double *der_pcg_Swo, double *krg, double *der_krg, double *mi_g, double *der_mi_g, double *Z,
	double *Bg, double *der_Bg, double *rho_g, double *der_rho_g,
	bool *acimapb, bool *zonaprodutora, int level,
	int *unW, int *ufW, int *uwbW, int *ueW, int *ucW, int *uwW, int *uebW, int *ubW, int *usW,
	int *unO, int *ufO, int *uwbO, int *ueO, int *ucO, int *uwO, int *uebO, int *ubO, int *usO,
	int *unG, int *ufG, int *uwbG, int *ueG, int *ucG, int *uwG, int *uebG, int *ubG, int *usG);