#pragma once

void coup_step_4(int coup_ts, int ite_coup,
	double BHP_r, double *BHP,
	double THP, double L, int nz, double *deltaz,
	double *Qw_total, double *Qo_total, double *Qg_total,
	double *PWB_ref, double *PWB, double *PWB_0, double *dPWB_D, double *PWB_coup,
	double *Qw_wellbore, double *Qo_wellbore, double *Qg_wellbore,
	double *Qw_wellbore_ref, double *Qo_wellbore_ref, double *Qg_wellbore_ref,
	int N_wellbore, double D_perf, double perf_dens,
	double teta, double rug, double gama, double teta_yue,
	double *po_coup, double *So_coup, double *Sw_coup, double *pb_coup,
	double *po_res, double *pw_res, double *pg_res,
	double *So_res, double *Sw_res, double *Sg_res, double *pb_res,
	double *po_res_0, double *So_res_0, double *Sw_res_0, double *pb_res_0, double *Sg_res_0, double *X_pb_0,
	double *phi_res_0, double *Bw_res_0, double *Rs_res_0, double *Bo_res_0, double *Bg_res_0,
	double *po_res_lv0, double *pb_res_lv0, double *So_res_lv0, double *Sw_res_lv0,
	double *varpo, double *varpb, double *varSo, double *varSw,
	double Qwsinj, double Qosinj, double Qgsinj,
	double *K_F, double *K_B, double *K_N, double *K_S, double *K_W, double *K_E,
	double *X_pb,bool *zonaprodutora,
	int res_timesteps, bool dontstopmenow,
	int p, int y, int u, int l, 
	int centerP, int centerSw, int centerSo, 
	int i, int k, int j, 
	int o, int m, int n, int dim,
	double *Kr, double *Kz, double pi, double T,
	int *IA, int *JA, double *A, double *B, double *X,
	double controladorP, double controladorSw, double controladorSo,
	double C1, double C2, double C3,
	double *Ageo_R, double *Ageo_PHI, double *Ageo_Z, double *Ageo_t,
	double *Adir_F, double *Adir_B, double *Adir_W, double *Adir_E, double *Adir_N, double *Adir_S,
	double *CHI_F, double *CHI_B, double *CHI_W, double *CHI_E, double *CHI_N, double *CHI_S,
	double *R_P, double *PHI_P, double *Z_P,
	double *R_PB, double *R_PF, double *delta_PHI, double *delta_Z, double deltat_res,
	double g, double Re, double Rw, 
	double rho_air, double rho_ws, double API, double yo, double rho_os, double yg, double rho_gs,
	double yn2, double yco2, double yh2s,
	double *lambda1_W, double *d_lambda1_po_W, double *d_lambda1_SW_W,
	double *lambda2_W, double *d_lambda2_po_W, double *d_lambda2_SW_W,
	double *lambda1_O, double *d_lambda1_po_O, double *d_lambda1_SW_O, double *d_lambda1_pb_O, double *d_lambda1_SO_O,
	double *lambda2_O, double *d_lambda2_po_O, double *d_lambda2_SW_O, double *d_lambda2_pb_O, double *d_lambda2_SO_O,
	double *lambda1_G, double *d_lambda1_po_G, double *d_lambda1_SW_G, double *d_lambda1_pb_G, double *d_lambda1_SO_G,
	double *lambda2_G, double *d_lambda2_po_G, double *d_lambda2_SW_G, double *d_lambda2_pb_G, double *d_lambda2_SO_G,
	double *lambda3_G, double *d_lambda3_po_G, double *d_lambda3_SW_G, double *d_lambda3_pb_G, double *d_lambda3_SO_G,
	double *lambda4_G, double *d_lambda4_po_G, double *d_lambda4_SW_G, double *d_lambda4_pb_G, double *d_lambda4_SO_G,
	double krwmax, double nw, double ng, double now, double nog, 
	double *phi_ref, double Cr,
	double Swc, double Sor, double Sgr, double epsilon,
	double pcwmin, double pcwmax, double pcgmin, double pcgmax,
	double *parametrosBw, double *parametrosmig, double *parametrosZ,
	double ppc, double Tpc, double Bwi, double Bpcw, double Bpcg, 
	int *iparm, double *dparm, int solver, int maxfct, int mnum, int phase,
	int error, int msglvl, int mtype, int nrhs, double ddum, int idum, void *pt,
	double *phi, double *der_phi, 
	double *rho_w, double *der_rho_w,
	double *mi_w, double *der_mi_w,
	double *Bw, double *der_Bw, 
	double *pcw, double *der_pcw,
	double *krw, double *der_krw,
	double *Rs, double *der_Rs_p, double *der_Rs_pb, 
	double *kro, double *der_kro_Sw, double *der_kro_So, 
	double *Bo, double *der_Bo_p, double *der_Bo_pb,
	double *mi_o, double *der_mi_o_p, double *der_mi_o_pb, 
	double *rho_o, double *der_rho_o_p, double *der_rho_o_pb,
	double *pcg, double *der_pcg_Sw, double *der_pcg_So,
	double *krg, double *der_krg, 
	double *mi_g, double *der_mi_g, double *Z, 
	double *Bg, double *der_Bg,
	double *rho_g, double *der_rho_g,
	int *unW, int *ufW, int *uwbW, int *ueW, int *ucW, int *uwW, int *uebW, int *ubW, int *usW,
	int *unO, int *ufO, int *uwbO, int *ueO, int *ucO, int *uwO, int *uebO, int *ubO, int *usO,
	int *unG, int *ufG, int *uwbG, int *ueG, int *ucG, int *uwG, int *uebG, int *ubG, int *usG,
	int *neg, bool *NpureW, int *contresint,
	char directory[], __int8 W_corre, __int8 WB_corre);

void coup_step_1(int coup_ts, int ite_coup,
	double BHP_r, double *BHP,
	double THP, double L, int nz, double *deltaz,
	double *Qw_total, double *Qo_total, double *Qg_total,
	double *PWB_ref, double *PWB, double *PWB_0, double *dPWB_D, double *PWB_coup,
	double *Qw_wellbore, double *Qo_wellbore, double *Qg_wellbore,
	double *Qw_wellbore_ref, double *Qo_wellbore_ref, double *Qg_wellbore_ref,
	int N_wellbore, double D_perf, double perf_dens,
	double teta, double rug, double gama, double teta_yue,
	double *po_coup, double *So_coup, double *Sw_coup, double *pb_coup,
	double *po_res, double *pw_res, double *pg_res, 
	double *So_res, double *Sw_res, double *Sg_res, double *pb_res,
	double *po_res_0, double *So_res_0, double *Sw_res_0, double *pb_res_0, double *Sg_res_0, double *X_pb_0,
	double *phi_res_0, double *Bw_res_0, double *Rs_res_0, double *Bo_res_0, double *Bg_res_0,
	double *po_res_lv0, double *pb_res_lv0, double *So_res_lv0, double *Sw_res_lv0,
	double *varpo, double *varpb, double *varSo, double *varSw,
	double Qwsinj, double Qosinj, double Qgsinj,
	double *K_F, double *K_B, double *K_N, double *K_S, double *K_W, double *K_E,
	double *X_pb,bool *zonaprodutora,
	int res_timesteps, bool dontstopmenow,
	int p, int y, int u, int l, 
	int centerP, int centerSw, int centerSo, 
	int i, int k, int j, 
	int o, int m, int n, int dim,
	double *Kr, double *Kz, double pi, double T,
	int *IA, int *JA, double *A, double *B, double *X,
	double controladorP, double controladorSw, double controladorSo,
	double C1, double C2, double C3,
	double *Ageo_R, double *Ageo_PHI, double *Ageo_Z, double *Ageo_t,
	double *Adir_F, double *Adir_B, double *Adir_W, double *Adir_E, double *Adir_N, double *Adir_S,
	double *CHI_F, double *CHI_B, double *CHI_W, double *CHI_E, double *CHI_N, double *CHI_S,
	double *R_P, double *PHI_P, double *Z_P,
	double *R_PB, double *R_PF, double *delta_PHI, double *delta_Z, double deltat_res,
	double g, double Re, double Rw,
	double rho_air, double rho_ws, double API, double yo, double rho_os, double yg, double rho_gs,
	double yn2, double yco2, double yh2s,
	double *lambda1_W, double *d_lambda1_po_W, double *d_lambda1_SW_W,
	double *lambda2_W, double *d_lambda2_po_W, double *d_lambda2_SW_W,
	double *lambda1_O, double *d_lambda1_po_O, double *d_lambda1_SW_O, double *d_lambda1_pb_O, double *d_lambda1_SO_O,
	double *lambda2_O, double *d_lambda2_po_O, double *d_lambda2_SW_O, double *d_lambda2_pb_O, double *d_lambda2_SO_O,
	double *lambda1_G, double *d_lambda1_po_G, double *d_lambda1_SW_G, double *d_lambda1_pb_G, double *d_lambda1_SO_G,
	double *lambda2_G, double *d_lambda2_po_G, double *d_lambda2_SW_G, double *d_lambda2_pb_G, double *d_lambda2_SO_G,
	double *lambda3_G, double *d_lambda3_po_G, double *d_lambda3_SW_G, double *d_lambda3_pb_G, double *d_lambda3_SO_G,
	double *lambda4_G, double *d_lambda4_po_G, double *d_lambda4_SW_G, double *d_lambda4_pb_G, double *d_lambda4_SO_G,
	double krwmax, double nw, double ng, double now, double nog, 
	double *phi_ref, double Cr,
	double Swc, double Sor, double Sgr, double epsilon,
	double pcwmin, double pcwmax, double pcgmin, double pcgmax,
	double *parametrosBw, double *parametrosmig, double *parametrosZ,
	double ppc, double Tpc, double Bwi, double Bpcw, double Bpcg, 
	int *iparm, double *dparm, int solver, int maxfct, int mnum, int phase,
	int error, int msglvl, int mtype, int nrhs, double ddum, int idum, void *pt,
	double *phi, double *der_phi, 
	double *rho_w, double *der_rho_w,
	double *mi_w, double *der_mi_w,
	double *Bw, double *der_Bw, 
	double *pcw, double *der_pcw,
	double *krw, double *der_krw,
	double *Rs, double *der_Rs_p, double *der_Rs_pb, 
	double *kro, double *der_kro_Sw, double *der_kro_So, 
	double *Bo, double *der_Bo_p, double *der_Bo_pb,
	double *mi_o, double *der_mi_o_p, double *der_mi_o_pb, 
	double *rho_o, double *der_rho_o_p, double *der_rho_o_pb,
	double *pcg, double *der_pcg_Sw, double *der_pcg_So,
	double *krg, double *der_krg, 
	double *mi_g, double *der_mi_g, double *Z, 
	double *Bg, double *der_Bg,
	double *rho_g, double *der_rho_g,
	int *unW, int *ufW, int *uwbW, int *ueW, int *ucW, int *uwW, int *uebW, int *ubW, int *usW,
	int *unO, int *ufO, int *uwbO, int *ueO, int *ucO, int *uwO, int *uebO, int *ubO, int *usO,
	int *unG, int *ufG, int *uwbG, int *ueG, int *ucG, int *uwG, int *uebG, int *ubG, int *usG,
	int *neg, bool *NpureW, int *contresint, 
	char directory[], __int8 W_corre, __int8 WB_corre);

void coup_step_2(int coup_ts, int ite_coup,
	double BHP_r, double *BHP,
	double THP, double L, int nz, double *deltaz,
	double *Qw_total, double *Qo_total, double *Qg_total,
	double *PWB_ref, double *PWB, double *PWB_0, double *dPWB_D, double *PWB_coup,
	double *Qw_wellbore, double *Qo_wellbore, double *Qg_wellbore,
	double *Qw_wellbore_ref, double *Qo_wellbore_ref, double *Qg_wellbore_ref,
	int N_wellbore, double D_perf, double perf_dens,
	double teta, double rug, double gama, double teta_yue,
	double *po_coup, double *So_coup, double *Sw_coup, double *pb_coup,
	double *po_res, double *pw_res, double *pg_res,
	double *So_res, double *Sw_res, double *Sg_res, double *pb_res,
	double *po_res_0, double *So_res_0, double *Sw_res_0, double *pb_res_0, double *Sg_res_0, double *X_pb_0,
	double *phi_res_0, double *Bw_res_0, double *Rs_res_0, double *Bo_res_0, double *Bg_res_0,
	double *po_res_lv0, double *pb_res_lv0, double *So_res_lv0, double *Sw_res_lv0,
	double *varpo, double *varpb, double *varSo, double *varSw,
	double Qwsinj, double Qosinj, double Qgsinj,
	double *K_F, double *K_B, double *K_N, double *K_S, double *K_W, double *K_E,
	double *X_pb,bool *zonaprodutora,
	int res_timesteps, bool dontstopmenow,
	int p, int y, int u, int l, 
	int centerP, int centerSw, int centerSo, 
	int i, int k, int j, 
	int o, int m, int n, int dim,
	double *Kr, double *Kz, double pi, double T,
	int *IA, int *JA, double *A, double *B, double *X,
	double controladorP, double controladorSw, double controladorSo,
	double C1, double C2, double C3,
	double *Ageo_R, double *Ageo_PHI, double *Ageo_Z, double *Ageo_t,
	double *Adir_F, double *Adir_B, double *Adir_W, double *Adir_E, double *Adir_N, double *Adir_S,
	double *CHI_F, double *CHI_B, double *CHI_W, double *CHI_E, double *CHI_N, double *CHI_S,
	double *R_P, double *PHI_P, double *Z_P,
	double *R_PB, double *R_PF, double *delta_PHI, double *delta_Z, double deltat_res,
	double g, double Re, double Rw, 
	double rho_air, double rho_ws, double API, double yo, double rho_os, double yg, double rho_gs,
	double yn2, double yco2, double yh2s,
	double *lambda1_W, double *d_lambda1_po_W, double *d_lambda1_SW_W,
	double *lambda2_W, double *d_lambda2_po_W, double *d_lambda2_SW_W,
	double *lambda1_O, double *d_lambda1_po_O, double *d_lambda1_SW_O, double *d_lambda1_pb_O, double *d_lambda1_SO_O,
	double *lambda2_O, double *d_lambda2_po_O, double *d_lambda2_SW_O, double *d_lambda2_pb_O, double *d_lambda2_SO_O,
	double *lambda1_G, double *d_lambda1_po_G, double *d_lambda1_SW_G, double *d_lambda1_pb_G, double *d_lambda1_SO_G,
	double *lambda2_G, double *d_lambda2_po_G, double *d_lambda2_SW_G, double *d_lambda2_pb_G, double *d_lambda2_SO_G,
	double *lambda3_G, double *d_lambda3_po_G, double *d_lambda3_SW_G, double *d_lambda3_pb_G, double *d_lambda3_SO_G,
	double *lambda4_G, double *d_lambda4_po_G, double *d_lambda4_SW_G, double *d_lambda4_pb_G, double *d_lambda4_SO_G,
	double krwmax, double nw, double ng, double now, double nog, 
	double *phi_ref, double Cr,
	double Swc, double Sor, double Sgr, double epsilon,
	double pcwmin, double pcwmax, double pcgmin, double pcgmax,
	double *parametrosBw, double *parametrosmig, double *parametrosZ,
	double ppc, double Tpc, double Bwi, double Bpcw, double Bpcg, 
	int *iparm, double *dparm, int solver, int maxfct, int mnum, int phase,
	int error, int msglvl, int mtype, int nrhs, double ddum, int idum, void *pt,
	double *phi, double *der_phi, 
	double *rho_w, double *der_rho_w,
	double *mi_w, double *der_mi_w,
	double *Bw, double *der_Bw, 
	double *pcw, double *der_pcw,
	double *krw, double *der_krw,
	double *Rs, double *der_Rs_p, double *der_Rs_pb, 
	double *kro, double *der_kro_Sw, double *der_kro_So, 
	double *Bo, double *der_Bo_p, double *der_Bo_pb,
	double *mi_o, double *der_mi_o_p, double *der_mi_o_pb, 
	double *rho_o, double *der_rho_o_p, double *der_rho_o_pb,
	double *pcg, double *der_pcg_Sw, double *der_pcg_So,
	double *krg, double *der_krg, 
	double *mi_g, double *der_mi_g, double *Z, 
	double *Bg, double *der_Bg,
	double *rho_g, double *der_rho_g,
	int *unW, int *ufW, int *uwbW, int *ueW, int *ucW, int *uwW, int *uebW, int *ubW, int *usW,
	int *unO, int *ufO, int *uwbO, int *ueO, int *ucO, int *uwO, int *uebO, int *ubO, int *usO,
	int *unG, int *ufG, int *uwbG, int *ueG, int *ucG, int *uwG, int *uebG, int *ubG, int *usG,
	int *neg, bool *NpureW, int *contresint,
	char directory[], __int8 W_corre, __int8 WB_corre);

	void coup_step_3(int coup_ts, int ite_coup,
	double BHP_r, double *BHP,
	double THP, double L, int nz, double *deltaz,
	double *Qw_total, double *Qo_total, double *Qg_total,
	double *PWB_ref, double *PWB, double *PWB_0, double *dPWB_D, double *PWB_coup,
	double *Qw_wellbore, double *Qo_wellbore, double *Qg_wellbore,
	double *Qw_wellbore_ref, double *Qo_wellbore_ref, double *Qg_wellbore_ref,
	int N_wellbore, double D_perf, double perf_dens,
	double teta, double rug, double gama, double teta_yue,
	double *po_coup, double *So_coup, double *Sw_coup, double *pb_coup,
	double *po_res, double *pw_res, double *pg_res,
	double *So_res, double *Sw_res, double *Sg_res, double *pb_res,
	double *po_res_0, double *So_res_0, double *Sw_res_0, double *pb_res_0, double *Sg_res_0, double *X_pb_0,
	double *phi_res_0, double *Bw_res_0, double *Rs_res_0, double *Bo_res_0, double *Bg_res_0,
	double *po_res_lv0, double *pb_res_lv0, double *So_res_lv0, double *Sw_res_lv0,
	double *varpo, double *varpb, double *varSo, double *varSw,
	double Qwsinj, double Qosinj, double Qgsinj,
	double *K_F, double *K_B, double *K_N, double *K_S, double *K_W, double *K_E,
	double *X_pb,bool *zonaprodutora,
	int res_timesteps, bool dontstopmenow,
	int p, int y, int u, int l, 
	int centerP, int centerSw, int centerSo, 
	int i, int k, int j, 
	int o, int m, int n, int dim,
	double *Kr, double *Kz, double pi, double T,
	int *IA, int *JA, double *A, double *B, double *X,
	double controladorP, double controladorSw, double controladorSo,
	double C1, double C2, double C3,
	double *Ageo_R, double *Ageo_PHI, double *Ageo_Z, double *Ageo_t,
	double *Adir_F, double *Adir_B, double *Adir_W, double *Adir_E, double *Adir_N, double *Adir_S,
	double *CHI_F, double *CHI_B, double *CHI_W, double *CHI_E, double *CHI_N, double *CHI_S,
	double *R_P, double *PHI_P, double *Z_P,
	double *R_PB, double *R_PF, double *delta_PHI, double *delta_Z, double deltat_res,
	double g, double Re, double Rw,
	double rho_air, double rho_ws, double API, double yo, double rho_os, double yg, double rho_gs,
	double yn2, double yco2, double yh2s,
	double *lambda1_W, double *d_lambda1_po_W, double *d_lambda1_SW_W,
	double *lambda2_W, double *d_lambda2_po_W, double *d_lambda2_SW_W,
	double *lambda1_O, double *d_lambda1_po_O, double *d_lambda1_SW_O, double *d_lambda1_pb_O, double *d_lambda1_SO_O,
	double *lambda2_O, double *d_lambda2_po_O, double *d_lambda2_SW_O, double *d_lambda2_pb_O, double *d_lambda2_SO_O,
	double *lambda1_G, double *d_lambda1_po_G, double *d_lambda1_SW_G, double *d_lambda1_pb_G, double *d_lambda1_SO_G,
	double *lambda2_G, double *d_lambda2_po_G, double *d_lambda2_SW_G, double *d_lambda2_pb_G, double *d_lambda2_SO_G,
	double *lambda3_G, double *d_lambda3_po_G, double *d_lambda3_SW_G, double *d_lambda3_pb_G, double *d_lambda3_SO_G,
	double *lambda4_G, double *d_lambda4_po_G, double *d_lambda4_SW_G, double *d_lambda4_pb_G, double *d_lambda4_SO_G,
	double krwmax, double nw, double ng, double now, double nog, 
	double *phi_ref, double Cr,
	double Swc, double Sor, double Sgr, double epsilon,
	double pcwmin, double pcwmax, double pcgmin, double pcgmax,
	double *parametrosBw, double *parametrosmig, double *parametrosZ,
	double ppc, double Tpc, double Bwi, double Bpcw, double Bpcg, 
	int *iparm, double *dparm, int solver, int maxfct, int mnum, int phase,
	int error, int msglvl, int mtype, int nrhs, double ddum, int idum, void *pt,
	double *phi, double *der_phi, 
	double *rho_w, double *der_rho_w,
	double *mi_w, double *der_mi_w,
	double *Bw, double *der_Bw, 
	double *pcw, double *der_pcw,
	double *krw, double *der_krw,
	double *Rs, double *der_Rs_p, double *der_Rs_pb, 
	double *kro, double *der_kro_Sw, double *der_kro_So, 
	double *Bo, double *der_Bo_p, double *der_Bo_pb,
	double *mi_o, double *der_mi_o_p, double *der_mi_o_pb, 
	double *rho_o, double *der_rho_o_p, double *der_rho_o_pb,
	double *pcg, double *der_pcg_Sw, double *der_pcg_So,
	double *krg, double *der_krg, 
	double *mi_g, double *der_mi_g, double *Z, 
	double *Bg, double *der_Bg,
	double *rho_g, double *der_rho_g,
	int *unW, int *ufW, int *uwbW, int *ueW, int *ucW, int *uwW, int *uebW, int *ubW, int *usW,
	int *unO, int *ufO, int *uwbO, int *ueO, int *ucO, int *uwO, int *uebO, int *ubO, int *usO,
	int *unG, int *ufG, int *uwbG, int *ueG, int *ucG, int *uwG, int *uebG, int *ubG, int *usG,
	int *neg, bool *NpureW, int *contresint,
	char directory[], __int8 W_corre, __int8 WB_corre);

void calc_guess_BHP(int ite_coup, int coup_ts,
	double BHP_r, double *BHP, double *BHP1, double *BHP2, double *indice1, double *indice2,
	double *deltaBHP, bool *controldelta, int *avoidshit, double P20, double THP);

void calc_THP(double prod_time, int coup_ts, bool *controlstartup, double *THP, double THP_i, double THP_adj);