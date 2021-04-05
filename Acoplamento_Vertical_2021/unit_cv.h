#pragma once

//entrando com a pressão em Pa retorna a pressão em Psia
double Pa_Psia(double P_Pa);

//entrando com a pressão em Psia retorna a pressão em Pa
double Psia_Pa(double P_Psia);

//entrando com a temperatura em Kelvin retorna a mesma em Far
double Kel_Far(double T_Kel);

//entrando com a temperatura em Farenheit retorna a mesma em Kelvin
double Far_Kel(double T_Far);

//Kelvin to Rankine
double Kel_Rank(double T_Kel);

//entrando com o comprimento em metros retorna o mesmo em pés
double m_ft(double L_m);

//entrando com o comprimento em pés retorna o mesmo em metros
double ft_m(double L_ft);

//entrando com o comprimento em polegadas retorna o mesmo em pés
double in_ft(double L_in);

//entrando com o comprimento em pés retorna o mesmo em metros
double ft_in(double L_ft);

//ft³ to m³
double ft3_m3(double V_ft3);

//m³ to ft³
double m3_ft3(double V_m3);

//entrando com o comprimento em pés cubicos retorna o mesmo em barris americanos
double  ft3_bbl(double V_ft3);

//entrando com o comprimento em barris americanos retorna o mesmo em pés cúbicos
double bbl_ft3(double V_bbl);

//ft³/s to bbl/d
double ft3s_bbld(double Q_ft3s);

//bbl/d to ft³/s
double bbld_ft3s(double Q_bbld);

//cP to Pa.s
double cP_Pas(double mu_cP);

//Pa.s to cP
double Pas_cP(double mu_Pas);

//bbl/d to m³/s
double bbld_m3s(double Q_bbld);

//m³/s to bbl/d 
double m3s_bbld(double Q_m3s);

//lb/ft³ to kg/m³
double lbft3_kgm3(double rho_lbft3);

//kg/m³ to lb/ft³ 
double kgm3_lbft3(double rho_kgm3);