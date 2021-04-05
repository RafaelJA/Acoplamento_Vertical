#pragma once

//entrando com a press�o em Pa retorna a press�o em Psia
double Pa_Psia(double P_Pa);

//entrando com a press�o em Psia retorna a press�o em Pa
double Psia_Pa(double P_Psia);

//entrando com a temperatura em Kelvin retorna a mesma em Far
double Kel_Far(double T_Kel);

//entrando com a temperatura em Farenheit retorna a mesma em Kelvin
double Far_Kel(double T_Far);

//Kelvin to Rankine
double Kel_Rank(double T_Kel);

//entrando com o comprimento em metros retorna o mesmo em p�s
double m_ft(double L_m);

//entrando com o comprimento em p�s retorna o mesmo em metros
double ft_m(double L_ft);

//entrando com o comprimento em polegadas retorna o mesmo em p�s
double in_ft(double L_in);

//entrando com o comprimento em p�s retorna o mesmo em metros
double ft_in(double L_ft);

//ft� to m�
double ft3_m3(double V_ft3);

//m� to ft�
double m3_ft3(double V_m3);

//entrando com o comprimento em p�s cubicos retorna o mesmo em barris americanos
double  ft3_bbl(double V_ft3);

//entrando com o comprimento em barris americanos retorna o mesmo em p�s c�bicos
double bbl_ft3(double V_bbl);

//ft�/s to bbl/d
double ft3s_bbld(double Q_ft3s);

//bbl/d to ft�/s
double bbld_ft3s(double Q_bbld);

//cP to Pa.s
double cP_Pas(double mu_cP);

//Pa.s to cP
double Pas_cP(double mu_Pas);

//bbl/d to m�/s
double bbld_m3s(double Q_bbld);

//m�/s to bbl/d 
double m3s_bbld(double Q_m3s);

//lb/ft� to kg/m�
double lbft3_kgm3(double rho_lbft3);

//kg/m� to lb/ft� 
double kgm3_lbft3(double rho_kgm3);