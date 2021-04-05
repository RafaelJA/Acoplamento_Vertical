//Contem todos os calculos para as propriedades do modelo Black Oil
//
#pragma once

//retorna o ponto de bolha para dada temperatura
double Bubble_Point(double T, double API, double SG_G, double GOR);

//calculo da tens�o superficial da mistura
double Tensuper(double rhol, double rhog);

//calculo da entalpia e da capacit�ncia t�rmica do liquido
void H_CP_Black(double T, double API, double *Hl, double *Cp);

//calculo da Rs com base no titulo m�ssico e do GOR
double Rs_Quality(double GOR, double x, double rhol_st, double rhog_st);

//calculo do titulo com base nas condi��es std
double Quality_Rs(double GOR, double Rs, double rhol_st, double rhog_st);

//termos de capacidade termica e entalpia
void water_H_CP(double P, double T, double *Hlw, double *Hlgw, double *Cplw, double *Cpgw);

//tens super agua e vapor
double tens_water(double P, double T);

double Z_factor(double P, double T, double SG_G);