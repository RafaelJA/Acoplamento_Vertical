double calc_q_I(double n, double A_I, double U_IS){
	
	//n perforation density
	//Cross-sectional area of each perforation
	//Velocities in the perforation	
	
	return n*A_I*U_IS;


}


void ouyang_multiphase_model_node(double *H_l, double *dpdx, double H_l_b, double rho_g, double rho_l, double U_sg, double U_sl, double q_I_l, double q_I_g, double mu_l, double mu_g, double D, double g,  double teta){

	
	double V_im=(q_I_l+q_I_g)/(M_PI*D);
	double U_m=U_sg+U_sl;
	double A=0.25l*M_PI*D*D;	
	double sigma=calc_surf_tens();	

	int pattern=flow_pattern_ouyang(H_l_b, mu_g, mu_l, rho_l, rho_g, U_m, g, teta, V_im, D, teta, U_sg);
	
	//SINGLE-PHASE FLOW
	if(pattern==0){
			
	}

	//BUBBLE FLOW
	else if(pattern==1){

		bubble_flow_slip_calc(*H_l, *dpdx, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, q_I_l, q_I_g, sigma, teta, g, D);
		
	}
	
	//STRATIFIED FLOW
	else if(pattern==2){
		
		stratified_flow_calc(*H_l, *dpdx, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, q_I_l, q_I_g, sigma, teta, g, D, A);

	}
	
	//ANNULAR-MIST FLOW
	else if(pattern==3){

		annularmist_flow_calc(*H_l, *dpdx, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, q_I_l, q_I_g, sigma, teta, g, D, A);


	}

	//INTERMITTENT FLOW
	else{

		intermittent_flow_calc(*H_l, *dpdx, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, q_I_l, q_I_g, sigma, teta, g, D, A);
	
		
	}	
		
	
}



double calc_hlD(double H_l){
	

	if((0.0l<H_l&& H_l<= 0.164l) || (0.228l<H_l&& H_l<= 0.337l) ){
		h_lD=0.7769409734*pow(H_l,0.6875921155);
	}
	else if((0.164l<H_l && H_l<=0.228l)|| (0.385l < H_l && H_l< 0.500l)){
		h_lD=-0.6896121749*H_l*H_l+1.2645448768*H_l+0.0286615347;
	} 
	else if(0.337l<H_l && H_l<=0.385l){
		h_lD=0.9441422481*H_l+0.0496745234;
	}
	else{
		h_lD=H_l;
	}

	double teta;
	double H_l_ref=0.5l;

	while(abs((H_l_ref-H_l)/H_l)>1e-10){
		teta=2.0l*acos(1.0l-2*h_lD);
		teta=2.0l*M_PI*H_l+sin(teta);	
		h_lD=0.5*(1.0l-cos(0.5l*teta));
		H_l_ref=0.5l*M_PI*(h_lD-sin(h_lD));
	}

}

double calc_A_G(double hlD, double D){
	double teta=2.0l*acos(1.0l-2*h_lD);
	return 0.25l*M_PI*D*D-0.125l*D*D*(teta-sin(teta));//A-A_L 
}

double calc_d_max_hinze(double U_m, double f_m, double D, double sigma, double rho_l){
	
	double k=1.14l;
	double kappa=2.0l*f_m*U_m*U_m*U_m/D; 

	return k*pow(sigma/rho_l,0.6l)*pow(kappa,-0.4l);


}

double calc_C_d_barry(double mu_g, double mu_l, double rho_g, double rho_l, double V_im, double d){

	double A=1.10535l;
	double alpha=2.5891;
	double beta=0.9879;
	double X=mu_g/mu_l;
	double P=rho_g/rho_l;	
	
	double Z=1.0l+((alpha-2.0l)*sqrt(X*P)+(beta-1.0l)*X*P)/(1+2*sqrt(X*P)+X*P);
		
	double AZ=A*Z;
	double C_1=(8.l/9.l)*(4.l+3.l*X)/(1.l+X);	
	double tau_i=sqrt(C_1-1.4*AZ);//1.4 is OK for the first step, not good, not bad, just OK
	double tau=100000000000;

	while(abs(tau_i-tau)>0.000000001l){
		
		tau=tau_i;
		tau_i=sqrt(C_1-tau*AZ);

	}
	tau=tau_i;
	
	double lambda=1.l/tau;
	double omega=1.l/(3.l*(1.l+X));

	double R_0=rho_l*V_im*d/mu_l;

	return 48.l/(R_0*(1.l+X))*(1.l+1.5l*X)*(sqrt(R_0)+AZ+tau*exp(-lambda*sqrt(R_0)))/(R_0/(1.l+X)+3.l*AZ+3.l*tau*exp(-omega*sqrt(R_0)));

}

bool is_dispersed_bubble_flow(double H_l, double mu_g, double mu_l, double rho_l, double rho_g, double U_m, double g, double teta, double V_im, double D, double sigma){
	
	if(H_l<0.48l){
		return 0;
	}

	double f_m=calc_f(rho_l*U_m*D/mu_l);
	
	double d=calc_d_max_hinze(U_m, f_m, D, sigma, rho_l);
	
	//Ouyang 2002
	double C_id=0.8l;
	double lambda=mu_g/mu_l;
	double A=(3.0l/8.0l)*(rho_l/(rho_l-rho_g))*U_m*U_m*f_m/(g*cos(teta));
	double B=(3.0l*lambda+2.0l)/(lambda+1.0l)*6.0l*C_id*mu_l*V_im/((rho_l-rho_g)*g*cos(teta));

	double d_cb=max(0.5l*(A+pow(A*A+4.0l*B,0.5l),0.5l*(A-pow(A*A+4.0l*B,0.5l));/////CHECAR AQUI


	if(d>d_cb){
		return 0;
	}
	//
	
	//Alves 2020
	double C_id=1.0l;
	
	double d_cb=(3.0/8.0l)*(rho_l/(rho_l-rho_g))/(g*cos(teta))*(f_m*U_m*U_m+2.0l*C_id*V_im*V_im*calc_C_d_barry(mu_g, mu_l, rho_g, rho_l, V_im, d));/////CHECAR


	if(d>d_cb){
		return 0;
	}
	//

	double d_cd=2.0l*pow((0.4l*sigma)/((rho_l-rho_g)*g),0.5l);
	
	if(d>d_cd){
		return 0;
	}

	return 1;

}

bool is_stratified_flow(double H_l, double D, double rho_l, double rho_g, double g, double teta, double V_im, double U_sg){

	double hlD=calc_hlD(H_l);
	double hl=hlD*D;
	double A_g=calc_A_G(hlD, D);
	double dAldhl=2.0l*sqrt(hl*(D-hl));
	double U_g=U_sg*0.25l*M_PI*D*D/A_g;

	double U_g_ref=(1.0l-hlD)*sqrt((A_g/dAldhl)*(rho_l*g*cos(teta)/(rho_g) + 0.5l*V_im*abs(V_im)/(D-hl)));
	

	if(U_g>U_g_ref){
		return 0;
	}
	
	return 1;

	
}

int flow_pattern_ouyang(double H_l, double mu_g, double mu_l, double rho_l, double rho_g, double U_m, double f_m, double g, double teta, double V_im, double D, double teta, double U_sg){

	if(H_l>0.9999999999999l){
		return 0;//SINGLE-PHASE FLOW	
	}

	if(is_dispersed_bubble_flow(H_l, mu_g, mu_l, rho_l, rho_g, U_m, g, teta, V_im, D, sigma)){
		return 1;//BUBBLE FLOW
	}

	if(is_stratified_flow(H_l, D, rho_l, rho_g, g, teta, V_im, U_sg)){
		return 2;//STRATIFIED FLOW
	}

	if(H_l<0.24l){
		return 3;//ANNULAR-MIST FLOW
	}
	
	return 4;//INTERMITTENT FLOW

}

double calc_Re_w(double q_I_l, double q_I_g, double rho_l, double rho_g, double mu_l, double mu_g){
	
	double C_I_l=q_I_l/(q_I_l+q_I_g);
	double rho_im=rho_l*C_I_l+rho_g*(1.0l-C_I_l);
	double mu_im=mu_l*C_I_l+mu_g*(1.0l-C_I_l);
	
	return rho_im*(q_I_l+q_I_g)/(M_PI*mu_im);
	

}

void bubble_flow_slip_calc(double *H_l, double *dpdx, double U_sl, double U_sg, double mu_l, double mu_g, double rho_l, double rho_g, double  q_I_l, double q_I_g, double sigma, double teta, double g, double D){	

	//SLIP CONSTANTS	
	double c_w=1.2l;
	double c_0=c_w-0.2l*sqrt(rho_g/rho_l);
	double U_0=1.53l*sqrt(sqrt(g*(rho_l-rho_g)*sigma)/rho_l)*sin(teta);
	//	
	
	double U_m=U_sg+U_sl;
	
	*H_l=1.0l-U_sg/(c_0*U_m+U_0);
	
	double mu_tp=mu_l**H_l+mu_g*(1.0l-*H_l);
	double rho_tp=rho_l**H_l+rho_g*(1.0l-*H_l);
	double U_tp=(rho_l*U_sl+rho_g*U_sg)/rho_tp;
	double q_I_tp=(rho_l*q_I_l+rho_g*q_I_g)/rho_tp;
	
	double Re_tp=rho_tp*U_tp*D/mu_tp;
	double Re_w=calc_Re_w(q_I_l, q_I_g, rho_l, rho_g, mu_l, mu_g);

	//calc_f_tp(Re_tp, Re_w)

	double tau_w=0.5l*f_tp*rho_tp*U_tp*U_tp;
		
	*dpdx=-(4.0l*tau_w/D+8.0l*rho_tp*U_tp*q_I_tp/(M_PI*D*D)+rho_tp*g*sin(teta));
	
}

void bubble_flow_noslip_calc(double *H_l, double *dpdx, double U_sl, double U_sg, double mu_l, double mu_g, double rho_l, double rho_g, double  q_I_l, double q_I_g, double sigma, double teta, double g, double D, double A){	

		
	double U_m=U_sg+U_sl;
	
	*H_l=U_sl/U_m;
	
	double mu_tp=mu_l**H_l+mu_g*(1.0l-*H_l);
	double rho_tp=rho_l**H_l+rho_g*(1.0l-*H_l);
	double U_tp=(rho_l*U_sl+rho_g*U_sg)/rho_tp;
	double q_I_tp=(rho_l*q_I_l+rho_g*q_I_g)/rho_tp;
	
	double Re_w=calc_Re_w(q_I_l, q_I_g, rho_l, rho_g, mu_l, mu_g);
	double Re_tp=rho_tp*U_tp*D/mu_tp;
	
	///calc_f_tp(Re_tp, Re_w)

	double tau_w=0.5l*f_tp*rho_tp*U_tp*U_tp;
	
	
	*dpdx=-(4.0l*tau_w/D+8.0l*rho_tp*U_tp*q_I_tp/(M_PI*D*D)+rho_tp*g*sin(teta));
	
}

void intermittent_flow_calc(double *H_l, double *dpdx, double U_sl, double U_sg, double mu_l, double mu_g, double rho_l, double rho_g, double  q_I_l, double q_I_g, double sigma, double teta, double g, double D, double A){
	
	//SLIP CONSTANTS	
	double c_w=1.2l;
	double c_0=c_w-0.2l*sqrt(rho_g/rho_l);	
	//	
	
	double U_m=U_sg+U_sl;
	
	//velocity of dispersed bubbles in liquid slug
	double U_b=1.53l*sqrt(sqrt(g*(rho_l-rho_g)*sigma)/rho_l)*sin(teta);
	double U_gdb=c_0*U_m+U_b;
	//

	//translational velocity of the elongated bubbles
	double Bo=(rho_l-rho_g)*g*D*D/sigma;//BOND, NUMBER OF BOND
	double U_DH_inf=(0.54l-1.76l*pow(Bo,-0.56l))*sqrt(g*D*(rho_l-rho_g)/rho_l);
	double U_DV_inf=0.345l*(1.l-exp(0.337l-0.1l*Bo))*sqrt(g*D*(rho_l-rho_g)/rho_l);
	double Re_inf=0.5l*rho_l*(U_DH_inf*cos(teta)+U_DV_inf*sin(teta))*D/mu_l;
	double f_m=min(1.0l, 0.316l*sqrt(Re_inf));
	double U_d=f_m*(U_DH_inf*cos(teta)+U_DV_inf*sin(teta));
	double U_t=c_0*U_m+U_d;
	//

	//liquid fraction in the liquid slug
	double H_ls=1.0l/(1.0l+pow(U_m/8.66l,1.39l));
	
	*H_l=H_ls+(U_gdb*(1.0l-H_ls)-U_sg)/U_t;

	double mu_tp=mu_l**H_l+mu_g*(1.0l-*H_l);
	double rho_tp=rho_l**H_l+rho_g*(1.0l-*H_l);
	double U_tp=(rho_l*U_sl+rho_g*U_sg)/rho_tp;
	double q_I_tp=(rho_l*q_I_l+rho_g*q_I_g)/rho_tp;	
	
	double Re_tp=rho_tp*U_tp*D/mu_tp;
	double Re_w=calc_Re_w(q_I_l, q_I_g, rho_l, rho_g, mu_l, mu_g);
	
	///calc_f_tp(Re_tp, Re_w)
	
	*dpdx=-2.0l*f_tp*rho_tp*U_tp*U_tp/D-2.0l*rho_tp*U_tp*q_I_tp/A-rho_tp*g*sin(teta);	
	

}

double calc_fi_stratified(double Fr_l, double Re_sl, double rho_l, double rho_g, double g, double D, double U_g){
	
	return (0.004l+0.5*Re_sl*10e-6)*pow(Fr_l,1.335l)*(rho_l*g*D/(rho_g*U_g*U_g));

}

void stratified_flow_calc(double *H_l, double *dpdx, double U_sl, double U_sg, double mu_l, double mu_g, double rho_l, double rho_g, double  q_I_l, double q_I_g, double sigma, double teta, double g, double D, double A){
	
	
	double Re_sl=rho_l*U_sl*D/mu_l;
	double Re_sg=rho_g*U_sg*D/mu_g;
	double Re_w=calc_Re_w(q_I_l, q_I_g, rho_l, rho_g, mu_l, mu_g);

	///calc_f_sl(Re_sl,Re_w)
	///calc_f_sg(Re_sg,Re_w)

	double dpdx_sg=-2.0l*f_sg*U_sg*U_sg*rho_g/D;

	double X2=f_sl*U_sl*U_sl*rho_l/(f_sg*U_sg*U_sg*rho_g);
	double Y=(rho_l-rho_g)*g*sin(teta)/dpdx_sg;

	double h_ld=riddermethodforlockhartmartinelli(1.0l, 0.0l, 1, n_threads, 1, D, A, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, Re_sl, Re_sg, Re_w, q_I_l, q_I_g, f_sl, f_sg, dpdx_sg, sigma, g, X2, Y);


	*H_l=0.5l*M_PI*(h_ld-sin(h_ld));//VERIFICAR

	double A_g=calc_A_G(h_ld, D);
	
	double S_i=D*sqrt(1.0l-pow(2.0l*h_ld-1.0l,2.0l));
	double S_g=M_PI*D-S_l;
	
	double D_g=4.0l*A_g/(S_g+S_i);

	double U_g=U_sg*A/A_g;
	double U_g=U_sl*A/(A-A_g);
	double U_i=U_g-U_l;

	double Re_g=rho_g*U_g*D_g/mu_g
	double Fr_l=U_l/sqrt(g**H_l);	
	
	///calc_f_wg(Re_g,Re_w)
	double tau_wg=0.5l*f_wg*rho_g*U_g*U_g;
		
	double f_i= calc_fi_stratified(Fr_l, Re_sl, rho_l, rho_g, g, D, U_g);
	double tau_i=0.5l*f_i*rho_g*U_i*abs(U_i);
	
	*dpdx=(-tau_i*S_i-tau_wg*S_g-rho_g*A_g*g*sin(teta)-2.0l*rho_g*U_g*q_I_g)/A_g;

	
}

double calc_Fe_annularmist(double mu_l, double U_sg, double rho_g, double sigma, double U_sl){
	
	double A=0.735l*pow(mu_l*mu_l*U_sg*U_sg*rho_g/(sigma*sigma*rho_l),0.074l)*pow(U_sg/U_sl,0.2l);
	
	return A/(A+1.0l);
	
}

double calc_fi_annularmist(double f_c, double Re_f, double sigma, double rho_c, double U_c, double D_c){
	
	return 0.24l*f_c*pow(Re_f,0.305l)*pow(sigma/(rho_c*U_c*U_c*D_c), 0.085l);

}


void annularmist_flow_calc(double *H_l, double *dpdx, double U_sl, double U_sg, double mu_l, double mu_g, double rho_l, double rho_g, double  q_I_l, double q_I_g, double sigma, double teta, double g, double D, double A){
	
	
	double Re_sl=rho_l*U_sl*D/mu_l;
	double Re_sg=rho_g*U_sg*D/mu_g;
	double Re_w=calc_Re_w(q_I_l, q_I_g, rho_l, rho_g, mu_l, mu_g);

	///calc_f_sl(Re_sl,Re_w)
	///calc_f_sg(Re_sg,Re_w)

	double dpdx_sg= -2.0l*f_sg*U_sg*U_sg*rho_g/D;

	double X2=f_sl*U_sl*U_sl*rho_l/(f_sg*U_sg*U_sg*rho_g);
	double Y=(rho_l-rho_g)*g*sin(teta)/dpdx_sg;

	double delta_ld=riddermethodforlockhartmartinelli(1.0l, 0.0l, 1, n_threads, 2, D, A, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, Re_sl, Re_sg, Re_w, q_I_l, q_I_g, f_sl, f_sg, dpdx_sg, sigma,  g, X2, Y);
	
	double A_c= A-M_PI*D*D*delta_ld*(1.0l-delta_ld);

	double Fe=calc_Fe_annularmist(mu_l, U_sg, rho_g, sigma, U_sl);	
		
	double H_c=A_c/A;
	double H_lc=Fe*U_sl/(U_sg+Fe*U_sl);

	*H_l=1.0l-H_c*(1.0l-H_lc);

	double S_i=M_PI*D*(1.0l-2.0l*V_D);
		
	double D_c=4.0l*A_c/S_i;
			
	double U_f=U_sl*(1.0l-Fe)*A/(1-A_c);
	double U_c=(U_sg+Fe*U_sl)*A/A_c;
	double U_i=U_c-U_f;

	double Fe=calc_Fe_annularmist(mu_l, U_sg, rho_g, sigma, U_sl);	
		
	double H_lc=Fe*U_sl/(U_sg+Fe*U_sl);
			
	double rho_c=rho_l*H_lc+rho_g*(1.0l-H_lc);
	
	double Re_c=rho_c*U_c*D_c/mu_c;
		
	//calc_f_c(Re_c)
	double f_i= calc_fi_annularmist(f_c, Re_f, sigma, rho_c, U_c, D_c);

	double tau_i=0.5l*rho_c*U_i*abs(U_i);
	
	*dpdx=(-tau_i*S_i-rho_c*A_c*g*sin(teta)-2.0l*rho_c*U_c*q_I_g/(1-Fe))/A_c;
	
}


double riddersmethodforlockhartmartinelli(double V_Dmax, double V_Dmin, int nivel, int n_threads, int pattern, double D, double A,
 double U_sl, double U_sg, double mu_l, double mu_g, double rho_l, double rho_g, double Re_sl, double Re_sg, double Re_w, 
 double  q_I_l, double q_I_g, double f_sl, double f_sg, double dpdx_sg, double sigma, double g, double X2, double Y){
	
	double *V_D_k = (double *)malloc(n_threads * sizeof(double));
	double *G_k = (double *)malloc(n_threads * sizeof(double));
	
	for (int k = 0; k < n_threads; k++) {
	
		V_D_k[k]=k*(V_Dmax-V_Dmin)/(n_threads-1)+V_Dmin;
		G_k[k]=	calc_G(pattern, V_D_k[k], D, A, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, Re_sl, Re_sg, Re_w, q_I_l, q_I_g, f_sl, f_sg, dpdx_sg, sigma,  g, X2, Y);	

	}
	
	int cont=0;
	int froot=-1;
	for (int k = 0; k < n_threads - 1; k++) {


		if (G_k[k] == 0.0l && k != n_threads - 1) {

			if (k == 0) return V_D_k[0];

			return riddermethodforlockhartmartinelli(V_D_k[k], V_Dmin, nivel + 1, n_threads, pattern, D, A, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, Re_sl, Re_sg, Re_w, q_I_l, q_I_g, f_sl, f_sg, dpdx_sg, sigma, g, X2, Y);

		
		}
	
		if ((abs(G_k[k + 1]) / G_k[k + 1]) != (abs(G_k[k]) / G_k[k]) && G_k[k]!=0.0l) {
			
			if(froot==-1){			
				froot=k;
			}
			
			cont++;
		}				
	
	}

	if (nivel>=5 && cont==0 && G_k[n_threads - 1] == 0.0l) return V_D_k[n_threads - 1];

	if(cont==1 && nivel<5)
	{
		if (froot + 1 == n_threads - 1) return riddermethodforlockhartmartinelli(V_D_k[froot + 1], V_Dmin, nivel + 1, 2*n_threads, pattern, D, A, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, Re_sl, Re_sg, Re_w, q_I_l, q_I_g, f_sl, f_sg, dpdx_sg, sigma, g, X2, Y);
		
		else return riddermethodforlockhartmartinelli(V_D_k[froot+1], V_Dmin, nivel+1, n_threads, pattern, D, A, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, Re_sl, Re_sg, Re_w, q_I_l, q_I_g, f_sl, f_sg, dpdx_sg, sigma,  g, X2, Y);

	}
	else if (cont == 0) {

		return riddermethodforlockhartmartinelli(V_Dmax, V_Dmin, nivel + 1, 2*n_threads, pattern, D, A, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, Re_sl, Re_sg, Re_w, q_I_l, q_I_g, f_sl, f_sg, dpdx_sg, sigma, g, X2, Y);
	
	}
	else{
		double V_D_l=V_D_k[froot];
		double V_D_h=V_D_k[froot+1];
		double G_l=G_k[froot];
		double G_h=G_k[froot+1];
		double inte=0;
		
		while(1){
		
			double V_D_m=0.5l*(V_D_h+V_D_l);
			double G_m=calc_G(pattern, V_D_m, D, A, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, Re_sl, Re_sg, Re_w, q_I_l, q_I_g, f_sl, f_sg, dpdx_sg, sigma,  g, X2, Y);

			if(G_m==0.0l) return V_D_m;

			double s=sqrt(G_m*G_m-G_l*G_h);
			double V_D_new=V_D_m+(V_D_m-V_D_l)*(G_l-G_h)*G_m/(abs(G_l-G_h)*s);
			double G_new=calc_G(pattern, V_D_new, D, A, U_sl, U_sg, mu_l, mu_g, rho_l, rho_g, Re_sl, Re_sg, Re_w, q_I_l, q_I_g, f_sl, f_sg, dpdx_sg, sigma,  g, X2, Y);
	
			if(G_new==0.0l) return V_D_new;
			
			if(abs(G_new)/G_new!=abs(G_m)/G_m){
				V_D_l=V_D_m;
				G_l=G_m;
				V_D_h=V_D_new;
				G_h=G_new;		
			}
			
			else if(abs(G_new)/G_new!=abs(G_l)/G_l)){
				V_D_h=V_D_new;	
				G_h=G_new;		
			}
			else if(abs(G_new)/G_new!=abs(G_h)/G_h)){
				V_D_l=V_D_new;
				G_l=G_new;			
			}


			if(abs(V_D_l-V_D_h))<0.00000001||inte>100) return 0.5l*(V_D_h+V_D_l);

			inte++;

		}
	}

}

double calc_G(int pattern, double V_D, double D, double A, double U_sl, double U_sg, double mu_l, double mu_g, double rho_l, double rho_g, double Re_sl, double Re_sg,  double Re_w, double  q_I_l, double q_I_g, double f_sl, double f_sg, double dpdx_sg, double sigma, double g, double X2, double Y){
	
	//Stratified Flow V_D=h_ld
	if(pattern==1){
		

		double H_l=0.5l*M_PI*(V_D-sin(V_D));//checar
		double A_g=calc_A_G(V_D, D);
		double A_l= A-A_g;

		double S_i=D*sqrt(1.0l-pow(2.0l*V_D-1.0l,2.0l));
		double S_l=D*(M_PI-acos(2.0l*V_D-1.0l));
		double S_g=M_PI*D-S_l;
	
		double D_l=4.0l*A_l/S_l;
		double D_g=4.0l*A_g/(S_g+S_i);

		double U_l=U_sl*A/A_l;
		double U_g=U_sg*A/A_g;
		double U_i=U_g-U_l;

		double Re_l=rho_l*U_l*D_l/mu_l;
		double Re_g=rho_g*U_g*D_g/mu_g;
		double Fr_l=U_l/sqrt(g*H_l);	
	
		///calc_f_wl(Re_l,Re_w)
		///calc_f_wg(Re_g,Re_w)		
		double f_i= calc_fi_stratified(Fr_l, Re_sl, rho_l, rho_g, g, D, U_g);

		double I=2.0l*(rho_l*U_l*q_I_l/A_l - rho_g*U_g*q_I_g/A_g)/(dpdx_sg);
		double F_1=(f_wg/f_sg)*(U_g/U_sg)*(U_g/U_sg)*(D*S_g/A_g + (f_i/f_wg)*U_i*abs(U_i)*D*S_i/(U_g*U_g)*(1.0l/A_l + 1.0l/A_g));
		double F_2=(f_wl/f_sl)*(U_l/U_sl)*(U_l/U_sl)*(D*S_l/A_l);

		return X2*F_2-F_1-4.0l*(Y+I);

	}
	

	//Annular Mist Flow V_D=delta_ld
	if(pattern==2){
	
		double A_f=M_PI*D*D*V_D*(1.0l-V_D);
		double A_c= A-A_f;

		double S_l=M_PI*D;	
		double S_i=M_PI*D*(1.0l-2.0l*V_D);
		
		double D_f=4.0l*A_f/S_l;
		double D_c=4.0l*A_c/S_i;
			
		double U_f=U_sl*(1.0l-Fe)*A/A_f;
		double U_c=(U_sg+Fe*U_sl)*A/A_c;
		double U_i=U_c-U_f;

		double Fe=calc_Fe_annularmist(mu_l, U_sg, rho_g, sigma, U_sl);	
		
		double H_c=A_c/A;
		double H_f=A_f/A;
		double H_lc=Fe*U_sl/(U_sg+Fe*U_sl);
			
		double rho_c=rho_l*H_lc+rho_g*(1.0l-H_lc);
		double mu_c=mu_l*H_lc+mu_g*(1.0l-H_lc);

		double Re_f=rho_l*U_f*D_f/mu_l;
		double Re_c=rho_c*U_c*D_c/mu_c;
		
		//calc_f_wl(Re_f,Re_w)
		//calc_f_c(Re_c)
		double f_i= calc_fi_annularmist(f_c, Re_f, sigma, rho_c, U_c, D_c);


		double I=(2.0l*(rho_l*U_l/A_l)*(q_I_l-Fe*q_I_g/(1.0l-Fe))-2.0l*rho_c*U_c*q_I_g/(A_c*(1.0l-Fe)))/dpdx_sg;
		double F_1=(f_i/f_sg)*(rho_c/rho_g)*(U_c/U_sg)*(U_c/U_sg)*D*S_i*(1.0l/A_l + 1.0l/ A_g);
		double F_2=(f_wl/f_sl)*(U_f/U_sl)*(U_f/U_sl)*(D*S_l/A_l);

		return X2*F_2-F_1-4.0l*(Y+I);
	}


}