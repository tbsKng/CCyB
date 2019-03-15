R_Ks = x(1);
Ls   = x(2);
L_ms = x(3)/n_m;
x_es = x(4);
x_ms = x(5);           		 
omega_bar_Hs=x(6);
omega_bar_Fs=x(7);
R_Ds=x(8);
Ds=x(9);


%%
%%%%%**********************************************************************   
% 1) A first set of endogenous variables has close form solution and 
% can be easily derived from a subset of equations
%%%%%**********************************************************************
R_DDs = (1/betta_s);  
R_Hs = 1-delta_H;     
q_Ks=1;               
q_Hs = 1;             
As=1;                 
g_Is= 0;   
g_I_1s=0 ; 
PIs=0;     
PHs=0;     
g_Hs= 0;   
g_H_1s=0 ; 

%%
sigma_m1_comp=sigma_m1;
sigma_e1_comp=sigma_e1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Structure of the file
%   (1) Guess values for the 9 unknown variables above
%   (2) Use model relationships (FOCs and mkt clearing conditions) to compute other model variables
%       (for explanations of model relationships used here, see comments to ThreeD_steadystate.m)
%   (3) Solve the system of equations consisting of the remaining 9 model relationships. 
%       Think of these equations as equilibrium 'residuals': the solver adjustes the initial guesses
%       for 9 variables above so that these 'residuals' become zero (a solution is thus found).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gamma_Hs=normcdf((log(omega_bar_Hs)-sigma_H^2/2)/sigma_H) + omega_bar_Hs*(1-normcdf((log(omega_bar_Hs)+sigma_H^2/2)/sigma_H));
Gamma_Fs=normcdf((log(omega_bar_Fs) - sigma_F^2/2)/sigma_F) + omega_bar_Fs*(1-normcdf((log(omega_bar_Fs)+sigma_F^2/2)/sigma_F));
G_Hs   = normcdf((log(omega_bar_Hs)-sigma_H^2/2)/sigma_H);
Gamma_H_1s = (normpdf((log(omega_bar_Hs)-sigma_H^2/2)/sigma_H)-omega_bar_Hs*normpdf((log(omega_bar_Hs)+sigma_H^2/2)/sigma_H))/(sigma_H*omega_bar_Hs) + (1-normcdf((log(omega_bar_Hs)+sigma_H^2/2)/sigma_H));
def_rate_Hs = normcdf((log(omega_bar_Hs)+sigma_H^2/2)/sigma_H)*100;
G_Fs   = normcdf((log(omega_bar_Fs)-sigma_F^2/2)/sigma_F);
Gamma_F_1s = (normpdf((log(omega_bar_Fs)-sigma_F^2/2)/sigma_F)-omega_bar_Fs*normpdf((log(omega_bar_Fs)+sigma_F^2/2)/sigma_F))/(sigma_F*omega_bar_Fs) + (1-normcdf((log(omega_bar_Fs)+sigma_F^2/2)/sigma_F));
def_rate_Fs = normcdf((log(omega_bar_Fs)+sigma_F^2/2)/sigma_F)*100;
%%-----------------------------------------------------------------
r_Ks = R_Ks - (1-delta_K);
YKs = r_Ks/alphaa;
Ks = ((YKs)/(Ls^(1-alphaa)))^(1/(alphaa-1)); 
s_K=(betta_s*(r_Ks+(1-delta_K))-1);
K_ss=(s_K/xi_K)^(1/(phi_K-1));
z_K=xi_K/phi_K*K_ss^phi_K;
PI_K=s_K*K_ss - z_K;
K_es=Ks-K_ss;
Is = delta_K*Ks;
%%-------------------------------------------------------------------------
%% ENTREPRENEURS
%%-------------------------------------------------------------------------
Ys = As*Ks^(alphaa)*Ls^(1-alphaa);
ws = (1-alphaa)*Ys/Ls;
omega_bar_es = x_es/R_Ks;
G_e_1s = normpdf((log(omega_bar_es)-sigma_e1_comp^2/2)/sigma_e1_comp )/(sigma_e1_comp *omega_bar_es);
G_es = normcdf((log(omega_bar_es)-sigma_e1_comp^2/2)/sigma_e1_comp );
Gamma_es   = normcdf((log(omega_bar_es) - sigma_e1_comp^2/2)/sigma_e1_comp ) + omega_bar_es*(1-normcdf((log(omega_bar_es)+sigma_e1_comp^2/2)/sigma_e1_comp ));
Gamma_e_1s = (normpdf((log(omega_bar_es)-sigma_e1_comp ^2/2)/sigma_e1_comp )-omega_bar_es*normpdf((log(omega_bar_es)+sigma_e1_comp ^2/2)/sigma_e1_comp ))/(sigma_e1_comp *omega_bar_es) + (1-normcdf((log(omega_bar_es)+sigma_e1_comp ^2/2)/sigma_e1_comp ));

chi_entrs=xi_entr*((1-Gamma_es)*R_Ks*K_es);


equity_entrs=(theta_entr*((1-Gamma_es)*R_Ks*K_es)+chi_entrs*(1-theta_entr)*entr_pop);

rho_es=(1-Gamma_es)*R_Ks*K_es/equity_entrs;
vi_entrs=betta_s*(1-theta_entr)*rho_es/(1-betta_s*rho_es*theta_entr);
b_es=q_Ks*K_es-equity_entrs;
R_Fs=omega_bar_es/(b_es/((r_Ks+(1-delta_K)*q_Ks)*K_es));




%% -------------------------------------------------------------------------
%%BANKS
%%-------------------------------------------------------------------------
def_rate_es = normcdf((log(omega_bar_es)+sigma_e1^2/2)/sigma_e1)*100;
% if a_e==0
% phi_F=(0.45*normcdf((norminv(def_rate_es/100*4)+sqrt((0.12*(2-(1-exp(-50*(def_rate_es/100*4)))/(1-exp(-50)))))*3.090232306167814)/(sqrt(1-(0.12*(2-(1-exp(-50*(def_rate_es/100*4)))/(1-exp(-50)))))))-(def_rate_es/100*4)*0.45)*(1-1.5*((0.11852-0.05478*log((def_rate_es/100*4)))^2))^(-1)*(1+(1-2.5)*((0.11852-0.05478*log((def_rate_es/100*4)))^2));
% else
phi_F=phi_Fs;
%end

rho_b_goods_Fs = (1-Gamma_Fs)*(Gamma_es-mu_e*G_es)*(r_Ks+(1-delta_K)*q_Ks)*K_es/(phi_F*b_es);
rho_b_goods_Hs = rho_b_goods_Fs;
vi_banks=betta_s*(1-theta_bank)*rho_b_goods_Fs/(1-betta_s*theta_bank*rho_b_goods_Fs);

lamb_b=betta_s*((1-theta_bank)+theta_bank*vi_banks);
rho_b_util_Fs=rho_b_goods_Fs*lamb_b;
rho_b_util_Hs=rho_b_goods_Hs*lamb_b;

R_tilde_Fs=rho_b_util_Fs*phi_F /( lamb_b*(1-Gamma_Fs));

%%-------------------------------------------------------------------------
%% HOUSEHOLDS
%%-------------------------------------------------------------------------
%%-------------------------------------------------------------------------
%% BORROWERS
%%-------------------------------------------------------------------------
omega_bar_ms = x_ms/R_Hs;
def_rate_ms = normcdf((log(omega_bar_ms)+sigma_m1^2/2)/sigma_m1)*100;
% if a_e==0
% phi_H=(0.45*normcdf((norminv(def_rate_ms/100*4)+sqrt(0.15)*3.090232306167814)/(sqrt(1-0.15)))-(def_rate_ms/100*4)*0.45);
% else
phi_H=phi_Hs;
%end
R_tilde_Hs=rho_b_util_Hs*phi_H /( lamb_b*(1-Gamma_Hs));


G_m_1s = normpdf((log(omega_bar_ms)-sigma_m1_comp^2/2)/sigma_m1_comp)/(sigma_m1_comp*omega_bar_ms);
Gamma_ms   = normcdf((log(omega_bar_ms)- sigma_m1_comp^2/2)/sigma_m1_comp) + omega_bar_ms*(1-normcdf((log(omega_bar_ms)+sigma_m1_comp^2/2)/sigma_m1_comp));
G_ms   = normcdf((log(omega_bar_ms)-sigma_m1_comp^2/2)/sigma_m1_comp);
Gamma_m_1s = (normpdf((log(omega_bar_ms)-sigma_m1_comp^2/2)/sigma_m1_comp)-omega_bar_ms*normpdf((log(omega_bar_ms)+sigma_m1_comp^2/2)/sigma_m1_comp))/(sigma_m1_comp*omega_bar_ms) + (1-normcdf((log(omega_bar_ms)+sigma_m1_comp^2/2)/sigma_m1_comp));
R_ms=R_tilde_Hs*(x_ms)/((Gamma_ms - mu_m*G_ms)*R_Hs);
UL_m_1s=varphi_m*L_ms^(eta);
Lambda_ms=UL_m_1s/ws;
UC_m_1s=Lambda_ms;
C_ms=(1/UC_m_1s)/(1-hab);
xi_ms=(betta_m*Lambda_ms*Gamma_m_1s)/((1-Gamma_Hs)*(betta_s*((1-theta_bank)+theta_bank*vi_banks))*(Gamma_m_1s - mu_m*G_m_1s));
ZZHms =v_m/((Lambda_ms-betta_m*Lambda_ms*((1-Gamma_ms)*R_Hs+Gamma_m_1s*omega_bar_ms*R_Hs)-xi_ms*(1-Gamma_Hs)*(betta_s*((1-theta_bank)+theta_bank*vi_banks))*((Gamma_ms - mu_m*G_ms)*R_Hs-(Gamma_m_1s- mu_m*G_m_1s)*omega_bar_ms*R_Hs)));

%-------------------------------------------------------------------------
%% SAVERS
%%-------------------------------------------------------------------------
L_ss=(Ls-L_ms*n_m)/n_s;
UL_s_1s=varphi_s*L_ss^(eta);
Lambda_ss=UL_s_1s/ws;
UC_s_1s = Lambda_ss;
C_ss=(1/UC_s_1s)/(1-hab);
Cs = C_ss*n_s + C_ms*n_m;
ZZHss =v_s/(Lambda_ss*(1-(1-delta_H)*betta_s));
q_Hs = 1;
Hs=(n_s*ZZHss+n_m*ZZHms);
IHs = delta_H*Hs;
H_ss=ZZHss;
H_ms=ZZHms;
b_ms=x_ms*(H_ms*q_Hs)/R_ms;
%%-------------------------------------------------------------------------

%% -------------------------------------------------------------------------
%%BANKS
%%-------------------------------------------------------------------------

equity_bank_Fs=phi_F*b_es;
equity_bank_Hs=phi_H*b_ms*n_m;
chi_banks=xi_bank*(rho_b_goods_Fs*equity_bank_Fs+rho_b_goods_Hs*equity_bank_Hs);
equity_banks=theta_bank*(rho_b_goods_Fs*equity_bank_Fs+rho_b_goods_Hs*equity_bank_Hs)+(1-theta_bank)*chi_banks*bank_pop;

Tr_Hs = ((omega_bar_Hs - Gamma_Hs + mu_H*G_Hs)*R_tilde_Hs*(n_m*H_ms*q_Hs*x_ms)/R_ms);
Tr_Fs = ((omega_bar_Fs - Gamma_Fs + mu_F*G_Fs)*R_tilde_Fs*(b_es));
Trs =Tr_Fs + Tr_Hs;
%%-------------------------------------------------------------------------
xi_es=-((1-Gamma_es)*R_Ks-Gamma_e_1s*R_Ks*K_es*((R_Fs*equity_entrs)/(q_Ks*K_es^2*R_Ks)))*(betta_s*(((1-theta_entr)+theta_entr*vi_entrs)))/((1-Gamma_Fs)*(betta_s*((1-theta_bank)+theta_bank*vi_banks))*(Gamma_es- mu_e*G_es) *R_Ks+(1-Gamma_Fs)*(betta_s*((1-theta_bank)+theta_bank*vi_banks))*(Gamma_e_1s- mu_e*G_e_1s)*((R_Fs*equity_entrs)/(q_Ks*K_es^2*R_Ks))*R_Ks*K_es - rho_b_util_Fs*phi_F);
omega_bar_ms= x_ms/R_Hs; 
PHI_m = (1-Gamma_ms)*R_Hs*q_Hs*H_ms;
av_defs= ((1-phi_F)*equity_bank_Fs*def_rate_Fs*4/phi_F + (1-phi_H)*equity_bank_Hs*def_rate_Hs*4/phi_H)/(Ds*n_s);


%%-------------------------------------------------------------------------
% All variables above are derives taking the initial 9 variables as given/known and the variables in 1.
% Now we find those variables by solving a system of 9 equations in 9 unknowns:
%%-------------------------------------------------------------------------
% System of 10 equations in 10 unknowns(to be set to zero in SS model solution)
z(1) = equity_banks   - (equity_bank_Fs+equity_bank_Hs)    ;                                     % equity of bankers
z(2) = n_m*C_ms + n_m*q_Hs*H_ms - n_m*ws*L_ms - n_m*b_ms - n_m*PHI_m+kappa_ins*Trs*a_b; % Budget constraint (borrowers)
z(3) = equity_banks + Ds - b_es - b_ms*n_m;
z(4) = -xi_es + Gamma_e_1s*(betta_s*(((1-theta_entr)+theta_entr*vi_entrs)))/((1-Gamma_Fs)*(betta_s*((1-theta_bank)+theta_bank*vi_banks))*(Gamma_e_1s- mu_e*G_e_1s));        % K FOC in entrepreneur's problem 
z(5) = Lambda_ms-betta_m*Lambda_ms*Gamma_m_1s*omega_bar_ms/b_ms*R_Hs*H_ms  -xi_ms*(rho_b_util_Hs*phi_H-(1-Gamma_Hs)*(betta_s*((1-theta_bank)+theta_bank*vi_banks))*(Gamma_m_1s - mu_m*G_m_1s)*omega_bar_ms/b_ms*R_Hs*H_ms); % FOC in household problem
z(6)=omega_bar_Hs-((1-phi_H)*R_Ds)/R_tilde_Hs;                              % Default cut off for mortgage bank
z(7)=omega_bar_Fs- ((1-phi_F)*R_Ds)/R_tilde_Fs;                             % Default cut off for corporate bank
z(8)=R_Ds-((1/betta_s)+(1-kappa_ins)*Trs/Ds)/((1-av_defs/400)+av_defs/400);
z(9)=Ds -(ws*L_ss+ (1-theta_entr)*(rho_es*equity_entrs-chi_entrs*entr_pop)+(1-theta_bank)*(rho_b_goods_Fs*equity_bank_Fs+rho_b_goods_Hs*equity_bank_Hs-chi_banks*bank_pop) -((1-kappa_ins)*Trs+kappa_ins*Trs*a_s) -(K_ss+s_K)+(r_Ks+(1-delta_K))*K_ss- delta_H*H_ss -C_ss + PIs+PHs+PI_K)/(-((1-av_defs/400)*R_Ds+av_defs/400*R_Ds-1));

% By solving this system of equations we get to know the 9 unknowns and then 
% replacing those into the definitions above we find the steady state values of all other variables.

