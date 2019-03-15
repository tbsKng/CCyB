//************************************************************
// This file solves the model as in "Optimal dynamic capital requirements"
// by  C. Mendicino,  K. Nikolov, J. Suarez and D. Supera, 
// forthcoming in Journal of Money, Credit, and Banking.
//************************************************************

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//1. Preamble 
//The preamble consists of the some declarations to setup the endogenous
//and exogenous variables, the parameters and assign values to these parameters.
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

//***************************************************
// DECLARATION ENDOGENOUS VARIABLES IN THE MODEL
//***************************************************
var b_e    // entrepreneurial debt
g_H        // housing investment adjustment cost
g_H_1      // first derivative of g_H with respect to I_H
IH         // housing investment 
PH         // Profit of housing capital producing firm
A          // productivity
b_m        // mortgage debt
C          // aggregate consumption
C_m        // consumption of borrowers (impatient households)
C_s        // consumption of savers (patient households)
D          // aggregate deposits
def_rate_e // corporate default rate
def_rate_m // mortgage default rate
def_rate_F // default rate of corporate banks
def_rate_H // default rate of mortgage banks
G_e        // share of entrepreneurial capital belonging to firms that default (BGG parameter)
G_e_1      // First derivative of G_e with respect to omega_e (BGG parameter)
G_F        // share of corporate loans belonging to corporate banks that default (BGG parameter)
G_H        // share of mortgage loans belonging to mortgage banks that default (BGG parameter)
g_I        // Capital investment adjustment cost 
g_I_1      // First derivative of g_I with respect to I_K
G_m        // share of housing belonging to households that default (BGG parameter)
G_m_1      // First derivative of G_m with respect to omega_m (BGG parameter)
Gamma_e    // Share of gross corporate revenues going to the bank (BGG parameter)
Gamma_e_1  // First derivative of Gamma_e with respect to omega_e(BGG parameter)
Gamma_F    // Share of gross corporate bank revenues going to depositors(BGG parameter)
Gamma_F_1  // First derivative of Gamma_F with respect to omega_F(BGG parameter)
Gamma_H    // Share of gross mortgage bank revenues going to depositors(BGG parameter)
Gamma_H_1  // First derivative of Gamma_H with respect to omega_H(BGG parameter)
Gamma_m    // Share of gross returns from housing going to the bank(BGG parameter)
Gamma_m_1  // First derivative of Gamma_m with respect to omega_m(BGG parameter)
H          // Aggregate housing supply
H_m        // Housing used by borrowers
H_s        // Housing used by savers
I          // Aggregate corporate investment
K          // Aggregate capital stock 
L          // Aggregate labour supply
L_m        // Borrowers' labour supply
L_s        // Savers' labour supply       
Lambda_m   // LM on the budget constraint in the household borrower's problem
Lambda_s   // LM on the budget constraint in the household saver's problem
omega_bar_e	// Idiosyncratic productivity shock below which the corporate borrower defaults
omega_bar_F // Idiosyncratic corporate bank loan return shock below which the corporate bank defaults
omega_bar_H	// Idiosyncratic mortgage bank loan return shock below which the mortgage bank defaults
omega_bar_m // Idiosyncratic housing return shock below which the mortgage borrower defaults
PI          // Profit of business capital producing firm 
q_H         // Housing price
q_K         // Capital price
R_D         // Deposit interest rate
R_F         // Corporate loan interest rate
R_H         // Aggregate financial return on housing (i.e. excluding imputed rents)
r_K         // Capital rental rate
R_K         // Capital rate of return
R_m         // Mortgage interest rate
R_tilde_F 	// Aggregate return on a diversified corporate loan portfolio (i.e. portfolio return after accounting for loan losses)
R_tilde_H 	// Aggregate return on a diversified housing loan portfolio (i.e. portfolio return after accounting for loan losses)
Tr          // Total deposit insurance payments to banks
Tr_F        // deposit insurance payments to corporate banks
Tr_H        // deposit insurance payments to mortgage banks
UC_m        // Borrowers' utility from consumption
UC_m_1      // Borrowers' marginal utility from consumption
UC_s        // Savers' utility from consumption
UC_s_1      // Savers' marginal utility from consumption
UH_m        // Borrowers' utility from housing services
UH_m_1      // Borrowers' marginal utility from housing services
UH_s        // Savers' utility from housing services
UH_s_1      // Savers' marginal utility from housing services
UL_m        // Borrowers' disutility from work
UL_m_1      // Borrowers' marginal disutility from work
UL_s        // Savers' disutility from work
UL_s_1      // Savers' marginal disutility from work
Util_m      // Total borrower utility
Util_s      // Total saver utility
w           // Wage rate
x_e         // Corporate leverage
x_m         // Household leverage
xi_e        // LM on bank's participation constraint in the entrepreneurs' problem
xi_m        // LM on bank's participation constraint in the household borrower's problem
Y           // Output
Y_obs       // All _obs variables are transformations to aid plotting. Transformations may differ across variables depending on what makes sense.
R_D_obs     // Please look at particular _obs variable definition in the code to see the exact transformation.
R_m_obs
R_H_obs
R_F_obs
H_m_obs
H_s_obs
b_m_obs
C_obs
C_m_obs
C_s_obs
D_obs
E_F_obs
I_obs
K_obs
L_obs
L_m_obs
L_s_obs
n_b_obs
n_e_obs
q_H_obs
q_K_obs
r_K_obs
R_K_obs
R_tilde_F_obs
R_tilde_H_obs
rho_F_obs
rho_H_obs
Tr_obs
Tr_F_obs
Tr_H_obs
w_obs
x_e_obs
x_m_obs
Vs           // Value function of household savers
Vm           // Value function of household borrowers
EJ           // Shock (housing preference)
ESe          // Shock (entrepreneur risk)
ESm          // Shock (housing risk)
ESH          // Shock (mortgage bank risk)
ESF          // Shock (corporate bank risk)
EdH          // Shock (housing depreciation)
EdK          // Shock (capital depreciation)
phi_F        // Minimum capital ratio for corporate banks
phi_H        // Minimum capital ratio for mortgage banks
b_e_obs
bsp_F        // Corporate loan return spread
bsp_H        // Household loan return spread
deltaH       // Housing depreciation rate
deltaK       // Capital depreciation rate
H_obs
IH_obs
av_def       // Average default for all banks
RDsp         // Deposit rate spread over the risk free rate
R_DD         // Deposit rate return net of bank default costs
W_m          // Household borrower net worth
W_m_obs
b_tot        // Total debt in the economy
b_tot_obs
A_news 
EJ_news 
ESm_news 
ESe_news 
ESF_news 
EdH_news 
EdK_news // News shocks
HW           // Housing wealth
HW_obs   
GDP_acc      // GDP as in the data: C + I + IH
GDP_acc_obs  
Spread_e     // Corporate loan spread
Spread_m     // Mortgage spead
write_off_m_new  // Write offs on mortgage loans
write_off_e_new // Write offs on corporate loans
K_e  // Capital entrepreneurs
K_s // Capital savers
PI_K //Profits capital management firms
s_K //Revenues (per unit of capital) capital management firms
z_K //Costs capital management firms
rho_e //return on equity of entrepreneurs
equity_entr //equity of entrepreneurs
chi_entr  //transfers from households to entrepreneurs
rho_b_goods_F //return on equity of corporate bank 
rho_b_goods_H //return on equity of mortgage bank 
rho_b_util_F //return on equity of corporate bank in utility terms
rho_b_util_H  //return on equity of mortgage bank in utility terms
equity_bank_F //equity of corporate bank
equity_bank_H //equity of mortgage bank
chi_bank //tranfers from households to bankers
equity_bank //equity of bankers
vi_entr //shadow value of equity entrepreneurs
vi_bank //shadow value of equity bankers
sav_cap //capital of savers/total capital
RW
R_DD_obs
phi_bar 
;

//***************************************************
// LIST OF PARAMETERS IN THE MODEL 
//***************************************************    
parameters n_s n_m sigma_epsiHd sigma_epsiHk sigma_epsiA sigma_epsiJ   sigma_epsiSe sigma_epsiSm sigma_epsiSF   pp  hab Cyphi_loan_H Cyphi_loan_F Cyphi_def_H Cyphi_def_F phi_Fs phi_Hs  alphaa delta_K delta_H betta_m betta_s mu_m mu_e mu_F mu_H sigma_e1 sigma_m1 sigma_F sigma_H varphi_s varphi_m v_s v_m chi_b chi_e eta   a_s a_b    
psi_i psi_h rhoA rhoJ  rhoSe rhoSm rhoSF  rhoSH   rhoHd rhoHk rw phi_K xi_K kappa_ins  
theta_bank
theta_entr
entr_pop
bank_pop
xi_entr
xi_bank
ind_totRisk
rho_CR
;

//***************************************************
// DECLARATION EXOGENOUS VARIABLES IN THE MODEL 
//***************************************************    
varexo epsiA  epsiJ  epsiSe 
epsiSm
 epsiSF 
%epsiSH 
epsiHd epsiHk
epsiA_news epsiJ_news epsiSe_news epsiSm_news epsiSF_news epsiHd_news epsiHk_news;

//************************************************************
// Assignment of parameter values
//************************************************************

load 3D_param

set_param_value('n_m',n_m);
set_param_value('n_s',n_s);
set_param_value('alphaa',alphaa);
set_param_value('delta_K',delta_K);
set_param_value('delta_H',delta_H);
set_param_value('betta_m',betta_m);
set_param_value('betta_s',betta_s);
set_param_value('mu_m',mu_m);
set_param_value('mu_e',mu_e);
set_param_value('mu_F',mu_F);
set_param_value('mu_H',mu_H);
set_param_value('sigma_e1',sigma_e1);
set_param_value('sigma_m1',sigma_m1);
set_param_value('sigma_F',sigma_F);
set_param_value('sigma_H',sigma_H);
set_param_value('varphi_s',varphi_s);
set_param_value('varphi_m',varphi_m);
set_param_value('v_s',v_s);
set_param_value('v_m',v_m);
set_param_value('chi_b',chi_b);
set_param_value('chi_e',chi_e);
set_param_value('eta',eta);
set_param_value('a_s',a_s);
set_param_value('a_b',a_b);
set_param_value('psi_i',psi_i);
set_param_value('psi_h',psi_h);
set_param_value('rhoA',rhoA);
set_param_value('rhoJ',rhoJ);
set_param_value('rhoSe',rhoSe);
set_param_value('rhoSm',rhoSm);
set_param_value('rhoSF',rhoSF);
set_param_value('rhoSH',rhoSH);
set_param_value('rhoHd',rhoHd);
set_param_value('rhoHk',rhoHk);
set_param_value('hab',hab);

set_param_value('phi_Fs',phi_Fs);
set_param_value('phi_Hs',phi_Hs);
set_param_value('Cyphi_def_H',Cyphi_def_H);
set_param_value('Cyphi_def_F',Cyphi_def_F);
set_param_value('Cyphi_loan_H',Cyphi_loan_H);
set_param_value('Cyphi_loan_F',Cyphi_loan_F);
%set_param_value('pp',pp);
set_param_value('sigma_epsiA',sigma_epsiA);
set_param_value('sigma_epsiJ',sigma_epsiJ);
set_param_value('sigma_epsiSe',sigma_epsiSe);
set_param_value('sigma_epsiSm',sigma_epsiSm);
set_param_value('sigma_epsiSF',sigma_epsiSF);
set_param_value('sigma_epsiHd',sigma_epsiHd);
set_param_value('sigma_epsiHk',sigma_epsiHk);

set_param_value('rw',rw);
set_param_value('phi_K',phi_K);
set_param_value('xi_K',xi_K);
set_param_value('kappa_ins',kappa_ins);

set_param_value('theta_bank',theta_bank);
set_param_value('theta_entr',theta_entr);
set_param_value('entr_pop',entr_pop);
set_param_value('bank_pop',bank_pop);
set_param_value('xi_entr',xi_entr);
set_param_value('xi_bank',xi_bank);

set_param_value('rho_CR',rho_CR);
set_param_value('ind_totRisk',ind_totRisk);
%rho_CR = 0.66;
%ind_totRisk = 1; 

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// DECLARATION OF THE MODEL
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//************************************************************
model;
//************************************************************
//************************************************************
// Households 
//************************************************************
//*******************************
// Utility functions
//*******************************
//
UC_s = log(C_s-hab*C_s(-1)); 
//
UC_m = log(C_m-hab*C_m(-1)); 
//
UL_s = varphi_s*L_s^(1+eta)/(1+eta); 
//
UL_m = varphi_m*L_m^(1+eta)/(1+eta); 
//
UH_s = EJ*v_s*log(H_s(0)); 
//
UH_m = EJ*v_m*log(H_m(0)); 
//
Util_s = UC_s - UL_s + UH_s;
//
Util_m = UC_m - UL_m + UH_m;
//
UC_s_1 = 1/(C_s-hab*C_s(-1));
//
UC_m_1 = 1/(C_m-hab*C_m(-1));
//
UL_s_1 = varphi_s*L_s^(eta);
//
UL_m_1 = varphi_m*L_m^(eta);
//
UH_s_1 = EJ*v_s/(H_s(0));
//
UH_m_1 = EJ*v_m/(H_m(0));
//*******************************
// WELFARE
//*******************************
// 
Vs=Util_s+betta_s*Vs(+1);
// 
Vm=Util_m+betta_m*Vm(+1);
//*******************************
// Savers
//*******************************
//  foc C_s
Lambda_s = UC_s_1;
//  foc L_s
UL_s_1 = w*Lambda_s;
//  foc D
Lambda_s = betta_s*Lambda_s(1)*(R_D-(1-kappa_ins)*Tr(+1)/D);
//  foc wrt rfr bond holding
Lambda_s = betta_s*Lambda_s(1)*(R_DD);
//  foc H_s
Lambda_s*(q_H) = UH_s_1(0) + betta_s*Lambda_s(1)*(1-deltaH(1))*q_H(1);
// foc K_s
Lambda_s*(q_K+s_K)=betta_s*Lambda_s(1)*(r_K(+1)+(1-deltaK(+1))*q_K(+1));
//  Budget Constraint Savers
C_s + q_K*(K_s+s_K) + q_H*(H_s-(1-deltaH)*H_s(-1))+ D= w*L_s 
+ (r_K+(1-deltaK)*q_K)*K_s(-1) - kappa_ins*Tr*a_s + PI + PH
+ PI_K+(1-av_def/400)*R_D(-1)*D(-1)+av_def/400*(R_D(-1)-(1-kappa_ins)*Tr/(av_def/400*D(-1)))*D(-1)
+ (1-theta_entr)*(rho_e*equity_entr(-1)-chi_entr*entr_pop)
+(1-theta_bank)*(rho_b_goods_F*equity_bank_F(-1)+rho_b_goods_H*equity_bank_H(-1)-chi_bank*bank_pop); 

//*******************************
// Borrowers
//*******************************
// Default cut off borrowers depends on leverage
omega_bar_m = x_m(-1)/R_H;
// Household leverage definition
x_m = R_m*b_m/(H_m*q_H);
// Housing rate of return
R_H = (1-deltaH)*q_H/q_H(-1);
// foc C_m    
Lambda_m = UC_m_1;
// foc L_m
UL_m_1 = w*Lambda_m;

//  foc mortgage interest rate (R_m)
(omega_bar_m(+1)/R_m)*R_H(+1)*q_H*H_m*(- betta_m*Lambda_m(1)*Gamma_m_1(1) + xi_m*(betta_s*Lambda_s(1)/Lambda_s(0)*((1-theta_bank)+theta_bank*vi_bank(+1)))*((1-Gamma_H(1))*(Gamma_m_1(1) - mu_m*G_m_1(1))))=0;

//  foc housing for borrowers (H_m)
UH_m_1(0) - Lambda_m*q_H + betta_m*Lambda_m(1)*((1-Gamma_m(1))*R_H(1)*q_H  - Gamma_m_1(+1)*(-omega_bar_m(+1)/H_m)*R_H(+1)*q_H*H_m) 
+xi_m*(betta_s*Lambda_s(1)/Lambda_s(0)*((1-theta_bank)+theta_bank*vi_bank(+1)))*((1-Gamma_H(1))*(Gamma_m(1)- mu_m*G_m(1))*R_H(1)*q_H+(1-Gamma_H(1))*(Gamma_m_1(1)- mu_m*G_m_1(1))*(-omega_bar_m(+1)/H_m)*R_H(+1)*q_H*H_m)= 0;                                                                                  

//  foc debt level for borrowers (b_m)
Lambda_m-betta_m*Lambda_m(+1)*Gamma_m_1(+1)*omega_bar_m(+1)/b_m*R_H(1)*q_H*H_m  
-xi_m*(rho_b_util_H(1)*phi_H-(betta_s*Lambda_s(1)/Lambda_s(0)*((1-theta_bank)+theta_bank*vi_bank(+1)))*(1-Gamma_H(1))*(Gamma_m_1(1)- mu_m*G_m_1(1))*(omega_bar_m(+1)/b_m)*R_H(+1)*q_H*H_m)=0;

//  Budget Constraint Borrowers
n_m*C_m + n_m*q_H*H_m - n_m*(1-Gamma_m)*R_H*q_H(-1)*H_m(-1) = n_m*w*L_m + n_m*b_m-kappa_ins*Tr*a_b;

//************************************************************
// Entrepreneurs
//************************************************************
//  Rate of return to capital
R_K = (r_K+(1-deltaK)*q_K)/q_K(-1);

//  Default cut off entrepreneurs depends on leverage
omega_bar_e = x_e(-1)/R_K;

//  Firm leverage definition
x_e = R_F*(q_K*K_e-equity_entr)/(q_K*K_e);

//  FOC entrepreneurial loan interest rate (R_F)
- (betta_s*Lambda_s(1)/Lambda_s*(((1-theta_entr)+theta_entr*vi_entr(+1))))*Gamma_e_1(1)*R_K(+1)*((q_K*K_e-equity_entr)/(q_K*K_e*R_K(+1))) + xi_e*(betta_s*Lambda_s(1)/Lambda_s(0)*((1-theta_bank)+theta_bank*vi_bank(+1)))*((1-Gamma_F(1))*(Gamma_e_1(1)- mu_e*G_e_1(1)))*R_K(+1)*((q_K*K_e-equity_entr)/(q_K*K_e*R_K(+1)))=0;
           
//  Foc capital (k)
(betta_s*Lambda_s(1)/Lambda_s*(((1-theta_entr)+theta_entr*vi_entr(+1))))*((1-Gamma_e(1))*R_K(1)-Gamma_e_1(1)*K_e*R_K(+1)*((R_F*equity_entr)/(q_K*K_e^2*R_K(+1)))) 
+ xi_e*((betta_s*Lambda_s(1)/Lambda_s(0)*((1-theta_bank)+theta_bank*vi_bank(+1)))*(1-Gamma_F(1))*(Gamma_e(1)- mu_e*G_e(1)) *R_K(1)+(betta_s*Lambda_s(1)/Lambda_s(0)*((1-theta_bank)+theta_bank*vi_bank(+1)))*(1-Gamma_F(1))*(Gamma_e_1(1)- mu_e*G_e_1(1))*((R_F*equity_entr)/(q_K*K_e^2*R_K(+1)))*R_K(+1)*K_e - rho_b_util_F(1)*phi_F);

//  Evolution of equity for entrepreneurs
equity_entr=(theta_entr*(rho_e*equity_entr(-1))+chi_entr*(1-theta_entr)*entr_pop);

//  Transfers of entering entrepreneurs from households
chi_entr=xi_entr*(rho_e*equity_entr(-1));

// Evolution of price of equity for entrepreneurs
vi_entr=betta_s*Lambda_s(1)/Lambda_s*((1-theta_entr)+theta_entr*vi_entr(+1))*rho_e(+1);

// Rate of return on entrepreneurs activity
#lamb_b=betta_s*Lambda_s(0)/Lambda_s(-1)*((1-theta_bank)+theta_bank*vi_bank(0));
rho_e=rho_b_util_F*phi_F(-1)*(1-Gamma_e)*R_K/(rho_b_util_F*phi_F(-1)-lamb_b*(1-Gamma_F(0))*(Gamma_e-mu_e*G_e)*R_K);

//************************************************************
// Bankers
//************************************************************
// Evolution of equity for bankers
equity_bank=(theta_bank*(rho_b_goods_F*(phi_F(-1)*(q_K(-1)*K_e(-1)-equity_entr(-1)))+rho_b_goods_H*(equity_bank(-1)-(phi_F(-1)*(q_K(-1)*K_e(-1)-equity_entr(-1)))))+(1-theta_bank)*chi_bank*bank_pop);

// Transfers to entering bankers from households ***
chi_bank=xi_bank*(rho_b_goods_F*(phi_F(-1)*(q_K(-1)*K_e(-1)-equity_entr(-1)))+rho_b_goods_H*(equity_bank(-1)-(phi_F(-1)*(q_K(-1)*K_e(-1)-equity_entr(-1)))));

// Evolution of price of equity for bankers
vi_bank=betta_s*Lambda_s(1)/Lambda_s*((1-theta_bank)+theta_bank*vi_bank(+1))*rho_b_goods_F(+1);

// Non-arbitrage condition for return on bankers equity
betta_s*Lambda_s(1)/Lambda_s*((1-theta_bank)+theta_bank*vi_bank(+1))*(rho_b_goods_F(+1)-rho_b_goods_H(+1))=0;

//************************************************************
// Banks for NFC
//************************************************************
// Default threshold for banks
omega_bar_F = (1-phi_F(-1))*R_D(-1)/R_tilde_F(0);

// Part Constr of lending to entrepreneurs renting capital
R_tilde_F = (Gamma_e-mu_e*G_e)*R_K*q_K(-1)*K_e(-1)/(q_K(-1)*K_e(-1)-equity_entr(-1));

//  Realized rate of return on equity for banks F (in terms of utility and goods)
rho_b_util_F = lamb_b*(1-Gamma_F)*R_tilde_F/phi_F(-1);
rho_b_goods_F=rho_b_util_F/lamb_b;

//  Definition of banks K equity
equity_bank_F=phi_F*b_e;
//************************************************************
// Banks for HH
//************************************************************
// Default threshold for banks
omega_bar_H = (1-phi_H(-1))*R_D(-1)/R_tilde_H(0);

// Part Constr of lending to HH
R_tilde_H = (Gamma_m - mu_m*G_m)*R_H*q_H(-1)*H_m(-1)/b_m(-1);

// Realized rate of return on equity for banks H (in terms of utility and goods)
rho_b_util_H = lamb_b*(1-Gamma_H)*R_tilde_H/phi_H(-1);
rho_b_goods_H=rho_b_util_H/lamb_b;

// Definition of banks H equity
equity_bank_H=phi_H*b_m*n_m;

//************************************************************
// Banks - market clearing
//************************************************************
//  Market clearing condition on bank assets and liabilities ***
equity_bank + D = b_e + b_m*n_m; 

//   Inside equity equilibrium condition 
equity_bank=equity_bank_F+equity_bank_H;

//************************************************************
// Consumption good production
//************************************************************
//  Output
Y = A*K(-1)^(alphaa)*L^(1-alphaa);
//  Capital rental rate
r_K = alphaa*Y/K(-1);
//  Wage rate
w = (1-alphaa)*Y/L;
//************************************************************
// Capital good production
//************************************************************
// Jermann adjustment costs (depend on I/K(-1)) - Jermann (1998) - FULL REFERENCE
//   Capital price under Jermann AC
q_K*(delta_K^(1/psi_i))*(I/K(-1))^(-1/psi_i)=1;
//  Capital stock evolution under Jermann AC
K = (1-deltaK)*K(-1) + ((delta_K^(1/psi_i))/(1-1/psi_i)*(I/K(-1))^(1-1/psi_i)+(delta_K-(delta_K^(1/psi_i))/(1-1/psi_i)*delta_K^(1-1/psi_i)))*K(-1);
%K = (1-deltaK)*K(-1) + ((delta_K^(1/psi_i))/(1-1/psi_i)*(I/K(-1))^(1-1/psi_i)+delta_K/(1-psi_i))*K(-1);
//  Flow profit of capital producers under Jermann AC
PI = q_K*((delta_K^(1/psi_i))/(1-1/psi_i)*(I/K(-1))^(1-1/psi_i)+(delta_K-(delta_K^(1/psi_i))/(1-1/psi_i)*delta_K^(1-1/psi_i)))*K(-1) - I; 
//************************************************************
//with CEE adjustment costs (depend on I/I(-1)) - Christiano, Eichenbaum and Evans (1995) - Journal of Political Economy
// Not used in current version of the model
g_I = psi_i/2*(I/I(-1)-1)^2;
// Not used in current version of the model (kept just in case)
g_I_1 = psi_i*(I/I(-1)-1);
// Capital stock evolution under CEE AC (not used in current version of the model)
%K = (1-deltaK)*K(-1) + I*(1-g_I );
//  Capital price under CEE AC (not used in current version of the model)
%q_K = 1 + g_I + I/I(-1)*g_I_1 - betta_s*Lambda_s(1)/Lambda_s*(I(1)/I)^2*g_I_1(1);
// Flow profit of capital producers under Jermann AC (not used in current version of the model)
%PI = q_K*I - (1 + g_I)*I; 
// Capital market clearing
K=K_s+K_e;
//************************************************************
// Capital management firms
//************************************************************
// Profits
PI_K=s_K*K_s - z_K;
//
s_K=xi_K*K_s^(phi_K-1);
//
z_K=xi_K/phi_K*K_s^phi_K;

//*******************************
// Goods market
//*******************************
//  Aggregate consumption definition
C = C_s*n_s + C_m*n_m;
//*******************************
// Labour market 
//*******************************
// Aggregate labour supply definition
L = L_s*n_s + L_m*n_m;
//*******************************
// Housing supply
//*******************************
// Housing price with Jermann Adj costs
q_H*(delta_H^(1/psi_h))*(IH/H(-1))^(-1/psi_h)=1;
// Housing stock evolution with Jermann AC
H = (1-deltaH)*H(-1) + ((delta_H^(1/psi_h))/(1-1/psi_h)*(IH/H(-1))^(1-1/psi_h)+(delta_H-(delta_H^(1/psi_h))/(1-1/psi_h)*delta_H^(1-1/psi_h)))*H(-1);
%H = (1-deltaH)*H(-1) + ((delta_H^(1/psi_h))/(1-1/psi_h)*(IH/H(-1))^(1-1/psi_h)+delta_H/(1-psi_h))*H(-1);
// Flow profits by housing producers with Jermann AC
PH= q_H*((delta_H^(1/psi_h))/(1-1/psi_h)*(IH/H(-1))^(1-1/psi_h)+(delta_H-(delta_H^(1/psi_h))/(1-1/psi_h)*delta_H^(1-1/psi_h)))*H(-1) - IH;
//*******************************
//with CEE Adj costs (not used)
g_H = psi_h/2*(IH/IH(-1)-1)^2;
// (not used)
g_H_1 = psi_h*(IH/IH(-1)-1);
//  Housing stock evolution with CEE AC (not used)
//H = (1-deltaH)*H(-1) + IH*(1-g_H);
// Housing price with CEE Adj costs (not used)
//q_H = 1 + g_H + IH/IH(-1)*g_H_1 - betta_s*Lambda_s(1)/Lambda_s*(IH(1)/IH)^2*g_H_1(1);
// (not used)
//PH= q_H*IH - (1 + g_H)*IH; 
// Aggregate housing market clearing (Supply (LHS) = Demand (RHS))
H = H_m*n_m + H_s*n_s; 
//*******************************
// Deposit insurance agency
//*******************************
// Deposit insurance transfers to corporate bank depositors whose banks have gone bankrupt
Tr_F = (omega_bar_F - Gamma_F + mu_F*G_F)*R_tilde_F*(q_K(-1)*K_e(-1)-equity_entr(-1));
%Tr_F = (omega_bar_F - (1-mu_F)*normcdf((log(omega_bar_F)-(ESF*sigma_F)^2/2)/(ESF*sigma_F)) -omega_bar_F*normcdf((log(omega_bar_F)-(ESF*sigma_F)^2/2)/(ESF*sigma_F))/)*R_tilde_F*(q_K(-1)*K_e(-1)-equity_entr(-1));
// Deposit insurance transfers to household bank depositors whose banks have gone bankrupt
Tr_H = (omega_bar_H - Gamma_H + mu_H*G_H)*R_tilde_H*(n_m*H_m(-1)*q_H(-1)*x_m(-1))/R_m(-1);
%Tr_H = (omega_bar_H - Gamma_H + mu_H*G_H)*R_tilde_H*(n_m*H_m(-1)*q_H(-1)*x_m(-1))/R_m(-1);
// Total deposit insurance transfers
Tr = Tr_F + Tr_H;

//************************************************************
// Distribution functions (BGG parameters)
//************************************************************
#sigma_m1_comp=ESm*sigma_m1;
#sigma_e1_comp=ESe*sigma_e1;
// 74
Gamma_m   = normcdf((log(omega_bar_m)- (sigma_m1_comp)^2/2)/(sigma_m1_comp)) + omega_bar_m*(1-normcdf((log(omega_bar_m)+(sigma_m1_comp)^2/2)/(sigma_m1_comp)));
// 
Gamma_e   = normcdf((log(omega_bar_e) - (sigma_e1_comp)^2/2)/(sigma_e1_comp)) + omega_bar_e*(1-normcdf((log(omega_bar_e)+(sigma_e1_comp) ^2/2)/(sigma_e1_comp)));
// 
Gamma_F   = normcdf((log(omega_bar_F) - (ESF*sigma_F)^2/2)/(ESF*sigma_F)) + omega_bar_F*(1-normcdf((log(omega_bar_F)+(ESF*sigma_F)^2/2)/(ESF*sigma_F)));
// 
Gamma_H   = normcdf((log(omega_bar_H)-(ESH*sigma_H)^2/2)/(ESH*sigma_H)) + omega_bar_H*(1-normcdf((log(omega_bar_H)+(ESH*sigma_H)^2/2)/(ESH*sigma_H)));
//
Gamma_m_1 = (normpdf((log(omega_bar_m)-(sigma_m1_comp)^2/2)/(sigma_m1_comp))-omega_bar_m*normpdf((log(omega_bar_m)+(sigma_m1_comp)^2/2)/(sigma_m1_comp)))/(sigma_m1_comp*omega_bar_m) + (1-normcdf((log(omega_bar_m)+(sigma_m1_comp)^2/2)/(sigma_m1_comp)));
// 
Gamma_e_1 = (normpdf((log(omega_bar_e)-(sigma_e1_comp) ^2/2)/(sigma_e1_comp))-omega_bar_e*normpdf((log(omega_bar_e)+(sigma_e1_comp) ^2/2)/(sigma_e1_comp)))/(sigma_e1_comp *omega_bar_e) + (1-normcdf((log(omega_bar_e)+(sigma_e1_comp)^2/2)/(sigma_e1_comp)));
// 
Gamma_F_1 = (normpdf((log(omega_bar_F)-(ESF*sigma_F)^2/2)/(ESF*sigma_F))-omega_bar_F*normpdf((log(omega_bar_F)+(ESF*sigma_F)^2/2)/(ESF*sigma_F)))/(ESF*sigma_F*omega_bar_F) + (1-normcdf((log(omega_bar_F)+(ESF*sigma_F)^2/2)/(ESF*sigma_F)));
// 
Gamma_H_1 = (normpdf((log(omega_bar_H)-(ESH*sigma_H)^2/2)/(ESH*sigma_H))-omega_bar_H*normpdf((log(omega_bar_H)+(ESH*sigma_H)^2/2)/(ESH*sigma_H)))/(ESH*sigma_H*omega_bar_H) + (1-normcdf((log(omega_bar_H)+(ESH*sigma_H)^2/2)/(ESH*sigma_H)));
//
G_m   = normcdf((log(omega_bar_m)-(sigma_m1_comp)^2/2)/(sigma_m1_comp));
//
G_e   = normcdf((log(omega_bar_e)-(sigma_e1_comp) ^2/2)/(sigma_e1_comp));
//
G_F   = normcdf((log(omega_bar_F)-(ESF*sigma_F)^2/2)/(ESF*sigma_F));
//
G_H   = normcdf((log(omega_bar_H)-(ESH*sigma_H)^2/2)/(ESH*sigma_H));
//
G_m_1 = normpdf((log(omega_bar_m)-(sigma_m1_comp)^2/2)/(sigma_m1_comp))/(sigma_m1_comp*omega_bar_m);
//
G_e_1 = normpdf((log(omega_bar_e)-(sigma_e1_comp) ^2/2)/(sigma_e1_comp))/(sigma_e1_comp *omega_bar_e);

//************************************************************
// Capital requirement
//************************************************************
// Capital requirement on mortgage banks
%phi_H = steady_state(phi_H)+Cyphi_H*(def_rate_m(+1)-steady_state(def_rate_m));
// Capital requirement on corporate banks
%phi_F = steady_state(phi_F)+Cyphi_F*(def_rate_e(+1)-steady_state(def_rate_e));


// 92 Capital requirement on mortgage banks
(phi_H-phi_F*RW) = 0;
%RW=(1-rho_rw)*(rw+shock_rw*ECR)+RW(-1)*rho_rw;

// 92 Housing risk weight (exogenous)
%RW=rw;
%phi_H  =  phi_F(-1)^(rho_CR)*(steady_state(phi_F)*(b_m/steady_state(b_m))^(Cyphi_F))^(1-rho_CR);
%phi_F  =  phi_F(-1)^(rho_CR)*(steady_state(phi_F)*(b_tot/steady_state(b_tot))^(Cyphi_F))^(1-rho_CR);
phi_H  =  steady_state(phi_H)+Cyphi_loan_H*(b_m-steady_state(b_m))+Cyphi_def_H*(def_rate_m(+1)-steady_state(def_rate_m));
phi_F = steady_state(phi_F)+Cyphi_loan_F*(b_e-steady_state(b_e))+Cyphi_def_F*(def_rate_e(+1)-steady_state(def_rate_e));
phi_bar = (phi_F*b_e+phi_H*b_m)/(b_e+b_m);
//************************************************************
// Depreciation shocks
//************************************************************
// Capital depreciation
deltaK = delta_K+EdK;
// Housing depreciation
deltaH = delta_H+EdH;

//*********************
//****Other definitions
//*********************
//Spreads
Spread_e=R_F-R_D(-1);
Spread_m=R_m-R_D(-1);
//Housing wealth
HW=H*q_H;
//write-offs
write_off_m_new=(def_rate_m/100)*(b_m(-1)-(1-mu_m)/(def_rate_m/100)*G_m*R_H*q_H(-1)*H_m(-1))/b_m(-1)*400;
write_off_e_new=(def_rate_e/100)*(b_e(-1)-(1-mu_e)/(def_rate_e/100)*G_e*R_K*q_K(-1)*K_e(-1))/b_e(-1)*400;
//Loans
b_tot=b_e+b_m*n_m;
b_e=(q_K*K_e-equity_entr);
//Borrowing households wealth
W_m = (1-Gamma_m)*R_H*q_H(-1)*H_m(-1) + w*L_m;
//Deposit spread
RDsp = 400*(R_D-R_DD);
//Average default of banks
av_def = ((1-phi_F(-1))*equity_bank_F(-1)*def_rate_F*4/phi_F(-1) + (1-phi_H(-1))*(equity_bank_H(-1))*def_rate_H*4/phi_H(-1))/(D(-1)*n_s);
//Return on Loans spread
bsp_H = 400*(R_tilde_H(1)-R_D);
bsp_F = 400*(R_tilde_F(1)-R_D);
//GDP DEFINITIONS
GDP_acc=C  + I + IH;
GDP_acc_obs = log(GDP_acc/steady_state(GDP_acc))*100;

//Capital of savers as a share of total capital
sav_cap=K_s/K;
//************************************************************
// Reporting variables in deviation from the Steady State
//************************************************************
def_rate_m = normcdf((log(omega_bar_m)+(sigma_m1_comp)^2/2)/(sigma_m1_comp))*100;
def_rate_e = normcdf((log(omega_bar_e)+(sigma_e1_comp)^2/2)/(sigma_e1_comp))*100;
def_rate_H = normcdf((log(omega_bar_H)+(ESH*sigma_H)^2/2)/(ESH*sigma_H))*100;
def_rate_F = normcdf((log(omega_bar_F)+(ESF*sigma_F)^2/2)/(ESF*sigma_F))*100;
Y_obs = log(Y/steady_state(Y))*100; 
R_D_obs = (R_D-steady_state(R_D))*400; 
R_DD_obs = (R_DD-steady_state(R_DD))*400; 
R_m_obs = (R_m-steady_state(R_m))*400;
R_H_obs = log(R_H/steady_state(R_H))*400;
R_F_obs = (R_F-steady_state(R_F))*400;
H_m_obs = log(H_m/steady_state(H_m))*100;
H_s_obs = log(H_s/steady_state(H_s))*100;
b_m_obs = log(b_m/steady_state(b_m))*100;
b_e_obs = log(b_e/steady_state(b_e))*100;
C_obs = log(C/steady_state(C))*100;
C_m_obs = log(C_m/steady_state(C_m))*100;
C_s_obs = log(C_s/steady_state(C_s))*100;
D_obs = log(D/steady_state(D))*100;
E_F_obs = log(equity_bank_F/steady_state(equity_bank_F))*100;
I_obs = log(I/steady_state(I))*100;
K_obs = log(K/steady_state(K))*100;
L_obs = log(L/steady_state(L))*100;
L_m_obs = log(L_m/steady_state(L_m))*100;
L_s_obs = log(L_s/steady_state(L_s))*100;
n_b_obs = log(equity_bank/steady_state(equity_bank))*100;
n_e_obs = log(equity_entr/steady_state(equity_entr))*100;
q_H_obs = log(q_H/steady_state(q_H))*100;
q_K_obs = log(q_K/steady_state(q_K))*100;
r_K_obs = log(r_K/steady_state(r_K))*400;
R_K_obs = (R_K-steady_state(R_K))*400;
R_tilde_F_obs = (R_tilde_F-steady_state(R_tilde_F))*400;
R_tilde_H_obs = (R_tilde_H-steady_state(R_tilde_H))*400;
rho_F_obs = (rho_b_goods_F-steady_state(rho_b_goods_F))*400;
rho_H_obs = (rho_b_goods_H-steady_state(rho_b_goods_H))*400;
Tr_obs = log(Tr/steady_state(Tr))*100;
Tr_F_obs = log(Tr_F/steady_state(Tr_F))*100;
Tr_H_obs = log(Tr_H/steady_state(Tr_H))*100;
w_obs = log(w/steady_state(w))*100;
x_e_obs = log(x_e/steady_state(x_e))*100;
x_m_obs = log(x_m/steady_state(x_m))*100;
H_obs = log(H/steady_state(H))*100;
IH_obs = log(IH/steady_state(IH))*100;
W_m_obs = log(W_m/steady_state(W_m))*100;
b_tot_obs=log(b_tot/steady_state(b_tot))*100;
HW_obs=log(HW/steady_state(HW))*100;

//************************************************************
// SHOCKS 
// Note that each of these shocks potentially has a component 
// that hits immediately and a 'news shock' component
//************************************************************
log(A)   = rhoA*log(A(-1)) - epsiA/100 - A_news(-4);
log(EJ)  = rhoJ*log(EJ(-1)) - epsiJ/100 - EJ_news(-4);
log(ESe) = rhoSe*log(ESe(-1)) - sigma_epsiSe*epsiSe/100 - ESe_news(-4) - sigma_epsiSF*epsiSF/100;
// 191
log(ESm) = rhoSm*log(ESm(-1)) - ind_totRisk*sigma_epsiSm*epsiSe/100 - sigma_epsiSm*epsiSm/100  - ESm_news(-4);
%log(ESm) = rhoSm*log(ESm(-1)) - sigma_epsiSm*epsiSm/100 - ESm_news(-4);
// 192
log(ESF) = rhoSF*log(ESF(-1)) - ind_totRisk*sigma_epsiSF*epsiSe/100 - ESF_news(-4) - sigma_epsiSF*epsiSF/100;
%log(ESF) = rhoSF*log(ESF(-1)) - sigma_epsiSF*epsiSF/100 - ESF_news(-4) ;
// 193
log(ESH) = rhoSH*log(ESH(-1)) - ind_totRisk*sigma_epsiSF*epsiSe/100  - sigma_epsiSm*epsiSm/100 - ESF_news(-4) - sigma_epsiSF*epsiSF/100;
%log(ESH) = rhoSH*log(ESH(-1)) - sigma_epsiSF*epsiSF/100 - ESF_news(-4);
(EdH)      = rhoHd*(EdH(-1)) + epsiHd/100+ EdH_news(-4);
(EdK)      = rhoHk*(EdK(-1)) + epsiHk/100+ EdK_news(-4);

A_news=epsiA_news/100;
EJ_news=epsiJ_news/100;
ESe_news=epsiSe_news/100;
ESm_news=epsiSm_news/100;
ESF_news=epsiSF_news/100;
EdH_news=epsiHd_news/100;
EdK_news=epsiHk_news/100;

end; 

//************************************************************
// Here specify the innovations
// Note that it is possible to shut down a shock by assigning it 
// a zero variance.
//************************************************************
shocks;

%var epsiA   = sigma_epsiA/2;
%var epsiJ   = sigma_epsiJ/2;
var epsiSe;  stderr 1;
var epsiSm;  stderr 1;
var epsiSF;  stderr 1;
%var epsiHd  = sigma_epsiHd/2;
%var epsiHk  = sigma_epsiHk/2;

%var epsiA_news   = sigma_epsiA/2;
%var epsiJ_news   = sigma_epsiJ/2;
%var epsiSe_news;  stderr 1;
%var epsiSm_news;  stderr 1;
%var epsiSF_news;  stderr 1;
%var epsiHd_news  = sigma_epsiHd/2;
%var epsiHk_news  = sigma_epsiHk/2;

end;


//*******************************
// SS
//*******************************
//resid;
//steady;
//check;
//*******************************
// SIMULATIONS - IRFS
//*******************************
// Stochastic simulations (order=1 for linear approximation, order=2 for 2nd order approximation)
stoch_simul(order=1,periods=200,nograph,noprint);
%   GDP_acc_obs C_obs I_obs q_K_obs Int_net_obs Pi_net_obs NFC_loans_GDP_obs b_tot_GDP_obs
%  ESe def_rate_e def_rate_m av_def q_H_obs R_H_obs R_F_obs IH_obs H_obs bsp_H R_m_obs r_K_obs
%   phi_F K_share_obs b_tot b_e_obs b_m_obs 
%EdK ESF ESm ECR EPm A_obs; //36
%;
/*
// Compute mean and standard deviations of variables of interest

mean_phi_F=oo_.mean(26);
mean_def_F=oo_.mean(27);
mean_def_H=oo_.mean(28);
mean_rho_F=oo_.mean(29);
mean_sav_cap=oo_.mean(30);

sd_gdp=sqrt(oo_.var(1,1));
sd_IK=sqrt(oo_.var(2,2));
sd_IH=sqrt(oo_.var(3,3));
sd_qH=sqrt(oo_.var(4,4));
sd_b_m=sqrt(oo_.var(5,5));
sd_b_e=sqrt(oo_.var(6,6));

sd_def_e=sqrt(oo_.var(14,14));
sd_def_m=sqrt(oo_.var(13,13));

sd_def_av=sqrt(oo_.var(9,9));

sd_rho_F=sqrt(oo_.var(10,10));  
sd_HW=sqrt(oo_.var(11,11));
sd_C=sqrt(oo_.var(12,12));

sd_spread_m=sqrt(oo_.var(15,15));
sd_spread_e=sqrt(oo_.var(16,16));

mean_def_av=(oo_.mean(9));

mean_write_off_e=(oo_.mean(14));
mean_write_off_m=(oo_.mean(13));
mean_gdp=oo_.mean(24);
mean_b_m=oo_.mean(17);
mean_b_e=oo_.mean(18);

mean_I=oo_.mean(19);
mean_IH=oo_.mean(20);
mean_H=oo_.mean(21);
mean_H_m=oo_.mean(22);
mean_H_s=oo_.mean(23);

mean_spread_m=(oo_.mean(15));
mean_spread_e=(oo_.mean(16));

mean_HWborr=n_m*mean_H_m/(n_m*mean_H_m+mean_H_s);

data_mean=[mean_write_off_e 
mean_write_off_m
mean_def_av
(n_m*mean_b_m)/mean_gdp
mean_b_e/mean_gdp
mean_IH/mean_gdp
mean_HWborr
mean_I/mean_gdp
(n_m*mean_H_m+mean_H_s)/mean_gdp
mean_spread_e
mean_spread_m
mean_phi_F
mean_def_F
mean_def_H
mean_rho_F
mean_sav_cap
];

// Save results in a .mat file
save data_momm sd_gdp sd_IK sd_IH sd_qH sd_b_m sd_b_e sd_def_e sd_def_m sd_def_av sd_rho_F sd_HW sd_C  sd_spread_m sd_spread_e data_mean
*/