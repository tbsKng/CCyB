%%%%%**********************************************************************
% This file solves for the steady state of the model (deterministic equilibrium in absence of shocks).
% The steady state file contains:
% 1. numerical solution for a subset of variables:
% Solve_SS.m is a file that reduces the problem to a system of equations
% and solves numerically for the previously selected unknowns. Equations to be solved are included in SS_eq.m file
% 2. derivations of the remaining variables in terms of the model parameters and the variables for which we solved in 1
%%%%%**********************************************************************
function [ys,check,x,data_SS,flaga] = ThreeD_steadystate(junk,ys)
format long
check=0;

%Load mat file with parameter values structure

load('3D_param')


input=[input_1
       input_2
       input_3
       input_4
       input_5
       input_6
       input_7
       input_8
       input_9];
  

%%%%%**********************************************************************   
% 1) A first set of endogenous variables has close form solution and 
% can be easily derived from a subset of equations
%%%%%**********************************************************************
R_DDs = (1/betta_s);  
rho_Fs = 0; 
rho_Hs = rho_Fs;      
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
%%%%%**********************************************************************
% 2. The remaining variables require the use of a solver (system of equations).
% We can then reduce the problem to a system of 9 equations in 9 unknowns.
% Note that it is not possible to use any model equation 2 times (with each model equation we have to derive the steady state value of one single variable!)
%%%%%**********************************************************************
% Solve_SS.m file reduce the problem to a system of 9 equations and solves numerically for the 9 previously selected unknowns.


options = optimset('TolFun',1e-28,'MaxFunEvals',2000,'Display','On');
[x fsfv flaga]=fsolve('Solve_SS',input,options,phi_Fs,phi_Hs,n_m,n_s,hab,delta_H,a_s,a_b,q_Ks,As,PIs,PHs,R_Hs,R_DDs,rho_Fs,rho_Hs,alphaa,delta_K,betta_m,betta_s,mu_m,mu_e,mu_F,mu_H,sigma_e1,sigma_m1,sigma_F,sigma_H,varphi_s,varphi_m,v_s,v_m,chi_b,chi_e,eta,phi_K,xi_K,kappa_ins,theta_bank,theta_entr,entr_pop,bank_pop,xi_entr,xi_bank);

SS_eq

% mortgage default rate
def_rate_ms = normcdf((log(omega_bar_ms)+sigma_m1_comp^2/2)/sigma_m1_comp)*100; 
% corporate default rate
def_rate_es = normcdf((log(omega_bar_es)+sigma_e1_comp^2/2)/sigma_e1_comp)*100; 
% Net rate of return on deposit equals the gross return minus transaction costs due to bank default
% Various utility measures 
UC_ss = log(C_ss*(1-hab));
UC_ms = log(C_ms*(1-hab));
%//3
UL_ss = varphi_s*L_ss^(1+eta)/(1+eta);
%//4
UL_ms = varphi_m*L_ms^(1+eta)/(1+eta);
%//5
UH_ss = v_s*log(H_ss);
%//6
UH_ms = v_m*log(H_ms);
%//7
Util_ss = UC_ss - UL_ss + UH_ss;
%//8
Util_ms = UC_ms - UL_ms + UH_ms;
UH_m_1s=v_m/H_ms;
UH_s_1s=v_s/H_ss;
% WELFARE
Welf_ss=Util_ss/(1-betta_s); % Savers welfare
Welf_ms=Util_ms/(1-betta_m); % Borrowers' welfare
%%-------------------------------------------------------------------------
%%-------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  HOUSE-KEEPING STUFF BELOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

av_def=av_defs;
b_e=b_es;
A = 1; EJ=1; 
ESe= 1; ESm= 1; ESH= 1; ESF= 1; 
EdH=0; EdK=0;
b_m = b_ms;
C_m = C_ms; C_s = C_ss; C = Cs;
D = Ds;
def_rate_e = def_rate_es;
def_rate_m = def_rate_ms;
def_rate_F = def_rate_Fs;
def_rate_H = def_rate_Hs;
G_e_1 = G_e_1s; G_e = G_es;
G_F = G_Fs; G_H = G_Hs;
g_I_1 = g_I_1s; g_I = g_Is;
g_H=g_Hs; g_H_1=g_H_1s;
G_m_1 = G_m_1s; G_m = G_ms;
Gamma_e = Gamma_es; Gamma_e_1 = Gamma_e_1s;
Gamma_F = Gamma_Fs; Gamma_F_1 = Gamma_F_1s;
Gamma_H = Gamma_Hs; Gamma_H_1 = Gamma_H_1s;
Gamma_m = Gamma_ms; Gamma_m_1 = Gamma_m_1s;
H_m = H_ms; H_s = H_ss; H = Hs; IH=IHs;
I = Is;
K = Ks;
L_m = L_ms; L_s = L_ss; L = Ls;
Lambda_m = Lambda_ms; Lambda_s = Lambda_ss;
omega_bar_e = omega_bar_es;
omega_bar_F = omega_bar_Fs;
omega_bar_H = omega_bar_Hs;
omega_bar_m = omega_bar_ms;
PI = PIs; PH=PHs;
q_H = q_Hs; q_K = q_Ks;
R_D = R_Ds; R_DD=R_DDs; R_F = R_Fs; R_H = R_Hs;
R_K = R_Ks; r_K = r_Ks;
R_m = R_ms;
R_tilde_F = R_tilde_Fs;
R_tilde_H = R_tilde_Hs;
Tr_F = Tr_Fs; Tr_H = Tr_Hs; Tr = Trs;
UC_m_1 = UC_m_1s;
UC_m = UC_ms;
UC_s_1 = UC_s_1s;
UC_s = UC_ss;
UH_m_1 = UH_m_1s;
UH_m = UH_ms;
UH_s_1 = UH_s_1s;
UH_s = UH_ss;
UL_m_1 = UL_m_1s;
UL_m = UL_ms;
UL_s_1 = UL_s_1s;
UL_s= UL_ss;
Util_m = Util_ms;
Util_s = Util_ss;

w = ws;
x_e = x_es; x_m = x_ms;
xi_e = xi_es; xi_m = xi_ms;
Y = Ys;

rho_e	=	rho_es	;
equity_entr	=	equity_entrs	;
chi_entr	=	chi_entrs	;
rho_b_goods_F	=	rho_b_goods_Fs	;
rho_b_goods_H	=	rho_b_goods_Hs	;
rho_b_util_F	=	rho_b_util_Fs	;
rho_b_util_H	=	rho_b_util_Hs	;
equity_bank_F	=	equity_bank_Fs	;
equity_bank_H	=	equity_bank_Hs	;
chi_bank	=	chi_banks	;
equity_bank	=	equity_banks	;
vi_entr	=	vi_entrs	;
vi_bank	=	vi_banks	;

Vs	=	Welf_ss	;
Vm	=	Welf_ms	;

deltaH=delta_H;
deltaK=delta_K;
%%-------------------------------------------------------------------------
% All _obs variables are in deviation from steady state. Hence all are zero
Y_obs = 0; 
R_D_obs = 0;
R_m_obs = 0;
R_H_obs = 0;
R_F_obs = 0;
H_m_obs = 0;
H_s_obs = 0;
b_m_obs = 0;
C_obs = 0;
C_m_obs = 0;
C_s_obs = 0;
D_obs = 0;
E_F_obs = 0;
I_obs = 0;
K_obs = 0;
L_obs = 0;
L_m_obs = 0;
L_s_obs = 0;
n_b_obs = 0;
n_e_obs = 0;
q_H_obs = 0;
q_K_obs = 0;
r_K_obs = 0;
R_K_obs = 0;
R_tilde_F_obs = 0;
R_tilde_H_obs = 0;
rho_F_obs = 0;
rho_H_obs = 0;
Tr_obs = 0;
Tr_F_obs = 0;
Tr_H_obs = 0;
w_obs = 0;
x_e_obs = 0;
x_m_obs = 0;
b_e_obs = 0;
%%-------------------------------------------------------------------------
% More reporting variables (definitions) computed below
%%-------------------------------------------------------------------------
A_news	=	0	;
EJ_news	=	0	;
ESm_news	=	0	;
ESe_news	=	0	;
ESF_news	=	0	;
EdH_news	=	0	;
EdK_news	=	0	;
HW	=	H	;
HW_obs	=	0	;

GDP_acc_obs	=	0	;



bsp_F=	400*(R_tilde_F-R_D);
bsp_H=	400*(R_tilde_H-R_D);

H = n_s*H_s + n_m*H_m;
H_obs = 0;
IH_obs = 0;
RDsp =400*(R_D-R_DD);
W_m = PHI_m + w*L_m;
W_m_obs = 0;
b_tot=b_e+b_m*n_m;
b_tot_obs=0;
GDP_acc=Cs  + Is + IHs;
write_off_m=(def_rate_m/100)*(b_m-(1-mu_m)/(def_rate_m/100)*G_m*R_H*q_H*H_m)/b_m*400;
write_off_e=(def_rate_e/100)*(b_e-(1-mu_e)/(def_rate_e/100)*G_e*R_K*q_K*K_es)/b_e*400;
Spread_e=R_F-R_D;
Spread_m=R_m-R_D;


write_off_m_new	=	write_off_m	;
write_off_e_new	=	write_off_e	;
K_e	=	K_es	;
K_s	=	K_ss	;
sav_cap=K_s/K;
RW = rw;
R_DD_obs = 0;
phi_bar = (phi_F*b_e+phi_H*b_m)/(b_e+b_m);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ENTIRE SS TO BE LOADED INTO DYNARE    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ys=[b_e
g_H
g_H_1
IH
PH
A
b_m
C
C_m
C_s
D
def_rate_e
def_rate_m
def_rate_F
def_rate_H
G_e
G_e_1
G_F
G_H
g_I
g_I_1
G_m
G_m_1
Gamma_e
Gamma_e_1
Gamma_F
Gamma_F_1
Gamma_H
Gamma_H_1
Gamma_m
Gamma_m_1
H
H_m
H_s
I
K
L
L_m
L_s
Lambda_m
Lambda_s
omega_bar_e
omega_bar_F
omega_bar_H
omega_bar_m
PI
q_H
q_K
R_D
R_F
R_H
r_K
R_K
R_m
R_tilde_F
R_tilde_H
Tr
Tr_F
Tr_H
UC_m
UC_m_1
UC_s
UC_s_1
UH_m
UH_m_1
UH_s
UH_s_1
UL_m
UL_m_1
UL_s
UL_s_1
Util_m
Util_s
w
x_e
x_m
xi_e
xi_m
Y
Y_obs
R_D_obs
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
Vs
Vm
EJ
ESe
ESm
ESH
ESF
EdH
EdK
phi_F
phi_H
b_e_obs
bsp_F
bsp_H
deltaH
deltaK
H_obs
IH_obs
av_def
RDsp
R_DD
W_m
W_m_obs
b_tot
b_tot_obs
A_news
EJ_news
ESm_news
ESe_news
ESF_news
EdH_news
EdK_news
HW
HW_obs
GDP_acc
GDP_acc_obs
Spread_e
Spread_m
write_off_m_new
write_off_e_new
K_e
K_s
PI_K
s_K
z_K
rho_e
equity_entr
chi_entr
rho_b_goods_F
rho_b_goods_H
rho_b_util_F
rho_b_util_H
equity_bank_F
equity_bank_H
chi_bank
equity_bank
vi_entr
vi_bank
sav_cap
RW
R_DD_obs
phi_bar];

size(ys);