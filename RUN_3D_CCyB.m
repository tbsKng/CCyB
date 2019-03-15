clear global
clear all
close all

% Load 3DM monetary calibration
%load calibrated_params_1st
%load 3D_param

% Load JMCB calibration
load calibrated_parameters

clear global
clear M_ oo_

hab = 0.6;
xi_K = 0.003;
rho_CR = 0;
save  3D_param

%Solve for the steady state of the model
[ys,check,x,data_SS,flaga] = ThreeD_steadystate;
SS_eq
phi_Fs_base=phi_F;
phi_Hs_base=phi_H;
%phi_Fs = 0.09;
%phi_Hs = rw*phi_Fs;

save  3D_param


% Parameters for deterministic simulation (not used at the moment)
%initval_file_ind    = 0; % loads the initial value file for the deterministic simulation
%trans_ind           = 0; % loads initial conditions and terminal conditions in case they are different ( relevant for transition from one steady state to another steady state)

%rho_CR              = 0.710300000000000; % persistency of CCyB rule
%xi_K                = 0.0005;

% Switch off CCyB for the benchmark simulation
Cyphi_loan_F        = 0;
Cyphi_loan_H        = 0;
Cyphi_def_F         = 0;
Cyphi_def_H         = 0;

ind_totRisk         = 1; % indicates wheter a shock should hit all default rates

% Determine size of risk shock

% [ys]=ThreeDAdj_RD_stoch_steadystate;
% av_defs         = ys(132);
% def_rate_es     = ys(12);
% def_rate_ms     = ys(13);
% def_rate_Fs     = ys(14);
% def_rate_Hs     = ys(15);
% omega_bar_es    = ys(42);
% omega_bar_ms    = ys(45);
% omega_bar_Fs     = ys(43);
% omega_bar_Hs     = ys(44);

% size_ESe = 0.5; % set to 50% of St.St. default rate
% size_ESF = 0.1; % set to 10pp default rate

% % Determine value of ESe, ESF, ESm, ESH
% ESe_tresh1=(norminv(def_rate_es/400*size_ESe)+sqrt((norminv(def_rate_es/400*size_ESe))^2-2*log(omega_bar_es)))/sigma_e1;
% ESe_tresh2=(norminv(def_rate_es/400*size_ESe)-sqrt((norminv(def_rate_es/400*size_ESe))^2-2*log(omega_bar_es)))/sigma_e1;

% ESm_tresh1=(norminv(def_rate_ms/400*size_ESe)+sqrt((norminv(def_rate_ms/400*size_ESe))^2-2*log(omega_bar_ms)))/sigma_m1;
% %ESm_tresh2=(norminv(def_rate_ms/400*0.5)-sqrt(norminv(def_rate_ms/400*0.5)^2-2*log(omega_bar_ms)))/sigma_m1;

% ESF_tresh1=(norminv(size_ESF/400)+sqrt((norminv(size_ESF/400))^2-2*log(omega_bar_Fs)))/sigma_F;
% %ESF_tresh2=(norminv(0.1)-sqrt(norminv(0.1)^2-2*log(omega_bar_Fs)))/sigma_F;

% ESH_tresh1=(norminv(size_ESF/400)+sqrt((norminv(size_ESF/400))^2-2*log(omega_bar_Hs)))/sigma_H;
%ESH_tresh2=(norminv(0.1)-sqrt(norminv(0.1)^2-2*log(omega_bar_Hs)))/sigma_H;

% Determine standard deviation of shock process
%sigma_epsiSe_news   = -log(ESe_tresh1)*100;
%sigma_epsiSm_news   = -log(ESm_tresh1)*100;
%sigma_epsiSF_news   = -log( ESF_tresh1)*100;

% Magnitude of the risk shocks
risk_shock_mag = 10;
%sigma_epsiSe_news   = -log(ESe_tresh1)*100;
%sigma_epsiSm_news   = -log(ESm_tresh1)*100;
%sigma_epsiSF_news   = -log( ESF_tresh1)*100;

sigma_epsiSe_news   = risk_shock_mag;
sigma_epsiSm_news   = risk_shock_mag;
%sigma_epsiSF_news   = risk_shock_mag;
%sigma_epsiSH_news   = risk_shock_mag;
sigma_epsiSF_news   = risk_shock_mag;
sigma_epsiSH_news   = risk_shock_mag;

sigma_epsiSe   = sigma_epsiSe_news;
sigma_epsiSm   = sigma_epsiSm_news;
sigma_epsiSF   = sigma_epsiSF_news;


rhoSe   = 0.9;
rhoSm   = 0.9;
rhoSF   = 0.9;
rhoSH   = 0.9;
ind_totRisk =1;

save('3D_param','phi_Fs','phi_Hs','ind_totRisk','rhoSe','rhoSm','rhoSF','rhoSH','rhoA','sigma_epsiA','sigma_epsiSe','sigma_epsiSm','sigma_epsiSF','sigma_epsiSe_news','sigma_epsiSm_news','sigma_epsiSF_news','-append')
    
%phi_Fs = 0.08;
%phi_Hs = phi_Fs*rw;
clear ys check xfs av_defs  def_rate_es def_rate_ms def_rate_Fs def_rate_Hs omega_bar_es omega_bar_Hs  omega_bar_Fs  omega_bar_ms
save 3D_param

variables ={'GDP_acc_obs','C_obs','I_obs','IH_obs','q_K_obs','q_H_obs','b_e_obs','b_m_obs','def_rate_e','def_rate_m','av_def','H_s_obs','H_m_obs','phi_F','phi_H'};
varnames=char('GDP','Consumption','Business Investment','Housing Investment','Price of Capital','House Prices','NFC Loans','Mortgage Loans','NFC Default Rate','HH Default Rate','Bank Default Rate','Housing Demand Savers','Housing Demand Borrowers','Capital Requirement Corp. Banks','Capital Requirement Mort. Banks');
var_Se 			= strcat(variables,'_epsiSe');
var_Sm 			= strcat(variables,'_epsiSm');
var_SF          = strcat(variables,'_epsiSF');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATION: 10% CR, no CCyB %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi_Fs_vec = [phi_Fs];

%phi_Fs_vec = [0.105 0.1 0.12 0.07];

for k = 1:length(phi_Fs_vec) %
    
    phi_Fs = phi_Fs_vec(k);
    phi_Hs = rw*phi_Fs;
    
    
    % 4pSimu simulates for all variables the responses to the shock (shock to all risk rates and bank risk)
    dynare ThreeDAdj_RD_stoch_4pSimu noclearall nostrict
    
    % Save the IRFs for the positive shock
    for v = 1:length(variables)
    	%data_A{k}(:,v) = eval(var_A{v});  
    	data_Se{k}(:,v)    = eval(var_Se{v});  
	    data_SeSF{k}(:,v)  = eval(var_SF{v});
        data_SmSH{k}(:,v)  = eval(var_Sm{v});          
    end
    % Substract the same positive shock from period 9 onwards

    %data_A_surprise{k}(9:end,:) = data_A{k}(9:end,:)-data_A{k}(1:end-8,:);
    %data_A_surprise{k}(:,end) = (data_A_surprise{k}(:,end)+phi_Fs)*100;
    %data_A{k}(:,end) = (data_A{k}(:,end)+phi_Fs)*100;
    %data_A_surprise{k}(:,1) = (data_A_surprise{k}(:,1)+phi_Hs)*100;
    %data_A{k}(:,1) = (data_A{k}(:,1)+phi_Hs)*100;
    
    data_Se_surprise{k} = data_Se{k};
    data_Se_surprise{k}(9:end,:) = data_Se{k}(9:end,:)-data_Se{k}(1:end-8,:);
    data_Se_surprise{k}(:,end-1) = (data_Se_surprise{k}(:,end-1)+phi_Fs)*100;
    data_Se{k}(:,end-1) = (data_Se{k}(:,end-1)+phi_Fs)*100;
    data_Se_surprise{k}(:,end) = (data_Se_surprise{k}(:,end)+phi_Hs)*100;
    data_Se{k}(:,end) = (data_Se{k}(:,end)+phi_Hs)*100;
    
    data_SeSF_surprise{k}  = data_SeSF{k};
    data_SmSH_surprise{k}  = data_SmSH{k};

    data_SeSF_surprise{k}(9:end,:) =  data_SeSF{k}(9:end,:)-data_SeSF{k}(1:end-8,:);
    data_SeSF_surprise{k}(:,end-1) = (data_SeSF_surprise{k}(:,end-1)+phi_Fs)*100;
    data_SeSF{k}(:,end-1) = (data_SeSF{k}(:,end-1)+phi_Fs)*100;
    data_SeSF_surprise{k}(:,end) = (data_SeSF_surprise{k}(:,end)+phi_Hs)*100;
    data_SeSF{k}(:,end) = (data_SeSF{k}(:,end)+phi_Hs)*100;
    
    data_SmSH_surprise{k}(9:end,:) =  data_SmSH{k}(9:end,:)-data_SmSH{k}(1:end-8,:);
    data_SmSH_surprise{k}(:,end-1) = (data_SmSH_surprise{k}(:,end-1)+phi_Fs)*100;
    data_SmSH{k}(:,end-1) = (data_SmSH{k}(:,end-1)+phi_Fs)*100;
    data_SmSH_surprise{k}(:,end) = (data_SmSH_surprise{k}(:,end)+phi_Hs)*100;
    data_SmSH{k}(:,end) = (data_SmSH{k}(:,end)+phi_Hs)*100;
    

%     data_SF_surprise{k} = data_SF{k};
%     data_SF_surprise{k}(9:end,:) = data_SF{k}(9:end,:)-data_SF{k}(1:end-8,:);
%     data_SF_surprise{k}(:,end-1) = (data_SF_surprise{k}(:,end-1)+phi_Fs)*100;
%     data_SF{k}(:,end-1) = (data_SF{k}(:,end-1)+phi_Fs)*100;
%     data_SF_surprise{k}(:,end) = (data_SF_surprise{k}(:,end)+phi_Hs)*100;
%     data_SF{k}(:,end) = (data_SF{k}(:,end)+phi_Hs)*100;
    
    
    %data_Sm{k}(:,v) = eval(var_Sm{v});

    % data_Sm_surprise{k} = data_Sm{k};
    % data_Sm_surprise{k}(9:end,:) =  data_Sm{k}(9:end,:)-data_Sm{k}(1:end-8,:);
    % data_Sm_surprise{k}(:,end) = (data_Sm_surprise{k}(:,end)+phi_Fs)*100;
    % data_Sm{k}(:,end) = (data_Sm{k}(:,end)+phi_Fs)*100;
    % data_Sm_surprise{k}(:,1) = (data_Sm_surprise{k}(:,1)+phi_Hs)*100;
    % data_Sm{k}(:,1) = (data_Sm{k}(:,1)+phi_Hs)*100;
end


% %%%%%%%%%%%%%%%%%%%%%%%%
% %% Countercyclical CR %%
% %%%%%%%%%%%%%%%%%%%%%%%%
%phi_Fs = 0.105;
%phi_Hs = phi_Fs*rw;
%save('3D_param','phi_Fs','phi_Hs','-append');

%Cyphi_F_vec         = [10 5];
%Cyphi_H_vec         = [10 5];
Cyphi_F_vec         = [0.25 0.3]; % first element sectoral response, second element response to total loans
Cyphi_H_vec         = [0.25 0.3]; % 0.25

% Create empty cells for pre-allocation
data_Se_CCB             = cell(7,1);
data_SF_CCB             = cell(7,1);
data_Sm_CCB             = cell(7,1);
data_SeSF_CCB           = cell(7,1);
data_SmSH_CCB           = cell(7,1);
data_Se_CCB_surprise    = cell(7,1);
data_SF_CCB_surprise    = cell(7,1);
data_Sm_CCB_surprise    = cell(7,1);
data_SeSF_CCB_surprise  = cell(7,1);
data_SmSH_CCB_surprise  = cell(7,1);


for k = 1:1
    
    
    % Sectoral Response
    if k == 1
        Cyphi_loan_F    = Cyphi_F_vec(k);
        Cyphi_loan_H    = 0;
        Cyphi_def_F     = 0;
        Cyphi_def_H     = Cyphi_H_vec(k);
        
        % Total Loan Response

        %rho_CR = 0.15;
        %rho_CR = 0.66;
        %rho_CR = 0.7103;

        save('3D_param','rho_CR','Cyphi_loan_F','Cyphi_loan_H','Cyphi_def_F','Cyphi_def_H','sigma_epsiA','-append')
        ind_totRisk =1;
        save('3D_param','ind_totRisk','-append');
              
        % All risk shock, response to morgage loans and NFC risk sensitivity
        dynare ThreeDAdj_RD_stoch_4pSimu noclearall nostrict
        for v = 1:length(variables)
            data_Se_CCB{k}(:,v)    = eval(var_Se{v});  
            data_SeSF_CCB{k}(:,v)  = eval(var_SF{v});
            data_SmSH_CCB{k}(:,v)  = eval(var_Sm{v});          
        end    

        data_Se_CCB_surprise{k} = data_Se_CCB{k};
        data_Se_CCB_surprise{k}(9:end,:) = data_Se_CCB{k}(9:end,:)-data_Se_CCB{k}(1:end-8,:);
        data_Se_CCB_surprise{k}(:,end-1) = (data_Se_CCB_surprise{k}(:,end-1)+phi_Fs)*100;
        data_Se_CCB{k}(:,end-1) = (data_Se_CCB{k}(:,end-1)+phi_Fs)*100;
        data_Se_CCB_surprise{k}(:,end) = (data_Se_CCB_surprise{k}(:,end)+phi_Hs)*100;
        data_Se_CCB{k}(:,end) = (data_Se_CCB{k}(:,end)+phi_Hs)*100;

        data_SmSH_CCB_surprise{k}  = data_SmSH_CCB{k};
        data_SmSH_CCB_surprise{k}(9:end,:) =  data_SmSH_CCB{k}(9:end,:)-data_SmSH_CCB{k}(1:end-8,:);        
        data_SmSH_CCB_surprise{k}(:,end-1) = (data_SmSH_CCB_surprise{k}(:,end-1)+phi_Fs)*100;
        data_SmSH_CCB{k}(:,end-1) = (data_SmSH_CCB{k}(:,end-1)+phi_Fs)*100;
        data_SmSH_CCB_surprise{k}(:,end) = (data_SmSH_CCB_surprise{k}(:,end)+phi_Hs)*100;
        data_SmSH_CCB{k}(:,end) = (data_SmSH_CCB{k}(:,end)+phi_Hs)*100;
               
        data_SeSF_CCB_surprise{k}  = data_SeSF_CCB{k};                
        data_SeSF_CCB_surprise{k}(9:end,:) =  data_SeSF_CCB{k}(9:end,:)-data_SeSF_CCB{k}(1:end-8,:);
        data_SeSF_CCB_surprise{k}(:,end-1) = (data_SeSF_CCB_surprise{k}(:,end-1)+phi_Fs)*100;
        data_SeSF_CCB{k}(:,end-1) = (data_SeSF_CCB{k}(:,end-1)+phi_Fs)*100;
        data_SeSF_CCB_surprise{k}(:,end) = (data_SeSF_CCB_surprise{k}(:,end)+phi_Hs)*100;
        data_SeSF_CCB{k}(:,end) = (data_SeSF_CCB{k}(:,end)+phi_Hs)*100;

        Cyphi_loan_H = 0;
        Cyphi_def_H = 0;
        Cyphi_def_F = 0;
        save('3D_param','Cyphi_loan_F','Cyphi_loan_H','Cyphi_def_F','Cyphi_def_H','-append')
        
        % Response to corporate loans
        dynare ThreeDAdj_RD_stoch_4pSimu noclearall nostrict
        for v = 1:length(variables)
            data_Se_CCB{2}(:,v)    = eval(var_Se{v});  
            data_SeSF_CCB{2}(:,v)  = eval(var_SF{v});       
        end

        data_Se_CCB_surprise{2} = data_Se_CCB{2};
        data_Se_CCB_surprise{2}(9:end,:) = data_Se_CCB{2}(9:end,:)-data_Se_CCB{2}(1:end-8,:);
        data_Se_CCB_surprise{2}(:,end-1) = (data_Se_CCB_surprise{2}(:,end-1)+phi_Fs)*100;
        data_Se_CCB{2}(:,end-1) = (data_Se_CCB{2}(:,end-1)+phi_Fs)*100;
        data_Se_CCB_surprise{2}(:,end) = (data_Se_CCB_surprise{2}(:,end)+phi_Hs)*100;
        data_Se_CCB{2}(:,end) = (data_Se_CCB{2}(:,end)+phi_Hs)*100;
    
        data_SeSF_CCB_surprise{2}  = data_SeSF_CCB{2};                
        data_SeSF_CCB_surprise{2}(9:end,:) =  data_SeSF_CCB{2}(9:end,:)-data_SeSF_CCB{2}(1:end-8,:);
        data_SeSF_CCB_surprise{2}(:,end-1) = (data_SeSF_CCB_surprise{2}(:,end-1)+phi_Fs)*100;
        data_SeSF_CCB{2}(:,end-1) = (data_SeSF_CCB{2}(:,end-1)+phi_Fs)*100;
        data_SeSF_CCB_surprise{2}(:,end) = (data_SeSF_CCB_surprise{2}(:,end)+phi_Hs)*100;
        data_SeSF_CCB{2}(:,end) = (data_SeSF_CCB{2}(:,end)+phi_Hs)*100;
        
        
        Cyphi_loan_H = Cyphi_H_vec(k);
        Cyphi_loan_F = 0;
        Cyphi_def_H = 0;
        Cyphi_def_F = 0;
        save('3D_param','Cyphi_loan_F','Cyphi_loan_H','Cyphi_def_F','Cyphi_def_H','-append')
        
        % Response to mortgage loans
        dynare ThreeDAdj_RD_stoch_4pSimu noclearall nostrict
        for v = 1:length(variables)
            data_SmSH_CCB{2}(:,v)  = eval(var_Sm{v});          
        end
        
        data_SmSH_CCB_surprise{2}  = data_SmSH_CCB{2};
        data_SmSH_CCB_surprise{2}(9:end,:) =  data_SmSH_CCB{2}(9:end,:)-data_SmSH_CCB{2}(1:end-8,:);        
        data_SmSH_CCB_surprise{2}(:,end-1) = (data_SmSH_CCB_surprise{2}(:,end-1)+phi_Fs)*100;
        data_SmSH_CCB{2}(:,end-1) = (data_SmSH_CCB{2}(:,end-1)+phi_Fs)*100;
        data_SmSH_CCB_surprise{2}(:,end) = (data_SmSH_CCB_surprise{2}(:,end)+phi_Hs)*100;
        data_SmSH_CCB{2}(:,end) = (data_SmSH_CCB{2}(:,end)+phi_Hs)*100;
               
        
              
        % data_SF_CCB_surprise{k} = data_SF_CCB{k};
        % data_SF_CCB_surprise{k}(9:end,:) = data_SF_CCB{k}(9:end,:)-data_SF_CCB{k}(1:end-8,:);
        % data_SF_CCB_surprise{k}(:,end-1) = (data_SF_CCB_surprise{k}(:,end-1)+phi_Fs)*100;
        % data_SF_CCB{k}(:,end-1) = (data_SF_CCB{k}(:,end-1)+phi_Fs)*100;
        % data_SF_CCB_surprise{k}(:,end) = (data_SF_CCB_surprise{k}(:,end)+phi_Hs)*100;
        % data_SF_CCB{k}(:,end) = (data_SF_CCB{k}(:,end)+phi_Hs)*100;
        
        
        % Response to both loans seperately
        Cyphi_loan_F    = Cyphi_F_vec(k);
        Cyphi_loan_H    = Cyphi_F_vec(k);
        Cyphi_def_F     = 0;
        Cyphi_def_H     = 0;

        save('3D_param','Cyphi_loan_F','Cyphi_loan_H','Cyphi_def_F','Cyphi_def_H','-append')
        dynare ThreeDAdj_RD_stoch_4pSimu noclearall nostrict
        for v = 1:length(variables)
                data_Se_CCB{6}(:,v)    = eval(var_Se{v});  
                data_SeSF_CCB{6}(:,v)  = eval(var_SF{v});
                data_SmSH_CCB{6}(:,v)  = eval(var_Sm{v});          
        end
    

        data_Se_CCB_surprise{6} = data_Se_CCB{6};
        data_Se_CCB_surprise{6}(9:end,:) = data_Se_CCB{6}(9:end,:)-data_Se_CCB{6}(1:end-8,:);
        data_Se_CCB_surprise{6}(:,end-1) = (data_Se_CCB_surprise{6}(:,end-1)+phi_Fs)*100;
        data_Se_CCB{6}(:,end-1) = (data_Se_CCB{6}(:,end-1)+phi_Fs)*100;
        data_Se_CCB_surprise{6}(:,end) = (data_Se_CCB_surprise{6}(:,end)+phi_Hs)*100;
        data_Se_CCB{6}(:,end) = (data_Se_CCB{6}(:,end)+phi_Hs)*100;

        data_SmSH_CCB_surprise{6}  = data_SmSH_CCB{6};
        data_SmSH_CCB_surprise{6}(9:end,:) =  data_SmSH_CCB{6}(9:end,:)-data_SmSH_CCB{6}(1:end-8,:);        
        data_SmSH_CCB_surprise{6}(:,end-1) = (data_SmSH_CCB_surprise{6}(:,end-1)+phi_Fs)*100;
        data_SmSH_CCB{6}(:,end-1) = (data_SmSH_CCB{6}(:,end-1)+phi_Fs)*100;
        data_SmSH_CCB_surprise{6}(:,end) = (data_SmSH_CCB_surprise{6}(:,end)+phi_Hs)*100;
        data_SmSH_CCB{6}(:,end) = (data_SmSH_CCB{6}(:,end)+phi_Hs)*100;
               
        data_SeSF_CCB_surprise{6}  = data_SeSF_CCB{6};                
        data_SeSF_CCB_surprise{6}(9:end,:) =  data_SeSF_CCB{6}(9:end,:)-data_SeSF_CCB{6}(1:end-8,:);
        data_SeSF_CCB_surprise{6}(:,end-1) = (data_SeSF_CCB_surprise{6}(:,end-1)+phi_Fs)*100;
        data_SeSF_CCB{6}(:,end-1) = (data_SeSF_CCB{6}(:,end-1)+phi_Fs)*100;
        data_SeSF_CCB_surprise{6}(:,end) = (data_SeSF_CCB_surprise{6}(:,end)+phi_Hs)*100;
        data_SeSF_CCB{6}(:,end) = (data_SeSF_CCB{6}(:,end)+phi_Hs)*100;

    end   
end

%% Plotting

RUN_CHART_Policy

%RUN_CHART_CR
%close all
