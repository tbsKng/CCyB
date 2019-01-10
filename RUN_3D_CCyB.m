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

save  3D_param

%Solve for the steady state of the model
[ys,check,x,data_SS,flaga] = ThreeD_steadystate;
SS_eq
phi_Fs_base=phi_F;
phi_Hs_base=phi_H;
phi_Fs = 0.09;
phi_Hs = rw*phi_Fs;

save  3D_param


% Parameters for deterministic simulation (not used at the moment)
%initval_file_ind    = 0; % loads the initial value file for the deterministic simulation
%trans_ind           = 0; % loads initial conditions and terminal conditions in case they are different ( relevant for transition from one steady state to another steady state)

rho_CR              = 0; % persistency of CCyB rule
%xi_K                = 0.0005;

% Switch off CCyB for the benchmark simulation
Cyphi_F             = 0;
Cyphi_H             = 0;

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

clear ys check xfs av_defs  def_rate_es def_rate_ms def_rate_Fs def_rate_Hs omega_bar_es omega_bar_Hs  omega_bar_Fs  omega_bar_ms
save 3D_param

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATION: 10% CR, no CCyB %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:1 %
    

    
    %sigma_epsiSe_news   = -log(ESe_tresh1)*100;
    %sigma_epsiSm_news   = -log(ESm_tresh1)*100;
    %sigma_epsiSF_news   = -log( ESF_tresh1)*100;

    sigma_epsiSe_news   = risk_shock_mag;
    sigma_epsiSm_news   = risk_shock_mag;
    %sigma_epsiSF_news   = risk_shock_mag;
    %sigma_epsiSH_news   = risk_shock_mag;
    sigma_epsiSF_news   = 80;
    sigma_epsiSH_news   = 80;
    
    sigma_epsiSe   = sigma_epsiSe_news;
    sigma_epsiSm   = sigma_epsiSm_news;
    sigma_epsiSF   = sigma_epsiSF_news;
    
    
    rhoSe   = 0.9;
    rhoSm   = 0.9;
    rhoSF   = 0.9;
    rhoSH   = 0.9;
    ind_totRisk =1;
    
    save('3D_param','ind_totRisk','rhoSe','rhoSm','rhoSF','rhoSH','rhoA','sigma_epsiA','sigma_epsiSe','sigma_epsiSm','sigma_epsiSF','sigma_epsiSe_news','sigma_epsiSm_news','sigma_epsiSF_news','-append')
    
    % 4pSimu simulates for all variables the responses to the shock (shock to all risk rates and bank risk)
    dynare ThreeDAdj_RD_stoch_4pSimu noclearall nostrict
    
    % Save the IRFs for the positive shock
    
    %data_A{k}  = [phi_H_epsiA GDP_acc_obs_epsiA C_obs_epsiA I_obs_epsiA IH_obs_epsiA q_K_obs_epsiA q_H_obs_epsiA b_e_obs_epsiA b_m_obs_epsiA def_rate_e_epsiA def_rate_m_epsiA av_def_epsiA RDsp_epsiA phi_F_epsiA ];
    data_Se{k} = [GDP_acc_obs_epsiSe C_obs_epsiSe I_obs_epsiSe IH_obs_epsiSe q_K_obs_epsiSe q_H_obs_epsiSe b_e_obs_epsiSe b_m_obs_epsiSe R_DD_obs_epsiSe def_rate_e_epsiSe def_rate_m_epsiSe av_def_epsiSe H_s_obs_epsiSe H_m_obs_epsiSe phi_F_epsiSe phi_H_epsiSe phi_bar_epsiSe];
    data_SF{k} = [GDP_acc_obs_epsiSF C_obs_epsiSF I_obs_epsiSF IH_obs_epsiSF q_K_obs_epsiSF q_H_obs_epsiSF b_e_obs_epsiSF b_m_obs_epsiSF R_DD_obs_epsiSF def_rate_e_epsiSF def_rate_m_epsiSF av_def_epsiSF H_s_obs_epsiSF H_m_obs_epsiSF phi_F_epsiSF phi_H_epsiSF phi_bar_epsiSF];
    
    %data_A_surprise{k} = data_A{k};
    
   
    % Substract the same positive shock from period 9 onwards

    %data_A_surprise{k}(9:end,:) = data_A{k}(9:end,:)-data_A{k}(1:end-8,:);
    %data_A_surprise{k}(:,end) = (data_A_surprise{k}(:,end)+phi_Fs)*100;
    %data_A{k}(:,end) = (data_A{k}(:,end)+phi_Fs)*100;
    %data_A_surprise{k}(:,1) = (data_A_surprise{k}(:,1)+phi_Hs)*100;
    %data_A{k}(:,1) = (data_A{k}(:,1)+phi_Hs)*100;
    
    data_Se_surprise{k} = data_Se{k};
    data_Se_surprise{k}(9:end,:) = data_Se{k}(9:end,:)-data_Se{k}(1:end-8,:);
    data_Se_surprise{k}(:,end-2) = (data_Se_surprise{k}(:,end-2)+phi_Fs)*100;
    data_Se{k}(:,end-2) = (data_Se{k}(:,end-2)+phi_Fs)*100;
    data_Se_surprise{k}(:,end-1) = (data_Se_surprise{k}(:,end-1)+phi_Hs)*100;
    data_Se{k}(:,end-1) = (data_Se{k}(:,end-1)+phi_Hs)*100;
    
    data_SF_surprise{k} = data_SF{k};
    data_SF_surprise{k}(9:end,:) = data_SF{k}(9:end,:)-data_SF{k}(1:end-8,:);
    data_SF_surprise{k}(:,end-2) = (data_SF_surprise{k}(:,end-2)+phi_Fs)*100;
    data_SF{k}(:,end-2) = (data_SF{k}(:,end-2)+phi_Fs)*100;
    data_SF_surprise{k}(:,end-1) = (data_SF_surprise{k}(:,end-1)+phi_Hs)*100;
    data_SF{k}(:,end-1) = (data_SF{k}(:,end-1)+phi_Hs)*100;
    
    
    % data_Sm{k} = [phi_H_epsiSm GDP_acc_obs_epsiSm C_obs_epsiSm I_obs_epsiSm IH_obs_epsiSm q_K_obs_epsiSm q_H_obs_epsiSm b_e_obs_epsiSm b_m_obs_epsiSm def_rate_e_epsiSm def_rate_m_epsiSm av_def_epsiSm RDsp_epsiSm phi_F_epsiSm ];
    
    % data_Sm_surprise{k} = data_Sm{k};
    % data_Sm_surprise{k}(9:end,:) =  data_Sm{k}(9:end,:)-data_Sm{k}(1:end-8,:);
    % data_Sm_surprise{k}(:,end) = (data_Sm_surprise{k}(:,end)+phi_Fs)*100;
    % data_Sm{k}(:,end) = (data_Sm{k}(:,end)+phi_Fs)*100;
    % data_Sm_surprise{k}(:,1) = (data_Sm_surprise{k}(:,1)+phi_Hs)*100;
    % data_Sm{k}(:,1) = (data_Sm{k}(:,1)+phi_Hs)*100;
    
    % Run the sectoral shocks NFC Risk + NFC Banks Risk
    % and HH Risk + Mortgage Bank Risk
    dynare ThreeDAdj_RD_stoch_Sectoral noclearall nostrict
    
    data_SeSF{k} =  [GDP_acc_obs_epsiSe C_obs_epsiSe I_obs_epsiSe IH_obs_epsiSe q_K_obs_epsiSe q_H_obs_epsiSe b_e_obs_epsiSe b_m_obs_epsiSe R_DD_obs_epsiSe def_rate_e_epsiSe def_rate_m_epsiSe av_def_epsiSe H_s_obs_epsiSe H_m_obs_epsiSe phi_F_epsiSe phi_H_epsiSe phi_bar_epsiSe];
    
    data_SmSH{k} =    [GDP_acc_obs_epsiSm C_obs_epsiSm I_obs_epsiSm IH_obs_epsiSm q_K_obs_epsiSm q_H_obs_epsiSm b_e_obs_epsiSm b_m_obs_epsiSm R_DD_obs_epsiSm def_rate_e_epsiSm def_rate_m_epsiSm av_def_epsiSm H_s_obs_epsiSm H_m_obs_epsiSm phi_F_epsiSm phi_H_epsiSm phi_bar_epsiSm];
    
    
    data_SeSF_surprise{k}  = data_SeSF{k};
    
    data_SmSH_surprise{k}  = data_SmSH{k};
    
    data_SeSF_surprise{k}(9:end,:) =  data_SeSF{k}(9:end,:)-data_SeSF{k}(1:end-8,:);
    data_SeSF_surprise{k}(:,end-2) = (data_SeSF_surprise{k}(:,end-2)+phi_Fs)*100;
    data_SeSF{k}(:,end-2) = (data_SeSF{k}(:,end-2)+phi_Fs)*100;
    data_SeSF_surprise{k}(:,end-1) = (data_SeSF_surprise{k}(:,end-1)+phi_Hs)*100;
    data_SeSF{k}(:,end-1) = (data_SeSF{k}(:,end-1)+phi_Hs)*100;
    
    data_SmSH_surprise{k}(9:end,:) =  data_SmSH{k}(9:end,:)-data_SmSH{k}(1:end-8,:);
    data_SmSH_surprise{k}(:,end-2) = (data_SmSH_surprise{k}(:,end-2)+phi_Fs)*100;
    data_SmSH{k}(:,end-2) = (data_SmSH{k}(:,end-2)+phi_Fs)*100;
    data_SmSH_surprise{k}(:,end-1) = (data_SmSH_surprise{k}(:,end-1)+phi_Hs)*100;
    data_SmSH{k}(:,end-1) = (data_SmSH{k}(:,end-1)+phi_Hs)*100;
    
end


% %%%%%%%%%%%%%%%%%%%%%%%%
% %% Countercyclical CR %%
% %%%%%%%%%%%%%%%%%%%%%%%%
%phi_Fs = 0.08;
%phi_Hs = phi_Fs*rw;
%save('3D_param','phi_Fs','phi_Hs','-append');

%Cyphi_F_vec         = [10 5];
%Cyphi_H_vec         = [10 5];
Cyphi_F_vec         = [0.1 0.3]; % first element sectoral response, second element response to total loans
Cyphi_H_vec         = [0.1 0.3];

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
        Cyphi_F     = Cyphi_F_vec(k);
        Cyphi_H     = Cyphi_H_vec(k);
        
        % Total Loan Response
    elseif k == 2
        Cyphi_F     = Cyphi_F_vec(k);
        %Cyphi_F     = 5;
        Cyphi_H     = 0;
        
        
    end
    
    %rho_CR = 0.15;
    %rho_CR = 0.66;
    %rho_CR = 0.7103;
    
    save('3D_param','rho_CR','Cyphi_F','Cyphi_H','sigma_epsiA','-append')
    ind_totRisk =1;
    save('3D_param','ind_totRisk','-append');
    
    if k == 1
        %Cyphi_F = 0.8;
        %Cyphi_H = 0.05;
        save('3D_param','Cyphi_F','Cyphi_H','-append')
        
        % All risk shock, response to morgage loans and NFC risk sensitivity
        dynare ThreeDAdj_RD_stoch_4pSimu_SecCCyB noclearall nostrict
        
        
    else
        % Response to total loans (if k in loop is k=2)
        dynare ThreeDAdj_RD_stoch_4pSimu noclearall nostrict
    end
    
    data_Se_CCB{k} = [GDP_acc_obs_epsiSe C_obs_epsiSe I_obs_epsiSe IH_obs_epsiSe q_K_obs_epsiSe q_H_obs_epsiSe b_e_obs_epsiSe b_m_obs_epsiSe R_DD_obs_epsiSe def_rate_e_epsiSe def_rate_m_epsiSe av_def_epsiSe H_s_obs_epsiSe H_m_obs_epsiSe phi_F_epsiSe phi_H_epsiSe phi_bar_epsiSe];
    %data_SF_CCB{k} = [GDP_acc_obs_epsiSF C_obs_epsiSF I_obs_epsiSF IH_obs_epsiSF q_K_obs_epsiSF q_H_obs_epsiSF b_e_obs_epsiSF b_m_obs_epsiSF R_DD_obs_epsiSF def_rate_e_epsiSF def_rate_m_epsiSF av_def_epsiSF H_s_obs_epsiSF H_m_obs_epsiSF phi_F_epsiSF phi_H_epsiSF];
    
    data_Se_CCB_surprise{k} = data_Se_CCB{k};
    data_Se_CCB_surprise{k}(9:end,:) = data_Se_CCB{k}(9:end,:)-data_Se_CCB{k}(1:end-8,:);
    data_Se_CCB_surprise{k}(:,end-2) = (data_Se_CCB_surprise{k}(:,end-2)+phi_Fs)*100;
    data_Se_CCB{k}(:,end-2) = (data_Se_CCB{k}(:,end-2)+phi_Fs)*100;
    data_Se_CCB_surprise{k}(:,end-1) = (data_Se_CCB_surprise{k}(:,end-1)+phi_Hs)*100;
    data_Se_CCB{k}(:,end-1) = (data_Se_CCB{k}(:,end-1)+phi_Hs)*100;
    
    % data_SF_CCB_surprise{k} = data_SF_CCB{k};
    % data_SF_CCB_surprise{k}(9:end,:) = data_SF_CCB{k}(9:end,:)-data_SF_CCB{k}(1:end-8,:);
    % data_SF_CCB_surprise{k}(:,end-2) = (data_SF_CCB_surprise{k}(:,end-2)+phi_Fs)*100;
    % data_SF_CCB{k}(:,end-2) = (data_SF_CCB{k}(:,end-2)+phi_Fs)*100;
    % data_SF_CCB_surprise{k}(:,end) = (data_SF_CCB_surprise{k}(:,end)+phi_Hs)*100;
    % data_SF_CCB{k}(:,end) = (data_SF_CCB{k}(:,end)+phi_Hs)*100;
    
    
    % Response to both loans
    %Cyphi_F = 0.1;
    %Cyphi_H = 0.1;
    save('3D_param','Cyphi_F','Cyphi_H','-append')
    dynare ThreeDAdj_RD_stoch_4pSimu_CCyB noclearall nostrict
    
    data_Se_CCB{6} = [GDP_acc_obs_epsiSe C_obs_epsiSe I_obs_epsiSe IH_obs_epsiSe q_K_obs_epsiSe q_H_obs_epsiSe b_e_obs_epsiSe b_m_obs_epsiSe R_DD_obs_epsiSe def_rate_e_epsiSe def_rate_m_epsiSe av_def_epsiSe H_s_obs_epsiSe H_m_obs_epsiSe phi_F_epsiSe phi_H_epsiSe phi_bar_epsiSe];
    
    data_Se_CCB_surprise{6} = data_Se_CCB{6};
    data_Se_CCB_surprise{6}(9:end,:) = data_Se_CCB{6}(9:end,:)-data_Se_CCB{6}(1:end-8,:);
    data_Se_CCB_surprise{6}(:,end-2) = (data_Se_CCB_surprise{6}(:,end-2)+phi_Fs)*100;
    data_Se_CCB{6}(:,end-2) = (data_Se_CCB{6}(:,end-2)+phi_Fs)*100;
    data_Se_CCB_surprise{6}(:,end-1) = (data_Se_CCB_surprise{6}(:,end-1)+phi_Hs)*100;
    data_Se_CCB{6}(:,end-1) = (data_Se_CCB{6}(:,end-1)+phi_Hs)*100;
    
    
    
    %    data_Sm_CCB{k} = [phi_H_epsiSm GDP_acc_obs_epsiSm C_obs_epsiSm I_obs_epsiSm IH_obs_epsiSm q_K_obs_epsiSm q_H_obs_epsiSm b_e_obs_epsiSm b_m_obs_epsiSm def_rate_e_epsiSm def_rate_m_epsiSm av_def_epsiSm RDsp_epsiSm phi_F_epsiSm ];
    
    %    data_Sm_CCB_surprise{k} = data_Sm_CCB{k};
    %    data_Sm_CCB_surprise{k}(9:end,:) =  data_Sm_CCB{k}(9:end,:)-data_Sm_CCB{k}(1:end-8,:);
    %    data_Sm_CCB_surprise{k}(:,end) = (data_Sm_CCB_surprise{k}(:,end)+phi_Fs)*100;
    %    data_Sm_CCB{k}(:,end) = (data_Sm_CCB{k}(:,end)+phi_Fs)*100;
    %    data_Sm_CCB_surprise{k}(:,1) = (data_Sm_CCB_surprise{k}(:,1)+phi_Hs)*100;
    %    data_Sm_CCB{k}(:,1) = (data_Sm_CCB{k}(:,1)+phi_Hs)*100;
    
    % Sectoral shocks, NFCs and HHs 
    if k == 1
        % HH sector risk, response to HH loans and NFC default sensitivity
        %Cyphi_F = 0.2;
        %Cyphi_H = 0.1;
        save('3D_param','Cyphi_F','Cyphi_H','-append')
        dynare ThreeDAdj_RD_stoch_Sectoral_HCCyB noclearall nostrict
        data_SmSH_CCB{k} =    [ GDP_acc_obs_epsiSm C_obs_epsiSm I_obs_epsiSm IH_obs_epsiSm q_K_obs_epsiSm q_H_obs_epsiSm b_e_obs_epsiSm b_m_obs_epsiSm R_DD_obs_epsiSm def_rate_e_epsiSm def_rate_m_epsiSm av_def_epsiSm H_s_obs_epsiSm H_m_obs_epsiSm phi_F_epsiSm phi_H_epsiSm phi_bar_epsiSm];
        
        % Response to HH loans only (switch off NFC default sensitivity)
        Cyphi_F = 0; 
        save('3D_param','Cyphi_F','-append')
        
        dynare ThreeDAdj_RD_stoch_Sectoral_HCCyB noclearall nostrict
        data_SmSH_CCB{5} =    [ GDP_acc_obs_epsiSm C_obs_epsiSm I_obs_epsiSm IH_obs_epsiSm q_K_obs_epsiSm q_H_obs_epsiSm b_e_obs_epsiSm b_m_obs_epsiSm R_DD_obs_epsiSm def_rate_e_epsiSm def_rate_m_epsiSm av_def_epsiSm H_s_obs_epsiSm H_m_obs_epsiSm phi_F_epsiSm phi_H_epsiSm phi_bar_epsiSm];
        
        
        data_SmSH_CCB_surprise{5}  = data_SmSH_CCB{5};
        data_SmSH_CCB_surprise{5}(9:end,:) =  data_SmSH_CCB{5}(9:end,:)-data_SmSH_CCB{5}(1:end-8,:);
        
        data_SmSH_CCB_surprise{5}(:,end-2) = (data_SmSH_CCB_surprise{5}(:,end-2)+phi_Fs)*100;
        data_SmSH_CCB{5}(:,end-2) = (data_SmSH_CCB{5}(:,end-2)+phi_Fs)*100;
        data_SmSH_CCB_surprise{5}(:,end-1) = (data_SmSH_CCB_surprise{5}(:,end-1)+phi_Hs)*100;
        data_SmSH_CCB{5}(:,end-1) = (data_SmSH_CCB{5}(:,end-1)+phi_Hs)*100;
        
       
        % NFC sector shock, response to NFC loans and sensitivity to HH default
        Cyphi_F = Cyphi_F_vec(k);
        Cyphi_H = Cyphi_H_vec(k);
        save('3D_param','Cyphi_F','Cyphi_H','-append')
        
        dynare ThreeDAdj_RD_stoch_Sectoral_FCCyB noclearall nostrict
        
        data_SeSF_CCB{k} =  [ GDP_acc_obs_epsiSe C_obs_epsiSe I_obs_epsiSe IH_obs_epsiSe q_K_obs_epsiSe q_H_obs_epsiSe b_e_obs_epsiSe b_m_obs_epsiSe R_DD_obs_epsiSe def_rate_e_epsiSe def_rate_m_epsiSe av_def_epsiSe H_s_obs_epsiSe H_m_obs_epsiSe phi_F_epsiSe phi_H_epsiSe phi_bar_epsiSe];
        
        % Response to NFC loans only (switch off HH default sensitivity)
        Cyphi_H = 0;
        save('3D_param','Cyphi_H','-append')
        
        dynare ThreeDAdj_RD_stoch_Sectoral_FCCyB noclearall nostrict
        
        data_SeSF_CCB{5} =  [ GDP_acc_obs_epsiSe C_obs_epsiSe I_obs_epsiSe IH_obs_epsiSe q_K_obs_epsiSe q_H_obs_epsiSe b_e_obs_epsiSe b_m_obs_epsiSe R_DD_obs_epsiSe def_rate_e_epsiSe def_rate_m_epsiSe av_def_epsiSe H_s_obs_epsiSe H_m_obs_epsiSe phi_F_epsiSe phi_H_epsiSe phi_bar_epsiSe];
        
        data_SeSF_CCB_surprise{5}  = data_SeSF_CCB{5};
        
        
        data_SeSF_CCB_surprise{5}(9:end,:) =  data_SeSF_CCB{5}(9:end,:)-data_SeSF_CCB{5}(1:end-8,:);
        data_SeSF_CCB_surprise{5}(:,end-2) = (data_SeSF_CCB_surprise{5}(:,end-2)+phi_Fs)*100;
        data_SeSF_CCB{5}(:,end-2) = (data_SeSF_CCB{5}(:,end-2)+phi_Fs)*100;
        data_SeSF_CCB_surprise{5}(:,end-1) = (data_SeSF_CCB_surprise{5}(:,end-1)+phi_Hs)*100;
        data_SeSF_CCB{5}(:,end-1) = (data_SeSF_CCB{5}(:,end-1)+phi_Hs)*100;
              

        % Response to both loans
        Cyphi_F = Cyphi_F_vec(k);
        Cyphi_H = Cyphi_H_vec(k);
        save('3D_param','Cyphi_F','Cyphi_H','-append')
        
        dynare ThreeDAdj_RD_stoch_Sectoral_CCyB noclearall nostrict
        data_SmSH_CCB{6} =    [ GDP_acc_obs_epsiSm C_obs_epsiSm I_obs_epsiSm IH_obs_epsiSm q_K_obs_epsiSm q_H_obs_epsiSm b_e_obs_epsiSm b_m_obs_epsiSm R_DD_obs_epsiSm def_rate_e_epsiSm def_rate_m_epsiSm av_def_epsiSm H_s_obs_epsiSm H_m_obs_epsiSm phi_F_epsiSm phi_H_epsiSm phi_bar_epsiSm];
        
        
        data_SmSH_CCB_surprise{6}  = data_SmSH_CCB{6};
        data_SmSH_CCB_surprise{6}(9:end,:) =  data_SmSH_CCB{6}(9:end,:)-data_SmSH_CCB{6}(1:end-8,:);
        
        data_SmSH_CCB_surprise{6}(:,end-2) = (data_SmSH_CCB_surprise{6}(:,end-2)+phi_Fs)*100;
        data_SmSH_CCB{6}(:,end-2) = (data_SmSH_CCB{6}(:,end-2)+phi_Fs)*100;
        data_SmSH_CCB_surprise{6}(:,end-1) = (data_SmSH_CCB_surprise{6}(:,end-1)+phi_Hs)*100;
        data_SmSH_CCB{6}(:,end-1) = (data_SmSH_CCB{6}(:,end-1)+phi_Hs)*100;
        
        data_SeSF_CCB{6} =  [ GDP_acc_obs_epsiSe C_obs_epsiSe I_obs_epsiSe IH_obs_epsiSe q_K_obs_epsiSe q_H_obs_epsiSe b_e_obs_epsiSe b_m_obs_epsiSe R_DD_obs_epsiSe def_rate_e_epsiSe def_rate_m_epsiSe av_def_epsiSe H_s_obs_epsiSe H_m_obs_epsiSe phi_F_epsiSe phi_H_epsiSe phi_bar_epsiSe];
        
        data_SeSF_CCB_surprise{6}  = data_SeSF_CCB{6};
        
        
        data_SeSF_CCB_surprise{6}(9:end,:) =  data_SeSF_CCB{6}(9:end,:)-data_SeSF_CCB{6}(1:end-8,:);
        data_SeSF_CCB_surprise{6}(:,end-2) = (data_SeSF_CCB_surprise{6}(:,end-2)+phi_Fs)*100;
        data_SeSF_CCB{6}(:,end-2) = (data_SeSF_CCB{6}(:,end-2)+phi_Fs)*100;
        data_SeSF_CCB_surprise{6}(:,end-1) = (data_SeSF_CCB_surprise{6}(:,end-1)+phi_Hs)*100;
        data_SeSF_CCB{6}(:,end-1) = (data_SeSF_CCB{6}(:,end-1)+phi_Hs)*100;

    else
        dynare ThreeDAdj_RD_stoch_Sectoral noclearall nostrict
        data_SeSF_CCB{k} =  [ GDP_acc_obs_epsiSe C_obs_epsiSe I_obs_epsiSe IH_obs_epsiSe q_K_obs_epsiSe q_H_obs_epsiSe b_e_obs_epsiSe b_m_obs_epsiSe R_DD_obs_epsiSe def_rate_e_epsiSe def_rate_m_epsiSe av_def_epsiSe H_s_obs_epsiSe H_m_obs_epsiSe phi_F_epsiSe phi_H_epsiSe phi_bar_epsiSe];
        
        data_SmSH_CCB{k} =    [ GDP_acc_obs_epsiSm C_obs_epsiSm I_obs_epsiSm IH_obs_epsiSm q_K_obs_epsiSm q_H_obs_epsiSm b_e_obs_epsiSm b_m_obs_epsiSm  R_DD_obs_epsiSm def_rate_e_epsiSm def_rate_m_epsiSm av_def_epsiSm H_s_obs_epsiSm H_m_obs_epsiSm phi_F_epsiSm phi_H_epsiSm phi_bar_epsiSm];  
    end
    
    
    data_SeSF_CCB_surprise{k}  = data_SeSF_CCB{k};
    
    data_SmSH_CCB_surprise{k}  = data_SmSH_CCB{k};
    
    
    
    data_SeSF_CCB_surprise{k}(9:end,:) =  data_SeSF_CCB{k}(9:end,:)-data_SeSF_CCB{k}(1:end-8,:);
    data_SeSF_CCB_surprise{k}(:,end-2) = (data_SeSF_CCB_surprise{k}(:,end-2)+phi_Fs)*100;
    data_SeSF_CCB{k}(:,end-2) = (data_SeSF_CCB{k}(:,end-2)+phi_Fs)*100;
    data_SeSF_CCB_surprise{k}(:,end-1) = (data_SeSF_CCB_surprise{k}(:,end-1)+phi_Hs)*100;
    data_SeSF_CCB{k}(:,end-1) = (data_SeSF_CCB{k}(:,end-1)+phi_Hs)*100;
    
    data_SmSH_CCB_surprise{k}(9:end,:) =  data_SmSH_CCB{k}(9:end,:)-data_SmSH_CCB{k}(1:end-8,:);
    data_SmSH_CCB_surprise{k}(:,end-2) = (data_SmSH_CCB_surprise{k}(:,end-2)+phi_Fs)*100;
    data_SmSH_CCB{k}(:,end-2) = (data_SmSH_CCB{k}(:,end-2)+phi_Fs)*100;
    data_SmSH_CCB_surprise{k}(:,end-1) = (data_SmSH_CCB_surprise{k}(:,end-1)+phi_Hs)*100;
    data_SmSH_CCB{k}(:,end-1) = (data_SmSH_CCB{k}(:,end-1)+phi_Hs)*100;
    
    
    
end

%% Plotting

varnames9=char('GDP','Consumption','Business Investment','Housing Investment','Price of Capital','House Prices','NFC Loans','Mortgage Loans','Real Rate','NFC Default Rate','HH Default Rate','Bank Default Rate','Housing Demand Savers','Housing Demand Borrowers','Capital Requirement Corp. Banks','Capital Requirement Mort. Banks','Average Capital Requirement');


RUN_CHART_Policy

%close all