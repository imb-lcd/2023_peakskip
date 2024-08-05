function [name_pars,SimPars] = get_original_value(nAllele,isModified,OME,TC)
%GET_ORIGINAL_VALURE Return the original values and names of the model
%parameters

% OME = 1; % omega
% TC = 1;  % time scaling

%% set up the original parameters
% Parameters
L0 = 0.6; % Time scale
c0 = 1.7; % Level scale

% p53 production and degradation rate
alpha_p = 0.2/L0/c0;  % units: [M/h]
% alpha_p = alpha_p*0.8
beta_p  = 1/L0/c0;    % units: [M/h] ?
beta_p = beta_p*8;

% Mdm2 mRNA transcription switch on and off
alpha_T = (3.5/L0/(c0^2))*10^1.35; % units: [1/M^2*1/h]
alpha_T = alpha_T*3.5;
beta_T  = 1.2/L0;     % units: [1/h]
% beta_T = beta_T*0.5;

% MDM2 mRNA production and degradation rate
alpha_m = 1/L0/c0;    % units: [M/h]
beta_m  = 1/L0;       % units: [1/h]

% Mdm2 protein production and degradation rate
alpha_M = 1/L0;       % units: [M/h]
beta_M  = 1/L0;       % units: [M/h]

% complex formation and dissociation
phi = 10^0.79/L0/c0;  % units: [1/M*1/h]
psi = 10^0.3/L0;      % units: [1/h]
% phi = phi.*2;
% psi = psi./2;

% q mRNA production and degradation
alpha_q = (0.002/L0/c0)*10^0.2; % units: [1/M*1/h] ?
% alpha_q = alpha_q*0.5;
beta_q  = 1/L0;       % units: [M/h]
% beta_q = beta_q*0.25;

% Q protein production and degradation
alpha_Q = 1/L0;       % units: [M/h]
% alpha_Q = alpha_Q*0.5;
beta_Q  = 1/L0;       % units: [M/h]
% beta_Q = beta_Q*0.5;

% Scale
sig = 6.02*0.2*10^2;
alpha_p = alpha_p*sig; 
% beta_p = beta_p/sig; ?
alpha_m = alpha_m*sig; 
alpha_T = alpha_T/(sig^2); 
phi = phi/sig;
% alpha_q = alpha_q/sig; ?

if isModified==0
    % SimPars = [alpha_p beta_p,... % p53 production and degradation
    %     alpha_T beta_T,... % Mdm2 mRNA transcription switch on and off
    %     alpha_m beta_m,... % MDM2 mRNA production and degradation
    %     alpha_M beta_M,... % Mdm2 protein production and degradation
    %     phi psi,...        % complex formation and dissociation
    %     alpha_q beta_q,... % q mRNA production and degradation
    %     alpha_Q beta_Q];   % Q protein production and degradation
    SimPars = [alpha_p*OME*TC       beta_p/OME*TC,... % p53 production and degradation
        alpha_T/OME^2*TC     beta_T*TC,... % Mdm2 mRNA transcription switch on and off
        alpha_m*OME*TC       beta_m*TC,... % MDM2 mRNA production and degradation
        alpha_M*TC beta_M*TC,... % Mdm2 protein production and degradation
        phi/OME*TC           psi*TC,...        % complex formation and dissociation
        alpha_q/OME*TC       beta_q*TC,... % q mRNA production and degradation
        alpha_Q*TC           beta_Q*TC];   % Q protein production and degradation
else
    SimPars = [alpha_p*OME*TC       beta_p/OME*TC,... % p53 production and degradation
        alpha_T/OME^2*TC     beta_T*TC,... % Mdm2 mRNA transcription switch on and off
        alpha_m*OME*TC       beta_m*TC,... % MDM2 mRNA production and degradation
        alpha_M*2/nAllele*TC beta_M*TC,... % Mdm2 protein production and degradation
        phi/OME*TC           psi*TC,...        % complex formation and dissociation
        alpha_q/OME*TC       beta_q*TC,... % q mRNA production and degradation
        alpha_Q*TC           beta_Q*TC];   % Q protein production and degradation
end

name_pars = ["alpha_p" "beta_p",... % p53 production and degradation
             "alpha_T" "beta_T",... % Mdm2 mRNA transcription switch on and off
             "alpha_m" "beta_m",... % MDM2 mRNA production and degradation
             "alpha_M" "beta_M",... % Mdm2 protein production and degradation
             "phi" "psi",...        % complex formation and dissociation
             "alpha_q" "beta_q",... % q mRNA production and degradation
             "alpha_Q" "beta_Q"];   % Q protein production and degradation
end

