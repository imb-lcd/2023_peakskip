%% Initialize
clear; close all; clc

%% Stochiometry
S = readtable('stoichiometry.txt');
S = table2array(S(1:11,2:23));

%% Parameters
L0 = 0.6; % Time scale
c0 = 1.7; % Level scale

% p53
alpha_p = 0.2/L0/c0;  % units: [M/h]
beta_p = 123/L0/c0; % units: [1/M*1/h]

% MDM2 mRNA
alpha_m = 1/L0/c0; % units: [M/h]
beta_m = 1/L0; % units: [1/h]

% Factor for KD treatments
gamma1 = 1;
gamma2 = 1;

% Mdm2 protein
alpha_M = 1/L0; % units: [1/h]
beta_M = 1/L0; % units: [1/h]

% q mRNA
alpha_q = 0.38/L0/c0; % units: [1/M*1/h]
beta_q = 1/L0; % units: [1/h]

% Q protein
alpha_Q = 1/L0; % units: [1/h]
beta_Q = 1/L0; % units: [1/h]

% Transcription switch
alpha_T = (3.5*10^1.35)/L0/(c0^2); % units: [1/M^2*1/h]
beta_T = 1.2/L0; % units: [1/h]

% Complex
phi = (10^0.79)/L0/c0; % units: [1/M*1/h]
psi = (10^0.3)/L0; % units: [1/h]

% Scale
sig = 6.02*0.2*10^2;
alpha_p = alpha_p*sig; beta_p = beta_p/sig; alpha_m = alpha_m*sig;...
    alpha_T = alpha_T/(sig^2); phi = phi/sig; alpha_q = alpha_q/sig;

% Collect parameters
SimPars = [alpha_p,beta_p,...
    alpha_m,beta_m,...
    gamma1,gamma2,...
    alpha_M,beta_M,...
    alpha_q,beta_q,...
    alpha_Q,beta_Q,...
    alpha_T,beta_T,...
    phi,psi];

%% Simulation

% Generate RNG for Gillespie algorithm
defaultStream = RandStream('mlfg6331_64');
RandStream.setGlobalStream(defaultStream)
defaultStream.Substream = 1;

% Initial condition =======================================================
X = [0, 0,  0,  0,  0,  0,  0,  0,  0,  0,	0];
%   [p, T1, m1, M1, T2, m2, M2, C1, C2, q,  Q];

% Simulation ==============================================================
ts = 0.5; click = 0; Tmax = 400; RT = 0;
SimVal = [];

while RT < Tmax
    
    % DNA damage-dependent initial Mdm2 suppresion ------------------------
    if RT < 4
        lambda_M = 10;
        lambda_C = 0.1;
    else
        lambda_M = 1;
        lambda_C = 1;
    end
    
    % Reactants -----------------------------------------------------------
    p = X(1); % p53
    T1 = X(2); T2 = X(5); % TX state
    m1 = X(3); m2 = X(6); % Mdm2 mRNA
    M1 = X(4); M2 = X(7); % Mdm2 protein
    C1 = X(8); C2 = X(9); % Complex
    q = X(10); Q = X(11); % q mRNA and Q protein
    
    % Rates ---------------------------------------------------------------
    R(1) = alpha_p;  % Unbound p53 production
    R(2) = alpha_T*(1-T1)*p^2;  R(8) = alpha_T*(1-T2)*p^2; % p53-depedent TX ON
    R(3) = beta_T*T1;           R(9) = beta_T*T2;    % TX OFF
    R(4) = alpha_m*T1;          R(10) = alpha_m*T2;  % MDM2 mRNA production
    R(5) = gamma1*beta_m*m1;    R(11) = gamma2*beta_m*m2;   % MDM2 mRNA degradation
    R(6) = alpha_M*m1;          R(12) = alpha_M*m2;  % Mdm2 protein production
    R(7) = beta_M*M1*lambda_M;  R(13) = beta_M*M2*lambda_M; % Mdm2 protein degradation
    R(14) = phi*p*M1*lambda_C;  R(16) = phi*p*M2*lambda_C;  % Complex formation
    R(15) = psi*C1;             R(17) = psi*C2;   % Complex dissociation
    R(18) = alpha_q*p^2;    % q mRNA production
    R(19) = beta_q*q;       % q mRNA degradation
    R(20) = alpha_Q*q;      % Q protein production
    R(21) = beta_Q*Q;       % Q protein degradation
    R(22) = beta_p*p*Q;     % Q-dependent p53 degradation
    
    % Gillespie algorithm -------------------------------------------------
    Rtot = sum(R);  % Total propensity
    r1 = rand();    % To choose timing of next reaction
    RT = RT - log(r1)/Rtot; % Timing of next reaction
    r2 = rand();    % To choose the next reaction
    idx = find(cumsum(R)/Rtot >= r2,1); % Index of next reaction
    
    X = X + S(:,idx)';  % Update reactants level after reaction
    X(X<0) = 0; % Prevent negatives
    
    % Output --------------------------------------------------------------
    if (RT > click*ts) % Record outputs at defined intervals only
        sim = zeros(1,8);
        sim(1) = click*ts;          % Time vector
        sim(2) = X(1)+X(8)+X(9);    % Total p53
        sim(3) = X(3);              % m1
        sim(4) = X(4)+X(8);         % Total M1
        sim(5) = X(6);              % m2
        sim(6) = X(7)+X(9);         % Total M2
        sim(7) = X(10);             % q
        sim(8) = X(11);             % Q
        SimVal = [SimVal;sim];
        click = click + 1;
    end
    
end

%% Signal processing
t = SimVal(:,1);    % Time vector
pTot = movmean(SimVal(:,2),3);  % Total p53
m1 = SimVal(:,3);
M1Tot = movmean(SimVal(:,4),3); % Total M1
m2 = SimVal(:,5);
M2Tot = movmean(SimVal(:,6),3); % Total M2
q = SimVal(:,7);
Q = SimVal(:,8);

pTot = detrend(pTot,3);
M1Tot = detrend(M1Tot,3);
M2Tot = detrend(M2Tot,3);

%% Plot output
close all
c=colororder;
fg=figure('Position',[2 32 958 964]);
tiledlayout(5,1,'Padding','compact','TileSpacing','compact')

% p53
nexttile(1)
plot(t,pTot,'Color',c(1,:),'LineWidth',2)
ylabel('p53')
ylim([-20 40])

% Mdm2 protein allele #1
nexttile(2)
plot(t,M1Tot,'Color',c(2,:),'LineWidth',2)
ylabel('Mdm2 allele 1')
ylim([-40 80])

% Mdm2 protein allele #2
nexttile(3)
plot(t,M2Tot,'Color',c(3,:),'LineWidth',2)
ylabel('Mdm2 allele 2')
ylim([-40 80])

% q mRNA
nexttile(4)
plot(t,q,'Color','k','LineWidth',2)
ylabel('q mRNA')

% Q protein
nexttile(5)
plot(t,Q,'Color','k','LineWidth',2)
ylabel('Q protein')

%% Get statistics

% Peaks
[~,plocs]=findpeaks(pTot,'MinPeakProminence',2,'MinPeakDistance',3);
[~,M1locs]=findpeaks(M1Tot,'MinPeakProminence',1.3,'MinPeakDistance',3);
[~,M2locs]=findpeaks(M2Tot,'MinPeakProminence',1.3,'MinPeakDistance',3);

% Valleys
ptrlocs=arrayfun(@(i) find(pTot(plocs(i):plocs(i+1))...
    ==min(pTot(plocs(i):plocs(i+1))))+plocs(i)-1,1:length(plocs)-1);
M1trlocs=arrayfun(@(i) find(M1Tot(M1locs(i):M1locs(i+1))...
    ==min(M1Tot(M1locs(i):M1locs(i+1))))+M1locs(i)-1,1:length(M1locs)-1);
M2trlocs=arrayfun(@(i) find(M2Tot(M2locs(i):M2locs(i+1))...
    ==min(M2Tot(M2locs(i):M2locs(i+1))))+M2locs(i)-1,1:length(M2locs)-1);

% Mark peaks
nexttile(2)
hold on
ylimit=get(gca,'YLim');
plot(t(plocs),0.9*ylimit(2)*ones(size(plocs)),'|','LineStyle','none'...
    ,'MarkerFaceColor',c(1,:),'MarkerEdgeColor',c(1,:),'MarkerSize',8 ...
    ,'LineWidth',1.5);
plot(t(M1locs),0.9*ylimit(2)*ones(size(M1locs)),'o','LineStyle','none'...
    ,'MarkerFaceColor',c(2,:),'MarkerEdgeColor',c(2,:),'MarkerSize',4);

nexttile(3)
hold on
ylimit=get(gca,'YLim');
plot(t(plocs),0.9*ylimit(2)*ones(size(plocs)),'|','LineStyle','none'...
    ,'MarkerFaceColor',c(1,:),'MarkerEdgeColor',c(1,:),'MarkerSize',8 ...
    ,'LineWidth',1.5);
plot(t(M2locs),0.9*ylimit(2)*ones(size(M2locs)),'o','LineStyle','none'...
    ,'MarkerFaceColor',c(3,:),'MarkerEdgeColor',c(3,:),'MarkerSize',4);

%% Skip detection
skip_detection_allele1_single
skip_detection_allele2_single
