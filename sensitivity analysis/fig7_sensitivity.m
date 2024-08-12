%% Initialize environment
clear; close all; clc

%% create table for the perturbation of individual parameters
[name_pars,~] = get_original_value(1,0,1,1); % get name of parameters
name_pars = strrep(name_pars,"_","");
name_pars = [name_pars+"_neg",name_pars+"_pos","ref"];

[X1,X2,X3,X4] = ndgrid(name_pars,1:6,[0 1],[0 1]);
tbl = table(X1(:),X2(:),X3(:),X4(:));
tbl.Properties.VariableNames = ["name_perturbation","nAlelles","isQ","isConstp53"];
tbl = tbl(:,[4 3 2 1]);

perturbation = cell(height(tbl),1);
for i = 1:height(tbl)
    OME = 1;
    TC = 1;
    [name_pars,SimPars] = get_original_value(tbl.nAlelles(i),tbl.isConstp53(i),OME,TC);
    name_pars = strrep(name_pars,"_","");

    name_perturbation = tbl.name_perturbation(i);

    idx = find(name_pars==extractBefore(name_perturbation,"_")); % index of the parameter to be perturbed

    if extractAfter(name_perturbation,"_")=="neg" % 0.5X
        factor = 0.5;
    elseif extractAfter(name_perturbation,"_")=="pos" % 2X
        factor = 2;
    else % 1X
        factor = 1; % factor of perturbation
    end

    SimPars(idx) = SimPars(idx).*factor;
    perturbation{i} = SimPars;
end
tbl.perturbation = perturbation;

clearvars -except tbl

idx = find(tbl.isConstp53==0 & tbl.isQ==1 & tbl.nAlelles==2);

ct = load("simulation_result_v2_ct.mat");
kd1 = load("simulation_result_v2_kd1.mat");

ct.tbl.Properties.VariableNames(6:8) = "ct_"+ct.tbl.Properties.VariableNames(6:8);
kd1.tbl.Properties.VariableNames(6:8) = "kd1_"+kd1.tbl.Properties.VariableNames(6:8);

tbl = [ct.tbl,kd1.tbl(:,6:8)];
tbl = tbl(idx,:);

clearvars -except tbl

%% summarize data
% s = RandStream('dsfmt19937','Seed',4088);
% id_sim = randperm(s,2000,500);

% the first moment (i.e., mean of amplitude and period)
tbl2 = tbl;
for i = 6:11
    v = tbl2.Properties.VariableNames(i);
    tbl2.(v{1}) = cellfun(@(s)mean(s),tbl2.(v{1}));
end

% the second moment (i.e., standard deviation of amplitude and period)
tbl5 = tbl;
for i = 6:11
    v = tbl5.Properties.VariableNames(i);
    tbl5.(v{1}) = cellfun(@(s)std(s),tbl5.(v{1}));
end

%% multiple allele, sensitivity of p53 amplitude/period, Q, control p53 or not
tbl3 = []; tbl6 = [];
[name_paras,~] = get_original_value(1,1,1,1);
for i = 1:length(name_paras)
    perturbation = strrep(name_paras(i),"_","");
    perturbation = [perturbation+"_neg","ref",perturbation+"_pos"];
    for j = 2 % number of alleles

        %%% no control p53, Q
        idx = tbl.isConstp53==0&tbl.isQ==1&ismember(tbl.name_perturbation,perturbation)&tbl.nAlelles==j;
        tmp = tbl2(idx,:);
        tbl3 = [tbl3;tmp.isConstp53(1),tmp.isQ(1),tmp.nAlelles(1),i,....
            {tmp.ct_p53period([1 3 2])'} {tmp.ct_skip1_prob([1 3 2])'},...
            {tmp.ct_p53amp([1 3 2])'},...
            {tmp.kd1_p53period([1 3 2])'} {tmp.kd1_skip1_prob([1 3 2])'},...
            {tmp.kd1_p53amp([1 3 2])'}];  

        tmp = tbl5(idx,:);
        tbl6 = [tbl6;tmp.isConstp53(1),tmp.isQ(1),tmp.nAlelles(1),i,....
            {tmp.ct_p53period([1 3 2])'} {tmp.ct_skip1_prob([1 3 2])'},...
            {tmp.ct_p53amp([1 3 2])'},...
            {tmp.kd1_p53period([1 3 2])'} {tmp.kd1_skip1_prob([1 3 2])'},...
            {tmp.kd1_p53amp([1 3 2])'}];  
    end
end
tbl3 = array2table(tbl3);
tbl3.Properties.VariableNames = ["isConstp53" "isQ" "nAlelles" "name_parameter",...
    "ct_p53period","ct_skip1_prob","ct_p53amp",...
    "kd1_p53period","kd1_skip1_prob","kd1_p53amp"];
tbl3.isConstp53 = cellfun(@(s)s,tbl3.isConstp53);
tbl3.isQ = cellfun(@(s)s,tbl3.isQ);
tbl3.nAlelles = cellfun(@(s)s,tbl3.nAlelles);
tbl3.name_parameter = cellfun(@(s)s,tbl3.name_parameter);
[G,TID] = findgroups(tbl3(:,1:4));
[~,I] = sort(G);
tbl3 = tbl3(I,:);

tbl6 = array2table(tbl6);
tbl6.Properties.VariableNames = ["isConstp53" "isQ" "nAlelles" "name_parameter",...
    "ct_p53period","ct_skip1_prob","ct_p53amp",...
    "kd1_p53period","kd1_skip1_prob","kd1_p53amp"];
tbl6.isConstp53 = cellfun(@(s)s,tbl6.isConstp53);
tbl6.isQ = cellfun(@(s)s,tbl6.isQ);
tbl6.nAlelles = cellfun(@(s)s,tbl6.nAlelles);
tbl6.name_parameter = cellfun(@(s)s,tbl6.name_parameter);
[G,TID] = findgroups(tbl6(:,1:4));
[~,I] = sort(G);
tbl6 = tbl6(I,:);

tbl4 = tbl3;
for i = 5:10
    v = tbl3.Properties.VariableNames(i);
    tbl4.(v{1}) = cellfun(@(s)sum(abs(diff(s))),tbl3.(v{1}));
end

tbl7 = tbl6;
for i = 5:10
    v = tbl6.Properties.VariableNames(i);
    tbl7.(v{1}) = cellfun(@(s)sum(abs(diff(s))),tbl6.(v{1}));
end

tbl8 = [groupsummary(tbl4,"isConstp53","none","mean",tbl4.Properties.VariableNames(5:10));
groupsummary(tbl7,"isConstp53","none","mean",tbl4.Properties.VariableNames(5:10))];

%% the first moment
label_para = ["{\itα}_{\itp}" "{\itβ}_{\itp}",...
    "{\itα}_{\itT}" "{\itβ}_{\itT}",...    
    "{\itα}_{\itm}" "{\itβ}_{\itm}",...
    "{\itα}_{\itM}" "{\itβ}_{\itM}",...
     "{\itφ}" "{\itΨ}" ,...
    "{\itα}_{\itq}" "{\itβ}_{\itq}" "{\itα}_{\itQ}" "{\itβ}_{\itQ}"];

close all
figure('Position',[2 542 1400/2 460/2])
y = [tbl4.kd1_p53amp,...
    tbl4.ct_p53amp];
b = bar(y,1);
b(2).FaceColor = [0.6 0.6 0.6];
b(1).EdgeColor = [1 1 1];
b(2).EdgeColor = [1 1 1];
xticks(1:14)
xticklabels(label_para)
legend(["1-allele KD" "Ctrl"])
set(gca,'Layer','top','TickDir','out','LineWidth',1,'fontsize',16,...
    'XTickLabelRotation',0)
ylabel(["Sensivity"])
title("Mean of amplitude")
axis padded
box off

% exportgraphics(gcf,'figS7A_ampMean.pdf','ContentType','vector')

figure('Position',[2+200 542 1400/2 460/2])
y = [tbl4.kd1_p53period,...
    tbl4.ct_p53period];
b = bar(y,1);
b(2).FaceColor = [0.6 0.6 0.6];
b(1).EdgeColor = [1 1 1];
b(2).EdgeColor = [1 1 1];
xticks(1:14)
xticklabels(label_para)
set(gca,'Layer','top','TickDir','out','LineWidth',1,'fontsize',16,...
    'XTickLabelRotation',0)
ylabel(["Sensivity"])
title("Mean of period")
axis padded
box off

% exportgraphics(gcf,'figS7A_perMean.pdf','ContentType','vector')

%% the second moment
label_para = ["{\itα}_{\itp}" "{\itβ}_{\itp}",...
    "{\itα}_{\itT}" "{\itβ}_{\itT}",...    
    "{\itα}_{\itm}" "{\itβ}_{\itm}",...
    "{\itα}_{\itM}" "{\itβ}_{\itM}",...
     "{\itφ}" "{\itΨ}" ,...
    "{\itα}_{\itq}" "{\itβ}_{\itq}" "{\itα}_{\itQ}" "{\itβ}_{\itQ}"];

close all
figure('Position',[2 542 1400/2 460/2])
y = [tbl7.kd1_p53amp,...
    tbl7.ct_p53amp];
b = bar(y,1);
b(2).FaceColor = [0.6 0.6 0.6];
b(1).EdgeColor = [1 1 1];
b(2).EdgeColor = [1 1 1];
xticks(1:14)
xticklabels(label_para)
legend(["1-allele KD" "Ctrl"])
set(gca,'Layer','top','TickDir','out','LineWidth',1,'fontsize',16,...
    'XTickLabelRotation',0)
ylabel(["Sensivity"])
title("S.D. of amplitude")
axis padded
box off

% exportgraphics(gcf,'figS7A_ampSD.pdf','ContentType','vector')

figure('Position',[2+200 542 1400/2 460/2])
y = [tbl7.kd1_p53period,...
    tbl7.ct_p53period];
b = bar(y,1);
b(2).FaceColor = [0.6 0.6 0.6];
b(1).EdgeColor = [1 1 1];
b(2).EdgeColor = [1 1 1];
xticks(1:14)
xticklabels(label_para)
set(gca,'Layer','top','TickDir','out','LineWidth',1,'fontsize',16,...
    'XTickLabelRotation',0)
ylabel(["Sensivity"])
title("S.D. of period")
axis padded
box off

% exportgraphics(gcf,'figS7A_perSD.pdf','ContentType','vector')

%%
close all
figure('Position',[931,50,900/3,630/2])
y = [tbl8.mean_kd1_p53amp,...
    tbl8.mean_ct_p53amp];
b = bar(y,1);
b(2).FaceColor = [0.6 0.6 0.6];
b(1).EdgeColor = [1 1 1];
b(2).EdgeColor = [1 1 1];
xticks(1:2)
xticklabels(["Mean" "S.D."])
set(gca,'Layer','top','TickDir','out','LineWidth',1,'fontsize',16,...
    'XTickLabelRotation',0)
ylabel("Average sensitivity")
title("p53 amplitude")
axis padded
box off

legend(["1-allele KD" "Ctrl"])

% exportgraphics(gcf,'fig7A_amp.pdf','ContentType','vector')

%%
close all
figure('Position',[931,50,900/3,630/2])
y = [tbl8.mean_kd1_p53period,...
    tbl8.mean_ct_p53period];
b = bar(y,1);
b(2).FaceColor = [0.6 0.6 0.6];
b(1).EdgeColor = [1 1 1];
b(2).EdgeColor = [1 1 1];
xticks(1:2)
xticklabels(["Mean" "S.D."])
set(gca,'Layer','top','TickDir','out','LineWidth',1,'fontsize',16,...
    'XTickLabelRotation',0)
ylabel("Average sensitivity")
title("p53 period")
axis padded
box off

% exportgraphics(gcf,'fig7A_per.pdf','ContentType','vector')