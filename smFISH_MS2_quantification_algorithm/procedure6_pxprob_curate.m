%% Initialize environment
clear; close all; clc

%% Quantify the signal intenstiy
ms2_data=[];
smFISH_data=[];
idx_data=[];
ksp_data=[];
img_data=[];

for pos=1:9
    % 1. get value_seg for all tracked cells
    fname_seg = sprintf('./Track data/seg_xy%d.tif',pos);
    load(sprintf('./Track data/centroid_xy%d.mat',pos));
    % use analyzed tracking file
    get_value_seg

    % 2. detect foci
    fname_label = sprintf('./Track data/seg_xy%d.tif',pos);
    fname_foci = sprintf('./Track data/MS2_xy%d_aligned.tif',pos);
    fname_pxprob...
        = sprintf('./Track data/MS2_xy%d_aligned_Probabilities.tiff',pos);
    fname_smFISH = sprintf('./Track data/Cy5_xy%d.tif',pos);
    pixel_filter5_pxprob_curate
    
    % collate data
    ms2_data=[ms2_data;ms2_array];
    smFISH_data=[smFISH_data;smFISH_array];
    idx_data=[idx_data;cen_idx];
    ksp_data=[ksp_data;ksp];
    img_data=[img_data;img_array];

end

% 3. Save
badQuant=or(isnan(ms2_data),isinf(ms2_data));
smFISH_data(badQuant)=[];
ms2_data(badQuant)=[];
idx_data(badQuant)=[];
ksp_data(badQuant)=[];
img_data(badQuant)=[];

%% Keep valid cells
validCells = log10(ksp_data)>-40;
ms2_data=ms2_data(validCells);
smFISH_data=smFISH_data(validCells);
idx_data=idx_data(validCells); 
ksp_data=ksp_data(validCells);
img_data=img_data(validCells);

%% Downsampling
rng(2)
[Y,E]=discretize(log(ms2_data),(-6.25:0.5:2));

binned_ms2=cell(1,max(Y));
binned_FISH=cell(1,max(Y));
binned_idx=cell(1,max(Y));
binned_img=cell(1,max(Y));

for i=1:max(Y)    
    binned_ms2{i}=ms2_data(Y==i); 
    binned_FISH{i}=smFISH_data(Y==i);
    binned_idx{i}=idx_data(Y==i,:);
    binned_img{i}=img_data(Y==i,:);
    
    if length(binned_ms2{i})>30
        idx=randperm(length(binned_ms2{i}),30);
        binned_ms2{i}=binned_ms2{i}(idx);
        binned_FISH{i}=binned_FISH{i}(idx);
        binned_idx{i}=binned_idx{i}(idx,:);
        binned_img{i}=binned_img{i}(idx,:);
    end
    
end

ms2_data_ds=[];
smFISH_data_ds=[];
idx_data_ds=[];
img_data_ds=[];

for i=1:length(binned_ms2)
    
    ms2_data_ds=[ms2_data_ds;binned_ms2{i}];
    smFISH_data_ds=[smFISH_data_ds;binned_FISH{i}];
    idx_data_ds=[idx_data_ds;binned_idx{i}];
    img_data_ds=[img_data_ds;binned_img{i}];
    
end

%% curate
badCells=[2,6,7,9,12,26,35,49,78,79,123,126];
ms2_data_ds(badCells)=[];
smFISH_data_ds(badCells)=[];
idx_data_ds(badCells)=[];
img_data_ds(badCells)=[];


