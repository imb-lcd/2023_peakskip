%% load nucleus segmentation
fname = fname_label;
[filepath,name,ext] = fileparts(fname);
info = imfinfo(fname);
seg = imread(fname, 'Info', info);

%% load foci channel
fname = fname_foci;
[filepath,name,ext] = fileparts(fname);
info = imfinfo(fname);
mCherry = imread(fname, 'Info', info);
mCherry(mCherry==0) = NaN;

%% load pixel probability
fname = fname_pxprob;
[filepath,name,ext] = fileparts(fname);
info = imfinfo(fname);
px_prob = imread(fname, 'Info', info);


%% load smFISH
fname = fname_smFISH;
[filepath,name,ext] = fileparts(fname);
info = imfinfo(fname);
smFISH = imread(fname, 'Info', info);

%% load seg
fname = fname_seg;
[filepath,name,ext] = fileparts(fname);
info = imfinfo(fname);
seg = imread(fname, 'Info', info);
% vn = value_seg;
nTrackedCell = length(value_seg);

%% filter pixels
v=(1:1:nTrackedCell)';

ms2_array=NaN(size(v));
smFISH_array=NaN(size(v));
img_array=cell(size(v));
ms2_mean=NaN(size(v));
ms2_median=NaN(size(v));
ksp=NaN(size(v));
cen_idx=[pos*ones(size(v)),NaN(size(v))];
SE=strel('disk',25);
% SE1=strel('disk',2);

close all
parfor u=1:length(v)
    
    tmp = seg==value_seg(u);
%     tmp = imerode(tmp,SE1);
    if any(tmp,'all') % if the segmented nucleus matches the tracked cell
        z = double(imgaussfilt(mCherry,1)).*tmp;
        p = double(imgaussfilt(px_prob,1)).*tmp;
        f = double(imgaussfilt(smFISH,1)).*tmp;
        
        % Find brightest pixels
        z1=z;
        z2=z1.*p;
        z3=reshape(z2,[],1);
        
        [max_px,max_idx]=sort(z3,'descend');
        max_px=max_px(1:3);
        
        % Calculate MS2 intensity
        z1(tmp==0)=NaN;
        background=median(z1,'all','omitnan');
        ms2_array(u)=mean(max_px)/background;
        
        % Find smFISH signal
        [r,c] = ind2sub(size(z),max_idx(1));
        mask = zeros(size(f));
        mask(r,c) = 1;
        mask=imdilate(mask,SE);
        f1 = f;
        f2 = f1.*mask; 
        f3=reshape(f2,[],1);        
        [max_px,max_idx]=sort(f3,'descend');
        
        % Calculate smFISH signal
        f1(tmp==0)=NaN;
        background=median(f1,'all','omitnan');
        smFISH_array(u)=mean(max_px(1:3))/background;
        
        % Get index
        cen_idx(u,2)=value_seg(u);
        
        % Get mean/median mCherry signal
        z1(tmp==0)=NaN;
        ms2_mean(u)=mean(z1,'all','omitnan');        
        ms2_median(u)=median(z1,'all','omitnan');
        
        if std(z1,[],'all','omitnan')>0
        [h,ksp(u),ksstat,cv]...
            =kstest((z1-mean(z1,'all','omitnan'))...
            /std(z1,[],'all','omitnan'));
        end
        
        stats=regionprops(tmp,'centroid');
        ctr=stats.Centroid;
        ctr=round(ctr);
        
        ctr(ctr<=40)=41;
        ctr(ctr>=1891)=1851;
        
        z4=z1(ctr(2)-40:ctr(2)+40,ctr(1)-40:ctr(1)+40);
        img_array{u}=z4;
    end
    
end