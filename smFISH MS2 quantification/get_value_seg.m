%% load nucleus segmentation
fname = fname_seg;
[filepath,name,ext] = fileparts(fname);
info = imfinfo(fname);
seg = imread(fname, 'Info', info);

%% generate foci dynamics for each tracked cell
idx = NaN(size(col));
for j = 1:size(col,1)
    value_seg = seg(fix(row(j)),fix(col(j)));
    if value_seg~=0
        idx(j) = value_seg;
    else
        idx(j) = nan;
    end
end

value_seg = NaN(size(col));
for j = 1:size(col,1)
    value_seg(j) = idx(j);
end
