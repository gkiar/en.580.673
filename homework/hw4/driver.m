%% The one script to rule them all

clear
close all

%% Reconstruction

fnames = dir('*STD.mat');
ims = {};
ims_rms = {};
plot = true;
for i=1:length(fnames)-1
    load(fnames(i).name)
    if i <= 3
        first = true;
    else
        first = false;
    end
    [ims{i}, ims_rms{i}] = Reconstruct(fdata, plot, first);
end