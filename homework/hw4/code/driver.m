%% The one script to rule them all

clear
close all

%% Reconstruction

fnames = dir('*STD.mat');
ims = {};
ims_rms = {};
plot = true;

for i = 13:length(fnames)
    load(fnames(i).name)
    if i <= 3=
        first = true;
    else
        first = false;
    end
    [ims{i}, ims_rms{i}] = Reconstruct(fdata, plot, first);
    
    namez = strsplit(fnames(i).name, '_');
    namez = namez(3);
    xlabel(namez)
    subplot(133);

    fg = roipoly();
    bg = roipoly();
    
    sloc = find(fg == 1);
    nloc = find(bg == 1);
    
    sig{i} = mean(ims_rms{i}(sloc));
    noi{i} = std(ims_rms{i}(nloc));
    
    SNR{i} = sig{i}/noi{i}
end
% panel plot
% title is the scan; and snr


%% thingymagiggy
load('mod_raw_017_STD.mat')
ims = {};
ims_rms = {};
plot = true;
first = false;

for i = 1:4
    [ims{i}, ims_rms{i}] = Reconstruct(fdata(:,:,:,:,1:i), plot, first);
    
    subplot(133);

    fg = roipoly();
    bg = roipoly();
    
    sloc = find(fg == 1);
    nloc = find(bg == 1);
    
    sig{i} = mean(ims_rms{i}(sloc));
    noi{i} = std(ims_rms{i}(nloc));
    
    SNR{i} = sig{i}/noi{i}
end