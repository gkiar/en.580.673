clear all
close all
load('whatisthis_8coil')

%% part 1
figure
for i=1:size(data,3)
    subplot(3,3,i)
    imagesc(log(abs(data(:,:,i))+1))
    colormap(gray); axis off;
end

figure
for i=1:size(data,3)
    subplot(3,3,i)
    imagesc(angle(data(:,:,i)))
    colormap(gray); axis off;
end

%% part 2

im = zeros(size(data));
for i=1:size(data,3)
    temp = fftshift(fft2(ifftshift(data(:,:,i))));
    im(:,:,1) = temp;
end

im_rss = zeros(size(im(:,:,1)));
for i=1:size(im,1)
    for j=1:size(im,2)
        im_rss(i,j) = sqrt(sum(abs(im(i,j,:)).^2));
    end
end

figure
imagesc(im_rss); axis off
colormap(gray)
caxis([(min(min(im_rss))/10) max(max(im_rss))])
