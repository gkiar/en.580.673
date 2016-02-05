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

im = fftshift(ifft2(ifftshift(data)));
im_rss = sqrt(sum(abs(im).^2, 3));

figure
subplot(2,1,1)
imagesc(abs(im_rss)); axis off; axis equal
colormap(gray)

subplot(2,1,2)
imagesc(angle(im(:,:,1))); axis off; axis equal
colormap(gray)
