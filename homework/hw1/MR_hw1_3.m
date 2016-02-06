clear all
close all
load('whatisthis_8coil')

%% part 1
figure(1)
for i=1:size(data,3)
    subplot(3,3,i)
    imagesc(log(abs(data(:,:,i))+1))
    colormap(gray); axis off;
end

figure(2)
for i=1:size(data,3)
    subplot(3,3,i)
    imagesc(angle(data(:,:,i)))
    colormap(gray); axis off;
end

%% part 2

im = fftshift(ifft2(ifftshift(data)));
im_rss = sqrt(sum(abs(im).^2, 3));

figure(3)
subplot(221)
imagesc(log(abs(im(:,:,1))+1)); axis off; axis equal
colormap(gray); title('Single coil magnitude image')

subplot(222)
imagesc(log(im_rss+1)); axis off; axis equal
colormap(gray); title('Reconstructed magnitude image')

subplot(223)
imagesc(angle(im(:,:,1))); axis off; axis equal
colormap(gray); title('Single coil phase image')

angi = abs(acos(sum(abs(data).*cos(angle(data)),3)./im_rss));
subplot(224)
imagesc(angi); axis off; axis equal
colormap(gray); title('Reconstructed phase image')