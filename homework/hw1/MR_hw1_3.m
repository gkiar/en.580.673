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

im = ifftshift(fft2(fftshift(data)));
im_rss = zeros(size(im(:,:,1)));
for i=1:size(im,1)
    for j=1:size(im,2)
        im_rss(i,j) = sqrt(sum(im(i,j,:).^2));
    end
end

figure
subplot(1,2,1)
imagesc(abs(im_rss)); axis off
colormap(gray)

subplot(1,2,2)
imagesc(angle(im_rss)); axis off
colormap(gray)
