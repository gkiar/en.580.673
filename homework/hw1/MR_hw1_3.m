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
IM = zeros(size(data(:,:,1)));
for i=1:size(data,1)
    for j=1:size(data,2)
        IM(i,j) = sqrt(sum(abs(data(i,j,:)).^2));
    end
end

% h = fspecial('gaussian', 100, 0.1);
% IM = imfilter(IM, h);

figure
subplot(1,3,1)
imagesc(log(IM + 1)); axis off

IM = ifftshift(IM);
im = ifft2(IM);
im = fftshift(im);

subplot(1,3,2)
imagesc(log(abs(im)+1)); axis off
subplot(1,3,3)
imagesc(angle(im)); axis off
colormap(gray) % it looks quite pinappular
