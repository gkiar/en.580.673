function [im, im_rss] = Reconstruct(coil_data, plot, first)

if nargin <2 
    plot = false;
end
if nargin <3
    first = false;
end

if first
    coil_data = sum(coil_data, 3);
    im = fftshift(ifft2(ifftshift(coil_data)));
    im_rss = RSS(im);
else
    coil_data = sum(coil_data, 5);
    im = fftshift(ifft2(ifftshift(coil_data)));
    im_rss = RSS(im);
end

if plot 
    figure
    colormap(gray)
    subplot(131)
    slize = ceil(size(im,3)/2);
    imagesc(log(abs(im(:,:,slize,1))+1)); axis off; axis image
    title('Single coil magnitude')

    subplot(133)
    imagesc(log(im_rss(:,:,slize)+1)); axis off; axis image
    title('RSS magnitude')

    subplot(132)
    imagesc(angle(im(:,:,slize))); axis off; axis image
    title('Single coil phase')
    
end

end