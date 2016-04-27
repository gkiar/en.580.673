load('mod_raw_025_NOI.mat')

figure

for i = 1:6
   subplot(2,3,i)
   dat = noise_data(:,1,1,i);
   ym = mean(imag(dat));
   ys = std(imag(dat));
   xm = mean(real(dat));
   xs = std(real(dat));
   
   for j = 1:length(dat)/4
        plot([0, real(dat(j))], [0, imag(dat(j))])
        hold on 
   end
   axis([-40 40 -40 40])
   xlabel(strcat('Mean: ', num2str(xm), '; STD: ', num2str(xs)))
   ylabel(strcat('Mean: ', num2str(ym), '; STD: ', num2str(ys)))
end


flat = [real(noise_data); imag(noise_data)];
figure
for i = 1:6
    subplot(2, 3, i)
    hist(flat(:,1,1,i), 10000)
    xlim([-30 30])
end

fift = fft(noise_data);
fiftflat = [real(fift); imag(fift)]/250;
figure
for i = 1:6
    subplot(2, 3, i)
    hist(fiftflat(:,1,1,i), 100)
end

magfitty = RSS(fift);
figure
subplot(121)
imagesc(magfitty)
subplot(122)
hist(magfitty, 100)
