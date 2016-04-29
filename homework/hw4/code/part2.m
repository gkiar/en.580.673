load('mod_raw_025_NOI.mat')

figure

for i = 1:6
   subplot(2,3,i)
   dat = noise_data(:,1,1,i);
   ym = mean(imag(dat));
   ys = std(imag(dat));
   xm = mean(real(dat));
   xs = std(real(dat));
   
   scatter(real(dat), imag(dat), '.')
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
    hist(fiftflat(:,1,1,i), 1000)
end

magfitty = RSS(fift);
figure
for i =1:6
    subplot(3, 3, i)
    hist(RSS(fift(:,:,1, i)), 10000)
    xlabel(strcat('Mean:',num2str(mean(RSS(fift(:,:,1,i)))),'; std:',num2str(std(RSS(fift(:,:,1,i))))))
end
subplot(3,3,9)
hist(magfitty, 10000)
xlabel(strcat('Mean:',num2str(mean(magfitty)),'; std:',num2str(std(magfitty))))