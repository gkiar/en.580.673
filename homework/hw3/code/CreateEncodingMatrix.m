
 clear all
 close all
 clc
 
%% Setup encoding/decoding matrices

N = 1024;


%a = ones(N,1);  % the signal (representing 5 equally spaced pixels in our FOV
a = (1:N).';

W = -1i*2*pi/N;

W = repmat(W,N,N);


n1 = (0:N-1)';
n1 = repmat(n1,1,N);

n2 = 0:N-1;
n2 = repmat(n2,N,1);

Wforward =            exp(   W.*n1.*n2);   % defined as in MATLAB's fft

  tic 
Winverse_calc = (1/N )*exp(-1*W.*n1.*n2);  % defined as in MATLAB's fft
time_calc = toc;

tic
Winverse_actual = inv(Wforward);
time_pinv = toc;
Winverse = Winverse_calc;

A1=fft(a);
A2=Wforward*a;

a1 = ifft(A1);
a2 = Winverse*A2;

% disp('Encoded vectors fft | Wforward');
% disp(num2str( [ A1 , A2 ]));
% fprintf('\n\n');
% 
% disp('Recovered vectors: original | ifft | Winverse');
% disp([ num2str( [a, a1 , a2 ] ) ]); 
% fprintf('\n\n');
% 
% disp('Encoding Matrix')
% disp(Wforward)
% fprintf('\n\n');
% 
% disp('Row-reduced Echelon form of the Encoding Matrix')
% disp(rref(Wforward))
% fprintf('\n\n');
% 
% disp('Rank of the Encoding Matrix')
% disp(rank(Wforward))
% fprintf('\n\n');
% 
% disp('Condition #')
% disp(cond(Wforward))
% fprintf('\n\n');
% 
% disp('Derived Inverse Matrix:');
% disp(Winverse_actual)
% fprintf('\n\n');
% 
% disp('Calculated Inverse Matrix:');
% disp(Winverse_calc)
% fprintf('\n\n');


% disp(' ... paused ...'); 
% pause
time_calc
time_pinv
%% Plot (after acquisition)

figure('Name', 'Net magnetization after DFT Encoding');
colororder = 'rbmg';
colororder = repmat(colororder, [1, ceil(N/length(colororder))+1]);
leg_string ={};
h_m=zeros(1,N);
h = zeros(1,N);
Mnet = zeros(1,N);

fs = 20;
	
subplot(1,2,1);
% plot the a vector
for j=1:length(a)
	plot( [0, real(a(j)) ], [ 0 imag(a(j)) ], ...
		[colororder(j), ':o'], 'linewidth', 1);
	hold on;
end

axis equal;
set(gca, 'fontsize', fs);
xlabel('x', 'fontsize', fs);
ylabel('y', 'fontsize', fs);
title(' M_{net} - the Measured Signal');
l = sum(a(:));
xlim([-l,l]); ylim([-l,l]);


for n=1:N
	M = Wforward(n,:) .* a';
	Mnet(n) = sum(M(:));
	
	subplot(1,2,1)
	h_m(n) = plot([0, real( Mnet(n) )], [ 0 imag( Mnet(n) ) ], ...
			'k-o', ...
			'linewidth', 2, 'markerfacecolor', colororder(n));
	
		for j = 1:length(M)
			h(j) = plot([0, real( M(j)) ], [ 0 imag(M(j)) ], ...
				[colororder(j),'-o'], ...
				'linewidth', 1, 'markerfacecolor', colororder(j));
		end
			
	leg_string = cat(1,leg_string, ['Wf row #', num2str(n)]);

	subplot(2,2,2);
	plot(n , abs(Mnet(n)), 'ko');
	
	if n==1
		xlim([0,N]);
		xlabel(' Measurement Index (n)');
		ylabel('|Mnet|');
		title('|Mnet|');
		hold on
	end;
		
	subplot(2,2,4);
	plot(n,  atan2(imag(Mnet(n)), real(Mnet(n)) ), 'ro');
	hold on;
	
	if n==1
		xlim([0,N]);
		xlabel(' Measurement Index (n)');
		ylabel('|Mnet|');
		title('phase Mnet');
		hold on
	end;
	
		
	shg
	if (N<=32)
		pause
	end
	delete(h)
end

a_rec = Winverse*Mnet.';

subplot(1,2,1)
axis equal;
set(gca, 'fontsize', fs);
xlabel('x', 'fontsize', fs);
ylabel('y', 'fontsize', fs);
title(' M_{net} - the Measured Magnetization');
l = sum(a(:));
xlim([-l,l]); ylim([-l,l]);

disp('Net signal measured: ')
disp(num2str(Mnet)); 
fprintf('\n\n');

disp('Coefficient vector | recovered')
disp(num2str([a , abs(a_rec)] ))
fprintf('\n\n');

disp('Rotation Angles induced by EncodingMatrix');
disp(atan2(imag(Wforward),  real(Wforward)) * 360 / (2*pi));
%%

figure('Name', ' FFT derived values')
subplot(2,1,1)
plot(abs(fft(a)), 'k.-')
xlabel('Frequency');
ylabel('|fft(a)|');

subplot(2,1,2)
plot(atan2(imag(fft(a)), real(fft(a))), 'r.-')
xlabel('Frequency');
ylabel('phase(fft(a))');



