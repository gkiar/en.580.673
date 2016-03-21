%Script to create DFT matrices

 clear all
 close all
 clc
 
%%

N = 3;


%a = ones(N,1);  % the signal (representing 5 equally spaced pixels in our FOV
a = (1:N).';
%a = circshift(a,[1 1])

W = -1i*2*pi/N;

W = repmat(W,N,N);


n1 = (0:N-1)';
n1 = repmat(n1,1,N);

n2 = 0:N-1;
n2 = repmat(n2,N,1);

Wforward = exp(W.*n1.*n2);            % defined as in MATLAB's fft
	
Winverse = (1/N )*exp(-1*W.*n1.*n2);  % defined as in MATLAB's fft

A1=fft(a);
A2=Wforward*a;

a1 = ifft(A1);
a2 = Winverse*A2;

disp('Encoded vectors fft | Wforward');
disp(num2str( [ A1 , A2 ]));
sprintf('\n\n');

disp('Recovered vectors: origina | ifft | Winverse');
disp(num2str( [a, a1 , a2 ]));
sprintf('\n\n');


%% Plot Experiment 2 (after acquisition)

figure('Name', 'Net magnetization after DFT Encoding');
colororder = 'rbmg';
colororder = repmat(colororder, [1, ceil(N/length(colororder))+1]);
leg_string ={};
h_m=zeros(1,N);
h = zeros(1,N);
Mnet = zeros(1,N);

fs = 20;
	
% plot the a vecto
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
	

	h_m(n) = plot([0, real( Mnet(n) )], [ 0 imag( Mnet(n) ) ], ...
			'k-o', ...
			'linewidth', 2, 'markerfacecolor', colororder(n));
	
		for j = 1:length(M)
			h(j) = plot([0, real( M(j)) ], [ 0 imag(M(j)) ], ...
				[colororder(j),'-o'], ...
				'linewidth', 1, 'markerfacecolor', colororder(j));
		end
		
		
	leg_string = cat(1,leg_string, ['Wf row #', num2str(n)]);

	shg
	pause
	delete(h)
end

a_rec = Winverse*Mnet.';

axis equal;
set(gca, 'fontsize', fs);
xlabel('x', 'fontsize', fs);
ylabel('y', 'fontsize', fs);
title(' M_{net} - the Measured Magnetization');
l = sum(a(:));
xlim([-l,l]); ylim([-l,l]);

disp('Net signal measured: ')
disp(num2str(Mnet)); 

disp('Coefficient vector | recovered')

disp(num2str([a , abs(a_rec)] ))



