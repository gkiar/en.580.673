function [Wx_encode, Wx_decode, Wy_encode, Wy_decode] = PlotIdealEncodingMatrices(Nx, Ny, A)

if nargin<2
	disp('Not enough input arguments');
	return
end

if nargin <3
	A = ones(Ny,Nx);
end


%% Encoding using the Fourier matrices -   
% confirm that the correct elements were calculated during experiment 5
  
% Demonstrate that the matrices we calculate using the measurements above
% are exctly the same ones that we calculate when we estimate the DFT
% matrices for Nx and Ny respectively
  
% calculate DFT matrix for Nx; apply Wx_encode (a.k.a. forward) to encode
% each row of A
Wx = -1i*2*pi/Nx;
Wx = repmat(Wx,Nx,Nx);

n1 = (0:Nx-1)';
n1 = repmat(n1,1,Nx);

n2 = 0:Nx-1;
n2 = repmat(n2,Nx,1);

Wx_encode =        exp(   Wx.*n1.*n2);   % defined as in MATLAB's fft
Wx_decode = (1/Nx)*exp(-1*Wx.*n1.*n2);   % defined as in MATLAB's fft

% calculate DFT matrix for Ny; apply Wy_encode (a.k.a. forward) to encode
% each column of A
Wy = -1i*2*pi/Ny;
Wy = repmat(Wy,Ny,Ny);

n1 = (0:Ny-1)';
n1 = repmat(n1,1,Ny);

n2 = 0:Ny-1;
n2 = repmat(n2,Ny,1);

Wy_encode =        exp(   Wy.*n1.*n2);   % defined as in MATLAB's fft
Wy_decode = (1/Ny)*exp(-1*Wy.*n1.*n2);  % defined as in MATLAB's fft

FxA  = fftshift( (Wx_encode*(ifftshift(A  ,2).')).',2);  % equivalent to applying and FFT in the x direction for every row
FxyA = fftshift( (Wy_encode* ifftshift(FxA,1)   )  ,1 );  % equivalent to appplying an FFT in the y direction for every col

% just for fun, recover:

FyA = fftshift( (Wx_decode*(ifftshift(FxyA  ,2).')).',2);  % equivalent to applying and FFT in the x direction for every row
Axy = fftshift( (Wy_decode* ifftshift(FyA,1)   )  ,1 );

% Equivalent expressions:
%    FxyA = fftshift(fft2(ifftshift(A)))
%  or
%    FxyA = fftshift( Wy_encode * ((Wx_encode*(ifftshift(A).')) .'))
%
% and to confirm recover A
%
%    Axy = fftshift(ifft2(ifftshift(FAxy)))
%  or
%    Axy = fftshift( Wy_decode * ((Wx_decode*(ifftshift(FAxy).')) .'))

% disp('A:');
% disp(num2str(A,'%2f.0'));
% disp('Axy (recovered):');
% disp(num2str(Axy,'%2f.0'));
% disp('Error in recovery:')
% disp(norm(A-Axy))


% Note order of iffshift, fftshift to insure that the center of the image
% returns to the center ...

%%  Visualize Encoding Matrices
	
VisualizeEncodingMatrix(Wx_encode, 'E5: Ideal Wx-Encode');
VisualizeEncodingMatrix(fftshift(Wx_encode), 'E5: Ideal FFTSHIFT Wx-Encode');
VisualizeEncodingMatrix(Nx*Wx_decode, 'E5: Ideal Wx-Decode');
VisualizeEncodingMatrix(Wy_encode, 'E5: Ideal Wy-Encode');
VisualizeEncodingMatrix(fftshift(Wy_encode), 'E5: Ideal FFTSHIFT Wy-Encode');
VisualizeEncodingMatrix(Ny*Wy_decode, 'E5: Ideal Wy-Decode');

	
% %% Plot the  Fourier Encoded Vectors - DFT Matrix
% % compare this to the magnetization as encoded with your EncodingMatrix
% % from E5
% 
% figure('Name', 'E6: Net magnetization after Fourier Encoding with DFT Matrices');
% colororder = 'rkbmg';
% fs=20;
% 
% subplot(1,2,1)
% for j=1:length(FxyA(:))
% 	plot( [0, real( FxyA(j)) ], [ 0 imag(FxyA(j)) ], ...
% 		[colororder(mod(j,length(colororder))+1),'-o'],...
% 		'linewidth', 1);
% 	hold on;
% 	
% end
% 
% axis equal;
% set(gca, 'fontsize', fs);
% xlabel('Real (DFT(A)) ', 'fontsize', fs);
% ylabel('Imag (DFT(A)) ', 'fontsize', fs);
% title(' M_{net} ');
% l = sum(A(:));
% xlim([-l,l]); ylim([-l,l]);
% 
% % without DC
% subplot(1,2,2)
% for j=1:length(FxyA(:))
% 	if j~=sub2ind(size(FxyA),ceil( (Ny+0.1)/2),ceil( (Nx+0.1)/2))
% 		
% 		plot( [0, real( FxyA(j)) ], [ 0 imag(FxyA(j)) ], ...
% 			[colororder(mod(j,length(colororder))+1),'-o'],...
% 			'linewidth', 1);
% 		hold on;
% 	end
% end
% 
% axis equal;
% set(gca, 'fontsize', fs);
% xlabel('Real (Encode(A)) ', 'fontsize', fs);
% ylabel('Imag (Encode(A)) ', 'fontsize', fs);
% title(' M_{net} - no DC');


