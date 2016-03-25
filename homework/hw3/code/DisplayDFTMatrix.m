% display DFT matrices with increasing N

f=figure;
for N=2:2:256

	
	%a = ones(N,1);  % the signal (representing 5 equally spaced pixels in our FOV
	a = (1:N).';  % the signal (representing 5 equally spaced pixels in our FOV
	%a = circshift(a,[1 1])
	
	W = -1i*2*pi/N;
	
	W = repmat(W,N,N);
	
	
	n1 = (0:N-1)';
	n1 = repmat(n1,1,N);
	
	n2 = 0:N-1;
	n2 = repmat(n2,N,1);
	
	Wforward = fftshift(exp(W.*n1.*n2));            % defined as in MATLAB's fft
	
	
	imagescn(cat(3,real(Wforward),imag(Wforward)), [], [], 10);
	
	set(f, 'name', ['DFT Matrix for N=', num2str(N)]);
	
	pause(0.1);

end	

f=figure;
N=256


%a = (1:N).';  % the signal (representing 5 equally spaced pixels in our FOV

W = -1i*2*pi/N;

W = repmat(W,N,N);


n1 = (0:N-1)';
n1 = repmat(n1,1,N);

n2 = 0:N-1;
n2 = repmat(n2,N,1);

Wforward = fftshift(exp(W.*n1.*n2));            % defined as in MATLAB's fft


imagescn(cat(3,real(Wforward),imag(Wforward)), [], [], 10);

set(f, 'name', ['DFT Matrix for N=', num2str(N)]);




