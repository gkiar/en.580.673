%% Script to demonstrate basic Frequency Encoding and Phase Encoding as used in MRI
%
% This is a demonstration of how we could use sequential linear gradients 
% in the magnetic field in MRI to separate groups of spins located at 
% different x and y coordinates. This demonstration should elucidate the
% following points:

% Assumptions:
%  - ignore relaxation
%  - encoding takes place after excitation (with an ideal 90)
%  - slice selection has already been done, and we are dealing with a 2D
%     plane of spin
%  - consider each pixel as a delta function with amplitude A
%


%% Clean the space

clc
clear all
close all


%% General: Define Variables used in experiments

dt     = 4;              % us, sampling interval between frequency encoding points
ndt    = 5;              % number of frequency encodings to read 
Gx     = 20;             % mT/m Frequency encoding gradient strength

Gy_max = 30;             % mT/m, max phase encoding gradient strength
ty     = 500;            % us, duration of phase encoding gradient
nPE    = 5;              % number of phase encoding steps

FOV    = 20;             % cm, Field-of-view, for both x and y coordinates

% Additional constants
gammabar  = 42.58;       % MHz/T
B0     = 0;              % T, field strength in the rotating frame
									
noise  = 0;              % flag for use of noise in measurements
sig    = 1*1e-1;           % standard deviation of noise; use %1e-5, -1,-2,-3,-4,-5
                         %  relative to Amp_max

Amp_max = 20;            % maximum amplitude for A vectors derived with randi


fs = 20;                 % fontsize for figures


%% General: Define Image matrix
% Underlying image data (spatial coefficients we will try and recover)
%  different rows correspond to different y-locations. Lots of choices
%  to test your code with both even and odd sampling patterns, as well as 
%  with large image matrices. 
% Choose one!

Achoices.A1_2     = [ 10 -2 ];    

Achoices.A1_3     = [ 10 -2  6];    

Achoices.A1_4     = [ 10 -2  6  6];    

Achoices.A1_128   = randi(Amp_max,1,128);    

Achoices.A2_1     = [ 10 ;...     
	                   5 ];    
                   
Achoices.A2_2     = [ 10 -2 ;...     
	                   5  8 ];        
                   
Achoices.A2_3     = [ 10 -2  6; ...  
	                   5  8  -3 ]; 

Achoices.A2_4     = [ 10 -2   6  6; ...    
	                   5  8  -3  4 ];        

Achoices.A3_1     = [ 10 ; ...  
	                   5 ; ... 
			           4 ];
				   
Achoices.A3_3     = [ 10 -2  6; ...  
	                   5  8  -3; ... 
			           4  -9  5];

Achoices.A4_1     = [ 10 ; ...  
	                   5 ; ... 
			           4 ; ...
			           3 ];
                   
Achoices.A4_4     = [ 10 -2   6  6; ...  
	                   5  8  -3  4; ... 
			           4  -9  5  7; ...
			           3   3 -2  7];

Achoices.A4_5i    = randi(Amp_max,4,5);
				   
Achoices.A10_10   = randi(Amp_max,10,10);

Achoices.A10_11   = randi(Amp_max,10,11);

Achoices.A11_10   = randi(Amp_max,11,10);

Achoices.A1_2i    = [10+3i -2-4i ];

Achoices.A1_3i    = [10+3i -2-4i  6-3i];

Achoices.A1_4i    = [10+3i  -2-4i  6-3i 6+2i];  

Achoices.A2_2i    = [10+3i -2-4i ;...
	                  5-1i  8+5i];   

Achoices.A2_4i    = [10+3i  -2-4i  6-3i 6+2i;...  
	                  5-1i   8+5i -3+5i 4+4i];
			          
				  
Achoices.A3_3i    = [10+3i -2-4i  6-3i; ...  
	                  5-1i  8+5i  -3+5i; ... 
			          4+2i  -9-3i  5+1i];
		 
Achoices.A4_4i    = [10+3i  -2-4i  6-3i 6+2i;...  
	                  5-1i   8+5i -3+5i 4+4i;... 
			          4+2i  -9-3i  5+1i 7-3i; ...
			          3-7i   3+2i -2-3i 7+3i];
		 
Achoices.A10_10i  = complex(randi(Amp_max,10,10),randi(Amp_max,10,10));

Achoices.A10_11i  = complex(randi(Amp_max,10,11),randi(Amp_max,10,11));

Achoices.A11_10i  = complex(randi(Amp_max,11,10),randi(Amp_max,11,10));

Achoices.Custom   = MakeAMatrix(5,5,Amp_max);


%% Choose Image Matrix

A = Achoices.A4_4;		% choose one of the examples provided				

% spatial coordinate matrix
Nx = size(A,2);
Ny = size(A,1);
Npixels = Nx * Ny;     % number of coefficients total

x = linspace(0,FOV,Nx+1)  - FOV/2;                    % cm, locations of pixel edges along x 
x = (x(1:end-1)+x(2:end))/2;                          % cm, pixel locations along x (delta fcns)

y = linspace(0,FOV,Ny+1).' - FOV/2;                   % cm, locations of pixel edges along y 
y = (y(1:end-1)+y(2:end))/2;                          % cm, pixel locations along y (delta fcns)

% noise parameters
if noise==1
	mu  = 0;
	noise1 = (mu + (sig*Amp_max)*randn(nPE,ndt));
	noise2 = (mu + (sig*Amp_max)*randn(Ny,Nx));
else
	noise1 = zeros(nPE,ndt);
	noise2 = zeros(Ny,Nx);
end
	
disp('Matrix we are trying to recover');
disp(A)


%% General: Display image - original

figure('Name', 'Image matrix to be recovered');
h_im = imagesc(x,y,abs(A));
h_ax = ancestor(h_im, 'axes');
colormap(gray);
axis image
clim=get(h_ax, 'clim');
hold on;
ii=1;
for xx = 1:length(x)
	for yy = 1:length(y)

		plot(x(xx),y(yy), 'ro', ...
			'markerfacecolor', 'r', 'markersize', 10);
		text(x(xx),y(yy), ...
			['A_{', num2str(ii), '}=', num2str(A(yy,xx), '%2.2f')], ...
			'Fontsize', 40/sqrt(length(A)), 'horizontalAlignment', 'Center', ...
			'VerticalAlignment', 'Bottom', 'color', 'r');
		ii=ii+1;
	end
end

set(h_ax, 'fontsize', fs);
xlabel('x', 'fontsize', fs);
ylabel('y', 'fontsize', fs);


%% General: Convert to common units

dt__s = dt * 1e-6;
ty__s = ty * 1e-6;
B0__T = B0;
Gx__mT_m = Gx;
Gx__T_m  = Gx__mT_m / 10^3;
Gy_max__mT_m = Gy_max;
Gy_max__T_m  = Gy_max__mT_m / 10^3;
gammabar_Hz_T = gammabar * 1e6;
gamma__rad_Ts = gammabar_Hz_T * 2 * pi ;
x__m = x / 100;
y__m = y / 100;
FOVx__m = FOV / 100;
FOVy__m = FOV / 100;


%% Experiment 1: No Encoding - Start 
% This experiment is skipped as it is demonstrated in the non-FE non PE
% case (Experiment 2). 


%% Experiment 2: Frequency Encoding - Start 
disp('Experiment 2: Frequency Encoding:');fprintf('\n');


%% Experiment 2: Acquisition
% In the rotating frame, the precession frequency (delta_omega) is given by
% the gyromagnetic ratio and the delta in the B-field due to the gradient
% in the x-direction
 
dw_rad__s = gamma__rad_Ts * (B0__T + x__m .* Gx__T_m);     %  rad/sec (or 2pi*Hz)

% Calculate the amount of precession (in rad) for each x-coordinate

dphi__rad = dw_rad__s * dt__s;

disp('E2: Rotation per gradient time step (deg) at each x location:');
disp(num2str(x));
disp(num2str(dphi__rad * 180/ (2*pi) ) );
fprintf('\n');

% delta-phi represents how much phase is accrued per dt when gradient Gx is on.
% We  use this to figure out how much each ensemble of spins precesses due
% to gradients.

% For each time point (t1:t_ndt), we sense a signal that is
% composed of a complex value correpondent to the net magnetization at that
% time. At  t1, no time has elapsed, and as such all magnetization
% is aligned along the y axis.  After that, we have precession dphi per dt

x__m2 = repmat(x__m, size(A,1),1);        % create spatial matrix for x coordinates
Mtotal = zeros(1,ndt);                     % initialize net complex magnetization observer
M = zeros(size(A,1),size(A,2),ndt);      % initialize intermediate individual magnetization vectors

dw_FE_rad__s = gamma__rad_Ts * (B0__T+ Gx__T_m*x__m2);  % precession frequency due to Gx

for n = 1:ndt
	phi_FE__rad = dw_FE_rad__s * dt__s * (n-1);	
	Mtotal(n) = sum( sum( A .* exp(1i*phi_FE__rad),1),2) + noise1(1,n);
	M(:,:,n) =            A .* exp(1i*phi_FE__rad);
end


%% Experiment 1: Plot Encoded Vectors (after acquisition)

figure('Name', 'Net magnetization no Spatial Localization or Encoding');
colororder = 'rkbmg';
colororder = repmat(colororder, [1, ceil(ndt/length(colororder))+1]);
leg_string ={};
h_m=0;


n = 1;

Mtmp = M(:,:,n).';
Mtmp = Mtmp(:);

for j=1:length(Mtmp)
	plot( [0, real( Mtmp(j)) ], [ 0 imag(Mtmp(j)) ], ...
		[colororder(n),':o'],...
		'linewidth', 1);
	hold on;
end

h_m(n) = plot([0, real( Mtotal(n)) ], [ 0 imag(Mtotal(n)) ], ...
	[colororder(n),'-o'], ...
	'linewidth', 4, 'markerfacecolor', colororder(n));
leg_string = cat(1,leg_string, ['t_{', num2str(n), '}']);
	
	
axis equal;
set(gca, 'fontsize', fs);
xlabel('x', 'fontsize', fs);
ylabel('y', 'fontsize', fs);
title(' M_{net} - the Measured Magnetization');
ll = sum(A(:));
xlim([-ll,ll]); ylim([-ll,ll]);

if nPE*ndt<20
	h_leg=legend(h_m, leg_string, 'location','eastoutside');
	set(h_leg, 'fontsize', fs);
end


%% Experiment 2: Plot Encoded Vectors (after acquisition)

figure('Name', 'Net magnetization after Frequency Encoding');
colororder = 'rkbmg';
colororder = repmat(colororder, [1, ceil(ndt/length(colororder))+1]);
leg_string ={};
h_m=zeros(1,ndt);
for n=1:ndt
	
	Mtmp = M(:,:,n).';
	Mtmp = Mtmp(:);
	
	for j=1:length(Mtmp)	
		plot( [0, real( Mtmp(j)) ], [ 0 imag(Mtmp(j)) ], ...
			[colororder(n),':o'],...
			'linewidth', 1);
		hold on;
	end
	h_m(n) = plot([0, real( Mtotal(n)) ], [ 0 imag(Mtotal(n)) ], ...
		[colororder(n),'-o'], ...
		'linewidth', 4, 'markerfacecolor', colororder(n));	
	leg_string = cat(1,leg_string, ['t_{', num2str(n), '}']);
	
	
end;

axis equal;
set(gca, 'fontsize', fs);
xlabel('x', 'fontsize', fs);
ylabel('y', 'fontsize', fs);
title(' M_{net} - the Measured Magnetization');
ll = sum(abs(A(:)));
xlim([-ll,ll]); ylim([-ll,ll]);

if nPE*ndt<20
	h_leg=legend(h_m, leg_string, 'location','eastoutside');
	set(h_leg, 'fontsize', fs);
end

% Note: you should see that delta functions that are colocalized in x cannot
% be distinguished. In other words, we can try and set up a system of
% equations but it should not work, regardless of how we approach the problem.


%% Experiment 2: Recover the Image

% We want to (attempt to) reverse calculate A from the data above so we
% have to set up a system of of the form:
%
%  Mnet       = EncodingMatrix   *   Avec
% (ndt x 1)    (ndt x Npixels)   * (Npixels x 1_
%
% where A_vec correspondes to the vectorized version of A(:). Each row of
% the encoding matrix corresponds to the phase values we imparted with the
% gradient above.
%
% Once we have the above setup, we can invert (or pseudoinvert) the
% EncodingMatrix to compute Avec
%
% Avec = pinv(EncodingMatrix) * Mnet

EncodingMatrix =  exp( -1i * ((0:ndt-1).'*dt__s) * dw_FE_rad__s(:).' );
I_EncodingMatrix = pinv(EncodingMatrix);

Avec = I_EncodingMatrix * Mtotal.';
Avec = reshape(Avec, size(A));

disp('E2: A'' as obtained by pinv:');
disp((Avec));

% What does this result mean? It means that we have an ill conditioned
% inversion problem! As can be seen by Figure 1, the vectors that are
% co-localized in x cannot be separated. By taking the pseudoinverse we find
% a solution that yields the average (A1+A3)/2, and(A2+A4)/2 in each
% tentative A. Note that it doesn't matter how large ndt gets! (Try it)

% One way to explore this idea is to look at the row-reduced (low-echelon)
% form of the encoding matrix:

RowReduced_EncodingMatrix = rref(EncodingMatrix);
Condition_EncodingMatrix  = cond(EncodingMatrix);
Rank_EncodingMatrix = rank (EncodingMatrix);
disp('E2: Row-reduced version of the Encoding Matrix');
disp(RowReduced_EncodingMatrix);
disp('E2: Condition number');
disp(Condition_EncodingMatrix);
disp('E2: Rank');
disp(Rank_EncodingMatrix);

% We should see that the shape of the matrix makes it impossible to invert!


%% Experiment 3: Frequency and Phase Encoding -Start 
disp('Experiment 3: Phase Encoding');fprintf('\n');

% We begin by applying a gradient in the y-direction, of one dt in duration.
% We vary the amplitude in fractions of a maxGy gradient amplitude(nPE=number of
% phase encodes), resulting in different phase accrual during that time period.


%% Experiment 3: Define and convert some intermediate values 

x__m2 = repmat(x__m, size(A,1),1);            % create spatial matrix for x gradient

Gy_max__mT_m = Gy_max;                        % max phase encoding gradient
y__m2  = repmat(y__m, 1, size(A,2));          % create spatial matrix for y gradient
PE     = linspace(-1,1,nPE).';                % column vector of fractional values for each PE step

Mtotal = zeros(nPE, ndt);                     % net complex magnetization observer
M = zeros(size(A,1),size(A,2),nPE,ndt);       % intermediate individual magnetization vectors

dw_FE_rad__s = gamma__rad_Ts * (B0__T+ x__m2 * Gx__T_m);
dw_PE_rad__s = gamma__rad_Ts * (B0__T+ y__m2 * Gy_max__T_m);

% Calculate the precesion angle for the maximum PE gradient for each
% y-coordinate

for m = 0:(nPE-1)
	% Calculate the phase encoded A (effectively A that has precessed 
	% about z) since the gradient occurs before frequency encoding, the
	% action of the gradient can be described as a single precession
	phi_rad__PE = dw_PE_rad__s * ty__s * PE(m+1);
	A_PE_preFE = exp(1i*phi_rad__PE).*A;
	
	% Now repeat the frequency encoding experiment for each PE. Note the
	% use of A_PE instead of A (a pre-phase warped version of A)
	
	for n = 0:(ndt-1)
		phi_FE__rad = dw_FE_rad__s * dt__s * n;
		Mtotal(m+1,n+1) = sum( sum( A_PE_preFE .* exp(1i*phi_FE__rad),1),2) + noise1(m+1,n+1);
		M(:,:,m+1,n+1) =            A_PE_preFE .* exp(1i*phi_FE__rad);
	end
	

end


%% Experiment 3: Plot Encoded Vectors (after acquisition)
figure('Name', 'Net magnetization after Frequency and Phase Encoding');
colororder = 'rkbmg';
colororder = repmat(colororder, [1, ceil(ndt*nPE/length(colororder))+1]);
symbolorder = 'ox^vs';
symbolorder = repmat(symbolorder, [1, ceil(ndt*nPE/length(symbolorder))+1]);
leg_string = cell(nPE,ndt);
h_mn = zeros(nPE, ndt);


for m = 1:nPE
	for n=1:ndt
	
		Mtmp = M(:,:,m,n).';
		Mtmp = Mtmp(:);
		
		for j=1:length(Mtmp)
			plot( [0, real( Mtmp(j)) ], [ 0 imag(Mtmp(j)) ], [colororder(n),':', symbolorder(m)]);
			hold on;
		end
		h_mn(m,n) = plot([0, real( Mtotal(m,n)) ], [ 0 imag(Mtotal(m,n)) ], [colororder(n),'-', symbolorder(m)], ...
			'linewidth', 2, 'markerfacecolor', colororder(n));
		leg_string{m,n} =  ['t_{', num2str(n), '}', '-nPE=', num2str(m)];
	end

end;

axis equal;
set(gca, 'fontsize', fs);
xlabel('x', 'fontsize', fs);
ylabel('y', 'fontsize', fs);
ll = sum(abs(A(:)));
xlim([-ll,ll]); ylim([-ll,ll]);
if nPE*ndt<20	
	legend(h_mn(:),leg_string{:}, 'location','eastoutside');
	set(h_leg, 'fontsize', fs);
end


%% Experiment 3: Recover the Image

% Setup the same system of equations:
%
%  Mnet             = EncodingMatrix          * Avec
% ( (ndt*nPE) x 1)    ((ndt*nPE) x Npixels)   * (Npixels x 1)
%
% where A_vec correspondes to the vectorized version of A(:). Each row of
% the encoding matrix corresponds to the phase values we imparted with the
% gradient above. Note that we have to expand the FE and PE ecnoding
% matrices to create a single encoding matrix, with one row per
% "combintation" of FE and PE. Hence, the EM should have size:
%   ndt*nPE x length(Avec)


%
% Avec = pinv(EncodingMatrix) * Mnet

EncodingMatrixFE =  exp( 1i * ((0:ndt-1).' * dt__s) * dw_FE_rad__s(:).' );
EncodingMatrixPE =  exp( 1i * (PE * ty__s) * dw_PE_rad__s(:).' );

EncodingMatrixPE = repmat( EncodingMatrixPE, ndt,1);

B = shiftdim(EncodingMatrixFE,-1);
C = repmat(B, [ nPE 1 1]);
EncodingMatrixFE = reshape(C,nPE*ndt,length(A(:)));

I_EncodingMatrix = pinv(EncodingMatrixFE .* EncodingMatrixPE);

RowReduced_EncodingMatrix = rref(EncodingMatrixFE .* EncodingMatrixPE);
Condition_EncodingMatrix  = cond(EncodingMatrixFE .* EncodingMatrixPE);
Rank_EncodingMatrix = rank (EncodingMatrixFE .* EncodingMatrixPE);
disp('E3: Row-reduced version of the Encoding Matrix');
disp(RowReduced_EncodingMatrix);
disp('E3: Condition number');
disp(Condition_EncodingMatrix);
disp('E3: Rank');
disp(Rank_EncodingMatrix);

Avec = I_EncodingMatrix * Mtotal(:);
Avec = reshape(Avec, size(A));

disp('E3: A as obtained by pinv (real):');
disp(real(Avec));
disp('E3: A as obtained by pinv (imag):');
disp(imag(Avec));


disp('E3: Original (real):');
disp(real(A));
disp('E3: Original (imag):');
disp(imag(A));


disp('E3: Error norm (A-Avec)');
disp(norm(A-Avec));


%% Experiment 3: Display image - recovered

figure('Name', 'Recovered Image');
h_im = imagesc(x,y,abs(Avec));
h_ax = ancestor(h_im, 'axes');
colormap(gray);
axis image

hold on;
ii=1;
for xx = 1:length(x)
	for yy = 1:length(y)

		plot(x(xx),y(yy), 'ro', ...
			'markerfacecolor', 'r', 'markersize', 10);
		text(x(xx),y(yy), ...
			['A_{', num2str(ii), '}=', num2str(real(Avec(yy,xx)), '%2.2f'),' + ', num2str(imag(Avec(yy,xx)), '%2.2f'),'i'    ], ...
			'Fontsize', 40/sqrt(length(A)), 'horizontalAlignment', 'Center', ...
			'VerticalAlignment', 'Bottom', 'color', 'r');
		ii=ii+1;
	end
end

set(h_ax, 'fontsize', fs, 'clim', clim);
xlabel('x', 'fontsize', fs);
ylabel('y', 'fontsize', fs);


%% Experiment 3: Display image - difference

diff_vec = A-Avec;
diff_vec = diff_vec(:);
err = abs(max(diff_vec));
err_ord = ceil(log10(err));

f=figure('Name', ['Difference Image (max error  x 10^', num2str(err_ord),')']);
h_im = imagesc(x,y,abs(A-Avec));
h_ax = ancestor(h_im, 'axes');
colormap(gray);
axis image

hold on;
ii=1;
for xx = 1:length(x)
	for yy = 1:length(y)

		plot(x(xx),y(yy), 'ro', ...
			'markerfacecolor', 'r', 'markersize', 10);
		text(x(xx),y(yy), ...
			['A_{', num2str(ii), '}=', num2str( abs((A(yy,xx) - Avec(yy,xx)) / 10^err_ord), 2)], ...
			'Fontsize', 40/sqrt(length(A)), 'horizontalAlignment', 'Center', ...
			'VerticalAlignment', 'Bottom', 'color', 'r');
		ii=ii+1;
	end
end

set(h_ax, 'fontsize', fs, 'clim', clim);
xlabel('x', 'fontsize', fs);
ylabel('y', 'fontsize', fs);


%% Experiment 4: Encoding with the Fourier Matrix
% We determine the phase accrual that each spin is to experience using 
% to sample the fourier coefficients
disp('Experiment 4: Fourier Encoding');fprintf('\n');

 % First: create spatial matrix for x & y
[x__m2,y__m2] = meshgrid(x__m,y__m);               
dx__m = FOVx__m/Nx; % spatial resolution, x
dy__m = FOVy__m/Ny; % spatial resolution, y

kx_max = 1/(2*dx__m); % max extent of k-space, x
ky_max = 1/(2*dy__m);  % max extent of k-space, y

dkx__1_m = 1/FOVx__m; % k-space resolution, x
dky__1_m = 1/FOVy__m; % k-space resolution, x

% Next: calculate k-space coordinates that need to be sampled (note that diff k
% should yield the same dkx as calculated before). kx and ky are vectors of
% that represent the k-coordinates of points to be sampled.

% Note: Adjust pixel coordinates to account for the fact that we are using 
% delta fcns to represent each  pixels. Thus, the FOV is one pixel width 
% smaller. Shift pixel coordinates by half a pixel width (if even # samples)
% This correction gets smaller with larger number of pixels

% shift so it's bottom heavy. We want our encoding matrix to look lke we
% expect it to look... i.e. 
if ~mod(Nx,2)  % even # pixels - sample DC component
	kx = -kx_max + (0:Nx-1)*dkx__1_m;
	x__m_shifted = x__m - dx__m/2;
else
	kx = -kx_max + (0:Nx-1)*dkx__1_m + dkx__1_m/2;
	x__m_shifted = x__m;
end
if Nx==1, kx = 0; end

if ~mod(Ny,2)  % even # pixels - sample DC component
	ky = -ky_max + (0:Ny-1)*dky__1_m;
	y__m_shifted = y__m - dy__m/2;
else
	ky = -ky_max + (0:Ny-1)*dky__1_m + dky__1_m/2;
	y__m_shifted = y__m;
end
if Ny==1, ky=0; end

% Construct encoding matrices directly from k-space coordinates
E_kx = zeros(Nx,Nx);
E_ky = zeros(Ny,Ny);

for kk = 1:Nx
    E_kx(kk,:) = exp(-2i*pi*kx(kk)*x__m_shifted);
end
for kk = 1:Ny
    E_ky(kk,:) = exp(-2i*pi*ky(kk)*y__m_shifted);
end

disp('E4: x coordinates: (m)'); disp(x__m)
disp('E4: y coordinates: (m)'); disp(y__m')
fprintf('\n');
disp('E4: kx coordinates: (cycles/m)'); disp(kx)
disp('E4: ky coordinates: (cycles/m)'); disp(ky)
fprintf('\n');


%% Experiment 4: Visualize encoding matrix using given code

if Nx<32
	VisualizeEncodingMatrix(E_kx, 'E4: E_{kx}');
end
if Ny<32
	VisualizeEncodingMatrix(E_ky, 'E4: E_{ky}');
end


%% Experiment 4: Display image space sampling and  the DFT matrix that has been chosen

DisplayImageAndKspace(FOVx__m, FOVy__m, x__m_shifted, y__m_shifted, dx__m, dy__m, kx, ky, kx_max, ky_max)


%% Experiment 4: Calculate Gx and Gy, and can convert to phase accruals, and 
% encode magnetization as we have done before. However, for this to be
% Fourier encoding, we will pre-wind the gradients in the readout
% direction.

% We have to determine the amplitude of the gradients given their 
% respective dts. We know FOVx = 1/dkx = 1/(gammabar*Gx*dt)
% Therefore, Gx = 1 / (gammabar * FOVx * dt) <1 / (Hz/T * m * s) = T/m >
% Note that dt__s is chosen by the user with other criteria in mind (SNR 
% and readout duration) and the gradient amplitude is calculated based on 
% that choice.

% First: calculate the amplitudes of the gradients given the assigned durations 
% of dt__s and ty__s

Gx_FE__T_m  = 1 / (gammabar * FOVx__m * dt__s); % positive gradient amplitude
Gy_dPE__T_m = 1 / (gammabar * FOVy__m * ty__s); % max positive gradient amplitude

disp('E5: Calculated FE delta in gradient amplitude:')
disp([num2str(Gx_FE__T_m * 1000) , ' mT/m']);

disp('E5: Calculated PE delta in gradient amplitude:')
disp([num2str(Gy_dPE__T_m * 1000) , ' mT/m']);

[x__m2_shifted, y__m2_shifted] = meshgrid(x__m_shifted, y__m_shifted);
dw_FE__rad_s = gamma__rad_Ts*(x__m2_shifted*Gx_FE__T_m); % dphi per dt, <note shift1!>
dphi_FE_prewind__rad = -dw_FE__rad_s*dt__s;
if ~mod(Nx,2) % even Nx 
% 	x__m2_shifted = FILL;
%	dphi_FE_prewind__rad = FILL;   % phase accrued by the echo time (kx=0)  
else
%	dw_FE__rad_s = FILL;           % dphi per dt
%	dphi_FE_prewind__rad = FILL;   % phase accrued by the echo tiem (kx=0)
end

dw_PE_max__rad_s = gamma__rad_Ts*(y__m2_shifted*Gy_dPE__T_m); % dphi per ty, with max Gy on <note shift!>

if ~mod(Ny,2) % even Ny
%	y__m2_shifted = FILL;
else
%	dw_PE_max__rad_s = FILL; % dphi per ty, with max Gy on
end

Mtotal = zeros(Ny, Nx);                       % net complex magnetization to outside observer
M = zeros(size(A,1),size(A,2),Ny,Nx);         % intermediate individual magnetization vectors

if ~mod(Ny,2) % column vector of fractional values for each PE step
	PEmult = (0:(Ny-1)) - Ny/2;
else
	PEmult = (0:(Ny-1)) - (Ny-1)/2;
end
if Ny==1, PEmult = 0; end

disp('E5: Calculated PE Multiplier steps:')
disp(PEmult);


%% Experiment 4: Simulate Acquisition (Encoding)
% Calculate the precesion angle for the maximum PE gradient for each
% y-coordinate

Ex = zeros(Nx,Nx);                      % Encoding matrix, x dimension
Ey = zeros(Ny,Ny);                      % Encoding matrix, y dimension

%Ex_prewind = FILL;   % vector that represents readout prewinder phase

for m = 0:(Ny-1) % for each PE
	% Phase encode. Since the PE occurs before FE, the
	% action of the gradient can be described as a single precession
	% about z by an angle determined both Gy_max * PEmult
	
	%phi_PE__rad = FILL;
	%A_PE_preFE = FILL;

	%Ey(:,m+1) = FILL;
	
	% Prewind the frequency encode
	%A_PE_preFE = FILL; % fcn of prev A_PE_preFE

	% Now repeat the frequency encoding experiment for each PE. Note the
	% use of A_PE instead of A (a pre-phase warped version of A)
	
	for n = 0:(Nx-1)
		%phi_FE__rad = FILL;
		
		%Ex(n+1,:)  = FILL;  % all rows the same
		
		Mtotal(m+1,n+1) = sum( sum( exp(-1i*phi_FE__rad) .* A_PE_preFE ,1),2) + noise2(m+1,n+1);
		M(:,:,m+1,n+1) =            exp(-1i*phi_FE__rad) .* A_PE_preFE;
    end
end

Ey = Ey.'; % tranpose for display 


%% Experiment 4: Visualize Encoding Matrices

if Nx<32
	VisualizeEncodingMatrix(Ex, 'E4: Ex-Encode');
	VisualizeEncodingMatrix(Ex_prewind, 'E4: Ex-Prewind');
	VisualizeEncodingMatrix(Ex .* repmat(Ex_prewind, [Nx, 1]), 'E4: Ex .* Ex-Prewind');
end
if Ny<32
	VisualizeEncodingMatrix(Ey, 'E4: Ey-Encode');
end


%% Experiment 5: Plot the Fourier Encoded Vectors - derived matrix

figure('Name', 'Net magnetization after Fourier Frequency and Phase Encoding');
colororder = 'rkbmg';
symbolorder = 'ox^vs';
leg_string = cell(Ny,Nx);
h_mn = zeros(Ny, Nx);


for m = 1:Ny
	for n=1:Nx
	
		Mtmp = M(:,:,m,n).';
		Mtmp = Mtmp(:);
		
		for j=1:length(Mtmp)
			plot( [0, real( Mtmp(j)) ], [ 0 imag(Mtmp(j)) ], ...
				[colororder(mod(n,length(colororder))+1),':', symbolorder(mod(m,length(symbolorder))+1)]);
			hold on;
		end
		h_mn(m,n) = plot([0, real( Mtotal(m,n)) ], [ 0 imag(Mtotal(m,n)) ], ...
			[colororder(mod(n,length(colororder))+1),'-', symbolorder(mod(m,length(symbolorder))+1)], ...
			'linewidth', 2, 'markerfacecolor', colororder(mod(n,length(colororder))+1));
		leg_string{m,n} =  ['t_{', num2str(n), '}', '-nPE=', num2str(m)];
	
	end

end;

axis equal;
set(gca, 'fontsize', fs);
xlabel('x', 'fontsize', fs);
ylabel('y', 'fontsize', fs);
ll = sum(abs(A(:)));
xlim([-ll,ll]); ylim([-ll,ll]);
if Ny*Nx<20	
	h_leg = legend(h_mn(:),leg_string{:}, 'location','eastoutside');
	set(h_leg, 'fontsize', fs);
end


%% Experiment 5: Recover the Image using Ideal Decoding Matrix

% Setup the same system of equations:
%
%  Mnet     = EncodingMatrix_y * EncodingMatrix_x * A
%  Ny x Nx    (Ny x Ny)          (Nx x Nx)        * (Ny x Nx)
%
% Each row of the encoding matrix corresponds to the phase values we imparted with the
% gradient above. Hence, we can now use the "decode" matrix or the ifft
% function to reconstruct our image! Note the equation above ignores
% transposes.
%       Arec = fftshift(ifft2(ifftshift(Mnet)))
% Look at CreateEncodingMatrix for examples of how to do this.

% Calculate IDFT matrices for Nx and Ny; 

%Wx = FILL;
%Wx = FILL; 

%n1 = FILL; 
%n1 = FILL; 

%n2 = FILL; 
%n2 = FILL; 

%Wx_decode  = FILL;   % scaling defined as in MATLAB's fft

%Wy = FILL; 
%Wy = FILL; 

%n1 = FILL; 
%n1 = FILL; 

%n2 = FILL; 
%n2 = FILL; 

%Wy_decode = FILL;   % sacling as defined as in MATLAB's fft

%Arec_DFT   = FILL; 

Arec_ifft2 = fftshift(ifft2(ifftshift(Mtotal)));

RowReduced_EncodingMatrix = rref(Wx_decode);
Condition_EncodingMatrix  = cond(Wx_decode);
Rank_EncodingMatrix = rank (Wx_decode);

disp('E5: X - Row-reduced version of the Encoding Matrix');
disp(RowReduced_EncodingMatrix);
disp('E5: X - Condition number');
disp(Condition_EncodingMatrix);
disp('E5: X - Rank');
disp(Rank_EncodingMatrix);

RowReduced_EncodingMatrix = rref(Wy_decode);
Condition_EncodingMatrix  = cond(Wy_decode);
Rank_EncodingMatrix = rank (Wy_decode);

disp('E5: Y - Row-reduced version of the Encoding Matrix');
disp(RowReduced_EncodingMatrix);
disp('E5: Y - Condition number');
disp(Condition_EncodingMatrix);
disp('E5: Y - Rank');
disp(Rank_EncodingMatrix);

disp('E5: A as obtained by DFT matrix:');
disp(num2str(Arec_DFT,'%1.1f  '));
disp('E5: A as obtained by ifft2:');
disp(num2str(Arec_ifft2,'%1.1f  '));
disp('E5: A- Original:');
disp(num2str(A,'%1.1f  ')); 
disp('E5: Error norm (A-Arec_DFT)');
disp(num2str(norm(A-Arec_DFT)));


%% Experiment 5: For comparison, Plot the ideal encoding matrices
% for image recovery with ifft2, these need to match exactly with Ex and Ey

if Nx<32 && Ny<32
	PlotIdealEncodingMatrices(Nx,Ny);
elseif Nx<32
	PlotIdealEncodingMatrices(Nx,1);
elseif Ny<32
	PlotIdealEncodingMatrices(1,Ny);
end



