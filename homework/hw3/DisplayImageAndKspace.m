function DisplayImageAndKspace(FOVx__m, FOVy__m, x__m, y__m, dx__m, dy__m, kx, ky, kx_max, ky_max)

% Inputs: 

%   FOVx__m,FOVy__m - scalar FOVs in x0 anb y-directions, meters;
%   x__m,y__m       - vector of x- and y-positions of each pixel (assumes grid)
%   dx__m,dy__m     - scalars of image space resolution in x and y
%   kx,ky           - vectors of kx positions to be sampled, 1/m
%   kx_max, ky_max  - scalars marking k-space boundaries (assumes symmetry)


[x__m2,y__m2] = meshgrid(x__m,y__m); 

Nx = length(x__m);
Ny = length(y__m);

figure('Name', 'Image and K-space sampling comparison');

subplot(1,2,1);

plot( [ -FOVx__m -FOVx__m  FOVx__m  FOVx__m -FOVx__m ]*0.5,... 
	  [ -FOVy__m  FOVy__m  FOVy__m -FOVy__m -FOVy__m ]*0.5,...
      'r-');
hold on;		
plot(x__m2, y__m2, 'ro', 'markerfacecolor', 'r');

for n=1:Nx
	plot(   [ -FOVx__m -FOVx__m]*0.5 + n*dx__m, ...
	        [ -FOVy__m  FOVy__m]*0.5, ...
			 'r:');
end

for n=1:Ny
	plot(   [ -FOVx__m  FOVx__m ]*0.5, ...
	        [ -FOVy__m -FOVy__m ]*0.5 + n*dy__m, ...
			 'r:');
end

title('Image space; pixel coordinates (m)');
xlabel(['x (m); Nx = ', num2str(Nx)]);
ylabel(['y (m); Ny = ', num2str(Ny)]);
axis equal
xlim ([ -FOVx__m, FOVx__m]*0.5*1.1);
ylim ([ -FOVy__m, FOVy__m]*0.5*1.1);
set(gca,'ydir', 'reverse');

subplot(1,2,2)
plot( [ -kx_max -kx_max kx_max  kx_max -kx_max], ...
	  [ -ky_max  ky_max ky_max -ky_max -ky_max],...
			 'k-');
hold on
[kx_2,ky_2] = meshgrid(kx,ky);                
plot(kx_2, ky_2, 'ko', 'markerfacecolor', 'k');
	title('kspace; sampling coordinates (cycles/mm)');
xlabel('k_{x} (cycles/m)');
ylabel('k_{y} (cycles/m)');
axis equal
xlim ([ -kx_max, kx_max]*1.1);
ylim ([ -ky_max, ky_max]*1.1);