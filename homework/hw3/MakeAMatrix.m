function A = MakeAMatrix(Ny,Nx,Amp)

if nargin<3
	Amp = 20;
end
if nargin<2
	Nx = Ny;
end
if nargin<1
	disp('Too few arguments...');
	return
end


A  = complex(randi(Amp,Ny,Nx),randi(Amp,Ny,Nx));