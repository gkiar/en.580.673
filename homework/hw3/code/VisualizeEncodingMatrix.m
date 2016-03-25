function f = VisualizeEncodingMatrix(E, Name)
% function to visualize encoding matrices 
% Herzka MRM Spring 2014

if nargin<2, Name = []; end

[Ny,Nx] = size(E);

A = ones(size(E));
Aencoded = A.*E;

f = figure;

if (Ny>32 || Nx > 32)

	subplot(3,2,1);
	imagesc(real(Aencoded));
	
	xlabel('N', 'fontsize', 20)
	ylabel('Measurement', 'fontsize', 20)
	xlim([-1 , Nx+1])
	ylim([-1 , Ny+1])
	title('Real','fontsize', 20);
	colormap(gray)
	axis image
	colorbar;
	
	subplot(3,2,2);
	imagesc(imag(Aencoded));
	xlabel('N', 'fontsize', 20)
	ylabel('Measurement', 'fontsize', 20)
	xlim([-1 , Nx+1])
	ylim([-1 , Ny+1])
	title('Imaginary','fontsize', 20);
	axis image
	colormap(gray);
	colorbar;
	
	h_s=subplot(3,2,6);
	imagesc(atan2(imag(Aencoded), real(Aencoded)));
	axis image
	colormap(gray);
	colorbar;
	set(h_s, 'position', [0.0785714 0.0785714 0.792857 0.559524], 'fontsize',20);
	xlabel('N', 'fontsize', 20)
	ylabel('Measurement', 'fontsize', 20)
	xlim([-1 , Nx+1])
	ylim([-1 , Ny+1])
	axis image
	title('Phase Angle');
	
	
else
		
	scale = 0.5;
	xc = linspace(0,2*pi,100);
	yc = scale*sin(xc);
	xc = scale*cos(xc);
	
	for x=0:Nx-1
		for y=0:Ny-1
			mArrow3( [x,y,0], [x+real(Aencoded(y+1,x+1)*scale), y+imag(Aencoded(y+1,x+1)*scale), 0],...
				'tipWidth', 0.08, 'stemWidth', 0.04);
			hold on
			plot( x+xc,y+yc, 'k:');
		end
	end
	
	axis equal
	xlabel('N', 'fontsize', 20)
	ylabel('Measurement', 'fontsize', 20)
	xlim([-1 , Nx+1])
	ylim([-1 , Ny+1])
	set(gca, 'ydir', 'reverse');
	
	if ~isempty(Name)
		title(['Encoding Matrix: ', Name], 'fontsize', 20)
	end
	set(gca, 'color', get(gcf, 'color'),...
		'ytick', 0:(Ny-1), ...
		'xtick', 0:(Nx-1), ...
		'fontsize', 20);
end
