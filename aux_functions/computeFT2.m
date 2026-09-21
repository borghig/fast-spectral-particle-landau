function F = computeFT2(f,x,y,M,L)
% FINUFFT type 1: scattered values -> a uniform 2-D Fourier-mode grid.
% FINUFFT uses 2*pi-periodic coordinates, hence the rescaling below.
tol = 1e-13;

xj = pi*(x(:)+L)/L;
yj = pi*(y(:)+L)/L;
opts.modeord = 1;               % Return [0,...,positive,negative] FFT order
% Store modes as F(ky,kx), consistently with fft2 and computeIFT2.
F = finufft2d1(xj,yj,f(:),-1,tol,M,M,opts).';
end
