function G = filter_GaussianG(M, epsi, vmax)
% Gaussian Fourier multiplier with variance epsi on [-vmax,vmax]^2.

kk   = ([0:M/2-1, -M/2:-1])'/(2*vmax);
[KX,KY] = ndgrid(kk,kk);

G   = exp(-2*pi^2*epsi*(KX.^2+KY.^2));

end
