function fI = computeIFT2(F, xq, yq, M, L)
% Evaluate a 2-D Fourier expansion at nonuniform points (FINUFFT type 2).
% F(row,col) = F(ky,kx), matching the original implementation.

tol = 1e-13;

% Map coordinates from [-L,L) to FINUFFT's 2*pi-periodic domain.
xq = pi*(xq(:)+L)/L;
yq = pi*(yq(:)+L)/L;

opts.modeord = 1;   % [0:M/2-1, -M/2:-1] ordering

% FINUFFT expects dimensions ordered as (kx,ky), hence F.'.
% The dot transpose is essential: it transposes without conjugating.
fI = finufft2d2(xq, yq, +1, tol, F.', opts)/M^2;

fI = real(fI);
end
