function f = compute_fKDE(w,v1,v2,v10,v20,prm)
% Reconstruct the filtered spectral particle density at evaluation points.

% spectral KDE
vmax = prm.vmax;
dvM  = 2*vmax/prm.M;
Ff = computeFT2(w, v1, v2, prm.M,vmax );
Ff = prm.G .* Ff / (dvM^2);
f = real(computeIFT2(Ff, v10, v20, prm.M, vmax));

% if we needed to evaluate at a grid, we reshape
sz = size(v10);
if numel(v10)~= sz(1)
    f = reshape(f, sz(1), sz(2));
end

end
