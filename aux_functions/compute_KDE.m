function f = compute_KDE(w,v1,v2,v10,v20,epsi)
% Reconstruct the Gaussian blob density at evaluation points (v10,v20).

Z1  = v10(:) - v1';
Z2  = v20(:) - v2';

PHI = psi_2d(Z1,Z2,epsi);
f = PHI*w;

% if we needed to evaluate at a grid, we reshape
sz = size(v10);
if numel(v10)~= sz(1)
    f = reshape(f, sz(1), sz(2));
end

end
