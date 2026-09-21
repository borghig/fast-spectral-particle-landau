function [dv1, dv2] = drift_blob(v, w, prm)
% Right-hand side of the blob ODE: d/dt [vv1; vv2] = RHS(vv1,vv2).
% Gaussian regularization and Maxwell-molecule collision tensor; C = 16.

vv1 = v(:,1);
vv2 = v(:,2);

% kernels (weights w are fixed)
[F1, F2] = tildeFNeps(w, vv1, vv2, prm.epsi);
gamma = 0;
[A11, A22, A12, A21] = gen_A(vv1, vv2, gamma);

% antisymmetric differences
DF1 = F1 - F1.';
DF2 = F2 - F2.';

% contracted tensors
AF1 = A11 .* DF1 + A12 .* DF2;
AF2 = A21 .* DF1 + A22 .* DF2;

% RHS (unscaled)
dv1 = -(AF1 * w);
dv2 = -(AF2 * w);

% Collision scaling: the operator is multiplied by 1/16.
if isfield(prm,'C'), C = prm.C; else, C = 16; end
dv1 = dv1 / C;
dv2 = dv2 / C;
end

%% Gaussian interaction and collision tensor

function [F1,F2] = tildeFNeps(w,v1,v2,epsi)
% type II, grid-free dynamics

feps = compute_KDE(w,v1,v2,v1,v2,epsi);
invfeps = 1./feps + 1./(feps');

% second, compute gradient
Z1 = v1 - v1';
Z2 = v2 - v2';
[gpsiZ1, gpsiZ2] = gpsi_2d(Z1,Z2,epsi);

F1 = (gpsiZ1.*invfeps)*w;
F2 = (gpsiZ2.*invfeps)*w;

end

function [A11,A22,A12,A21]= gen_A(v1,v2, gamma)

Z1 = v1 - v1';
Z2 = v2 - v2';
Znorm = sqrt(Z1.^2 + Z2.^2);

% A matrix
A11 = Znorm.^gamma.*(Znorm.^2 - Z1.^2);
A22 = Znorm.^gamma.*(Znorm.^2 - Z2.^2);
A12 = Znorm.^gamma.*(-Z1.*Z2);
A21 = A12;

end
