function [U1part,U2part] = drift_spectral_particle(v, w, prm)
% Reconstruct density and Landau flux, then evaluate J/f at particles.

M = prm.M;
vmax_pad = prm.vmax + prm.padd;
dvM_pad = 2*vmax_pad/M;

Ff = computeFT2(w,v(:,1),v(:,2),M,vmax_pad);
Ff = prm.G.*Ff/(dvM_pad^2);
fM = real(ifft2(Ff));

[~,~,FJ1,FJ2] = Landau_Jflux(vmax_pad, fM, prm);

J1 = real(computeIFT2(FJ1,v(:,1),v(:,2),M,vmax_pad));
J2 = real(computeIFT2(FJ2,v(:,1),v(:,2),M,vmax_pad));

f = real(computeIFT2(Ff,v(:,1),v(:,2),M,vmax_pad));

f = max(f(:),prm.save_division);
U1part = J1(:)./f;
U2part = J2(:)./f;
end
