function [y1,y2] = gpsi_2d(v1,v2,epsi)

y1 = psi_2d(v1,v2,epsi).*(-v1/epsi);
y2 = psi_2d(v1,v2,epsi).*(-v2/epsi);
end
