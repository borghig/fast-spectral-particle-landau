function y = psi_2d(v1,v2,epsi)

% Returns a matrix of size v1 and v1 with values of psi evaluated at
% (v1,v2).

d = 2;
y = 1/(2*pi*epsi)^(d/2)*exp(-(v1.^2+v2.^2)/(2*epsi));
end
