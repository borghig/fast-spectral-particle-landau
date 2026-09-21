function M4 = solution_BKW_2D_mom4(t)
% Exact fourth moment for the 2D BKW solution:
%
%   M4(t) = int_{R^2} |v|^4 f(t,v) dv
%
% For the BKW solution with
%   K(t) = 1 - exp(-t/8)/2,
% one obtains
%   M4(t) = 8*K(t)*(2 - K(t)) = 8 - 2*exp(-t/4).
%
% For t < 0, this matches solution_BKW_2D by returning the Maxwellian value.

K = 1 - 0.5*exp(-t/8);

M4 = 8*K.*(2 - K);

% Match convention in solution_BKW_2D
M4(t < 0) = 8;

end
