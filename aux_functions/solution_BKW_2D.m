function f = solution_BKW_2D(t, v1, v2)

Kt = @(t) 1 - exp(-t/8)/2;
f_exact = @(t,v1,v2) (1/(2*pi*Kt(t))) * exp(-(v1.^2 + v2.^2)/(2*Kt(t))) ...
    .* ((2*Kt(t) - 1)/Kt(t) + (1 - Kt(t)) * (v1.^2 + v2.^2)/(2*Kt(t)^2));

f = f_exact(t, v1, v2);

if t<0
    Kt = 1;
    f = (1/(2*pi*Kt)) * exp(-(v1.^2 + v2.^2)/(2*Kt)) ...
    .* ((2*Kt - 1)/Kt + (1 - Kt) * (v1.^2 + v2.^2)/(2*Kt^2));
end

end
