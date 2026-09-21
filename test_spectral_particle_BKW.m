% Spectral particle demonstration against the exact BKW solution.

%% Settings
M = 64;      % Fourier modes per axis: 32, 48, 64, or 128
N = 100^2;    % Total particle count; must be a perfect square
tmax = 8;    % Final time;

%% Solve and visualize
addpath(fileparts(mfilename('fullpath')));
setup_paths();
prm = initialization(M, N, 'BKW');
validateattributes(tmax, {'numeric'}, {'scalar','real','finite','positive'});
n_intervals = round(tmax/prm.dt_error);

% Include the final diagnostic and protect the step count against roundoff.
prm.tmax = tmax + 10*eps(tmax);
prm.times_error = (0:n_intervals)*prm.dt_error;

[out, prm] = method_spectral_particle(prm);
out = compute_error(out, prm);
fprintf('BKW final relative L2 error: %.3e\n', out.rel_L2_error(end));
