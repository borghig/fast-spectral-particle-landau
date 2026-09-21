function experiment_MvsN_fixed_time(recompute)
% BKW error at t = 8 (Figure 3b).
% With no argument, reuse saved runs; pass true to recompute this experiment.
if nargin < 1, recompute = false; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
setup_paths();
results_dir = fullfile(root,'figures generation','results');
if ~isfolder(results_dir), mkdir(results_dir); end

%% Settings
t = 8;
dt = 1e-3;
n_error = 200; % error grid: n_error x n_error
Ms_left = [32, 48, 64, 128];
Nsides_left = [100, 200];
Ms_right = [48, 64];
Nsides_right = [40, 60, 80, 100, 120]; % N = Nside^2
blob_file = fullfile(results_dir, 'experiment1_BKW_blob.mat');
experiment_blob_only();
blob_Nsides = [40, 60, 80];
blob_dt = 1e-2;
rerun = recompute;
filename = fullfile(results_dir, 'experiment_MvsN_fixed_time.mat');

assert(t > 0 && dt > 0 && ...
    abs(t/dt - round(t/dt)) < 1e-10, 'Choose t as a multiple of dt.');

settings = struct('t', t, 'dt', dt, 'testname', 'BKW', ...
    'timeint', 'RK4', 'epsi', 0, 'G', 1, 'n_error', n_error);
results = struct('M', {}, 'Nside', {}, 'error', {});
if isfile(filename) && ~rerun
    saved = load_paper_data(filename, 'settings', 'results');
    if isequal(saved.settings, settings), results = saved.results; end
end

%% Load or compute: combinations needed by either panel
[ML, NL] = ndgrid(Ms_left, Nsides_left);
[MR, NR] = ndgrid(Ms_right, Nsides_right);
pairs = unique([ML(:), NL(:); MR(:), NR(:)], 'rows');
for k = 1:size(pairs, 1)
    M = pairs(k, 1); Nside = pairs(k, 2);
    if any([results.M] == M & [results.Nside] == Nside), continue; end
    fprintf('M = %d, N = %d^2, t = %g\n', M, Nside, t);
    prm = initialization(M, Nside^2, settings.testname);
    prm = set_error_grid(prm, n_error);
    prm.plot = false;
    prm.epsi = settings.epsi;
    prm.G = settings.G;
    prm.timeint = settings.timeint;
    prm.dt = dt;
    prm.tmax = t + 10*eps(t); % protect floor(tmax/dt) from roundoff
    prm.dt_error = t;
    prm.times_error = [0, t];
    [out, prm] = method_spectral_particle(prm);
    out = compute_error(out, prm);
    err = out.rel_L2_error(end);
    assert(isfinite(err) && err > 0, 'Invalid error for M=%d, N=%d^2.', M, Nside);
    results(end+1) = struct('M', M, 'Nside', Nside, 'error', err);
    save(filename, 'settings', 'results');
end

%% Plot and export at the fixed time
Eleft = zeros(size(ML)); Eright = zeros(size(MR));
for k = 1:numel(ML)
    j = find([results.M] == ML(k) & [results.Nside] == NL(k), 1);
    Eleft(k) = results(j).error;
end
for k = 1:numel(MR)
    j = find([results.M] == MR(k) & [results.Nside] == NR(k), 1);
    Eright(k) = results(j).error;
end

% Re-evaluate saved Blob particles on the same error grid (no simulations).
blob_data = load_paper_data(blob_file, 'results');
Rb = blob_data.results;
Rb = Rb(~cellfun(@isempty, {Rb.out}));
blob_N = []; blob_error = [];
for Nside = sort(blob_Nsides)
    j = find(arrayfun(@(r) r.Nside == Nside && ...
        r.prm.dt == blob_dt && strcmp(r.prm.testname, settings.testname), Rb), ...
        1, 'last');
    assert(~isempty(j), 'Missing saved Blob run: N=%d^2, dt=%g.', Nside, blob_dt);
    [delta, it] = min(abs(Rb(j).out.times_error - t));
    assert(delta < 1e-10*max(1, abs(t)), ...
        'No saved Blob error at t=%g for N=%d^2.', t, Nside);
    blob_N(end+1) = Rb(j).N;
    pb = set_error_grid(Rb(j).prm, n_error);
    v = squeeze(Rb(j).out.vs_error(it,:,:));
    x = pb.v1_error(:); y = pb.v2_error(:);
    f = zeros(size(x));
    % Evaluate in batches to limit the memory used by compute_KDE.
    for first = 1:1000:numel(x)
        q = first:min(first+999, numel(x));
        f(q) = compute_KDE(pb.w, v(:,1), v(:,2), x(q), y(q), pb.epsi);
    end
    f_exact = pb.f_exact(t, x, y);
    blob_error(end+1) = norm(f-f_exact(:), 2)/norm(f_exact(:), 2);
end

fig = figure('Position', [120, 220, 850, 330]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
loglog(Ms_left, Eleft, '-o', 'LineWidth', 1.2);
xticks(Ms_left);
xlabel('$M$ (num. of modes per dim.)', 'Interpreter', 'latex');
legend(compose('$N = %d^2$', Nsides_left), ...
    'Interpreter', 'latex', 'Location', 'best');
grid on; box on; axis padded;
ylabel('rel. $L^2$ error ($t = 8$)', 'Interpreter', 'latex');

nexttile;
h_spectral = loglog(Nsides_right.^2, Eright.', '-o', 'LineWidth', 1.2);
hold on;
default_colors = get(groot, 'factoryAxesColorOrder');
set(h_spectral(Ms_right == 48), 'Color', default_colors(2,:), 'Marker', 's');
set(h_spectral(Ms_right == 64), 'Color', default_colors(3,:), 'Marker', '^');
loglog(blob_N, blob_error, '-.d', 'LineWidth', 1.2, ...
    'Color', default_colors(4,:));
% Reference slopes in the total particle count N; adjust ref_scale to move them.
ref_scale = 40;
Nref = logspace(log10(min(Nsides_right.^2)), log10(max(Nsides_right.^2)), 100);
[~, iNref] = min(Nsides_right);
eref = ref_scale * min(Eright(:, iNref));
loglog(Nref, 2*eref*(Nref/Nref(1)).^(-1), 'k--', 'LineWidth', 1.1);
loglog(Nref, eref*(Nref/Nref(1)).^(-4), 'k:', 'LineWidth', 1.1);
xticks(Nsides_right.^2);
xticklabels(compose('%d^{2}', Nsides_right));
xlabel('$N$ (num. of particles)', 'Interpreter', 'latex');
legend([compose('$M = %d$', Ms_right(:)); {'Blob'; 'Slope $N^{-1}$'; 'Slope $N^{-4}$'}], ...
    'Interpreter', 'latex', 'Location', 'northeast');
grid on; box on; axis padded;
ylabel('rel. $L^2$ error', 'Interpreter', 'latex');
exportgraphics(fig, fullfile(results_dir, 'experiment_MvsN_fixed_time.pdf'), ...
    'ContentType', 'vector');

end

function prm = set_error_grid(prm, n_error)
prm.n_error = n_error;
prm.v_error = linspace(-prm.vmax, prm.vmax, n_error);
[prm.v1_error, prm.v2_error] = meshgrid(prm.v_error);
prm.dv_error = abs(prm.v_error(2)-prm.v_error(1));
end
