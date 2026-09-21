function experiment_conservation_entropy(recompute)
% BKW conservation, fourth moment, and entropy (Figure 5).
% With no argument, reuse saved runs; pass true to recompute this experiment.
if nargin < 1, recompute = false; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
setup_paths();
results_dir = fullfile(root,'figures generation','results');
if ~isfolder(results_dir), mkdir(results_dir); end

%% Settings

testname = 'BKW';
% Selected spectral particle configurations: [M, Nside, dt]
% Delete or add rows here to change the comparison.
configs = [ ...
    48,  80, 1e-2; ...
    64, 100, 1e-3; ...
];

dt_error = 1e-1;

% Blob baseline
blob_Nside = 80;
blob_dt = 1e-2;

spectral_file = fullfile(results_dir, 'experiment_MvsN.mat');
blob_file = fullfile(results_dir, 'experiment1_BKW_blob.mat');
experiment_blob_only();
%% Load or compute
if isfile(spectral_file)
    saved_data = load_paper_data(spectral_file, 'run_info', 'results');

    if ~isfield(saved_data, 'run_info') || ~isfield(saved_data, 'results')
        error('The existing file does not contain run_info and results: %s', ...
            spectral_file);
    end

    run_info = saved_data.run_info;
    results = saved_data.results;
else
    run_info = struct();
    run_info.name = 'Spectral particle experiment';
    run_info.testname = testname;
    run_info.method = 'spectral_particle';
    run_info.date = string(datetime('now'));

    results = struct( ...
        'method',  {}, ...
        'M',       {}, ...
        'Nside',   {}, ...
        'N',       {}, ...
        'dt',      {}, ...
        'runtime', {}, ...
        'prm',     {}, ...
        'out',     {});
end

required_diagnostics = {'momentum', 'energy', 'mom4', 'log_entropy'};
selected_idx = zeros(size(configs,1),1);
% Reuse complete configurations and run only the missing ones
for icfg = 1:size(configs,1)
    M = configs(icfg,1);
    Nside = configs(icfg,2);
    N = Nside^2;
    dt = configs(icfg,3);

    idx_existing = [];

    for j = 1:numel(results)
        if isempty(results(j).out) || isempty(results(j).M) || ...
                isempty(results(j).Nside) || isempty(results(j).dt)
            continue;
        end

        has_diagnostics = all(isfield(results(j).out, ...
            required_diagnostics));

        if strcmp(results(j).method, 'spectral_particle') && ...
                results(j).M == M && ...
                results(j).Nside == Nside && ...
                results(j).dt == dt && has_diagnostics
            idx_existing = j;
            break;
        end
    end

    if ~isempty(idx_existing) && ~recompute
        selected_idx(icfg) = idx_existing;

        fprintf(['Reusing M=%d, N=%d^2, dt=%.0e from ' ...
            'experiment_MvsN.mat.\n'], M, Nside, dt);
        continue;
    end

    fprintf('\nRunning M=%d, N=%d^2, dt=%.0e\n', M, Nside, dt);

    prm = initialization(M, N, testname);
    prm.plot = false;
    prm.epsi = 0;
    prm.G = 1;
    prm.timeint = 'RK4';
    prm.dt = dt;
    prm.dt_error = dt_error;
    prm.times_error = 0:prm.dt_error:prm.tmax;

    t_start = tic;
    [out, prm] = method_spectral_particle(prm);
    runtime = toc(t_start);

    out = compute_error(out, prm);

    % Reuse an empty preallocated entry, or append a new one.
    k = idx_existing;
    if isempty(k), k = find(arrayfun(@(r) isempty(r.out), results), 1); end
    if isempty(k)
        k = numel(results) + 1;
    end

    results(k).method = 'spectral_particle';
    results(k).M = M;
    results(k).Nside = Nside;
    results(k).N = N;
    results(k).dt = dt;
    results(k).runtime = runtime;
    results(k).prm = prm;
    results(k).out = out;

    selected_idx(icfg) = k;

    run_info.diagnostic_configs = configs;
    run_info.last_updated = string(datetime('now'));

    save(spectral_file, 'run_info', 'results', '-v7.3');

    fprintf('Saved M=%d, N=%d^2, dt=%.0e. Runtime: %.2f seconds.\n', ...
        M, Nside, dt, runtime);
end

run_info.diagnostic_configs = configs;
run_info.last_updated = string(datetime('now'));
save(spectral_file, 'run_info', 'results', '-v7.3');

selected_results = results(selected_idx);
% Load the Blob baseline
if ~isfile(blob_file)
    error('Blob results file not found: %s', blob_file);
end

blob_data = load_paper_data(blob_file, 'results');
blob_results = blob_data.results;
blob_results = blob_results(~cellfun(@isempty, {blob_results.out}));

idx_blob = find( ...
    [blob_results.Nside] == blob_Nside & ...
    [blob_results.dt] == blob_dt, ...
    1, 'last');

if isempty(idx_blob)
    error('Missing Blob baseline: N=%d^2, dt=%.0e.', ...
        blob_Nside, blob_dt);
end

blob_result = blob_results(idx_blob);

if ~all(isfield(blob_result.out, required_diagnostics))
    error('The saved Blob baseline does not contain all diagnostics.');
end
%% Plot and export
set(groot, 'DefaultLineLineWidth', 1.1);
set(groot, 'DefaultAxesFontSize', 9);

nConfigs = size(configs,1);
plot_colors = lines(nConfigs);

% Keep the Blob baseline visually fixed across every diagnostic tile.
matlab_default_colors = lines(7);
blob_color = matlab_default_colors(4,:);

labels = cell(nConfigs,1);
for icfg = 1:nConfigs
    dtExponent = round(log10(configs(icfg,3)));

    labels{icfg} = sprintf( ...
        '$M=%d$, $N=%d^2$, $\\Delta t=10^{%d}$', ...
        configs(icfg,1), configs(icfg,2), dtExponent);
end
blob_label = sprintf('Blob, $N=%d^2$, $\\Delta t=10^{%d}$', ...
    blob_Nside, round(log10(blob_dt)));

f_diag = figure('Position', [100, 200, 700, 450]);
tl_diag = tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

diagnostic_fields = {'momentum', 'energy', 'mom4'};
diagnostic_titles = { ...
    'Momentum error', ...
    'Energy error', ...
    'Fourth-moment error'};

for ifield = 1:numel(diagnostic_fields)
    ax = nexttile(tl_diag);
    field = diagnostic_fields{ifield};

    for icfg = 1:nConfigs
        t = selected_results(icfg).out.times_error;
        y = selected_results(icfg).out.(field);

        semilogy(ax, t, max(y, 1e-16), ...
            'Color', plot_colors(icfg,:), ...
            'DisplayName', labels{icfg});
        hold(ax, 'on');
    end

    t_blob = blob_result.out.times_error;
    y_blob = blob_result.out.(field);

    semilogy(ax, t_blob, max(y_blob, 1e-16), '-.', ...
        'Color', blob_color, ...
        'DisplayName', blob_label);

    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, '$t$', 'Interpreter', 'latex');
    title(ax, diagnostic_titles{ifield});
end

% Entropy functional H(f) = integral f log(f) dv.
ax_entropy = nexttile(tl_diag);

for icfg = 1:nConfigs
    t = selected_results(icfg).out.times_error;
    H = selected_results(icfg).out.log_entropy;

    plot(ax_entropy, t, H, ...
        'Color', plot_colors(icfg,:), ...
        'DisplayName', labels{icfg});
    hold(ax_entropy, 'on');
end

plot(ax_entropy, blob_result.out.times_error, ...
    blob_result.out.log_entropy, '-.', ...
    'Color', blob_color, ...
    'DisplayName', blob_label);

% Exact entropy decay of the two-dimensional BKW solution.
tmax_reference = max(blob_result.out.times_error);
for icfg = 1:nConfigs
    tmax_reference = max(tmax_reference, ...
        max(selected_results(icfg).out.times_error));
end

t_reference = linspace(0, tmax_reference, 500);
H_reference = entropy_BKW_2D(t_reference);

plot(ax_entropy, t_reference, H_reference, 'k:', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'exact');

grid(ax_entropy, 'on');
box(ax_entropy, 'on');
xlabel(ax_entropy, '$t$', 'Interpreter', 'latex');
%    'Interpreter', 'latex');
title(ax_entropy, 'Entropy decay');
legend(ax_entropy, 'Location', 'best', 'Interpreter', 'latex');

exportgraphics(f_diag, ...
    fullfile(results_dir, 'experiment_conservation_entropy.pdf'), ...
    'ContentType', 'vector');

end

function H = entropy_BKW_2D(t)
%ENTROPY_BKW_2D Exact H(f) = integral f*log(f) for the 2D BKW solution.

input_size = size(t);
t = t(:);

K = 1 - 0.5*exp(-t/8);
a = (2*K - 1)./K;
b = (1 - K)./K;

H = zeros(size(t));
gamma_euler = 0.5772156649015329;

% Limiting value at t = 0, where K = 1/2 and a = 0.
idx_initial = a <= 1e-12;
H(idx_initial) = -log(pi) - 1 - gamma_euler;

% Maxwellian equilibrium limit K = 1.
idx_equilibrium = b <= 1e-12 & ~idx_initial;
H(idx_equilibrium) = -log(2*pi) - 1;

% Intermediate BKW states.
idx = ~idx_initial & ~idx_equilibrium;
q = a(idx)./b(idx);

H(idx) = ...
    -log(2*pi*K(idx)) ...
    -1./K(idx) ...
    +log(a(idx)) ...
    +b(idx).*(1 + exp(q).*expint(q));

H = reshape(H, input_size);
end
