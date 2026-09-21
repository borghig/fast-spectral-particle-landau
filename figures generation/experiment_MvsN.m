function experiment_MvsN(recompute)
% BKW error over time and computational cost (Figures 3a and 4).
% With no argument, reuse saved runs; pass true to recompute this experiment.
if nargin < 1, recompute = false; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
setup_paths();
results_dir = fullfile(root,'figures generation','results');
if ~isfolder(results_dir), mkdir(results_dir); end

%% Settings

Ms = [32, 48, 64];
Nsides = {
    [30, 40, 60], ...
    [40, 60, 80], ...
    [80, 100, 120]
};
dts = [1e-3, 1e-3, 1e-3];
testname = 'BKW';

filename = fullfile(results_dir, ...
    'experiment_MvsN.mat');

%% Load or compute
if isfile(filename)
    saved_data = load_paper_data(filename, 'run_info', 'results');

    if ~isfield(saved_data, 'run_info') || ~isfield(saved_data, 'results')
        error('The existing file does not contain run_info and results: %s', ...
            filename);
    end

    run_info = saved_data.run_info;
    results = saved_data.results;

    fprintf('Loaded %d completed runs from: %s\n', ...
        nnz(~cellfun(@isempty, {results.out})), filename);
else
    run_info = struct();

    nRuns = sum(cellfun(@numel, Nsides));
    results = struct( ...
        'method',  cell(nRuns,1), ...
        'M',       cell(nRuns,1), ...
        'Nside',   cell(nRuns,1), ...
        'N',       cell(nRuns,1), ...
        'dt',      cell(nRuns,1), ...
        'runtime', cell(nRuns,1), ...
        'prm',     cell(nRuns,1), ...
        'out',     cell(nRuns,1));
end

% Keep the metadata synchronized with the current requested experiment.
run_info.name = 'Spectral particle experiment';
run_info.testname = testname;
run_info.Ms = Ms;
run_info.Nsides = Nsides;
run_info.dts = dts;
run_info.epsi = 0;
run_info.timeint = 'RK4';
run_info.method = 'spectral_particle';
if ~isfield(run_info, 'date')
    run_info.date = string(datetime('now'));
end
run_info.last_updated = string(datetime('now'));

save(filename, 'run_info', 'results', '-v7.3');

for iM = 1:numel(Ms)
    M = Ms(iM);
    dt = dts(iM);

    for Nside = Nsides{iM}
        N = Nside^2;

        % A run is complete only if M, Nside and dt all match and its
        % output structure is nonempty.
        idx_existing = [];
        for j = 1:numel(results)
            if isempty(results(j).out) || isempty(results(j).M) || ...
                    isempty(results(j).Nside)
                continue;
            end

            saved_dt = results(j).dt;
            if isempty(saved_dt) && ~isempty(results(j).prm) && ...
                    isfield(results(j).prm, 'dt')
                saved_dt = results(j).prm.dt;
            end

            if results(j).M == M && results(j).Nside == Nside && ...
                    ~isempty(saved_dt) && saved_dt == dt
                idx_existing = j;
                break;
            end
        end

        if ~isempty(idx_existing) && ~recompute
            fprintf(['Skipping spectral particles: M = %d, N = %d^2, ' ...
                'dt = %.0e (already saved).\n'], M, Nside, dt);
            continue;
        end

        fprintf(['\nRunning spectral particles: M = %d, N = %d^2, ' ...
            'dt = %.0e\n'], M, Nside, dt);

        prm = initialization(M, N, testname);
        prm.plot = false;
        prm.epsi = 0;
        prm.G = 1;
        prm.timeint = 'RK4';
        prm.dt = dt;
        prm.dt_error = 1e-1;
        prm.times_error = 0:prm.dt_error:prm.tmax;

        tStart = tic;
        [out, prm] = method_spectral_particle(prm);
        runtime = toc(tStart);
        out = compute_error(out, prm);

        % Reuse a preallocated empty entry, or append if none is available.
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

        run_info.last_updated = string(datetime('now'));
        save(filename, 'run_info', 'results', '-v7.3');

        fprintf('Saved M = %d, N = %d^2, dt = %.0e.\n', ...
            M, Nside, dt);
    end
end

fprintf('\nFinished all experiments.\nSaved to: %s\n', filename);

%% Plot and export

filename = fullfile(results_dir, ...
    'experiment_MvsN.mat');
saved = load_paper_data(filename);
run_info = saved.run_info;
results = saved.results;

set(groot, 'DefaultLineLineWidth', 1.1);

Ms = run_info.Ms;
Nsides = run_info.Nsides;
allNsides = unique([Nsides{:}]);
colors = lines(numel(allNsides));
errfield = 'rel_L2_error';

figure('Position', [100, 200, 700, 250]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for iM = 1:numel(Ms)
    M = Ms(iM);
    nexttile;

    for Nside = Nsides{iM}
        idx = find([results.M] == M & [results.Nside] == Nside & ...
            [results.dt] == 1e-3, 1);
        iColor = find(allNsides == Nside, 1);

        semilogy(results(idx).out.times_error, ...
            results(idx).out.(errfield), ...
            'Color', colors(iColor,:), ...
            'DisplayName', sprintf('$N = %d^2$', Nside));
        hold on;
        axis padded
    end

    grid on;
    box on;
    xlabel('$t$', 'Interpreter', 'latex');
    title(sprintf('$M = %d$', M), 'Interpreter', 'latex');
    legend('Location', 'best', 'Interpreter', 'latex');

    if iM == 1
        ylabel('rel. $L^2$ error', 'Interpreter', 'latex');
        % ylim([7e-4,2e-2]);
    end
end

% commonYLim = [min(AX(2).YLim(1), AX(3).YLim(1)), ...
% set(AX(2:3), 'YLim', commonYLim);

exportgraphics(gcf, fullfile(results_dir, 'experiment_MvsN.pdf'), ...
    'ContentType', 'vector');

%% Cost, accuracy and complexity comparison

spectral_file = fullfile(results_dir, 'experiment_MvsN.mat');
blob_file = fullfile(results_dir, 'experiment1_BKW_blob.mat');
experiment_blob_only();

if ~isfile(spectral_file)
    error('Spectral particle results file not found: %s', spectral_file);
end
if ~isfile(blob_file)
    error('Blob results file not found: %s', blob_file);
end

spectral_data = load_paper_data(spectral_file, 'run_info', 'results');
blob_data = load_paper_data(blob_file, 'results');

Rs = spectral_data.results;
Rb = blob_data.results;

% Keep completed runs only.
Rs = Rs(~cellfun(@isempty, {Rs.out}));
Rb = Rb(~cellfun(@isempty, {Rb.out}));

Ms_cost = spectral_data.run_info.Ms;
Nsides_cost = spectral_data.run_info.Nsides;
dts_cost = spectral_data.run_info.dts;
errfield_cost = 'rel_L2_error';

% Gather the currently requested spectral particle configurations.
spectral_curves = repmat(struct( ...
    'M', [], 'N', [], 'runtime', [], 'error', []), ...
    numel(Ms_cost), 1);

for iM = 1:numel(Ms_cost)
    M_cost = Ms_cost(iM);

    spectral_curves(iM).M = M_cost;

    for Nside_cost = Nsides_cost{iM}
        idx = find( ...
            [Rs.M] == M_cost & ...
            [Rs.Nside] == Nside_cost & ...
            [Rs.dt] == dts_cost(iM), ...
            1, 'last');

        if isempty(idx)
            warning('Missing spectral run: M=%d, N=%d^2, dt=%.0e.', ...
                M_cost, Nside_cost, dts_cost(iM));
            continue;
        end

        err = Rs(idx).out.(errfield_cost);

        spectral_curves(iM).N(end+1) = Rs(idx).N;
        spectral_curves(iM).runtime(end+1) = Rs(idx).runtime;
        spectral_curves(iM).error(end+1) = err(end);
    end
end

% Load the requested Blob runs. These three runs use dt = 1e-2.
blob_Nsides = [40, 60, 80];
blob_dt = 1e-2;
blob_N = [];
blob_runtime = [];
blob_error = [];

for Nside_blob = blob_Nsides
    idx = find( ...
        [Rb.Nside] == Nside_blob & ...
        [Rb.dt] == blob_dt, ...
        1, 'last');

    if isempty(idx)
        warning('Missing Blob run: N=%d^2, dt=%.0e.', ...
            Nside_blob, blob_dt);
        continue;
    end

    err = Rb(idx).out.(errfield_cost);

    blob_N(end+1) = Rb(idx).N;
    blob_runtime(end+1) = Rb(idx).runtime;
    blob_error(end+1) = err(end);
end

method_colors = lines(numel(Ms_cost) + 1);
method_markers = {'o', 's', '^', 'd'};

f_cost = figure('Position', [120, 220, 780, 300]);
tl_cost = tiledlayout(1, 2, ...
    'Padding', 'compact');
% Left: runtime versus final relative L2 error
nexttile(tl_cost);

for iM = 1:numel(Ms_cost)
    if isempty(spectral_curves(iM).N)
        continue;
    end

    loglog(spectral_curves(iM).runtime, ...
        spectral_curves(iM).error, ...
        ['-' method_markers{iM}], ...
        'Color', method_colors(iM,:), ...
        'DisplayName', sprintf('$M=%d$', ...
            spectral_curves(iM).M));
    hold on;
end

if ~isempty(blob_N)
    loglog(blob_runtime, blob_error, ...
        ['-.' method_markers{end}], ...
        'Color', method_colors(end,:), ...
        'DisplayName', 'Blob');
end

grid on;
box on;
axis padded;
xlabel('Computational time (s)','Interpreter', 'latex');
ylabel('rel. $L^2$ error ($t = 8$)', 'Interpreter', 'latex');
legend('Location', 'best', 'Interpreter', 'latex');
% Right: total particle count versus runtime
nexttile(tl_cost);

all_N = [];
all_runtime = [];

for iM = 1:numel(Ms_cost)
    if isempty(spectral_curves(iM).N)
        continue;
    end

    Nplot = spectral_curves(iM).N(:);
    tplot = spectral_curves(iM).runtime(:);

    % Sort by particle count
    [Nplot, order] = sort(Nplot);
    tplot = tplot(order);

    % Compute empirical slope: runtime approximately C*N^slope
    valid_fit = Nplot > 0 & tplot > 0;

    if nnz(valid_fit) >= 2
        p = polyfit(log(Nplot(valid_fit)), ...
            log(tplot(valid_fit)), 1);

        computed_slope = p(1);

        curve_label = sprintf( ...
            '$M=%d$ (slope %.2f)', ...
            spectral_curves(iM).M, computed_slope);
    else
        curve_label = sprintf( ...
            '$M=%d$', spectral_curves(iM).M);
    end

    loglog(Nplot, tplot, ...
        ['-' method_markers{iM}], ...
        'Color', method_colors(iM,:), ...
        'DisplayName', curve_label);
    hold on;

    all_N = [all_N, Nplot.'];
    all_runtime = [all_runtime, tplot.'];
end

if ~isempty(blob_N)
    loglog(blob_N, blob_runtime, ...
        ['-.' method_markers{end}], ...
        'Color', method_colors(end,:), ...
        'DisplayName', 'Blob');

    all_N = [all_N, blob_N];
    all_runtime = [all_runtime, blob_runtime];
end

% Reference slopes with N denoting the total number of particles
valid = all_N > 0 & all_runtime > 0;

if any(valid)
    Nmin = min(all_N(valid));
    Nmax = max(all_N(valid));

    Nref = logspace(log10(Nmin), log10(Nmax), 100);
    tref = min(all_runtime(valid))/3;

    loglog(Nref, 1*tref*(Nref/Nmin), 'k--', ...
        'DisplayName', 'Slope $N$');

    loglog(Nref, 2*tref*(Nref/Nmin).^2, 'k:', ...
        'DisplayName', 'Slope $N^2$');
end

grid on;
box on;
axis padded;

xlabel('$N$ (num. of particles)', 'Interpreter', 'latex');
ylabel('Computational time (s)','Interpreter', 'latex');

legend('Location', 'best', 'Interpreter', 'latex', 'NumColumns',1);

exportgraphics(f_cost, ...
    fullfile(results_dir, 'experiment_MvsN_cost_accuracy.pdf'), ...
    'ContentType', 'vector');
end
