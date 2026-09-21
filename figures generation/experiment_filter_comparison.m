function experiment_filter_comparison(recompute)
% Gaussian filtering at M = 48 and M = 128 (Figure 6).
% With no argument, reuse saved runs; pass true to recompute this experiment.
if nargin < 1, recompute = false; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
setup_paths();
results_dir = fullfile(root,'figures generation','results');
if ~isfolder(results_dir), mkdir(results_dir); end

%% Settings

Ms = [48,128];
experiment_blob_only();
for M = Ms
    Nside = 80;
    N = Nside^2;
    dt = 1e-2;
    testname = 'BKW';

    epsi_factors = [0.1, 0.3, 0.64];
    n_spectral_cases = 1 + numel(epsi_factors);
    n_cases = n_spectral_cases + 1;

    labels = cell(1, n_cases);
    labels{1} = 'no filtering';
    for j = 1:numel(epsi_factors)
        labels{j + 1} = sprintf('Gaussian, \\epsilon = %.2g h^2', ...
            epsi_factors(j));
    end
    labels{end} = 'Blob method';

    filename = fullfile(results_dir, sprintf('experiment_filter_M%d_N80.mat', M));
    blob_filename = fullfile(results_dir, 'experiment1_BKW_blob.mat');

    use_cache = false;
    if isfile(filename) && ~recompute
        saved = load_paper_data(filename);
        use_cache = saved.M == M && saved.Nside == Nside && saved.dt == dt ...
            && isequal(saved.epsi_factors, epsi_factors) ...
            && numel(saved.results) == n_cases ...
            && all(arrayfun(@(r) ~isempty(r.out) && ...
                isfield(r.out,'rel_L2_error'), saved.results));
    end

    %% Load or compute
    if ~use_cache
        results = struct( ...
            'label', labels, ...
            'epsi_factor', cell(1, n_cases), ...
            'epsi', cell(1, n_cases), ...
            'runtime', cell(1, n_cases), ...
            'prm', cell(1, n_cases), ...
            'out', cell(1, n_cases));

        save(filename, 'M', 'Nside', 'N', 'dt', 'testname', ...
            'epsi_factors', 'results', '-v7.3');

        for i = 1:n_spectral_cases
            prm = initialization(M, N, testname, 5);
            prm.plot = false;
            prm.dt = dt;
            prm.timeint = 'euler';
            prm.dt_error = 10*prm.dt;
            prm.times_error = 0:prm.dt_error:prm.tmax;

            h = abs(prm.v(1,2) - prm.v(2,2));

            if i == 1
                epsi_factor = 0;
                prm.epsi = 0;
                prm.G = ones(M);
            else
                epsi_factor = epsi_factors(i - 1);
                prm.epsi = epsi_factor*h^2;
                prm.G = filter_GaussianG(M, prm.epsi, prm.vmax);
            end

            fprintf('\nRunning %s\n', labels{i});

            tStart = tic;
            [out, prm] = method_spectral_particle(prm);
            runtime = toc(tStart);
            out = compute_error(out, prm);

            results(i).epsi_factor = epsi_factor;
            results(i).epsi = prm.epsi;
            results(i).runtime = runtime;
            results(i).prm = prm;
            results(i).out = out;

            save(filename, 'M', 'Nside', 'N', 'dt', 'testname', ...
                'epsi_factors', 'results', '-v7.3');
        end

        % Load the matching Blob result instead of running method_blob again.
        if ~isfile(blob_filename)
            error('Blob results file not found: %s', blob_filename);
        end

        blob_data = load_paper_data(blob_filename, 'results');
        blob_results = blob_data.results;

        % Ignore unfinished entries if the Blob experiment was interrupted.
        is_complete = ~cellfun(@isempty, {blob_results.out});
        blob_results = blob_results(is_complete);

        idx_blob = find( ...
            strcmp({blob_results.method}, 'blob') & ...
            ([blob_results.Nside] == Nside) & ...
            ([blob_results.dt] == dt));

        if isempty(idx_blob)
            error(['No saved Blob run was found for N = %d^2 and dt = %.0e ' ...
                'in %s.'], Nside, dt, blob_filename);
        elseif numel(idx_blob) > 1
            error(['More than one saved Blob run was found for N = %d^2 and ' ...
                'dt = %.0e.'], Nside, dt);
        end

        blob_run = blob_results(idx_blob);
        i_blob = n_cases;

        fprintf('\nLoading %s: N = %d^2, dt = %.0e\n', ...
            labels{i_blob}, Nside, dt);

        results(i_blob).epsi_factor = NaN;
        results(i_blob).epsi = blob_run.prm.epsi;
        results(i_blob).runtime = blob_run.runtime;
        results(i_blob).prm = blob_run.prm;
        results(i_blob).out = blob_run.out;

        save(filename, 'M', 'Nside', 'N', 'dt', 'testname', ...
            'epsi_factors', 'results', '-v7.3');

        fprintf('\nFinished all experiments.\nSaved to: %s\n', filename);

    end

    %% Plot and export

    set(groot, 'DefaultLineLineWidth', 1.1);
    set(groot, 'DefaultAxesFontSize', 9);

    filename = fullfile(results_dir, sprintf('experiment_filter_M%d_N80.mat', M));
    saved = load_paper_data(filename);
    M = saved.M;
    Nside = saved.Nside;
    results = saved.results;

    standard_colors = [ ...
        0.0000, 0.4470, 0.7410; ...
        0.8500, 0.3250, 0.0980; ...
        0.9290, 0.6940, 0.1250; ...
        0.4940, 0.1840, 0.5560; ...
        0.4660, 0.6740, 0.1880; ...
        0.3010, 0.7450, 0.9330; ...
        0.6350, 0.0780, 0.1840];

    % Keep MATLAB's fourth standard color (purple) for Blob only.
    plot_colors = [standard_colors([1, 2, 3, 5], :); ...
        standard_colors(4, :)];
    styles = repmat({'-'}, 1, numel(results));
    styles{end} = '-.';
    errfield = 'rel_L2_error';

    figure('Position', [100, 200, 700, 230]);

    % Preserve the original 1:2 width ratio between kernels and errors.
    tl = tiledlayout(1, 5);

    % Left: one-dimensional kernel slices at v2 = 0

    ax1 = nexttile(tl, 1, [1, 2]);

    i_blob = numel(results);
    n_spectral_cases = i_blob - 1;
    L = results(1).prm.vmax;
    epsiBlob = results(i_blob).prm.epsi;

    % Use an odd number of points so the slice is evaluated exactly at v1 = 0.
    x = linspace(-1.5, 1.5, 4001);

    % Fourier modes used by the spectral particle method.
    k = [0:M/2-1, -M/2:-1];
    xi = k/(2*L);
    E = exp(2i*pi*x(:)*xi);

    % Exact one-dimensional slices at v2 = 0 of the trigonometric particle
    % kernels used in the unfiltered and three Gaussian-filtered simulations.
    spectralKernelSlices = cell(1, n_spectral_cases);

    for i = 1:n_spectral_cases
        Gused = results(i).prm.G;

        if isvector(Gused)
            if numel(Gused) ~= M
                error('The stored filter for run %d has %d rather than %d modes.', ...
                    i, numel(Gused), M);
            end

            % Interpret a vector as G(k1,k2) = g(k1)g(k2).
            g = Gused(:);
            Gslice = g*sum(g);
        elseif isequal(size(Gused), [M, M])
            Gslice = sum(Gused, 2);
        else
            error(['The stored filter for run %d has size %s; expected an ' ...
                'M-vector or an M-by-M array.'], i, mat2str(size(Gused)));
        end

        spectralKernelSlices{i} = real(E*Gslice)/(2*L)^2;
    end

    % Check the unfiltered kernel normalization at the origin.
    [~, zeroXIndex] = min(abs(x));
    actualUnfilteredPeak = spectralKernelSlices{1}(zeroXIndex);
    expectedUnfilteredPeak = M^2/(2*L)^2;

    fprintf(['Kernel-slice check: M = %d, L = %.16g, ' ...
        'unfiltered peak = %.16g, expected = %.16g.\n'], ...
        M, L, actualUnfilteredPeak, expectedUnfilteredPeak);

    if abs(actualUnfilteredPeak - expectedUnfilteredPeak) ...
            > 1e-10*expectedUnfilteredPeak
        error(['Unexpected unfiltered kernel peak: obtained %.16g, ' ...
            'expected %.16g. Check the stored Fourier multiplier.'], ...
            actualUnfilteredPeak, expectedUnfilteredPeak);
    end

    % Slice at v2 = 0 of the continuous two-dimensional Blob Gaussian.
    blobGaussianSlice = exp(-x.^2/(2*epsiBlob))/(2*pi*epsiBlob);

    kernels = [spectralKernelSlices, {blobGaussianSlice}];
    kernelLabels = {results.label};
    kernelLabels{1} = 'Unfiltered Fourier';
    kernelLabels{end} = 'Blob Gaussian';

    fprintf('Kernel-slice peaks at v1 = v2 = 0:\n');
    for i = 1:numel(kernels)
        fprintf('  %s: %.16g\n', kernelLabels{i}, kernels{i}(zeroXIndex));
    end

    for i = 1:numel(kernels)
        plot(ax1, x, kernels{i}, ...
            'LineStyle', styles{i}, ...
            'Color', plot_colors(i,:), ...
            'DisplayName', kernelLabels{i});
        hold(ax1, 'on');
    end

    grid(ax1, 'on');
    box(ax1, 'on');
    xlim(ax1, [-1.5, 1.5]);
    xlabel(ax1, '$v$', 'Interpreter', 'latex');
    title(ax1, '1D kernel slice');

    % Right: relative L2 errors

    ax2 = nexttile(tl, 3, [1, 3]);

    for i = 1:numel(results)
        semilogy(ax2, results(i).out.times_error, ...
            results(i).out.(errfield), ...
            'LineStyle', styles{i}, ...
            'Color', plot_colors(i,:), ...
            'DisplayName', results(i).label);
        hold(ax2, 'on');
    end

    grid(ax2, 'on');
    box(ax2, 'on');
    xlabel(ax2, '$t$', 'Interpreter', 'latex');
    ylabel(ax2, 'rel. $L^2$ error', 'Interpreter', 'latex');
    title(ax2, sprintf('$M=%d$, $N=%d^2$', M, Nside), ...
        'Interpreter', 'latex');
    legend(ax2, 'Location', 'best');

    figname = fullfile(results_dir, sprintf('experiment_filter_M%d_N%d.pdf', M, Nside));

    exportgraphics(gcf,figname);
end
end
