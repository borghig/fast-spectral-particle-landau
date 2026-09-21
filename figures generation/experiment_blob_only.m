function experiment_blob_only(recompute)
% Generate the blob baseline used in the BKW paper figures.
% With no argument, reuse saved runs; pass true to recompute this experiment.
if nargin < 1, recompute = false; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
setup_paths();
results_dir = fullfile(root,'figures generation','results');
if ~isfolder(results_dir), mkdir(results_dir); end

%% Settings

Nsides = [40, 60, 80 ];
dts    = [1e-2, 1e-2, 1e-2];

testname = 'BKW';

% The Blob method does not depend on M. This value only defines the grid
% used by compute_error for density-based errors and diagnostics.
M_blob = 64;
% Save file
filename = fullfile(results_dir, ...
    "experiment1_" + testname + "_blob.mat");

if isfile(filename) && ~recompute
    saved = load_paper_data(filename, 'results');
    complete = true;
    for j = 1:numel(Nsides)
        complete = complete && any(arrayfun(@(r) ...
            ~isempty(r.out) && r.Nside == Nsides(j) && r.dt == dts(j), ...
            saved.results));
    end
    if complete, return; end
end

%% Compute and save
run_info = struct();
run_info.name     = 'Experiment 1: Blob method only';
run_info.testname = testname;
run_info.Nsides   = Nsides;
run_info.dts      = dts;
run_info.M_blob   = M_blob;
run_info.methods  = {'blob'};
run_info.date     = string(datetime('now'));

% Store the same fields as the spectral particle experiments.
nRuns = numel(Nsides);
results = struct( ...
    'method',  cell(nRuns,1), ...
    'M',       cell(nRuns,1), ...
    'Nside',   cell(nRuns,1), ...
    'N',       cell(nRuns,1), ...
    'dt',      cell(nRuns,1), ...
    'runtime', cell(nRuns,1), ...
    'prm',     cell(nRuns,1), ...
    'out',     cell(nRuns,1) ...
);

% Create the output file immediately, then checkpoint it after every run.
save(filename, 'run_info', 'results', '-v7.3');
% Blob runs
for iN = 1:numel(Nsides)
    Nside = Nsides(iN);
    N = Nside^2;
    dt = dts(iN);

    fprintf('\nRunning blob: N = %d^2, dt = %.0e\n', Nside, dt);

    vmax_part = 5;
    prm = initialization(M_blob, N, testname,vmax_part);
    prm.plot = false;
    prm.dt = dt;
    prm.dt_error = 10*prm.dt;
    prm.tmax        = 8;
    prm.times_error = 0:prm.dt_error:prm.tmax;

    t_start = tic;
    [out, prm] = method_blob(prm);
    runtime = toc(t_start);

    out = compute_error(out, prm);

    results(iN).method  = 'blob';
    results(iN).M       = NaN;
    results(iN).Nside   = Nside;
    results(iN).N       = N;
    results(iN).dt      = dt;
    results(iN).runtime = runtime;
    results(iN).prm     = prm;
    results(iN).out     = out;

    save(filename, 'run_info', 'results', '-v7.3');

    fprintf('Saved blob, N = %d^2. Runtime: %.2f seconds.\n', ...
        Nside, runtime);
end

fprintf('\nFinished all Blob experiments.\n');
fprintf('Saved to: %s\n', filename);
end
