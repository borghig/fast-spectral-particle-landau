function experiment_bump_on_tail(recompute)
% Bump-on-tail density and marginal evolution (Figure 8).
% With no argument, reuse saved runs; pass true to recompute this experiment.
if nargin < 1, recompute = false; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
setup_paths();
results_dir = fullfile(root,'figures generation','results');
if ~isfolder(results_dir), mkdir(results_dir); end

%% Settings

testname = 'BumpOnTail';

M = 64;
Nside = 140;
N = Nside^2;
dt = 1e-3;
tmax = 8;
plot_times = [0,1,3,8];

filename = fullfile(results_dir,'experiment_BumpOnTail.mat');

%% Load or compute
if ~isfile(filename) || recompute
    prm = initialization(M,N,testname,6);
    prm.plot = false;
    prm.epsi = 0;
    prm.G = ones(M);
    prm.timeint = 'RK4';
    prm.dt = dt;
    prm.tmax = tmax;
    prm.dt_error = 1;
    prm.times_error = 0:prm.dt_error:prm.tmax;

    % Finer grid for density reconstruction and visualization.
    prm.n_error = 201;
    prm.v_error = linspace(-prm.vmax,prm.vmax,prm.n_error);
    [prm.v1_error,prm.v2_error] = meshgrid(prm.v_error);
    prm.dv_error = prm.v_error(2)-prm.v_error(1);

    fprintf('\nBump on tail: M=%d, N=%d^2, dt=%.0e\n', ...
        M,Nside,dt);

    t_start = tic;
    [out,prm] = method_spectral_particle(prm);
    runtime = toc(t_start);

    densities = zeros(numel(plot_times),prm.n_error,prm.n_error);
    marginals = zeros(numel(plot_times),prm.n_error);

    for i = 1:numel(plot_times)
        [~,idx] = min(abs(out.times_error-plot_times(i)));
        v = squeeze(out.vs_error(idx,:,:));
        f = compute_fKDE(prm.w,v(:,1),v(:,2), ...
            prm.v1_error,prm.v2_error,prm);

        densities(i,:,:) = f;
        marginals(i,:) = sum(f,1)*prm.dv_error;
    end

    f_equilibrium = prm.f_equilibrium( ...
        prm.v1_error,prm.v2_error);
    marginal_equilibrium = sum(f_equilibrium,1)*prm.dv_error;

    save(filename,'M','Nside','N','dt','tmax','plot_times', ...
        'runtime','prm','out','densities','marginals', ...
        'f_equilibrium','marginal_equilibrium','-v7.3');

end

%% Plot and export

filename = fullfile(results_dir,'experiment_BumpOnTail.mat');
saved = load_paper_data(filename);
prm = saved.prm;
plot_times = saved.plot_times;
densities = saved.densities;
marginals = saved.marginals;
marginal_equilibrium = saved.marginal_equilibrium;

set(groot,'DefaultLineLineWidth',1.1);
set(groot,'DefaultAxesFontSize',9);

cmax = max(densities,[],'all');
ymax = max([marginals(:);marginal_equilibrium(:)]);

fig_density = figure('Position',[100,200,800,220]);
tl_density = tiledlayout(1,4, ...
    'TileSpacing','compact','Padding','compact');

axlim = 5;

for i = 1:numel(plot_times)
    ax = nexttile(tl_density,i);
    imagesc(ax,prm.v_error,prm.v_error, ...
        squeeze(densities(i,:,:)));
    set(ax,'YDir','normal');
    xlim(ax,[-axlim,axlim]);
    ylim(ax,[-axlim,axlim]);
    axis(ax,'square');
    clim(ax,[0,cmax]);

    xlabel(ax,'$v_1$','Interpreter','latex');
    if i == 1
        ylabel(ax,'$v_2$','Interpreter','latex');
    elseif i ==4
        colorbar(ax);
    end
    title(ax,sprintf('$t = %g$',plot_times(i)), ...
        'Interpreter','latex');
end

exportgraphics(fig_density, ...
    fullfile(results_dir, 'experiment_bump_on_tail_density.pdf'), ...
    'ContentType','vector');

fig_marginal = figure('Position',[100,200,800,220]);
tl_marginal = tiledlayout(1,4, ...
    'TileSpacing','compact','Padding','compact');

for i = 1:numel(plot_times)
    ax = nexttile(tl_marginal,i);
    plot(ax,prm.v_error,marginals(i,:),'HandleVisibility','off');
    hold(ax,'on');

    plot(ax,prm.v_error,marginal_equilibrium,'k:', ...
        'DisplayName','Maxwellian');

    if i == 1
        legend(ax,'Location','best');
    end

    grid(ax,'on');
    box(ax,'on');
    xlim(ax,[-axlim,axlim]);
    ylim(ax,[0,1.05*ymax]);
    xlabel(ax,'$v_1$','Interpreter','latex');
    if i == 1
        ylabel(ax,'density marginal','Interpreter','latex');
    end
    title(ax,sprintf('$t = %g$',plot_times(i)), ...
        'Interpreter','latex');
end

exportgraphics(fig_marginal, ...
    fullfile(results_dir, 'experiment_bump_on_tail_marginal.pdf'), ...
    'ContentType','vector');
end
