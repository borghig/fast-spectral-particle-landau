function out = compute_error(out,prm)
% Reconstruct density, evaluate benchmark diagnostics, and optionally animate.

% Trubnikov supplies an exact second-moment law, not the full density used
% by the other benchmarks.  Compute its diagnostic directly from stored
% particle locations.  With the adapted C=1/16 benchmark, solver time and
% physical time coincide.
if strcmpi(prm.testname,'Trubnikov')
    out = compute_trubnikov_temperature(out,prm);
    return
end

Nt = numel(out.times_error);
out.rel_L2_error = zeros(Nt,1);
out.rel_Linf_error = zeros(Nt,1);
out.momentum = zeros(Nt,1);
out.energy = zeros(Nt,1);
out.mom4 = zeros(Nt,1);
out.log_entropy = zeros(Nt,1);

for kerr = 1:Nt
    f_ref = prm.f_exact(prm.times_error(kerr), ...
        prm.v1_error,prm.v2_error);

    switch out.method
        case 'Spectral particle'
            v = squeeze(out.vs_error(kerr,:,:));
            f = compute_fKDE(prm.w,v(:,1),v(:,2), ...
                prm.v1_error,prm.v2_error,prm);
            out.fs_error(kerr,:,:) = f;
        case 'Blob'
            v = squeeze(out.vs_error(kerr,:,:));
            f = compute_KDE(prm.w,v(:,1),v(:,2), ...
                prm.v1_error,prm.v2_error,prm.epsi);
            out.fs_error(kerr,:,:) = f;
    end

    df = f-f_ref;
    out.rel_L2_error(kerr) = norm(df(:),2)/norm(f_ref(:),2);
    out.rel_Linf_error(kerr) = norm(df(:),inf)/norm(f_ref(:),inf);

    w = prm.w(:);
    mom1 = sum(w.*v(:,1));
    mom2 = sum(w.*v(:,2));
    speed2 = v(:,1).^2+v(:,2).^2;
    out.momentum(kerr) = hypot( ...
        mom1-prm.u(1),mom2-prm.u(2));
    out.energy(kerr) = abs(sum(w.*speed2)-prm.T);
    M4 = solution_BKW_2D_mom4(out.times_error(kerr));
    out.mom4(kerr) = abs(sum(w.*speed2.^2)-M4);

    dV = prm.dv_error^2;
    f_safe = max(f(:),realmin);
    out.log_entropy(kerr) = sum(f_safe.*log(f_safe))*dV;
end

if prm.plot
    [~,ix0] = min(abs(prm.v_error));
    xsec = prm.v_error(ix0);
    t = out.times_error(1);
    f = squeeze(out.fs_error(1,:,:));
    v = squeeze(out.vs_error(1,:,:));

    figure('Position',[200,200,800,600]);
    tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    title(tl,sprintf('%s method: M = %d, N = %d', ...
        out.method,prm.M,prm.N))

    ax_f = nexttile(tl,1);
    img_f = imagesc(ax_f,prm.v_error,prm.v_error,f);
    set(ax_f,'YDir','normal')
    axis(ax_f,'equal','tight')
    colorbar(ax_f)
    xlabel(ax_f,'v_1')
    ylabel(ax_f,'v_2')
    title(ax_f,sprintf('%s: density f(v,t), t = %.3g',out.method,t));

    ax_sc = nexttile(tl,2);
    sc = scatter(ax_sc,v(:,1),v(:,2));
    xlabel(ax_sc,'v_1')
    ylabel(ax_sc,'v_2')
    title(ax_sc,'Particles');
    xlim(ax_sc,[-prm.vmax,prm.vmax]);
    ylim(ax_sc,[-prm.vmax,prm.vmax]);

    ax_sec = nexttile(tl,3);
    h_sec = plot(ax_sec,prm.v_error,f(:,ix0),'LineWidth',1.5);
    hold(ax_sec,'on')
    f_exact_slice = prm.f_exact(t,prm.v1_error(:,ix0),prm.v2_error(:,ix0));
    h_sec_exact = plot(ax_sec,prm.v_error,f_exact_slice,'--','LineWidth',1.5);
    grid(ax_sec,'on')
    xlabel(ax_sec,'v_2')
    ylabel(ax_sec,sprintf('f(%.2g,v_2,t)',xsec))
    title(ax_sec,sprintf('Section at v_1 ~= 0, t = %.3g',t));
    legend(ax_sec,{out.method,'Exact BKW'},'Location','best')

    ax_diag = nexttile(tl,4);
    semilogy(ax_diag,out.times_error,out.rel_L2_error,'LineWidth',1.5)
    hold(ax_diag,'on')
    semilogy(ax_diag,out.times_error,out.momentum,'LineWidth',1.5)
    semilogy(ax_diag,out.times_error,out.energy,'LineWidth',1.5)
    grid(ax_diag,'on')
    xlabel(ax_diag,'time')
    ylabel(ax_diag,'error')
    legend(ax_diag,{'Relative L^2','Momentum','Energy'},'Location','best')
    title(ax_diag,'Density and moment errors');

    for kerr = 1:Nt
        t = out.times_error(kerr);
        f = squeeze(out.fs_error(kerr,:,:));
        v = squeeze(out.vs_error(kerr,:,:));
        img_f.CData = f;
        h_sec.YData = f(:,ix0);
        h_sec_exact.YData = prm.f_exact(t,prm.v1_error(:,ix0),prm.v2_error(:,ix0));
        sc.XData = v(:,1);
        sc.YData = v(:,2);
        title(ax_f,sprintf('Density f(v,t), t = %.3g',t));
        title(ax_sec,sprintf('Section at v_1 ~= 0, t = %.3g',t));
        drawnow limitrate
        pause(0.05);
    end
end

end

function out = compute_trubnikov_temperature(out,prm)
% Directional temperatures and error in the Trubnikov moment law.

physical_time = out.times_error(:);
out.times_error = physical_time;

Nt = numel(physical_time);
w = prm.w(:);
mass = sum(w);
w = w/mass;

out.temperature_x = zeros(Nt,1);
out.temperature_y = zeros(Nt,1);
out.temperature_mean = zeros(Nt,1);
out.temperature_difference = zeros(Nt,1);
out.temperature_anisotropy = zeros(Nt,1);
out.temperature_anisotropy_exact = ...
    exp(-physical_time/prm.tau_trubnikov);
out.temperature_error = zeros(Nt,1);

out.rel_L2_error = zeros(Nt,1);
out.rel_Linf_error = zeros(Nt,1);
out.momentum = zeros(Nt,1);
out.energy = zeros(Nt,1);
out.mom4 = nan(Nt,1);
out.log_entropy = zeros(Nt,1);

for kerr = 1:Nt
    v = squeeze(out.vs_error(kerr,:,:));

    switch out.method
        case 'Spectral particle'
            f = compute_fKDE(prm.w,v(:,1),v(:,2), ...
                prm.v1_error,prm.v2_error,prm);
        case 'Blob'
            f = compute_KDE(prm.w,v(:,1),v(:,2), ...
                prm.v1_error,prm.v2_error,prm.epsi);
    end
    out.fs_error(kerr,:,:) = f;

    f_ref = prm.f_exact(physical_time(kerr), ...
        prm.v1_error,prm.v2_error);
    df = f-f_ref;
    out.rel_L2_error(kerr) = norm(df(:),2)/norm(f_ref(:),2);
    out.rel_Linf_error(kerr) = norm(df(:),inf)/norm(f_ref(:),inf);

    ux = sum(w.*v(:,1));
    uy = sum(w.*v(:,2));
    Tx = sum(w.*(v(:,1)-ux).^2);
    Ty = sum(w.*(v(:,2)-uy).^2);
    deltaT = Tx-Ty;

    out.temperature_x(kerr) = Tx;
    out.temperature_y(kerr) = Ty;
    out.temperature_mean(kerr) = 0.5*(Tx+Ty);
    out.temperature_difference(kerr) = deltaT;
    out.temperature_anisotropy(kerr) = deltaT/prm.deltaT0;
    out.temperature_error(kerr) = abs( ...
        out.temperature_anisotropy(kerr) ...
        -out.temperature_anisotropy_exact(kerr));

    out.momentum(kerr) = hypot(ux-prm.u(1),uy-prm.u(2));
    out.energy(kerr) = abs(Tx+Ty-prm.T);

    f_safe = max(f(:),realmin);
    out.log_entropy(kerr) = ...
        sum(f_safe.*log(f_safe))*prm.dv_error^2;
end

out.particle_mass = mass;

end
