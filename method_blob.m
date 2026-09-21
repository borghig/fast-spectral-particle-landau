function [out,prm] = method_blob(prm)
% Type-II blob baseline with explicit Euler time integration.

% --- Initial condition ---
v = prm.v;  % particles location
w = prm.w;  % particles weights

% --- Stored diagnostics ---
Nt           = floor(prm.tmax / prm.dt);
out.times_error   = prm.times_error;
out.fs_error = zeros(numel(out.times_error),prm.n_error, prm.n_error);
out.vs_error = zeros(numel(out.times_error),prm.N,2);
out.method = 'Blob';

kerr = 0;

% --- Time stepping ---
for kt = 0:Nt
    t = kt*prm.dt;

    % diagnostics error cadence
    if abs(t - kerr*prm.dt_error) < prm.dt/4
        kerr = kerr + 1;

        % Save particles; compute_error reconstructs density after timing.
        out.vs_error(kerr,:,:)   = v;

        fprintf('\b\b\b\b\b\b%05.2f%%', t/prm.tmax*100);

    end

    % --- Time integrator ---
    [dv1, dv2] = drift_blob(v, w, prm);

    v = v + prm.dt*[dv1, dv2];
end

end
