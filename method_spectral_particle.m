function [out,prm] = method_spectral_particle(prm)
% Spectral particle method for the homogeneous Landau equation.
% Use initialization to set parameters; compute_error reconstructs
% densities and diagnostics from the saved particle positions.

% --- Initial condition ---
v = prm.v;  % particles location
w = prm.w;  % particles weights

% --- Stored diagnostics ---
Nt           = floor(prm.tmax / prm.dt);
out.times_error   = prm.times_error;
out.fs_error = zeros(numel(out.times_error),prm.n_error, prm.n_error);
out.vs_error = zeros(numel(out.times_error),prm.N,2);
out.method = 'Spectral particle';

kerr = 0;

fprintf('.....')

% --- Time stepping ---
for kt = 0:Nt
    t = kt*prm.dt;

    % diagnostics error cadence
    if abs(t - kerr*prm.dt_error) < prm.dt/4
        kerr = kerr + 1;

        % save data
        out.vs_error(kerr,:,:)   = v;

        fprintf('\b\b\b\b\b\b%05.2f%%', t/prm.tmax*100);

    end

    % --- Time integrator ---

    switch prm.timeint
        case 'euler'
            [dv1, dv2] = drift_spectral_particle(v, w, prm);
            v = v - prm.dt*[dv1, dv2];
        case "RK4"

            % k1
            [dv1, dv2] = drift_spectral_particle(v, w, prm); k1 = -[dv1, dv2];

            % k2
            [dv1, dv2] = drift_spectral_particle(v + 0.5*prm.dt*k1, w, prm); k2 = -[dv1, dv2];

            % k3
            [dv1, dv2] = drift_spectral_particle(v + 0.5*prm.dt*k2, w, prm); k3 = -[dv1, dv2];

            % k4
            [dv1, dv2] = drift_spectral_particle(v + prm.dt*k3, w, prm); k4 = -[dv1, dv2];

            % RK4 update
            v = v + prm.dt/6 * (k1 + 2*k2 + 2*k3 + k4);

        otherwise
            error('Unknown time integrator: %s', prm.timeint);
    end

    % projection to periodic domain
    v = mod(v + prm.vmax, 2*prm.vmax) - prm.vmax;

end

end
