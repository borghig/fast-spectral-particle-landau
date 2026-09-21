function prm = initialization(M,N,testname,vmax_part)
% Collect simulation parameters and a few handy derived quantities.

% Defaults give the BKW example with Gaussian filtering and Euler stepping.
if nargin < 1, M = 64; end
if nargin < 2, N = 80^2; end
if nargin < 3, testname = 'BKW'; end
if nargin < 4, vmax_part = 5; end
validateattributes(N, {'numeric'}, {'scalar','integer','positive'});
if N < 4 || round(sqrt(N))^2 ~= N
    error('N must be a perfect square with at least two particles per axis.');
end
validateattributes(M, {'numeric'}, {'scalar','integer','positive'});

% --- Domain and Time
prm.d     = 2;     % dimension

prm.vmax  = 10;
prm.dt    = 1e-2;
prm.tmax  = 8;
prm.timeint = 'euler'; % or 'RK4'

prm.plot = true;

% --- Particles and density
prm.N     = N;  % particle number
prm.vmax_part = vmax_part;

% --- Fourier initialization
prm.M     = M;    % Fourier modes (per axis)
root = fileparts(fileparts(mfilename('fullpath')));
mode_file = fullfile(root, 'modes', sprintf('mymodes_%d.mat', prm.M));
if ~isfile(mode_file)
    error('Available mode counts are M = 32, 48, 64, 128.');
end
load(mode_file, ...
    'F','I11','I22','I21','I12')
prm.gamma = 0;
fourier_scale = (prm.vmax/pi)^(prm.gamma + 2 + prm.d);
prm.F   = F  * fourier_scale;
prm.I11 = I11* fourier_scale;
prm.I12 = I12* fourier_scale;
prm.I21 = I21* fourier_scale;
prm.I22 = I22* fourier_scale;
prm.padd = 0;

prm.dvM  = 2*prm.vmax/prm.M;   % Fourier grid spacing (physical)
prm.vM  = -prm.vmax + prm.dvM*(0:prm.M-1)';
[v1_M, v2_M] = meshgrid(prm.vM);
prm.v1_M = v1_M;
prm.v2_M = v2_M;

% ---- initial data at t = 0 ----

switch testname
    case 'BKW'
        [prm.v,prm.w,prm.f_M,prm.u,prm.T,t0] = initial_data_BKW(prm);
        dv = abs(prm.v(1,2) - prm.v(2,2));
        prm.epsi = 0; %0.64 * dv^2;
        prm.G = 1;
        prm.testname = 'BKW';
        prm.f_exact = @(t,v1,v2)solution_BKW_2D(t0+t, v1,v2);
    case {'Trubnikov','TRUBNIKOV','trubnikov'}
        prm.gamma = 0;
        prm.rho = 1;
        prm.C_collision = 1/16;
        prm.Tx0 = 0.75;
        prm.Ty0 = 0.5;
        prm.Teq = 0.5*(prm.Tx0+prm.Ty0);
        prm.deltaT0 = prm.Tx0-prm.Ty0;
        prm.tau_trubnikov = 1/(8*prm.C_collision*prm.rho);

        [prm.v,prm.w,prm.f_M,prm.u,prm.T,t0] = ...
            initial_data_Trubnikov(prm);
        dv = abs(prm.v(1,2) - prm.v(2,2));
        prm.epsi = 0.64 * dv^2;
        prm.G = 1;
        prm.testname = 'Trubnikov';
        prm.f_exact = @(t,v1,v2) solution_Trubnikov_2D( ...
            t0+t,v1,v2,prm.Tx0,prm.Ty0,prm.tau_trubnikov);
        prm.deltaT_exact = @(t) ...
            prm.deltaT0*exp(-t/prm.tau_trubnikov);
        prm.Tx_exact = @(t) ...
            prm.Teq + 0.5*prm.deltaT_exact(t);
        prm.Ty_exact = @(t) ...
            prm.Teq - 0.5*prm.deltaT_exact(t);

        prm.tmax = 8;
    case {'BumpOnTail','BOT','bump-on-tail','Bump-on-tail'}
        prm.gamma = 0;
        prm.rho = 1;
        prm.bump_mass = 0.1;
        prm.bump_center = [3,0];
        prm.T_core = 1;
        prm.T_bump = 0.1;

        [prm.v,prm.w,prm.f_M,prm.u,prm.T] = ...
            initial_data_BumpOnTail(prm);
        dv = abs(prm.v(1,2) - prm.v(2,2));
        prm.epsi = 0.64 * dv^2;
        prm.testname = 'BumpOnTail';
        prm.T_equilibrium = ...
            (prm.T-sum(prm.u.^2))/prm.d;
        prm.f_initial = @(v1,v2) bump_on_tail_2D( ...
            v1,v2,prm.bump_mass,prm.bump_center, ...
            prm.T_core,prm.T_bump);
        prm.f_equilibrium = @(v1,v2) maxwellian_2D( ...
            v1,v2,prm.rho,prm.u,prm.T_equilibrium);
        prm.tmax = 8;
    otherwise
        error('Unknown test name: %s', testname);
end

% --- Derived quantities ---
prm.dvM  = 2*prm.vmax/prm.M;
prm.G = filter_GaussianG(prm.M, prm.epsi, prm.vmax);
prm.save_division = 1e-13;

% ---- Error/Diagnostic computation ----
prm.n_error = 120;
prm.v_error = linspace(-prm.vmax,prm.vmax,prm.n_error);
[v1_error,v2_error] = meshgrid(prm.v_error);
prm.v1_error = v1_error;
prm.v2_error = v2_error;
prm.dv_error = abs(prm.v_error(2)-prm.v_error(1));

prm.dt_error = 1e-1;
prm.times_error = 0:prm.dt_error:prm.tmax;

end

%% Initial distributions

function [v_part,w_part,f_M,u,T,t0] = initial_data_BKW(prm)

t0 = 0;
v = linspace(-prm.vmax_part,prm.vmax_part,sqrt(prm.N));
[v1,v2] = meshgrid(v);
v_part = [v1(:),v2(:)];
w_part = solution_BKW_2D(t0,v_part(:,1),v_part(:,2));
w_part = w_part/sum(w_part);
f_M = solution_BKW_2D(t0,prm.v1_M,prm.v2_M);
u = [0 0];
T = 2;

end

function [v_part,w_part,f_M,u,T,t0] = initial_data_Trubnikov(prm)

t0 = 0;

v = linspace(-prm.vmax_part,prm.vmax_part,sqrt(prm.N));
[v1,v2] = meshgrid(v);
v_part = [v1(:),v2(:)];

w_part = solution_Trubnikov_2D(t0,v_part(:,1),v_part(:,2), ...
    prm.Tx0,prm.Ty0,prm.tau_trubnikov);
w_part = w_part/sum(w_part);

f_M = solution_Trubnikov_2D(t0,prm.v1_M,prm.v2_M, ...
    prm.Tx0,prm.Ty0,prm.tau_trubnikov);

u = [0 0];
T = prm.Tx0+prm.Ty0;

end

function f = solution_Trubnikov_2D(t,v1,v2,Tx0,Ty0,tau)
% Gaussian with the exact Trubnikov directional-temperature evolution.
% It supplies a density in the same format as solution_BKW_2D for plots;
% the exact benchmark quantity remains the temperature difference.

Teq = 0.5*(Tx0+Ty0);
deltaT = (Tx0-Ty0)*exp(-t/tau);
Tx = Teq+0.5*deltaT;
Ty = Teq-0.5*deltaT;

f = exp(-v1.^2./(2*Tx)-v2.^2./(2*Ty)) ...
    ./(2*pi*sqrt(Tx*Ty));

end

function [v_part,w_part,f_M,u,T] = initial_data_BumpOnTail(prm)

v = linspace(-prm.vmax_part,prm.vmax_part,sqrt(prm.N));
[v1,v2] = meshgrid(v);
v_part = [v1(:),v2(:)];

w_part = bump_on_tail_2D(v_part(:,1),v_part(:,2), ...
    prm.bump_mass,prm.bump_center,prm.T_core,prm.T_bump);
w_part = w_part/sum(w_part);

f_M = bump_on_tail_2D(prm.v1_M,prm.v2_M, ...
    prm.bump_mass,prm.bump_center,prm.T_core,prm.T_bump);

% Exact momentum and energy of the Gaussian mixture.  These determine the
% theoretical Maxwellian equilibrium independently of the particle grid.
u = prm.bump_mass*prm.bump_center;
T = (1-prm.bump_mass)*prm.d*prm.T_core ...
    +prm.bump_mass*(sum(prm.bump_center.^2) ...
    +prm.d*prm.T_bump);

end

function f = bump_on_tail_2D(v1,v2,alpha,b,Tcore,Tbump)

f = (1-alpha)*maxwellian_2D(v1,v2,1,[0,0],Tcore) ...
    +alpha*maxwellian_2D(v1,v2,1,b,Tbump);

end

function f = maxwellian_2D(v1,v2,rho,u,T)

f = rho/(2*pi*T)*exp( ...
    -((v1-u(1)).^2+(v2-u(2)).^2)/(2*T));

end
