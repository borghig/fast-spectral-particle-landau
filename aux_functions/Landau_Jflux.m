function [J1, J2, FJ1, FJ2] = Landau_Jflux(vmax, f, prm)
% Maxwell-molecule Landau flux on the M-by-M Fourier grid, scaled by 1/16.

M  = prm.M;        % output resolution

if ~isequal(size(f), [M, M])
    f = reshape(f, M, M);
end

% Use the preloaded  modes
F   = prm.F;
I11 = prm.I11;
I12 = prm.I12;
I21 = prm.I21;
I22 = prm.I22;

[F, I11, I22, I12, I21] = enforce_symmetry(F, I11, I22, I12, I21);

k = [0:M/2-1, -M/2:-1];
[kk1, kk2] = meshgrid(k, k);

% Fourier coefficients on the uniform grid
Ff = fft2(f);
fp = real(ifft2(Ff));

CFL = 16*(vmax/pi);

% First flux component
J1a = ifft2(Ff .* (F - I11)) .* ifft2(Ff .* kk1);
J1b = -ifft2(Ff .* I21) .* ifft2(Ff .* kk2);
J1c = -ifft2(Ff .* (F .* kk1 - I11 .* kk1 - I12 .* kk2)) .* fp;

J1p = real(1i * (J1a + J1b + J1c)) / CFL;

% Second flux component
J2a = ifft2(Ff .* (F - I22)) .* ifft2(Ff .* kk2);
J2b = -ifft2(Ff .* I12) .* ifft2(Ff .* kk1);
J2c = -ifft2(Ff .* (F .* kk2 - I21 .* kk1 - I22 .* kk2)) .* fp;

J2p = real(1i * (J2a + J2b + J2c)) / CFL;

% Flux coefficients and grid values
FJ1 = fft2(J1p);
FJ2 = fft2(J2p);

J1 = real(ifft2(FJ1));
J2 = real(ifft2(FJ2));

end

function [F, I11, I22, I12, I21] = enforce_symmetry(F, I11, I22, I12, I21)
% Enforce reality, coordinate parity, and tensor symmetry of the kernel modes.

    % Basic checks
    M = size(F,1);
    if size(F,2) ~= M
        error('F must be square.');
    end
    if any(size(I11) ~= [M M]) || any(size(I22) ~= [M M]) || ...
       any(size(I12) ~= [M M]) || any(size(I21) ~= [M M])
        error('All input arrays must be MxM and the same size as F.');
    end
    if mod(M,2) ~= 0
        error('M must be even (Fourier grid [0:M/2-1,-M/2:-1]).');
    end

    %----------------------------------------------------------------------
    % 1) Project everything to real (strip tiny imaginary parts)
    %----------------------------------------------------------------------
    F   = real(F);
    I11 = real(I11);
    I22 = real(I22);
    I12 = real(I12);
    I21 = real(I21);

    %----------------------------------------------------------------------
    % 2) Enforce symmetry/antisymmetry in m1 (row index)
    %----------------------------------------------------------------------
    for m = 1:(M/2 - 1)
        ipos = m + 1;          % row index for +m
        ineg = M - m + 1;      % row index for -m  (since m<0 => i = M+m+1)

        % Even in m1: F, I11, I22
        row_pos = F(ipos, :);
        row_neg = F(ineg, :);
        meanrow = 0.5 * (row_pos + row_neg);
        F(ipos, :) = meanrow;
        F(ineg, :) = meanrow;

        row_pos = I11(ipos, :);
        row_neg = I11(ineg, :);
        meanrow = 0.5 * (row_pos + row_neg);
        I11(ipos, :) = meanrow;
        I11(ineg, :) = meanrow;

        row_pos = I22(ipos, :);
        row_neg = I22(ineg, :);
        meanrow = 0.5 * (row_pos + row_neg);
        I22(ipos, :) = meanrow;
        I22(ineg, :) = meanrow;

        % Odd in m1: I12, I21
        row_pos = I12(ipos, :);
        row_neg = I12(ineg, :);
        baserow = 0.5 * (row_pos - row_neg);
        I12(ipos, :) =  baserow;
        I12(ineg, :) = -baserow;

        row_pos = I21(ipos, :);
        row_neg = I21(ineg, :);
        baserow = 0.5 * (row_pos - row_neg);
        I21(ipos, :) =  baserow;
        I21(ineg, :) = -baserow;
    end

    % Enforce axes values for I12,I21 to be 0 (m1 = 0 and m1 = -M/2 rows)
    I12(1, :)       = 0;   % m1 = 0
    I12(M/2 + 1, :) = 0;   % m1 = -M/2 (Nyquist)
    I21(1, :)       = 0;
    I21(M/2 + 1, :) = 0;

    %----------------------------------------------------------------------
    % 3) Enforce symmetry/antisymmetry in m2 (column index)
    %----------------------------------------------------------------------
    for m = 1:(M/2 - 1)
        jpos = m + 1;          % column index for +m
        jneg = M - m + 1;      % column index for -m

        % Even in m2: F, I11, I22
        col_pos = F(:, jpos);
        col_neg = F(:, jneg);
        meancol = 0.5 * (col_pos + col_neg);
        F(:, jpos) = meancol;
        F(:, jneg) = meancol;

        col_pos = I11(:, jpos);
        col_neg = I11(:, jneg);
        meancol = 0.5 * (col_pos + col_neg);
        I11(:, jpos) = meancol;
        I11(:, jneg) = meancol;

        col_pos = I22(:, jpos);
        col_neg = I22(:, jneg);
        meancol = 0.5 * (col_pos + col_neg);
        I22(:, jpos) = meancol;
        I22(:, jneg) = meancol;

        % Odd in m2: I12, I21
        col_pos = I12(:, jpos);
        col_neg = I12(:, jneg);
        basecol = 0.5 * (col_pos - col_neg);
        I12(:, jpos) =  basecol;
        I12(:, jneg) = -basecol;

        col_pos = I21(:, jpos);
        col_neg = I21(:, jneg);
        basecol = 0.5 * (col_pos - col_neg);
        I21(:, jpos) =  basecol;
        I21(:, jneg) = -basecol;
    end

    % Enforce axes values for I12,I21 to be 0 (m2 = 0 and m2 = -M/2 columns)
    I12(:, 1)       = 0;   % m2 = 0
    I12(:, M/2 + 1) = 0;   % m2 = -M/2
    I21(:, 1)       = 0;
    I21(:, M/2 + 1) = 0;

    %----------------------------------------------------------------------
    % 4) Enforce symmetry under swap (m1 <-> m2)
    %
    %   F(m1,m2)   = F(m2,m1)
    %   I12(m1,m2) = I12(m2,m1)
    %   I22(m1,m2) = I11(m2,m1)
    %----------------------------------------------------------------------
    % scalar kernel
    F = 0.5 * (F + F.');

    % off-diagonal tensor component
    I12 = 0.5 * (I12 + I12.');
    I21 = I12;

    % diagonal tensor components: make them transpose of each other
    tmp = 0.5 * (I11 + I22.');
    I11 = tmp;
    I22 = tmp.';

    %----------------------------------------------------------------------
    % 5) Final clean-up: enforce reality and I21 = I12
    %----------------------------------------------------------------------
    F   = real(F);
    I11 = real(I11);
    I22 = real(I22);
    I12 = real(I12);
    I21 = I12;
end
