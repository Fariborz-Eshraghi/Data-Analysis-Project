%% SMA Twin-Beam TMD – Modal and FRF Analysis (English Comments)

% -------------------------------------------------------------------------
%  This script evaluates the natural frequencies of a twin NiTi SMA beam
%  (two cantilevered beams in parallel) versus temperature.  It also
%  estimates the start (As) and finish (Af) temperatures of the
%  austenitic transformation and converts them to Joule-heating currents.
%  A simple FRF is generated at key temperatures for verification.
%
%  Geometry and material data follow the paper’s specimen, except the beam
%  length is reduced from 0.140 m to 0.137 m (3 mm shorter).
% -------------------------------------------------------------------------

clear; close all; clc;

%% ------------------------------------------------------------------------
% 1. Geometry and material parameters
% -------------------------------------------------------------------------
L  = 0.137;    % single beam length       [m]  (137 mm)
D  = 0.004;    % outer diameter           [m]  (4 mm)
t  = 0.0005;   % wall thickness           [m]  (0.5 mm)

I  = pi*(D^4 - (D-2*t)^4)/64;                    % 2nd moment of area
A_cs  = pi*(D^2 - (D-2*t)^2)/4;                  % cross-section area
A_surf = pi*D*L;                                 % lateral surface area

m1 = 0.0176;   % tip mass of each beam     [kg]
m2 = m1;       % TMD mass (= tip mass)

% SMA (NiTi) modulus data
E_M = 32e9;    % martensite modulus        [Pa]
E_A = 39e9;    % austenite modulus         [Pa]

% Transformation temperatures (literature / paper)
Ms = 30; Mf = 40;   % martensite start / finish [°C]
As = 50; Af = 70;   % austenite  start / finish [°C]

% Electrical & convective data (for current estimate)
rho_e  = 8e-7;                % resistivity        [Ω·m]
R_beam = rho_e*L/A_cs;        % electrical R       [Ω]
h_conv = 12;                  % natural convection [W/(m²·K)]
T_amb  = 20;                  % ambient temp.      [°C]

%% ------------------------------------------------------------------------
% 2. Modal analysis versus temperature
% -------------------------------------------------------------------------
T_vec = 20:0.5:80;                      % 20 °C → 80 °C  (0.5 °C step)
nT    = numel(T_vec);

f1 = zeros(1,nT);                       % first natural frequency [Hz]
f2 = zeros(1,nT);                       % second natural freq.    [Hz]

for i = 1:nT
    T = T_vec(i);

    % Austenite fraction xi(T) – linear interpolation
    if     T <= As, xi = 0;
    elseif T >= Af, xi = 1;
    else           , xi = (T - As)/(Af - As);
    end

    % Temperature-dependent modulus
    E_T = E_M + xi*(E_A - E_M);

    % Equivalent beam tip stiffness  k = 3EI/L³
    k_beam = 3*E_T*I/L^3;

    % 2-DOF mass & stiffness matrices
    M = diag([m1, m2]);
    K = [ k_beam,   -k_beam;
         -k_beam,    k_beam];

    % Eigenvalue problem (undamped)
    [~, D_eig] = eig(K, M);
    w_n = sqrt(real(diag(D_eig)));      % natural ang. freq. [rad/s]

    f1(i) = w_n(1)/(2*pi);
    f2(i) = w_n(2)/(2*pi);
end

% Plot natural frequencies
figure;
plot(T_vec, f1, 'b-', T_vec, f2, 'r-', 'LineWidth', 1.5);
xlabel('Temperature  [°C]'); ylabel('Natural frequency  [Hz]');
legend('Mode 1','Mode 2','Location','northwest');
title('Natural frequencies vs. temperature');

% Numerical slope of f1(T)
df1 = diff(f1)./diff(T_vec);

%% ------------------------------------------------------------------------
% 3. Estimate As and Af from slope curve
% -------------------------------------------------------------------------
[maxSlope, idxMax] = max(df1);

% As:   first point where slope reaches 10 % of max
idxStart = find(df1 > 0.10*maxSlope, 1, 'first');
if isempty(idxStart)
    As_est = 50;                % fallback to paper value
else
    As_est = T_vec(idxStart);
end

% Af:   point right after max slope
Af_est = T_vec(idxMax+1);

fprintf('Estimated  As ≈ %.2f °C\n', As_est);
fprintf('Estimated  Af ≈ %.2f °C\n', Af_est);

%% ------------------------------------------------------------------------
% 4. Convert As_est & Af_est to steady-state currents (Joule heating model)
%    T_ss = T_amb + I² R / (h A)  ➜  I = sqrt( (T_ss - T_amb) h A / R )
% -------------------------------------------------------------------------
I_as = sqrt( (As_est - T_amb)*h_conv*A_surf / R_beam );
I_af = sqrt( (Af_est - T_amb)*h_conv*A_surf / R_beam );

fprintf('Estimated  I_as ≈ %.2f  A\n', I_as);
fprintf('Estimated  I_af ≈ %.2f  A\n', I_af);

%% ------------------------------------------------------------------------
% 5. FRF of mode 1 at key temperatures (for visual check)
% -------------------------------------------------------------------------
T_key   = [20, As_est, Af_est, 80];           % cold, As, Af, hot
omega   = linspace(2*pi*0.5, 2*pi*60, 1000);  % 0.5 Hz → 60 Hz
colors  = lines(numel(T_key));

figure; hold on;
for j = 1:numel(T_key)
    T = T_key(j);

    % fraction & modulus
    if     T <= As, xi = 0;
    elseif T >= Af, xi = 1;
    else           , xi = (T - As)/(Af - As);
    end
    E_T     = E_M + xi*(E_A - E_M);
    k_beam  = 3*E_T*I/L^3;

    % system matrices
    M = diag([m1, m2]);
    K = [ k_beam,   -k_beam;
         -k_beam,    k_beam];
    C = 1e-4 * K;                       % simple Rayleigh damping

    % FRF H11 (base force → displacement of mass 1)
    H11 = zeros(size(omega));
    for k = 1:numel(omega)
        s  = 1i*omega(k);
        x  = (-omega(k)^2*M + s*C + K) \ [1; 0];
        H11(k) = abs(x(1));
    end

    plot(omega/(2*pi), H11, 'Color', colors(j,:), ...
        'DisplayName', sprintf('T = %.1f °C', T));
end
xlabel('Frequency  [Hz]'); ylabel('|H_{11}|  [m/N]');
title('FRF of mode 1 at key temperatures');
legend show; grid on;
xlim([0 60]);
