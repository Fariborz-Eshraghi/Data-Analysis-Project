%%  SMA Twin-Beam TMD – 2-DOF analytical model (host-structure stiffness added)
% ------------------------------------------------------------------------------
%  * twin cantilevered NiTi beams (length L = 0.137 m, diameter D = 4 mm)
%  * each beam carries the tip mass m1 = 17.6 g; the masses are tied together
%  * the upper mass is connected to the supporting floor by a stiffness
%    k_struct : this removes the rigid-body mode that produced 0 Hz before
%  * temperature-dependent Young’s modulus E(T) => beam tip stiffness k_beam(T)
%  * outputs:
%       – natural frequencies f1(T), f2(T)
%       – estimated As   (start of austenite)
%       – estimated Af   (finish of austenite)
%       – currents I_as, I_af needed to reach those temperatures
%       – FRF of DOF-1 at four key temperatures
% ------------------------------------------------------------------------------

clear; close all; clc;

%% ---------------------------------------------------------------------------
% 1.  Geometry, masses and constant material data
% ---------------------------------------------------------------------------
L  = 0.137;          % beam length (m)
D  = 0.004;          % outer diameter (m)
t  = 0.0005;         % wall thickness (m)

I  = pi*(D^4-(D-2*t)^4)/64;                 % 2nd moment of area (m^4)
A_cs  = pi*(D^2-(D-2*t)^2)/4;               % cross-sectional area (m^2)
A_surf = pi*D*L;                            % lateral area for convection (m^2)

m1 = 0.0176;        % tip mass of each beam (kg)
m2 = m1;            % TMD lumped mass (kg)

% host-structure vertical stiffness (added DOF constraint)
k_struct = 1.0e4;   % N/m   ← tune if you know the floor’s own frequency
c_struct = 0.01*k_struct;   % 1 % critical damping on that DOF (optional)

%% ---------------------------------------------------------------------------
% 2.  NiTi material: modulus and transformation temperatures
% ---------------------------------------------------------------------------
E_M = 32e9;   % martensite modulus (Pa)
E_A = 39e9;   % austenite  modulus (Pa)

Ms = 30; Mf = 40;   % martensite start / finish (°C)
As = 50; Af = 70;   % austenite  start / finish (°C)

%% ---------------------------------------------------------------------------
% 3.  Electrical / convective data (for current estimate)
% ---------------------------------------------------------------------------
rho_e  = 8e-7;                 % resistivity  (Ω·m)
R_beam = rho_e*L/A_cs;         % beam resistance (Ω)
h_conv = 12;                   % natural-convection h (W/m²K)
T_amb  = 20;                   % ambient temperature (°C)

%% ---------------------------------------------------------------------------
% 4.  Modal analysis versus temperature
% ---------------------------------------------------------------------------
T_vec = 20:0.5:80;                 % analysed temperatures (°C)
nT    = numel(T_vec);
f1 = zeros(1,nT);
f2 = zeros(1,nT);

for i = 1:nT
    T = T_vec(i);

    % linear austenite fraction
    xi = max(0, min(1, (T-As)/(Af-As)));

    % temperature-dependent modulus
    E_T = E_M + xi*(E_A - E_M);

    % equivalent cantilever tip stiffness  k = 3EI/L^3
    k_beam = 3*E_T*I/L^3;

    % 2-DOF mass & stiffness matrices (host stiffness added)
    M = diag([m1, m2]);
    K = [ k_struct + k_beam ,   -k_beam ;
          -k_beam           ,    k_beam ];

    % undamped natural frequencies
    wn = sqrt(eig(K,M));                % rad/s
    f1(i) = wn(1)/(2*pi);               % mode 1  (host)
    f2(i) = wn(2)/(2*pi);               % mode 2  (beam + mass)
end

figure;
plot(T_vec, f1, 'b-', T_vec, f2, 'r-', 'LineWidth',1.4);
xlabel('Temperature  [°C]'); ylabel('Natural frequency  [Hz]');
legend('Mode 1','Mode 2','Location','northwest');
title('Natural frequencies vs. temperature');

% numerical slope of f1(T)
df1 = diff(f1)./diff(T_vec);

%% ---------------------------------------------------------------------------
% 5.  Estimate As (start) and Af (finish) temperatures
% ---------------------------------------------------------------------------
[maxSlope, idxMax] = max(df1);          % Af → peak slope
Af_est = T_vec(idxMax+1);

% As → first point where slope reaches 10 % of its maximum
idxStart = find(df1 > 0.10*maxSlope, 1, 'first');
if isempty(idxStart)
    As_est = 50;                        % fallback to paper value
else
    As_est = T_vec(idxStart);
end

fprintf('\nEstimated  As ≈ %.2f  °C\n', As_est);
fprintf('Estimated  Af ≈ %.2f  °C\n', Af_est);

%% ---------------------------------------------------------------------------
% 6.  Convert those temperatures to steady-state currents
%      T_ss = T_amb + I² R / (h A)  ⇒  I = sqrt((T-Tamb) h A / R)
% ---------------------------------------------------------------------------
I_as = sqrt( (As_est - T_amb)*h_conv*A_surf / R_beam );
I_af = sqrt( (Af_est - T_amb)*h_conv*A_surf / R_beam );

fprintf('Required current at As  ≈ %.2f  A\n', I_as);
fprintf('Required current at Af  ≈ %.2f  A\n\n', I_af);

%% ---------------------------------------------------------------------------
% 7.  FRF of DOF-1 at key temperatures
% ---------------------------------------------------------------------------
T_key  = [20, As_est, Af_est, 80];          % cold, As, Af, hot
omega  = linspace(2*pi*0.5, 2*pi*60, 1000); % 0.5 Hz → 60 Hz
colors = lines(numel(T_key));

figure; hold on;
for j = 1:numel(T_key)
    T = T_key(j);
    xi = max(0, min(1, (T-As)/(Af-As)));
    E_T = E_M + xi*(E_A - E_M);
    k_beam = 3*E_T*I/L^3;

    M = diag([m1, m2]);
    K = [ k_struct + k_beam ,  -k_beam ;
          -k_beam           ,   k_beam ];

    % add small damping on each DOF
    C = diag([c_struct, 0]) + 1e-4*K;

    H11 = zeros(size(omega));
    for k = 1:numel(omega)
        w = omega(k);
        H = (-w^2*M + 1i*w*C + K)\[1;0];
        H11(k) = abs(H(1));
    end
    plot(omega/(2*pi), H11, 'Color', colors(j,:), ...
        'DisplayName', sprintf('T = %.1f °C', T));
end
xlabel('Frequency  [Hz]'); ylabel('|H_{11}|  [m/N]');
title('FRF of DOF-1 at key temperatures');
grid on; legend show; xlim([0 60]);
