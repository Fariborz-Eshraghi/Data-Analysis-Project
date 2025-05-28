clc;
clear all;
close all;

%% --------- Geometric and mechanical properties of the beam (from paper) ----------
L       = 0.137;          % beam length (m)
R       = 2e-3;           % outer radius (m)
r       = 1.5e-3;         % inner radius (m)
D       = 2 * R;          % outer diameter (m)
t       = R - r;          % wall thickness (m)
J       = (pi/4)*(R^4 - r^4);    % second moment of area (m^4)
theta_m = 6450;           % material density (kg/m^3)
mu      = theta_m * pi * (R^2 - r^2);  % mass per unit length (kg/m)
ms      = 0.0176;         % tip mass (kg)

%% --------- Thermal and electrical properties of SMA ---------
E_a    = 52.7e9;      % Young's modulus of austenite (Pa)
E_m    = 32.3e9;      % Young's modulus of martensite (Pa)
rho_a  = 100e-8;      % resistivity of austenite (Ω·m)
rho_m  = 90e-8;       % resistivity of martensite (Ω·m)
zeta_a = 0.009;       % damping ratio in austenite
zeta_m = 0.0122;      % damping ratio in martensite
T0     = 22;          % ambient temperature (°C)
As     = 55; Af = 65; % austenite start/finish temperatures (°C)

h      = 15;            % convective heat-transfer coefficient (W/m²·K)
A_cs   = pi*(D^2 - (D-2*t)^2)/4;  % beam cross-sectional area (m^2)
A_ext  = pi*D*2*L;      % external surface area for convection (m^2)
R_L    = 2*L / A_cs;    % length-to-area ratio for resistance

%% --------- Currents to be evaluated ----------
I_vals = [0, 6.1, 6.2 ,6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7];

%% --------- Frequency vector for plotting ----------
f_range = linspace(15,35,1000);
Omega   = 2*pi*f_range;

%% ===== FIGURE 1: Force–Displacement FRF (Mode 1 & Mode 2) =====
figure(1); clf; hold on;
for I = I_vals
    %--- compute martensite fraction xi via temperature iteration
    xi_guess = 1;
    for iter = 1:500
        rho_eff = rho_a + xi_guess*(rho_m - rho_a);
        R_elec  = rho_eff * R_L;
        T       = T0 + (R_elec * I^2)/(h*A_ext);
        if      T < As,    xi = 1;
        elseif  T > Af,    xi = 0;
        else               xi = (Af - T)/(Af - As);
        end
        xi_guess = xi;
    end

    %--- temperature-dependent properties
    E    = E_a + xi*(E_m - E_a);
    zeta = zeta_a + xi*(zeta_m - zeta_a);

    %% --- solve characteristic equation for first two modes ---
    df    = 1e-4;
    f_det = 0:df:500;      % extended to 500 Hz
    w_det = 2*pi*f_det;
    dets  = zeros(size(w_det));
    c0    = (mu/(E*J))^0.25;

    for k=1:numel(w_det)
        wk    = w_det(k);
        gamma = sqrt(wk)*c0;
        tau   = E*J*gamma^3;
        lam   = wk^2*ms;
        Hk    = [ ...
            0,1,0,1; ...
            1,0,1,0; ...
            -sin(gamma*L), -cos(gamma*L), sinh(gamma*L), cosh(gamma*L); ...
            -tau*cos(gamma*L)+lam*sin(gamma*L), ...
             tau*sin(gamma*L)+lam*cos(gamma*L), ...
             tau*cosh(gamma*L)+lam*sinh(gamma*L), ...
             tau*sinh(gamma*L)+lam*cosh(gamma*L)];
        dets(k) = det(Hk);
    end

    %--- find all peaks, then pick the two lowest-frequency modes ---
    [~, locs_all] = findpeaks(-abs(dets));
    if numel(locs_all) < 2
        error('Only one peak found; second mode not identified. Adjust findpeaks parameters.');
    end
    freqs_all    = locs_all * df;
    freqs_sorted = sort(freqs_all);
    f_nat1 = freqs_sorted(1);
    f_nat2 = freqs_sorted(2);
    omega1 = 2*pi*f_nat1;
    omega2 = 2*pi*f_nat2;
    fprintf('I=%.1fA: f1=%.2fHz, f2=%.2fHz, T=%.2f°C, xi=%.2f\n', I,f_nat1,f_nat2,T,xi);

    %--- loop over modes to build & plot FRFs ---
    for modeIdx = 1:2
        if modeIdx==1
            omega = omega1; lineStyle = '-';
        else
            omega = omega2; lineStyle = '--';
        end

        %--- compute mode shape coefficients ---
        gamma  = (omega^2*mu/(E*J))^0.25;
        tau    = E*J*gamma^3;
        lambda = omega^2*ms;
        Lg     = gamma*L;
        A_red  = [ ...
            0,1,0; ...
            -cos(Lg), sinh(Lg), cosh(Lg); ...
            tau*sin(Lg)+lambda*cos(Lg), ...
            tau*cosh(Lg)+lambda*sinh(Lg), ...
            tau*sinh(Lg)+lambda*cosh(Lg)];
        bvec   = [-1; sin(Lg); -(-tau*cos(Lg)+lambda*sin(Lg))];
        B_rest = A_red \ bvec;
        B      = [1; B_rest];

        %--- define φ(x) and its 3rd derivative at x=0 ---
        phi = @(x) real( ...
            B(1)*sin(gamma*x) + B(2)*cos(gamma*x) + ...
            B(3)*sinh(gamma*x) + B(4)*cosh(gamma*x));
        dx    = 1e-6;
        d2phi = @(x) (phi(x+dx)-2*phi(x)+phi(x-dx))/dx^2;
        d3phi = @(x) (d2phi(x+dx)-d2phi(x-dx))/(2*dx);
        phi3  = d3phi(0);

        %--- modal mass integral ---
        phiL   = phi(L);
        intPhi = integral(phi,0,L);
        m_mod  = mu*intPhi + ms*phiL;

        %--- force–displacement FRF ---
        numerator  = Omega.^2 .* phi3 .* m_mod;
        denom_disp = m_mod .* (-Omega.^2 + 2j*zeta*omega.*Omega + omega^2);
        G_SY       = -2*E*J * (numerator ./ denom_disp);

        %--- plot |G_SY| ---
        plot(f_range, abs(G_SY), ...
             'LineStyle',lineStyle,'LineWidth',1.5, ...
             'DisplayName',sprintf('I=%.1fA, mode %d',I,modeIdx));
    end
end

xlabel('Frequency (Hz)');
ylabel('|G_{SY}(j\omega)| [N/m]');
title('Force–Displacement FRF (Mode 1 & Mode 2)');
legend('Location','best');
grid on;
xlim([15 35]);
