clc;
clear all;
close all;

%% --------- Geometric and mechanical properties of the beam (from paper) ----------
L = 0.137;           % beam length (m)
R = 2e-3;            % outer radius (m)
r = 1.5e-3;          % inner radius (m)
D = 2 * R;           % outer diameter (m)
t = R - r;           % wall thickness (m)
J = (pi/4)*(R^4 - r^4);    % second moment of area (m^4)
theta_m = 6450;             % material density (kg/m^3)
mu = theta_m * pi * (R^2 - r^2);  % mass per unit length (kg/m)
ms = 0.0176;                 % tip mass (kg)

%% --------- Thermal and electrical properties of SMA ---------
E_a = 52.7e9;      % Young's modulus of austenite (Pa)
E_m = 32.3e9;      % Young's modulus of martensite (Pa)
rho_a = 100e-8;    % resistivity of austenite (Ω·m)
rho_m = 90e-8;     % resistivity of martensite (Ω·m)
zeta_a = 0.009;    % damping ratio in austenite
zeta_m = 0.0122;   % damping ratio in martensite
T0 = 22;           % ambient temperature (°C)
As = 55; Af = 65;  % austenite start/finish temperatures (°C)

h = 15;                                 % convective heat-transfer coefficient (W/m²·K)
A_cs = pi * (D^2 - (D - 2*t)^2) / 4;    % beam cross-sectional area (m^2)
A_ext = pi * D * 2 * L;                 % external surface area for convection (m^2)
R_L = 2 * L / A_cs;                     % length-to-area ratio for resistance

%% --------- Currents to be evaluated ----------
I_vals = [0, 6.3, 6.5, 7];

%% --------- Frequency sweep for GSY calculation ----------
f_range = linspace(15, 35, 1000);
Omega = 2 * pi * f_range;

figure; hold on;

for I = I_vals

    % --- initial guess: fully martensitic (ξ = 1)
    xi_guess = 1;
    for iter = 1:3
        rho_eff = rho_a + xi_guess * (rho_m - rho_a);
        R_elec  = rho_eff * R_L;
        T       = T0 + (R_elec * I^2) / (h * A_ext);
        if T < As
            xi = 1;
        elseif T > Af
            xi = 0;
        else
            xi = (Af - T) / (Af - As);
        end
        xi_guess = xi;  % update for next iteration
    end

    % --- temperature-dependent material properties
    E    = E_a + xi * (E_m - E_a);
    zeta = zeta_a + xi * (zeta_m - zeta_a);

    % -------- characteristic equation solution and first natural frequency --------
    df    = 1e-4;
    f_det = 0:df:100;
    w_det = 2 * pi * f_det;
    dets  = zeros(size(w_det));
    c0    = (mu/(E * J))^0.25;

    for k = 1:numel(w_det)
        wk    = w_det(k);
        gamma = sqrt(wk) * c0;
        tau   = E * J * gamma^3;
        lam   = wk^2 * ms;
        Hk = [ 0, 1, 0, 1;
               1, 0, 1, 0;
              -sin(gamma*L), -cos(gamma*L), sinh(gamma*L), cosh(gamma*L);
              -tau*cos(gamma*L)+lam*sin(gamma*L), ...
               tau*sin(gamma*L)+lam*cos(gamma*L), ...
               tau*cosh(gamma*L)+lam*sinh(gamma*L), ...
               tau*sinh(gamma*L)+lam*cosh(gamma*L) ];
        dets(k) = det(Hk);
    end

    [~, locs] = findpeaks(-abs(dets), 'MinPeakProminence',0.1);
    f_nat     = locs(1) * df;
    omega1    = 2 * pi * f_nat;
    fprintf("I = %.1f A | T = %.2f °C | ξ = %.2f | E = %.2f GPa | f₁ = %.2f Hz\n", ...
        I, T, xi, E/1e9, f_nat);

    % ---------- compute mode shape coefficients ----------
    gamma  = (omega1^2 * mu / (E * J))^0.25;
    tau    = E * J * gamma^3;
    lambda = omega1^2 * ms;

    Lg = gamma * L;
    A_reduced = [
        0, 1, 0;
        -cos(Lg), sinh(Lg), cosh(Lg);
        tau*sin(Lg)+lambda*cos(Lg), tau*cosh(Lg)+lambda*sinh(Lg), tau*sinh(Lg)+lambda*cosh(Lg)
    ];
    b = [
        -1;
        sin(Lg);
        -(-tau*cos(Lg)+lambda*sin(Lg))
    ];
    B_rest = A_reduced \ b;
    B      = [1; B_rest];

    % ---------- define mode shape φ(x) ----------
    phi = @(x) real( ...
        B(1)*sin(gamma*x) + ...
        B(2)*cos(gamma*x) + ...
        B(3)*sinh(gamma*x) + ...
        B(4)*cosh(gamma*x));

    % ---------- derivatives at the root and integrals ----------
    dx   = 1e-6;
    dphi = @(x) (phi(x+dx) - phi(x-dx)) / (2*dx);
    d2phi = @(x) (phi(x+dx) - 2*phi(x) + phi(x-dx)) / dx^2;
    d3phi = @(x) (d2phi(x+dx) - d2phi(x-dx)) / (2*dx);
    phi3   = d3phi(0);
    phiL   = phi(L);
    int_phi = integral(phi, 0, L);
    m1      = mu * integral(@(x) phi(x).^2, 0, L) + ms * phiL^2;

    % ---------- compute GSY ----------
    numerator   = Omega.^2 .* phi3 .* (mu * int_phi + ms * phiL);
    denominator = m1 * (-Omega.^2 + 2j * zeta * omega1 * Omega + omega1^2);
    GSY = -2 * E * J * numerator ./ denominator;

    % ---------- plot ----------
    plot(f_range, abs(GSY), 'LineWidth', 2, 'DisplayName', ['I = ' num2str(I) ' A']);
end

xlabel('Frequency (Hz)');
ylabel('|G_{SY}(j\Omega)| [N/m]');
title('GSY – Mode 1 Only (Thermo-electromechanical SMA Model)');
legend show;
grid on;
xlim([18 28]);

%% Appendix: martensite fraction
function xi = xi_T(T)
    As = 55; Af = 65;
    xi = zeros(size(T));
    xi(T <= As)    = 1;
    mid = T>As & T<Af;
    xi(mid) = (Af - T(mid))/(Af - As);
end
