% ┌─────────────────────────────────────────────────────────────────────────┐ %
% │                           MPC Helicopter Main                           │ %
% └─────────────────────────────────────────────────────────────────────────┘ %
clear all
clc


% Helicopter model parameters ──────────────────────────────────────────────────
bx = 2;
by = 2;
bz = 18;
bpsi = 110;
kx = -0.5;
ky = -0.5;
kpsi = -5;
ki = 2;
parameters = [bx; by; bz; bpsi; kx; ky; kpsi; ki];
g = 9.81;
Ts = 0.1;
x_constraints = [
    -10, 10;
    -10, 10;
    -10, 10;
    -30, 30;
    -30, 30;
    -20, 20;
     0, 2*pi;
    -250, 250;
];
u_constraints = [
    -10, 10;
    -10, 10;
    -10, 10;
    0, 10;
];
model = Helicopter(parameters, Ts, x_constraints, u_constraints);


% Guide points for reference trajectory ────────────────────────────────────────
nf = 4; % Number of flat-outputs
Tend = 10;

N_intervals = 24;
N_guide = N_intervals + 1; % Number of guide points (last == first)
T_guide = linspace(0, Tend, N_guide);
Z_guide = zeros(N_guide, nf);

theta = 0;
delta = 2*pi/N_intervals;
radius = 5;
for i = 1:N_guide
    Z_guide(i, :) = [radius*cos(theta), radius*sin(theta), 0, 0.5*pi + theta];
    theta = theta + delta;
end 


% Reference trajectory filling ─────────────────────────────────────────────────
N_basis = 2; % Number of basis functions
order = 0;   % Maximum differentiation order

% Basis functions
syms x
basis = sym('x', [1 N_basis]);
for k = 1:N_basis
    basis(k) = x^(k-1);
end

% Small m(t) matrix function
function m = m(t, basis, order)
    syms x
    m = zeros(order+1, length(basis));
    for i = 0:order
        for j = 1:length(basis)
            derivative = diff(basis(j), x, i);
            m(i+1, j) = subs(derivative, x, t);
        end
    end
end

% Trajectory filling
Z_ref = [];
x_ref = [];
u_ref = [];
% alpha_tensor = zeros(nf, N_basis, length(T_guide) - 1);
for i = 1:length(T_guide) - 1

    % Time interval
    t0 = T_guide(i);
    t1 = T_guide(i+1);

    % Matrix M
    m0 = m(t0, basis, order);
    m1 = m(t1, basis, order);
    M0 = kron(eye(nf), m0);
    M1 = kron(eye(nf), m1);
    M = [M0; M1];

    % Guide points in the interval
    z_bar = [Z_guide(i, 1:nf)'; Z_guide(i+1, 1:nf)'];

    % Solve the system and reshape
    alpha = M\z_bar;
    alpha = reshape(alpha, [], nf)';
    % alpha_tensor(:, :, i) = alpha;

    % Reference function z(t) and its derivatives
    z = alpha*basis.';
    dz = diff(z, x);
    ddz = diff(dz, x);

    % Generate missing points
    N_filling = ceil((t1 - t0)/Ts);
    T_filling = linspace(t0, t1, N_filling);
    for j = 1:length(T_filling) - 1

        % Evaluate z(t) and derivatives
        z_t = double(subs(z, x, T_filling(j)))';
        z_t(nf) = wrapTo2Pi(z_t(nf)); % wrap the angle state
        dz_t = double(subs(dz, x, T_filling(j)))';
        ddz_t = double(subs(ddz, x, T_filling(j)))';

        % Compute remaining states and inputs
        dxb = cos(z_t(4))*dz_t(1) + sin(z_t(4))*dz_t(2);
        dyb = -sin(z_t(4))*dz_t(1) + cos(z_t(4))*dz_t(2);
        dzb = dz_t(3);
        dpsi = dz_t(4);

        x_t = [z_t(1), z_t(2), z_t(3), dxb, dyb, dzb, z_t(4), dpsi]; 

        ux = (cos(z_t(4))*(ddz_t(1) - kx*dz_t(1)) + sin(z_t(4))*(ddz_t(2) - kx*dz_t(2)))/bx;
        uy = (cos(z_t(4))*(ddz_t(2) - ky*dz_t(2)) + sin(z_t(4))*(-ddz_t(1) + ky*dz_t(1)))/by;
        uz = (ddz_t(3) - g)/bz;
        upsi = (ddz_t(4) - kpsi*dz_t(4))/bpsi;

        u_t = [ux, uy, uz, upsi];

        % Store the reference flat-output, states and inputs
        Z_ref = [Z_ref; z_t];
        x_ref = [x_ref; x_t];
        u_ref = [u_ref; u_t];        
    end
end


% Plot Reference / Guide ───────────────────────────────────────────────────────
figure(1);
arrow_length = 0.03;

% Guide points
guide_points = scatter(Z_guide(:, 1), Z_guide(:, 2), 15, 'filled', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', '#808080');
hold on;
% for i = 1:length(Z_guide)
%     x_arrow = arrow_length * cos(Z_guide(i, 4));
%     y_arrow = arrow_length * sin(Z_guide(i, 4));
%     quiver(Z_guide(i, 1), Z_guide(i, 2), x_arrow, y_arrow, 'AutoScale', 'off', 'Color', '#808080');
% end
% hold on;

% Reference points
reference_points = scatter(Z_ref(:, 1), Z_ref(:, 2), 5, 'filled', 'MarkerFaceColor', '#808080');
hold on;
for i = 1:length(Z_ref)
    x_arrow = arrow_length * cos(Z_ref(i, 4));
    y_arrow = arrow_length * sin(Z_ref(i, 4));
    quiver(Z_ref(i, 1), Z_ref(i, 2), x_arrow, y_arrow, 'AutoScale', 'off', 'Color', '#808080');
end
hold on;

% Labels
legend([guide_points, reference_points],{'Guide Points', 'Reference trajectory'}, 'Location', 'northwest');
xlabel('x1'); ylabel('x2');
grid on;
axis equal;
hold on;


% MPC ──────────────────────────────────────────────────────────────────────────

% MPC parameters
x0 = zeros(model.n,1);
N  = 5;                             % prediction horizon
% Q  = eye(model.n);                  % state cost
% R  = eye(model.m);                  % input cost
Q = 10*diag([50, 50, 5, 10, 3, 3, 1, 2]);
R = diag([2, 2, 2, 2]);
preview = 1;                        % MPC preview flag
formulation = 0;                    % MPC formulation flag
debug = 0;                          % MPC debug flag
Tend = 9;
t = 0:Ts:Tend;
Nsteps = length(t) - (N+1);         % number of MPC optimization steps

% Optimization 
mpc = MPC(model, x0, Tend, N, Q, R, x_ref, u_ref, preview, formulation, debug);
[x_mpc, u] = mpc.optimize();

% Plot
%%%%%%%%%%
% <<<<<<<<
%%%%%%%%%%
% Wait for figure
pause(1);

% Real trajectory
for i = 1:Nsteps
    x_line = plot(x_mpc(1:i, 1), x_mpc(1:i, 2), 'blue', 'LineWidth', 1);
    x_line.Color(4) = 0.5; % line transparency 50%
    hold on;
    x_points = scatter(x_mpc(1:i, 1), x_mpc(1:i, 2), 5, 'blue', 'filled');
    hold on;
    quiver(x_mpc(1:i, 1), x_mpc(1:i, 2), arrow_length * cos(x_mpc(1:i, 7)), arrow_length * sin(x_mpc(1:i, 7)), 'AutoScale', 'off', 'Color', 'blue');
    hold on;
    
    % legend([ref_points, x_points],{'Reference trajectory', 'Real trajectory'}, 'Location', 'northwest');
    legend([guide_points, reference_points, x_points],{'Guide Points', 'Reference trajectory', 'Real trajectory'}, 'Location', 'northwest');
    hold on;

    pause(0.05);
    if i < Nsteps
        delete(x_line);
    end
end