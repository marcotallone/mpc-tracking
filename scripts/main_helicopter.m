% ┌─────────────────────────────────────────────────────────────────────────┐ %
% │                           MPC Helicopter Main                           │ %
% └─────────────────────────────────────────────────────────────────────────┘ %
clear all
clc


% Helicopter model parameters ──────────────────────────────────────────────────
bx = 2;
by = 2;
bz = 18;
bpsi = 111;
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
    -3, 3;
    -3, 3;
    -2, 2;
     0, 2*pi;
    -25, 25;
    % -5, 5;
    % -5, 5;
];
u_constraints = [
    -1, 1;
    -1, 1;
    -1, 1;
      0, 1;
];
model = Helicopter(parameters, Ts, x_constraints, u_constraints);


% Guide points for reference trajectory ────────────────────────────────────────
nf = 4; % Number of flat-outputs

% Generation parameters
N_intervals = 24;
N_guide = N_intervals + 1; % Number of guide points (last == first)
N_points_filling = 5; % Number of points between guide points
Tend = N_intervals * N_points_filling * Ts;

% Guide
T_guide = linspace(0, Tend, N_guide);
Z_guide = zeros(N_guide, nf);

theta = 0;
delta = 2*pi/N_intervals;
radius = 0.5;
for i = 1:N_guide
    Z_guide(i, :) = [radius*cos(theta), radius*sin(theta), 0, 0.5*pi + theta];
    theta = theta + delta;
end 


% Reference trajectory filling ─────────────────────────────────────────────────
N_basis = 2; % Number of basis functions
order = 1;   % Maximum differentiation order

% Guide vectors for first and second derivatives
dZ_guide = [zeros(1, nf)];
% ddZ_guide = [zeros(1, nf)];

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
for i = 1:2 %length(T_guide) - 1

    % Time interval
    t0 = T_guide(i);
    t1 = T_guide(i+1);

    % Matrix M
    m0 = m(t0, basis, order);
    m1 = m(t1, basis, 0);
    M0 = kron(eye(nf), m0);
    M1 = kron(eye(nf), m1);
    M = [M0; M1];

    % Guide points in the interval
    z_bar = [
        Z_guide(i, 1:nf)';
        dZ_guide(i, 1:nf)';
        % ddZ_guide(i, 1:nf)';
        Z_guide(i+1, 1:nf)';
    ];

    % Solve the system and reshape
    alpha = M\z_bar;
    alpha = reshape(alpha, [], nf)';
    % alpha_tensor(:, :, i) = alpha;

    % Reference function z(t) and its derivatives
    z = alpha*basis.';
    dz = diff(z, x);
    ddz = diff(dz, x);

    % Append to guide derivatives values for next interval
    dZ_guide = [dZ_guide; double(subs(dz, x, t1))'];
    % ddZ_guide = [ddZ_guide; double(subs(ddz, x, t1))'];

    % Generate missing points
    N_filling = ceil((t1 - t0)/Ts) + 1;
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
        % dxint = 0.0;    % augmented states reference always zero
        % dyint = 0.0;    % augmented states reference always zero

        x_t = [z_t(1), z_t(2), z_t(3), dxb, dyb, dzb, z_t(4), dpsi]; %, dxint, dyint];

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
arrow_length = 0.01;

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


% Simulate system to test trajectory ───────────────────────────────────────────
x0 = x_ref(1, :)';
% x = model.simulate(x0, u_ref, Tend);
for i = 1:length(u_ref)
    x_sim = model.simulate(x0, u_ref(i, :)', Ts);
    x0 = x_sim;

    % Plot
    x_line = plot(x_sim(1), x_sim(2), 'green', 'LineWidth', 1);
    x_line.Color(4) = 0.5; % line transparency 50%
    hold on;
    x_points = scatter(x_sim(1), x_sim(2), 5, 'green', 'filled');
    hold on;
    quiver(x_sim(1), x_sim(2), arrow_length * cos(x_sim(7)), arrow_length * sin(x_sim(7)), 'AutoScale', 'off', 'Color', 'green');
    hold on;
    target = scatter(Z_ref(i, 1), Z_ref(i, 2), 20, 'filled', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'red');
    hold on;
    legend([guide_points, reference_points, x_points, target],{'Guide Points', 'Reference trajectory', 'Real trajectory', 'Target'}, 'Location', 'northwest');
    hold on;

    pause(0.05);
    if i < length(u_ref)
        delete(x_line);
        delete(target);
    end
end


% MPC ──────────────────────────────────────────────────────────────────────────

% Multiply references for multiple laps
n_laps = 1;
x_ref = repmat(x_ref, n_laps, 1);
u_ref = repmat(u_ref, n_laps, 1);
Tend = Tend*n_laps;

% MPC parameters
x0 = zeros(model.n,1);
% x0(9) = x0(1) - x_ref(1,1);
% x0(10) = x0(2) - x_ref(1,2);
% x0 = [0.5000; 0; 0; 0.2588; 0.0341; 0; 1.5708; 0.5236; 0; 0];
N  = 18;
Q = diag([50, 50, 5, 10, 3, 3, 1, 2]);  % 15, 15]);
R = diag([2, 2, 2, 2]);
preview = 1;                        % MPC preview flag
formulation = 0;                    % MPC formulation flag
debug = 0;                          % MPC debug flag
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
    target = scatter(Z_ref(i, 1), Z_ref(i, 2), 20, 'filled', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'red');
    hold on;
    legend([guide_points, reference_points, x_points, target],{'Guide Points', 'Reference trajectory', 'Real trajectory', 'Target'}, 'Location', 'northwest');
    hold on;

    pause(0.05);
    if i < Nsteps
        delete(x_line);
        delete(target);
    end
end