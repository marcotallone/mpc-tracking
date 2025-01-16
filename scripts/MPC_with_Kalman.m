clear all
clc

% Model parameters
r  = 0.1;                           % wheel radius
L  = 0.5;                           % distance between wheels
Ts = 0.1;                           % sampling time
x_constraints = [
    -2 2;                           % x position constraints
    -2 2;                           % y position constraints
     0 2*pi;                        % heading angle constraints
];
u_constraints = [-30, 30];          % angular velocity constraints

model_n = 3;
model_p = 2;
Q_tilde = 0.75*1e-3*eye(model_n);        % Process noise covariance
R_tilde = 1e-2*eye(model_p);        % Measurement noise covariance
P0 = eye(model_n);                  % Initial state covariance

model = Unicycle(r, L, Ts, x_constraints, u_constraints, P0, Q_tilde, R_tilde);

% Simulation parameters
x0 = [0; 0; pi/4;];                 % initial state
Tend = 90*Ts;                       % simulation time
t = 0:Ts:Tend;                      % vector of time steps

% Reference trajectory
x_ref = zeros(length(t), model.n);
u_ref = zeros(length(t), model.m);
for i = 1:length(t)
    x_ref(i, :) = [cos(t(i)), sin(t(i)), wrapTo2Pi(0.5 * pi + t(i))];
    u_ref(i, :) = [(2-L)/(2*r), (2+L)/(2*r)];
end

% Noise components generation for Kalman filter
% w = mvnrnd(zeros(model.n, 1), Q_tilde, length(0:Ts:Tend));
% v = mvnrnd(zeros(model.p, 1), R_tilde, length(0:Ts:Tend));

% MPC parameters
N  = 5;                             % prediction horizon
Q  = 1000*eye(model.n);             % state cost
R  = eye(model.m);                  % input cost
preview = 1;                        % MPC preview flag
formulation = 0;                    % MPC formulation flag
noise = 1;                          % MPC noise flag
debug = 0;                          % MPC debug flag
Nsteps = length(t) - (N+1);         % number of MPC optimization steps

% Optimization
mpc = noiseMPC(model, x0, Tend, N, Q, R, x_ref, u_ref, preview, formulation, noise, debug);
[x, u] = mpc.optimize();

% Plot

% Main trajectory plot
figure(1);

% Reference trajectory
ref_points = scatter(x_ref(:, 1), x_ref(:, 2), 5, 'filled', 'MarkerFaceColor', '#808080');
hold on;
arrow_length = 0.03;
for i = 1:length(t)
    x_arrow = arrow_length * cos(x_ref(i, 3));
    y_arrow = arrow_length * sin(x_ref(i, 3));
    quiver(x_ref(i, 1), x_ref(i, 2), x_arrow, y_arrow, 'AutoScale', 'off', 'Color', '#808080');
end
legend(ref_points,{'Reference trajectory'}, 'Location', 'northwest');

% Labels
title('Trajectory Tracking with MPC (Non-Linear Unicycle System)');
xlabel('x1'); ylabel('x2');
grid on;
axis equal;
hold on;

%%%%%%%%%%
% <<<<<<<<
%%%%%%%%%%
% Wait for figure
pause(1);

% Real trajectory
for i = 1:Nsteps
    x_line = plot(x(1:i, 1), x(1:i, 2), 'blue', 'LineWidth', 1);
    x_line.Color(4) = 0.5; % line transparency 50%
    hold on;
    x_points = scatter(x(1:i, 1), x(1:i, 2), 5, 'blue', 'filled');
    hold on;
    quiver(x(1:i, 1), x(1:i, 2), arrow_length * cos(x(1:i, 3)), arrow_length * sin(x(1:i, 3)), 'AutoScale', 'off', 'Color', 'blue');
    legend([ref_points, x_points],{'Reference trajectory', 'Real trajectory'}, 'Location', 'northwest');
    hold on;

    pause(0.05);
    if i < Nsteps
        delete(x_line);
    end
end