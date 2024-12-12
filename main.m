% MPC Reference Tracking with Unicycle Model
clear all
clc

% Model parameters -------------------------------------------------------------
r  = 0.1;                           % wheel radius
L  = 0.5;                           % distance between wheels
Ts = 0.1;                           % sampling time
omega_limits = [-30, 30];           % angular velocity limits
x_limits = [-2, 2];                 % x position limits
y_limits = [-2, 2];                 % y position limits
theta_limits = [-2*pi, 2*pi];           % heading angle limits

% Initialize unicycle model object
unicycle = Unicycle(r, L, Ts, omega_limits, x_limits, y_limits, theta_limits);
n = unicycle.n;                     % number of states
m = unicycle.m;                     % number of inputs

% Simulation parameters and MPC parameters -------------------------------------
x0 = [0; 0; 2*pi;];               % initial state
N  = 3;                            % prediction horizon
t  = 0:Ts:(30*Ts);                  % time vector
T  = length(t)-(N+1);               % simulation time (number of steps)
Q  = 1000*eye(n);                   % state cost
% Q = [
%     1000 0 0;
%     0 1000 0;
%     0 0 0.0001];
R  = eye(m);                        % input cost
preview = 1;                        % MPC preview flag
formulation = 0;                    % formulation to use for MPC problem

% Quadprog optimization options
options = optimoptions('quadprog', 'Display', 'none');  % silent
% options = optimoptions('quadprog', 'Display', 'iter');  % verbose

% Reference trajectory ---------------------------------------------------------
x_ref = zeros(length(t), n);
u_ref = zeros(length(t), m);
for i = 1:length(t)
    x_ref(i, :) = [cos(t(i)), sin(t(i)), wrapTo2Pi(0.5 * pi + t(i))];
    u_ref(i, :) = [(2-L)/(2*r), (2+L)/(2*r)];
end

% MPC Reference Tracking -------------------------------------------------------

% Initialize the MPC optimization problem and solve it
mpc = MPC(Q, R, N, unicycle, x_ref, u_ref, x0, Ts, t, T, preview, formulation, options);
[x, u] = mpc.optimize();

% Plot -------------------------------------------------------------------------

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

% Real trajectory
for i = 1:T
    x_line = plot(x(1:i, 1), x(1:i, 2), 'blue', 'LineWidth', 1);
    x_line.Color(4) = 0.5; % line transparency 50%
    hold on;
    x_points = scatter(x(1:i, 1), x(1:i, 2), 5, 'blue', 'filled');
    hold on;
    quiver(x(1:i, 1), x(1:i, 2), arrow_length * cos(x(1:i, 3)), arrow_length * sin(x(1:i, 3)), 'AutoScale', 'off', 'Color', 'blue');
    legend([ref_points, x_points],{'Reference trajectory', 'Real trajectory'}, 'Location', 'northwest');
    hold on;

    pause(0.05);
    if i < T
        delete(x_line);
    end
end