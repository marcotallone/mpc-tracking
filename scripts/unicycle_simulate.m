% ┌─────────────────────────────────────────────────────────────────────────┐ %
% │                            MPC Tracking Main                            │ %
% └─────────────────────────────────────────────────────────────────────────┘ %
clear all
clc

% Unicycle model parameters ────────────────────────────────────────────────────
r  = 0.05;                           % wheel radius
L  = 0.2;                           % distance between wheels
Ts = 0.1;                           % sampling time
x_constraints = 100*[
    -2 2;                           % x position constraints
    -2 2;                           % y position constraints
    -pi 3*pi;                        % heading angle constraints
];
u_constraints = 100*[-30, 30];          % angular velocity constraints
states = 3;
outputs = 2;
Q_tilde = 0.75*1e-3*eye(states);    % Process noise covariance
R_tilde = 1e-2*eye(outputs);        % Measurement noise covariance
P0 = eye(states);                   % Initial state covariance

model = Unicycle(r, L, Ts, x_constraints, u_constraints, P0, Q_tilde, R_tilde);

% Reference Trajectory Generation ──────────────────────────────────────────────

% % Circle trajectory
% N_guide = 100;
% radius = 0.5;
% shape = "circle";
% [x_ref, u_ref, Tend] = model.generate_trajectory(N_guide, shape, radius);

% % Leminscate trajectory
% N_guide = 100;
% a = 0.5;
% shape = "leminscate";
% [x_ref, u_ref, Tend] = model.generate_trajectory(N_guide, shape, a);

% Batman trajectory
N_guide = 200;
shape = "batman";
[x_ref, u_ref, Tend] = model.generate_trajectory(N_guide, shape);

% % Arbitrary trajectory
% N_intervals = 24;
% N_guide = N_intervals + 1;
% N_points_filling = 5; % Number of points between guide points
% Tend = N_intervals * N_points_filling * Ts;
% T_guide = linspace(0, Tend, N_guide);
% Z_guide = zeros(N_guide, outputs);
% theta = 0;
% delta = 2*pi/N_intervals;
% radius = 0.5;
% for i = 1:N_guide
%     Z_guide(i, :) = [radius*cos(theta), radius*sin(theta), 0.5*pi + theta];
%     theta = theta + delta;
% end 
% shape = "arbitrary";
% [x_ref, u_ref, Tend] = model.generate_trajectory(N_guide, shape, {N_points_filling, Z_guide});

% % Multiply periodic  references for multiple laps
% n_laps = 2;
% x_ref = repmat(x_ref, n_laps, 1);
% u_ref = repmat(u_ref, n_laps, 1);
% Tend = Tend*n_laps;

% Simulate the system ──────────────────────────────────────────────────────────
x0 = x_ref(1, :)';                  % initial state
Nsteps = length(x_ref);             % number of steps
x = zeros(Nsteps, states);          % state vector
x(1, :) = x0';                      % initial state

for i = 1:Nsteps-1
    % Control input
    u = u_ref(i, :)';
    
    % Simulate the system
    x(i+1, :) = model.simulate(x(i, :)', u, Ts);
end

% Plot

% Main trajectory plot
figure(1);

% Reference trajectory
ref_points = scatter(x_ref(:, 1), x_ref(:, 2), 5, 'filled', 'MarkerFaceColor', '#808080');
hold on;
arrow_length = 0.01;
for i = 1:length(x_ref)
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
