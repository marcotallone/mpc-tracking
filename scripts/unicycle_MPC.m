% MPC tracking of a unicycle model
clear all
clc


% Unicycle model parameters ────────────────────────────────────────────────────
r  = 0.03;                          % wheel radius
L  = 0.3;                           % distance between wheels
Ts = 0.1;                           % sampling time
x_constraints = [
    -2 2;                           % x position constraints
    -2 2;                           % y position constraints
    -pi 3*pi;                       % heading angle constraints
];
u_constraints = [-50, 50];          % angular velocity constraints
states = 3;                         % number of states
outputs = 2;                        % number of outputs
Q_tilde = 0.75*1e-3*eye(states);    % Process noise covariance
R_tilde = 1e-2*eye(outputs);        % Measurement noise covariance
P0 = eye(states);                   % Initial state covariance

model = Unicycle(r, L, Ts, x_constraints, u_constraints, P0, Q_tilde, R_tilde);


% Reference Trajectory Generation ──────────────────────────────────────────────
% (comment / uncomment the desired trajectory)


% Circle trajectory
N_guide = 100;
radius = 0.5;
shape = "circle";
[x_ref, u_ref, Tend] = model.generate_trajectory(N_guide, shape, radius);
max_start = 0.05;


% % Lemniscate trajectory
% N_guide = 100;
% a = 1;
% shape = "lemniscate";
% [x_ref, u_ref, Tend] = model.generate_trajectory(N_guide, shape, a);
% max_start = 0.05;


% % Arbitrary trajectory
% N_guide = 5;
% Z_guide = [
%     1, 0;
%     0, 1;
%     -1, 0;
%     0, -1;
%     1, 0;
% ];
% N_points_filling = 25;
% shape = "arbitrary";
% N_basis = 2;
% order = 0;
% [x_ref, u_ref, Tend] = model.generate_trajectory(N_guide, shape, {N_points_filling, N_basis, order, Z_guide});
% max_start = 0.05;


% % Multiply periodic references for multiple laps
% n_laps = 2;
% x_ref = repmat(x_ref, n_laps, 1);
% u_ref = repmat(u_ref, n_laps, 1);
% Tend = Tend*n_laps;


% MPC ──────────────────────────────────────────────────────────────────────────

% Initial condition
% x0 = zeros(model.n,1);            % origin initial state
% x0 = x_ref(1, :)';                % first reference initial state
% x0 = [1; 0; pi/4];                % custom initial state (x, y, theta)
x0 = x_ref(1, :)' + max_start*rand(states, 1);  % random initial condition

% MPC parameters
N  = 10;                            % prediction horizon
Q  = 1e3*eye(model.n);              % state cost
R  = eye(model.m);                  % input cost
preview = 1;                        % MPC preview flag
formulation = 0;                    % MPC formulation flag
noise = 0;                          % MPC noise flag
debug = 0;                          % MPC debug flag

% Simulation time and steps
x_ref = [x_ref; x_ref(1:N+1, :)];   % add N steps to complete a full loop
u_ref = [u_ref; u_ref(1:N+1, :)];   % add N steps to complete a full loop
Tend = Tend + (N+1)*Ts;
t = 0:Ts:Tend;                      % vector of time steps
Nsteps = length(t) - (N+1);         % number of MPC optimization steps

% Optimization
mpc = MPC(model, x0, Tend, N, Q, R, x_ref, u_ref, preview, formulation, noise, debug);
[x, u] = mpc.optimize();


% Plot ─────────────────────────────────────────────────────────────────────────

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
xlabel('x'); ylabel('y');
grid on;
axis equal;
hold on;

% % Set plot limits
% xlim([-0.6, 0.6]);
% ylim([-0.6, 0.6]);
% hold on;

% ────────────────────
% Wait for figure here
pause(1);

% Real trajectory
for i = 1:Nsteps
    x_line = plot(x(1:i, 1), x(1:i, 2), 'blue', 'LineWidth', 1);
    x_line.Color(4) = 0.5; % line transparency 50%
    hold on;
    x_points = scatter(x(1:i, 1), x(1:i, 2), 5, 'blue', 'filled');
    hold on;
    quiver(x(1:i, 1), x(1:i, 2), arrow_length * cos(x(1:i, 3)), arrow_length * sin(x(1:i, 3)), 'AutoScale', 'off', 'Color', 'blue');
    % target = scatter(x_ref(i, 1), x_ref(i, 2), 20, 'filled', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'red');
    target = scatter(x_ref(i, 1), x_ref(i, 2), 20, 'red', 'filled');
    hold on;
    legend([ref_points, x_points, target],{'Reference trajectory', 'Real trajectory', 'Target'}, 'Location', 'northwest');
    hold on;

    pause(0.05);
    if i < Nsteps
        delete(x_line);
        delete(target);
    end
end


% % GIF ──────────────────────────────────────────────────────────────────────────
% % Main trajectory plot
% figure(1);
% filename = 'images/unicycle_output.gif'; % Output GIF filename

% % Reference trajectory
% ref_points = scatter(x_ref(:, 1), x_ref(:, 2), 5, 'filled', 'MarkerFaceColor', '#808080');
% hold on;
% arrow_length = 0.02;
% for i = 1:length(x_ref)
%     x_arrow = arrow_length * cos(x_ref(i, 3));
%     y_arrow = arrow_length * sin(x_ref(i, 3));
%     quiver(x_ref(i, 1), x_ref(i, 2), x_arrow, y_arrow, 'AutoScale', 'off', 'Color', '#808080');
% end
% % legend(ref_points,{'Reference trajectory'}, 'Location', 'northwest');

% % Labels
% title('Trajectory Tracking with MPC (Non-Linear Unicycle System)');
% xlabel('x'); ylabel('y');
% grid on;
% axis equal;
% hold on;
% axis tight; % Adjust axis limits to fit the data tightly
% hold on;

% % Set plot limits
% % xlim([-1.5, 1.5]);
% xlim([-0.6, 0.6]);
% ylim([-0.6, 0.6]);
% hold on;

% % Adjust figure to fit tightly around the plot
% set(gca, 'LooseInset', get(gca, 'TightInset'));

% % Capture initial frame for GIF
% frame = getframe(gca);
% img = frame2im(frame);
% [imind, cm] = rgb2ind(img, 256); 
% imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);

% % Real trajectory animation and GIF capture
% for i = 1:Nsteps
%     x_line = plot(x(1:i, 1), x(1:i, 2), 'blue', 'LineWidth', 1);
%     x_line.Color(4) = 0.5; % line transparency 50%
%     hold on;
%     x_points = scatter(x(1:i, 1), x(1:i, 2), 5, 'blue', 'filled');
%     hold on;
%     quiver(x(1:i, 1), x(1:i, 2), arrow_length * cos(x(1:i, 3)), arrow_length * sin(x(1:i, 3)), 'AutoScale', 'off', 'Color', 'blue');
%     hold on;
%     target = scatter(x_ref(i, 1), x_ref(i, 2), 20, 'filled', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'red');
%     hold on;
%     legend([ref_points, x_points, target],{'Reference trajectory', 'Real trajectory', 'Target'}, 'Location', 'northwest');
%     hold on;
%     % Capture frame for GIF
%     frame = getframe(gca);
%     img = frame2im(frame);
%     [imind, cm] = rgb2ind(img, 256);
%     imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);

%     % pause(0.05);
%     if i < Nsteps
%         delete(x_line);
%         delete(target);
%     end
% end