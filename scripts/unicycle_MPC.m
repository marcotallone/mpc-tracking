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


% % Circle trajectory
% N_guide = 100;
% radius = 0.5;
% shape = "circle";
% [x_ref, u_ref, Tend] = model.generate_trajectory(N_guide, shape, radius);
% max_start = 0.05;


% Leminscate trajectory
N_guide = 100;
a = 1;
shape = "leminscate";
[x_ref, u_ref, Tend] = model.generate_trajectory(N_guide, shape, a);
max_start = 0.05;


% % Arbitrary trajectory
% N_guide = 9;
% Z_guide = [
%     1, 1;
%     2.14, 1.85;
%     3.4, 1.91;
%     4.66, 1.41;
%     6, 1;
%     7.04, 1.43;
%     8, 2;
%     8.92, 2.63;
%     9.96, 2.37;
% ];
% N_points_filling = 3;
% shape = "arbitrary";
% N_basis = 2;
% order = 0;
% [x_ref, u_ref, Tend] = model.generate_trajectory(N_guide, shape, {N_points_filling, N_basis, order, Z_guide});


% % Multiply periodic references for multiple laps
% n_laps = 2;
% x_ref = repmat(x_ref, n_laps, 1);
% u_ref = repmat(u_ref, n_laps, 1);
% Tend = Tend*n_laps;


% MPC ──────────────────────────────────────────────────────────────────────────

% MPC parameters
% x0 = zeros(model.n,1);            % origin initial state
% x0 = x_ref(1, :)';                % first reference initial state
% x0 = [1; 0; pi/4];                  % custom initial state (x, y, theta)


N  = 10;                            % prediction horizon
Q  = 1e3*eye(model.n);              % state cost
R  = eye(model.m);                  % input cost
preview = 1;                        % MPC preview flag
formulation = 0;                    % MPC formulation flag
noise = 0;                          % MPC noise flag
debug = 0;                          % MPC debug flag
t = 0:Ts:Tend;                      % vector of time steps
Nsteps = length(t) - (N+1);         % number of MPC optimization steps

for index = 1:10

    % Set initial condition to a random point around x_ref(1, :) inside the ball of radius max_start
    x0 = x_ref(1, :)' + max_start*randn(3, 1);


    % Optimization
    mpc = MPC(model, x0, Tend, N, Q, R, x_ref, u_ref, preview, formulation, noise, debug);
    [x, u] = mpc.optimize();

end


% % Plot
% 
% % Main trajectory plot
% figure(1);
% 
% % Reference trajectory
% ref_points = scatter(x_ref(:, 1), x_ref(:, 2), 5, 'filled', 'MarkerFaceColor', '#808080');
% hold on;
% arrow_length = 0.01;
% for i = 1:length(x_ref)
%     x_arrow = arrow_length * cos(x_ref(i, 3));
%     y_arrow = arrow_length * sin(x_ref(i, 3));
%     quiver(x_ref(i, 1), x_ref(i, 2), x_arrow, y_arrow, 'AutoScale', 'off', 'Color', '#808080');
% end
% legend(ref_points,{'Reference trajectory'}, 'Location', 'northwest');
% 
% % Labels
% title('Trajectory Tracking with MPC (Non-Linear Unicycle System)');
% xlabel('x1'); ylabel('x2');
% grid on;
% axis equal;
% hold on;
% 
% %%%%%%%%%%
% % <<<<<<<<
% %%%%%%%%%%
% % Wait for figure
% pause(1);
% 
% % Real trajectory
% for i = 1:Nsteps
%     x_line = plot(x(1:i, 1), x(1:i, 2), 'blue', 'LineWidth', 1);
%     x_line.Color(4) = 0.5; % line transparency 50%
%     hold on;
%     x_points = scatter(x(1:i, 1), x(1:i, 2), 5, 'blue', 'filled');
%     hold on;
%     quiver(x(1:i, 1), x(1:i, 2), arrow_length * cos(x(1:i, 3)), arrow_length * sin(x(1:i, 3)), 'AutoScale', 'off', 'Color', 'blue');
%     legend([ref_points, x_points],{'Reference trajectory', 'Real trajectory'}, 'Location', 'northwest');
%     hold on;
% 
%     pause(0.05);
%     if i < Nsteps
%         delete(x_line);
%     end
% end




% GIF
% 
% % Plot ─────────────────────────────────────────────────────────────────────────
% % Main trajectory plot
% figure(1);
% filename = 'images/unicycle_MPC.gif'; % Output GIF filename
% 
% % Reference trajectory
% ref_points = scatter(x_ref(:, 1), x_ref(:, 2), 5, 'filled', 'MarkerFaceColor', '#808080');
% hold on;
% arrow_length = 0.03;
% for i = 1:length(x_ref)
%     x_arrow = arrow_length * cos(x_ref(i, 3));
%     y_arrow = arrow_length * sin(x_ref(i, 3));
%     quiver(x_ref(i, 1), x_ref(i, 2), x_arrow, y_arrow, 'AutoScale', 'off', 'Color', '#808080');
% end
% legend(ref_points,{'Reference trajectory'}, 'Location', 'northwest');
% 
% % Labels
% title('Trajectory Simulation for the Unicycle Model');
% xlabel('x'); ylabel('y');
% grid on;
% axis equal;
% hold on;
% 
% % Capture initial frame for GIF
% frame = getframe(gcf);
% img = frame2im(frame);
% [imind, cm] = rgb2ind(img, 256); 
% imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.05);
% 
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
%     legend([ref_points, ref_points, x_points, target],{'Guide Points', 'Reference trajectory', 'Real trajectory', 'Target'}, 'Location', 'northwest');
%     hold on;
%     % Capture frame for GIF
%     frame = getframe(gcf);
%     img = frame2im(frame);
%     [imind, cm] = rgb2ind(img, 256);
%     imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
% 
%     pause(0.05);
%     if i < Nsteps
%         delete(x_line);
%         delete(target);
%     end
% end
