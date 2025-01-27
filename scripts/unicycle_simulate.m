% Trajectory generation and simulation for the unicycle model
clear all
clc

% Unicycle model parameters ────────────────────────────────────────────────────
r  = 0.03;                          % wheel radius
L  = 0.3;                           % distance between wheels
Ts = 0.1;                           % sampling time
x_constraints = 100*[
    -2 2;                           % x position constraints
    -2 2;                           % y position constraints
    -pi 3*pi;                       % heading angle constraints
];
u_constraints = 100*[-30, 30];      % angular velocity constraints
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


% % Leminscate trajectory
% N_guide = 100;
% a = 1;
% shape = "leminscate";
% [x_ref, u_ref, Tend] = model.generate_trajectory(N_guide, shape, a);


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

% Loop first reference as last reference
x_ref = [x_ref; x_ref(1, :)];
u_ref = [u_ref; u_ref(1, :)];
Tend = Tend + Ts;


% Simulate the system ──────────────────────────────────────────────────────────
Nsteps = length(x_ref);             % number of steps
x = zeros(Nsteps, states);          % state vector
x(1, :) = x_ref(1, :);              % initial state
MSE_x = zeros(Nsteps-1, 1);         % mean squared error for state
for i = 1:Nsteps-1                  % Simulate the system
    x(i+1, :) = model.simulate(x(i, :)', u_ref(i, :)', Ts);
    MSE_x(i) = norm(x(i, :) - x_ref(i, :))^2;
end

max_MSE = max(MSE_x);               % maximum mean squared error for state
disp(['Max Mean Squared Error for State: ', num2str(max_MSE)]);
MSE_x = mean(MSE_x);                % mean squared error for state
disp(['Mean Squared Error for State: ', num2str(MSE_x)]);


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
title('Trajectory Simulation for the Unicycle Model');
xlabel('x'); ylabel('y');
grid on;
axis equal;
hold on;

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
    legend([ref_points, x_points],{'Reference trajectory', 'Real trajectory'}, 'Location', 'northwest');
    hold on;

    pause(0.05);
    if i < Nsteps
        delete(x_line);
    end
end


% % GIF ──────────────────────────────────────────────────────────────────────────
% % Main trajectory plot
% figure(1);
% filename = 'images/unicycle_simulation.gif'; % Output GIF filename
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
%     legend([ref_points, x_points], {'Reference trajectory', 'Real trajectory'}, 'Location', 'northwest');
%     hold on;
% 
%     % Capture frame for GIF
%     frame = getframe(gcf);
%     img = frame2im(frame);
%     [imind, cm] = rgb2ind(img, 256);
%     imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
% 
%     pause(0.05);
%     if i < Nsteps
%         delete(x_line);
%     end
% end
