% Trajectory generation and simulation for the helicopter model
clear all
clc

% Helicopter model parameters ──────────────────────────────────────────────────
bx = 2; by = 2; bz = 18; bpsi = 111;
kx = -0.5; ky = -0.5; kpsi = -5; ki = 2;
parameters = [bx; by; bz; bpsi; kx; ky; kpsi; ki];
g = 9.81;
Ts = 0.1;                           % sampling time
x_constraints = 4*[                 % state constraints
    -10, 10;
    -10, 10;
    -10, 10;
    -10, 10;
    -10, 10;
    -10, 10;
    -pi, 3*pi;
    -25, 25;
];
u_constraints = 100*[               % input constraints
    -1, 1;
    -1, 1;
    -1, 1;
    -1, 1;
];

states = 8;                         % number of states
outputs = 4;                        % number of outputs
Q_tilde = 0.75*1e-3*eye(states);    % Process noise covariance
R_tilde = 1e-2*eye(outputs);        % Measurement noise covariance
P0 = eye(states);                   % Initial state covariance

model = Helicopter(parameters, Ts, x_constraints, u_constraints, P0, Q_tilde, R_tilde);


% Reference Trajectory Generation ──────────────────────────────────────────────
% (comment / uncomment the desired trajectory)


% Circle trajectory
N_guide = 100;
radius = 0.5;
shape = "circle";
[x_ref, u_ref, Tend] = model.generate_trajectory(N_guide, shape, radius);


% % Lemniscate trajectory
% N_guide = 100;
% a = 1;
% shape = "lemniscate";
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


% Simulate the system ──────────────────────────────────────────────────────────
Nsteps = length(x_ref);             % number of steps
x = zeros(Nsteps, states);          % state vector
x(1, :) = x_ref(1, :);              % initial state
for i = 1:Nsteps-1                  % Simulate the system
    x(i+1, :) = model.simulate(x(i, :)', u_ref(i, :)', Ts);
end


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
title('Trajectory Simulation for the Helicopter Model');
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


% GIF ──────────────────────────────────────────────────────────────────────────
% % Main trajectory plot
% figure(1);
% filename = 'images/helicopter_simulation.gif'; % Output GIF filename

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

% % Labels
% title('Trajectory Simulation for the Helicopter Model');
% xlabel('x'); ylabel('y');
% grid on;
% axis equal;
% hold on;

% % Capture initial frame for GIF
% frame = getframe(gcf);
% img = frame2im(frame);
% [imind, cm] = rgb2ind(img, 256); 
% imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.05);

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

%     % Capture frame for GIF
%     frame = getframe(gcf);
%     img = frame2im(frame);
%     [imind, cm] = rgb2ind(img, 256);
%     imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);

%     pause(0.05);
%     if i < Nsteps
%         delete(x_line);
%     end
% end