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
model = Unicycle(r, L, Ts, x_constraints, u_constraints);

% Simulation parameters
x0 = [1; 0; pi/2;];                 % initial state
Tend = 90*Ts;                       % simulation time
t = 0:Ts:Tend;                      % vector of time steps

% Reference trajectory
x_ref = zeros(length(t), model.n);
u_ref = zeros(length(t), model.m);
for i = 1:length(t)
    x_ref(i, :) = [cos(t(i)), sin(t(i)), wrapTo2Pi(0.5 * pi + t(i))];
    u_ref(i, :) = [(2-L)/(2*r), (2+L)/(2*r)];
end

% Covariance matrices for noise

Qtilde = 1e-4*eye(model.n);           % Process noise covariance
Rtilde = 1e-4*eye(model.p);           % Measurement noise covariance
P0 = eye(model.n);               % Initial state covariance

% Noise components genration for Kalman filter
% model state n = 3, model output measurements = 3
function [w, v] = generateNoise(model, Qtilde, Rtilde, Tend, Ts)
    w = mvnrnd(zeros(model.n, 1), Qtilde, length(0:Ts:Tend))';
    v = mvnrnd(zeros(model.p, 1), Rtilde, length(0:Ts:Tend))';
end

% Initial conditions
[w, v] = generateNoise(model, Qtilde, Rtilde, Tend, Ts);
w = w'; v = v';


% Kalman filter
model.initEKF(x0, P0, Qtilde, Rtilde);


% Simulation loop 
x = [x0'];
x_sim = [(x0 + w(1,:)')'];
x_hat = [(x0 + w(1,:)')'];

for i = 1:length(t)

    % Simulate the real system
    x_real = model.simulate(x(end, :)', u_ref(i, :)', Ts);

    % Simulate with noise
    x_simulated = x_real + w(i, :)';

    % EKF step
    y = model.output(x_simulated, u_ref(i, :)') + v(i, :)';
    x_kalman = model.stepEKF(y, u_ref(i, :)');

    % Update arrays
    x = [x; x_real'];
    x_sim = [x_sim; x_simulated'];
    x_hat = [x_hat; x_kalman'];

end


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
Nsteps = length(t);
for i = 1:Nsteps

    % % Real trajectory

    % x_line = plot(x(1:i, 1), x(1:i, 2), 'blue', 'LineWidth', 1);
    % x_line.Color(4) = 0.5; % line transparency 50%
    % hold on;
    % x_points = scatter(x(1:i, 1), x(1:i, 2), 5, 'blue', 'filled');
    % hold on;
    % quiver(x(1:i, 1), x(1:i, 2), arrow_length * cos(x(1:i, 3)), arrow_length * sin(x(1:i, 3)), 'AutoScale', 'off', 'Color', 'blue');
    % legend([ref_points, x_points],{'Reference trajectory', 'Real trajectory'}, 'Location', 'northwest');
    % hold on;

    % pause(0.05);
    % if i < Nsteps
    %     delete(x_line);
    % end

    % Simulated noisy trajectory
    x_sim_line = plot(x_sim(1:i, 1), x_sim(1:i, 2), 'blue', 'LineWidth', 1);
    x_sim_line.Color(4) = 0.5; % line transparency 50%
    hold on;
    x_sim_points = scatter(x_sim(1:i, 1), x_sim(1:i, 2), 5, 'blue', 'filled');
    hold on;
    quiver(x_sim(1:i, 1), x_sim(1:i, 2), arrow_length * cos(x_sim(1:i, 3)), arrow_length * sin(x_sim(1:i, 3)), 'AutoScale', 'off', 'Color', 'blue');
    legend([ref_points, x_sim_points],{'Reference trajectory', 'Simulated noisy trajectory'}, 'Location', 'northwest');
    hold on;

    % Kalman filter estimated trajectory
    x_hat_line = plot(x_hat(1:i, 1), x_hat(1:i, 2), 'red', 'LineWidth', 1);
    x_hat_line.Color(4) = 0.5; % line transparency 50%
    hold on;
    x_hat_points = scatter(x_hat(1:i, 1), x_hat(1:i, 2), 5, 'red', 'filled'); 
    hold on;
    quiver(x_hat(1:i, 1), x_hat(1:i, 2), arrow_length * cos(x_hat(1:i, 3)), arrow_length * sin(x_hat(1:i, 3)), 'AutoScale', 'off', 'Color', 'red');
    legend([ref_points, x_sim_points, x_hat_points],{'Reference trajectory', 'Simulated noisy trajectory', 'Kalman filter estimated trajectory'}, 'Location', 'northwest');
    hold on;

    pause(0.05);
    if i < Nsteps
        delete(x_sim_line);
        delete(x_hat_line);
    end

end