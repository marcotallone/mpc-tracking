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
Q_tilde = 0.25*1e-3*eye(model_n);        % Process noise covariance
R_tilde = 1e-5*eye(model_p);        % Measurement noise covariance
P0 = eye(model_n);                  % Initial state covariance

model = Unicycle(r, L, Ts, x_constraints, u_constraints, P0, Q_tilde, R_tilde);

% Simulation parameters
x0 = [1; 0; pi/2;];                 % initial state
Tend = 63*Ts;                       % simulation time
t = 0:Ts:Tend;                      % vector of time steps

% Reference trajectory
x_ref = zeros(length(t), model.n);
u_ref = zeros(length(t), model.m);
for i = 1:length(t)
    x_ref(i, :) = [cos(t(i)), sin(t(i)), wrapTo2Pi(0.5 * pi + t(i))];
    u_ref(i, :) = [(2-L)/(2*r), (2+L)/(2*r)];
end

% Noise components generation for Kalman filter
w = mvnrnd(zeros(model.n, 1), Q_tilde, length(0:Ts:Tend));
v = mvnrnd(zeros(model.p, 1), R_tilde, length(0:Ts:Tend));


% Simulation arrays

% Simulation without noise
x = [x0'];

% Simulation real system with noise
x_real = [x0'];

% Simulation Kalman filter with noise
% x_hat = [ ( x0 + mvnrnd(zeros(model.n, 1), P0) )' ];
% x_hat = [ ( x0 + (randn(1, model.n) * chol(P0))' )' ];
x_hat = [ ( x0 + w(1,:)' )' ];

% Simulation without Kalman filter
% x_sim = [x_hat(1, :)];


% Simulation loop 
for i = 2:length(t)

    % Simulate the system without noise
    x = [x; (model.simulate(x(end, :)', u_ref(i, :)', Ts))'];

    % Simulate the real system
    update_real = model.simulate(x_real(end, :)', u_ref(i, :)', Ts) + w(i, :)';
    x_real = [x_real; update_real'];

    % Extended Kalman filter step
    y = model.output(update_real, u_ref(i, :)') + v(i, :)';
    x_hat = [x_hat; (model.EKF_step(x_hat(end, :)', u_ref(i, :)', y))'];

    % Simulate the system without Kalman filter
    % x_sim = [x_sim; (model.simulate(x_hat(i-1, :)', u_ref(i, :)', Ts))'];
    % x_sim = [x_sim; (model.simulate(x_hat(i-1, :)', u_ref(i, :)', Ts))'];

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
xlabel('x'); ylabel('y');
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

    % Simulation without noise (black dotted line only)
    x_line = plot(x(1:i, 1), x(1:i, 2), 'black', 'LineWidth', 1, 'LineStyle', '--');
    x_line.Color(4) = 0.5; % line transparency 50%
    hold on;
    legend([ref_points, x_line],{'Reference trajectory', 'No Noise'}, 'Location', 'northwest');
    hold on;

    % Simulation real system with noise
    x_real_line = plot(x_real(1:i, 1), x_real(1:i, 2), 'blue', 'LineWidth', 1);
    x_real_line.Color(4) = 0.5; % line transparency 50%
    hold on;
    x_real_points = scatter(x_real(1:i, 1), x_real(1:i, 2), 5, 'blue', 'filled');
    hold on;
    quiver(x_real(1:i, 1), x_real(1:i, 2), arrow_length * cos(x_real(1:i, 3)), arrow_length * sin(x_real(1:i, 3)), 'AutoScale', 'off', 'Color', 'blue');
    legend([ref_points, x_line, x_real_points],{'Reference trajectory', 'No Noise', 'With noise'}, 'Location', 'northwest');
    hold on;

    % Kalman filter estimated trajectory
    x_hat_line = plot(x_hat(1:i, 1), x_hat(1:i, 2), 'red', 'LineWidth', 1);
    x_hat_line.Color(4) = 0.5; % line transparency 50%
    hold on;
    x_hat_points = scatter(x_hat(1:i, 1), x_hat(1:i, 2), 5, 'red', 'filled');
    hold on;
    quiver(x_hat(1:i, 1), x_hat(1:i, 2), arrow_length * cos(x_hat(1:i, 3)), arrow_length * sin(x_hat(1:i, 3)), 'AutoScale', 'off', 'Color', 'red');
    legend([ref_points, x_line, x_real_points, x_hat_points],{'Reference trajectory', 'No Noise', 'With noise', 'Kalman filter'}, 'Location', 'northwest');
    hold on;

    % Simulation without Kalman filter
    % x_no_kal_line = plot(x_sim(1:i, 1), x_sim(1:i, 2), 'green', 'LineWidth', 1);
    % x_no_kal_line.Color(4) = 0.5; % line transparency 50%
    % hold on;
    % x_no_kal_points = scatter(x_sim(1:i, 1), x_sim(1:i, 2), 5, 'green', 'filled');
    % hold on;
    % quiver(x_sim(1:i, 1), x_sim(1:i, 2), arrow_length * cos(x_sim(1:i, 3)), arrow_length * sin(x_sim(1:i, 3)), 'AutoScale', 'off', 'Color', 'green');
    % legend([ref_points, x_line, x_real_points, x_hat_points, x_no_kal_points],{'Reference trajectory', 'No Noise', 'With noise', 'Kalman filter', 'No Kalman filter'}, 'Location', 'northwest');
    % hold on;

    pause(0.05);
    if i < Nsteps
        delete(x_line);
        delete(x_real_line);
        delete(x_hat_line);
        % delete(x_no_kal_line);
    end

end


% % Plot the state and the estimated state for x, y, theta over time
figure(2);

% x
subplot(3, 1, 1);
plot(t, x(:, 1), 'black', 'LineWidth', 1, 'LineStyle', '--');
hold on;
plot(t, x_real(:, 1), 'blue', 'LineWidth', 1);
hold on;
plot(t, x_hat(:, 1), 'red', 'LineWidth', 1);
hold on;
% plot(t, x_sim(:, 1), 'green', 'LineWidth', 1);
% hold on;
xlabel('Time [s]');
ylabel('x');
legend('No Noise', 'With noise', 'Kalman filter');
grid on;

% y
subplot(3, 1, 2);
plot(t, x(:, 2), 'black', 'LineWidth', 1, 'LineStyle', '--');
hold on;
plot(t, x_real(:, 2), 'blue', 'LineWidth', 1);
hold on;
plot(t, x_hat(:, 2), 'red', 'LineWidth', 1);
hold on;
% plot(t, x_sim(:, 2), 'green', 'LineWidth', 1);
% hold on;
xlabel('Time [s]');
ylabel('y');
legend('No Noise', 'With noise', 'Kalman filter');
grid on;

% theta
subplot(3, 1, 3);
plot(t, x(:, 3), 'black', 'LineWidth', 1, 'LineStyle', '--');
hold on;
plot(t, x_real(:, 3), 'blue', 'LineWidth', 1);
hold on;
plot(t, x_hat(:, 3), 'red', 'LineWidth', 1);
hold on;
% plot(t, x_sim(:, 3), 'green', 'LineWidth', 1);
% hold on;
xlabel('Time [s]');
ylabel('theta');
legend('No Noise', 'With noise', 'Kalman filter');
grid on;
