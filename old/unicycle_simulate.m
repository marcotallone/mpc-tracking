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

% % Batman trajectory
% N_guide = 100;
% shape = "batman";
% [x_ref, u_ref, Tend] = model.generate_trajectory(N_guide, shape);

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



% Test 1
N_guide = 9;
N_points_filling = 3;
% Tend = (N_guide -1) * N_points_filling * Ts;
% T_guide = linspace(0, Tend, N_guide);
% for i = 1:N_guide
%     Z_guide(i, :) = [T_guide(i), T_guide(i)^2];
% end 
Z_guide = [
    1, 1;
    2.14, 1.85;
    3.4, 1.91;
    4.66, 1.41;
    6, 1;
    7.04, 1.43;
    8, 2;
    8.92, 2.63;
    9.96, 2.37;
];


% Test 2
% 
% N1 = 96;  % Number of points in [-8, -1]
% N2 = 8;   % Number of points in (-1,1)
% N3 = 96;  % Number of points in [1,8]
% N_guide = N1 + N2 + N3;
% 
% x1 = linspace(-8, -2, N1);   % Uniform in [-8, -1]
% x2 = linspace(-2, 2, N2+2);  % Include -1 and 1, remove duplicates later
% x3 = linspace(2, 8, N3);     % Uniform in [1, 8]
% 
% l = unique([x1, x2(2:end-1), x3]);  % Merge and remove duplicate -1 and 1
% 
% Z_guide = zeros(N_guide, 2);
% syms t real;
% 
% x = (abs(t)/t) * (0.3*abs(t) + 0.2*abs(abs(t)-1) + 2.2*abs(abs(t)-2) ...
%     - 2.7*abs(abs(t)-3) - 3*abs(abs(t)-5) + 3*abs(abs(t)-7) ...
%     + 5*sin((pi/4)*(abs(abs(t)-3) - abs(abs(t)-4) + 1)) ...
%     + (5/4)*(abs(abs(t)-4) - abs(abs(t)-5) - 1)^3 ...
%     - 5.3*cos(((pi/2) + asin(47/53)) * ((abs(abs(t)-7) - abs(abs(t)-8) - 1)/2)) ...
%     + 2.8);
% 
% y = (3/2)*abs(abs(t)-1) - (3/2)*abs(abs(t)-2) - (29/4)*abs(abs(t)-4) ...
%     + (29/4)*abs(abs(t)-5) + (7/16)*(abs(abs(t)-2) - abs(abs(t)-3) - 1)^4 ...
%     + 4.5*sin((pi/4)*(abs(abs(t)-3) - abs(abs(t)-4) - 1)) ...
%     - 3*(sqrt(2)/5) * (abs(abs(abs(t)-5) - abs(abs(t)-7)))^(5/2) ...
%     + 6.4*sin(((pi/2) + asin(47/53)) * ((abs(abs(t)-7) - abs(abs(t)-8) + 1)/2) + asin(56/64)) ...
%     + 4.95;
% 
% for i = 1:N_guide
%     Z_guide(i, :) = [subs(x, t, l(i)), subs(y, t, l(i))];
% end
% 
% N_points_filling = 2;

shape = "arbitrary";
N_basis = 3;
order = 1;
[x_ref, u_ref, Tend] = model.generate_trajectory(N_guide, shape, {N_points_filling, N_basis, order, Z_guide});

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
