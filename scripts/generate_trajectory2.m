clear all
clc


% Helicopter model parameters ──────────────────────────────────────────────────
bx = 2;
by = 2.1;
bz = 18;
bpsi = 111;
kx = -0.5;
ky = -0.4;
kpsi = -5;
ki = 2;
parameters = [bx; by; bz; bpsi; kx; ky; kpsi; ki];
g = 9.81;
Ts = 0.1;
x_constraints = [
    -6, 6;
    -6, 6;
    -6, 6;
    -3, 3;
    -3, 3;
    -2, 2;
     0, 2*pi;
    -25, 25;
];
u_constraints = [
    -1, 1;
    -1, 1;
    -1, 1;
      0, 1;
];
model = Helicopter(parameters, Ts, x_constraints, u_constraints);


% Guide points for flat-outputs ────────────────────────────────────────────────

% Number of flat-outputs
nf = 4;

% Generation parameters
N_intervals = 96;
N_guide = N_intervals + 1;
Tend = N_intervals * Ts;
theta = 0;
delta = 2*pi/N_intervals;
m_theta = delta/Ts; % angulat coefficient of theta(t) = m_theta * t
radius = 0.5;

% Guide points and derivatives
T_guide = linspace(0, Tend, N_guide);
Z_guide = zeros(N_guide, nf);
dZ_guide = zeros(N_guide, nf);
ddZ_guide = zeros(N_guide, nf);
for i = 1:N_guide
    Z_guide(i, :) = [radius*cos(theta), radius*sin(theta), 0, 0.5*pi + theta];
    dZ_guide(i, :) = [-radius*m_theta*sin(theta), radius*m_theta*cos(theta), 0, m_theta];
    ddZ_guide(i, :) = [-radius*(m_theta^2)*cos(theta), -radius*(m_theta^2)*sin(theta), 0, 0];
    theta = theta + delta;
end 

% Analytical definition
syms t real;
ftheta = @(t) m_theta*t;
z = [radius*cos(ftheta(t)), radius*sin(ftheta(t)), 0, 0.5*pi + ftheta(t)];
dz = diff(z, t);
ddz = diff(dz, t);

% Reference trajectory filling with remaining state components ─────────────────

% Anlytical definition of remaining states and inputs
dxb = cos(z(4))*dz(1) + sin(z(4))*dz(2);
dyb = -sin(z(4))*dz(1) + cos(z(4))*dz(2);
dzb = dz(3);
dpsi = dz(4);
ux = (cos(z(4))*(ddz(1) - kx*dz(1)) + sin(z(4))*(ddz(2) - kx*dz(2))/bx);
uy = (cos(z(4))*(ddz(2) - ky*dz(2)) + sin(z(4))*(-ddz(1) + ky*dz(1))/by);
uz = (ddz(3) + g)/bz;
upsi = (ddz(4) - kpsi*dz(4))/bpsi;

x = [z(1), z(2), z(3), dxb, dyb, dzb, z(4), dpsi];
u = [ux, uy, uz, upsi];

% Reference trajectory
x_ref = [];
u_ref = [];
for i = 1:N_guide - 1
    x_t = double(subs(x, t, T_guide(i)));
    u_t = double(subs(u, t, T_guide(i)));

    x_ref = [x_ref; x_t];
    u_ref = [u_ref; u_t];
end


% Plotting ─────────────────────────────────────────────────────────────────────

% Plot z(t) and guide points
figure(2);

ts = linspace(0, Tend, 1000);
zs = zeros(length(ts), nf);
for i = 1:length(ts)
    zs(i, :) = double(subs(z, t, ts(i)));
end
    
subplot(4, 1, 1);
plot(ts, zs(:, 1), 'blue');
hold on;
scatter(T_guide, Z_guide(:, 1), 15, 'filled', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'blue');
hold on;
subplot(4, 1, 2);
plot(ts, zs(:, 2), 'red');
hold on;
scatter(T_guide, Z_guide(:, 2), 15, 'filled', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'red');
hold on;
subplot(4, 1, 3);
plot(ts, zs(:, 3), 'magenta');
hold on;
scatter(T_guide, Z_guide(:, 3), 15, 'filled', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'magenta');
hold on;
subplot(4, 1, 4);
plot(ts, zs(:, 4), 'green');
hold on;
scatter(T_guide, Z_guide(:, 4), 15, 'filled', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'green');
hold on;

% Plot inputs
figure(3);
us = zeros(length(ts), 4);
for i = 1:length(ts)
    us(i, :) = double(subs(u, t, ts(i)));
end

subplot(4, 1, 1);
plot(ts, us(:, 1), 'blue');
hold on;
scatter(T_guide(1:end-1), u_ref(:, 1), 15, 'filled', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'blue');
hold on;
subplot(4, 1, 2);
plot(ts, us(:, 2), 'red');
hold on;
scatter(T_guide(1:end-1), u_ref(:, 2), 15, 'filled', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'red');
hold on;
subplot(4, 1, 3);
plot(ts, us(:, 3), 'magenta');
hold on;
scatter(T_guide(1:end-1), u_ref(:, 3), 15, 'filled', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'magenta');
hold on;
subplot(4, 1, 4);
plot(ts, us(:, 4), 'green');
hold on;
scatter(T_guide(1:end-1), u_ref(:, 4), 15, 'filled', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'green');
hold on;


% % Plot Reference / Guide ───────────────────────────────────────────────────────
figure(1);
arrow_length = 0.01;

% Guide points
guide_points = scatter(Z_guide(:, 1), Z_guide(:, 2), 15, 'filled', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', '#808080');
hold on;
for i = 1:length(Z_guide)
    x_arrow = arrow_length * cos(Z_guide(i, 4));
    y_arrow = arrow_length * sin(Z_guide(i, 4));
    quiver(Z_guide(i, 1), Z_guide(i, 2), x_arrow, y_arrow, 'AutoScale', 'off', 'Color', '#808080');
end
hold on;

% Labels
legend(guide_points,{'Guide Points'}, 'Location', 'northwest');
xlabel('x1'); ylabel('x2');
grid on;
axis equal;
hold on;



% Simulate system to test trajectory ───────────────────────────────────────────
% x0 = x_ref(1, :)';
% for i = 1:N_guide - 1
% 
%     x_sim = model.simulate(x0, u_ref(i, :)', Ts);
%     % x0 = x_ref(i + 1, :)'; % reposition to next reference
%     x0 = x_sim'; % continuous simulation
% 
%     x_line = plot(x_sim(1), x_sim(2), 'blue', 'LineWidth', 1);
%     x_line.Color(4) = 0.5; % line transparency 50%
%     hold on;
%     x_points = scatter(x_sim(1), x_sim(2), 5, 'blue', 'filled');
%     hold on;
%     quiver(x_sim(1), x_sim(2), arrow_length * cos(x_sim(7)), arrow_length * sin(x_sim(7)), 'AutoScale', 'off', 'Color', 'blue');
%     hold on;
%     legend([guide_points, x_points],{'Guide Points', 'Real trajectory'}, 'Location', 'northwest');
%     hold on;
% 
%     pause(0.05);
%     if i < length(u_ref)
%         delete(x_line);
%     end
% end



% MPC ──────────────────────────────────────────────────────────────────────────

% Multiply references for multiple laps
n_laps = 2;
Z_guide = repmat(Z_guide(1:end-1, :), n_laps, 1);
x_ref = repmat(x_ref, n_laps, 1);
u_ref = repmat(u_ref, n_laps, 1);
Tend = Tend*n_laps;

% MPC parameters
x0 = zeros(model.n,1);
N  = 18;
Q = diag([50, 50, 5, 10, 3, 3, 1, 2]);
R = diag([2, 2, 2, 2]);
preview = 1;
formulation = 0;
debug = 0;
times = 0:Ts:Tend;
Nsteps = length(times) - (N+1);

% Optimization 
mpc = MPC(model, x0, Tend, N, Q, R, x_ref, u_ref, preview, formulation, debug);
[x_mpc, u_mpc] = mpc.optimize();

%%%%%%%%%%
% <<<<<<<<
%%%%%%%%%%

% Plot
pause(1);
for i = 1:Nsteps
    x_line = plot(x_mpc(1:i, 1), x_mpc(1:i, 2), 'blue', 'LineWidth', 1);
    x_line.Color(4) = 0.5; % line transparency 50%
    hold on;
    x_points = scatter(x_mpc(1:i, 1), x_mpc(1:i, 2), 5, 'blue', 'filled');
    hold on;
    quiver(x_mpc(1:i, 1), x_mpc(1:i, 2), arrow_length * cos(x_mpc(1:i, 7)), arrow_length * sin(x_mpc(1:i, 7)), 'AutoScale', 'off', 'Color', 'blue');
    hold on;
    target = scatter(Z_guide(i, 1), Z_guide(i, 2), 20, 'filled', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'red');
    hold on;
    legend([guide_points, x_points, target],{'Guide Points', 'Real trajectory', 'Target'}, 'Location', 'northwest');
    hold on;

    pause(0.05);
    if i < Nsteps
        delete(x_line);
        delete(target);
    end
end