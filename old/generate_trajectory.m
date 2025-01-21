clear all
clc


% Helicopter model parameters ──────────────────────────────────────────────────
bx = 2;
by = 2;
bz = 18;
bpsi = 111;
kx = -0.5;
ky = -0.5;
kpsi = -5;
ki = 2;
parameters = [bx; by; bz; bpsi; kx; ky; kpsi; ki];
g = 9.81;
Ts = 0.1;
x_constraints = [
    -10, 10;
    -10, 10;
    -10, 10;
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


% Guide points for reference trajectory ────────────────────────────────────────
nf = 4; % Number of flat-outputs

% Generation parameters
N_intervals = 24;
N_guide = N_intervals + 1; % Number of guide points (last == first)
N_points_filling = 5; % Number of points between guide points
Tend = N_intervals * N_points_filling * Ts;

% Guide
T_guide = linspace(0, Tend, N_guide);
Z_guide = zeros(N_guide, nf);

theta = 0;
delta = 2*pi/N_intervals;
radius = 0.5;
for i = 1:N_guide
    Z_guide(i, :) = [radius*cos(theta), radius*sin(theta), 0, 0.5*pi + theta];
    theta = theta + delta;
end 


% % Reference trajectory filling ─────────────────────────────────────────────────
% N_basis = 2; % Number of basis functions
% order = 2;   % Maximum differentiation order

% % Guide vectors for first and second derivatives
% dZ_guide = [zeros(1, nf)];
% ddZ_guide = [zeros(1, nf)];

% % Basis functions
% syms x
% basis = sym('x', [1 N_basis]);
% for k = 1:N_basis
%     basis(k) = x^(k-1);
% end

% % Small m(t) matrix function
% function m = m(t, basis, order)
%     syms x
%     m = zeros(order+1, length(basis));
%     for i = 0:order
%         for j = 1:length(basis)
%             derivative = diff(basis(j), x, i);
%             m(i+1, j) = subs(derivative, x, t);
%         end
%     end
% end

% % Trajectory filling
% Z_ref = [];
% x_ref = [];
% u_ref = [];

% x_simul2 = [];
% x_simul3 = [];

% alpha_tensor = zeros(nf, N_basis, length(T_guide) - 1);
% for i = 1:length(T_guide) - 1

%     % Time interval
%     t0 = T_guide(i);
%     t1 = T_guide(i+1);

%     % Matrix M
%     m0 = m(t0, basis, 0); % first go, only order 0 (no derivatives)
%     m1 = m(t1, basis, 0);
%     M0 = kron(eye(nf), m0);
%     M1 = kron(eye(nf), m1);
%     M = [M0; M1];

%     % Guide points in the interval
%     z_bar = [
%         Z_guide(i, 1:nf)';
%         Z_guide(i+1, 1:nf)';
%     ];

%     % Solve the system and reshape coefficients
%     alpha = M\z_bar;
%     alpha = reshape(alpha, [], nf)';
%     alpha_tensor(:, :, i) = alpha;

%     % Reference function z(t) and its derivatives
%     z = alpha*basis.';
%     dz = diff(z, x);
%     ddz = diff(dz, x);

%     % Generate missing points
%     N_filling = ceil((t1 - t0)/Ts) + 1;
%     T_filling = linspace(t0, t1, N_filling);
%     for j = 1:length(T_filling) - 1

%         % Evaluate z(t) and derivatives
%         z_t = double(subs(z, x, T_filling(j)))';
%         z_t(nf) = wrapTo2Pi(z_t(nf)); % wrap the angle state
%         dz_t = double(subs(dz, x, T_filling(j)))';
%         ddz_t = double(subs(ddz, x, T_filling(j)))';

%         % Compute remaining states and inputs
%         dxb = cos(z_t(4))*dz_t(1) + sin(z_t(4))*dz_t(2);
%         dyb = -sin(z_t(4))*dz_t(1) + cos(z_t(4))*dz_t(2);
%         dzb = dz_t(3);
%         dpsi = dz_t(4);

%         x_t = [z_t(1), z_t(2), z_t(3), dxb, dyb, dzb, z_t(4), dpsi];

%         ux = (cos(z_t(4))*(ddz_t(1) - kx*dz_t(1)) + sin(z_t(4))*(ddz_t(2) - kx*dz_t(2)))/bx;
%         uy = (cos(z_t(4))*(ddz_t(2) - ky*dz_t(2)) + sin(z_t(4))*(-ddz_t(1) + ky*dz_t(1)))/by;
%         uz = (ddz_t(3) + g)/bz;
%         upsi = (ddz_t(4) - kpsi*dz_t(4))/bpsi;

%         u_t = [ux, uy, uz, upsi];

%         % Store the reference flat-output, states and inputs
%         Z_ref = [Z_ref; z_t];
%         x_ref = [x_ref; x_t];
%         u_ref = [u_ref; u_t];

%         % Simulate the real system
%         if j == 1
%             x02 = x_t';
%             if i == 1
%                 x03 = x_t';
%             end
%         end
%         x_sim2 = model.simulate(x02, u_t', Ts);
%         x02 = x_sim2;
%         x_simul2 = [x_simul2; x_sim2'];
%         x_sim3 = model.simulate(x03, u_t', Ts);
%         x03 = x_sim3;
%         x_simul3 = [x_simul3; x_sim3'];

%     end

%     % % NOOOOOOOOOOOOOO
%     % % Simulate the system of last point
%     % x_sim = model.simulate(x_t', u_t', Ts);
%     % 
%     % % Reverse-obtain flat-outputs derivatives from state
%     % dz1 = cos(x_sim(7))*x_sim(4) - sin(x_sim(7))*x_sim(5);
%     % dz2 = sin(x_sim(7))*x_sim(4) + cos(x_sim(7))*x_sim(5);
%     % dz3 = x_sim(6);
%     % dz4 = x_sim(8);
%     % 
%     % u_ = bx*u_t(1);
%     % v_ = by*u_t(2);
%     % c_ = cos(x_sim(7));
%     % s_ = sin(x_sim(7));
%     % a_ = -kx*dz1;
%     % b_ = -kx*dz2;
%     % f_ = -ky*dz2;
%     % g_ = ky*dz1;
%     % alpha_ = u_ + c_*a_ + s_*b_;
%     % beta_ = v_ + c_*f_ - s_*g_;
%     % 
%     % % Solve the linear system:
%     % A_ = [c_, s_; -s_, c_];
%     % B_ = [alpha_; beta_];
%     % ddz_ = A_\B_;
%     % ddz1 = ddz_(1);
%     % ddz2 = ddz_(2);
%     % 
%     % ddz3 = bz*u_t(3) + g;
%     % ddz4 = bpsi*u_t(4) + kpsi*dz4;
%     % 
%     % % Store the reference derivatives
%     % dZ_guide = [dZ_guide; [dz1, dz2, dz3, dz4]];
%     % ddZ_guide = [ddZ_guide; [ddz1, ddz2, ddz3, ddz4]];
%     % 
%     % % Re-compute z(t) with new derivative information
%     % m0 = m(t0, basis, order);
%     % m1 = m(t1, basis, order);
%     % M0 = kron(eye(nf), m0);
%     % M1 = kron(eye(nf), m1);
%     % M = [M0; M1];
%     % z_bar = [
%     %     Z_guide(i, 1:nf)';
%     %     dZ_guide(i, 1:nf)';
%     %     ddZ_guide(i, 1:nf)';
%     %     Z_guide(i+1, 1:nf)';
%     %     dZ_guide(i+1, 1:nf)';
%     %     ddZ_guide(i+1, 1:nf)';
%     % ];
%     % alpha = M\z_bar;
%     % alpha = reshape(alpha, [], nf)';
%     % alpha_tensor(:, :, i) = alpha;
%     % z = alpha*basis.';
%     % dz = diff(z, x);
%     % ddz = diff(dz, x);
%     % for j = 1:length(T_filling) - 1
%     % 
%     %     % Evaluate z(t) and derivatives
%     %     z_t = double(subs(z, x, T_filling(j)))';
%     %     z_t(nf) = wrapTo2Pi(z_t(nf)); % wrap the angle state
%     %     dz_t = double(subs(dz, x, T_filling(j)))';
%     %     ddz_t = double(subs(ddz, x, T_filling(j)))';
%     % 
%     %     % Compute remaining states and inputs
%     %     dxb = cos(z_t(4))*dz_t(1) + sin(z_t(4))*dz_t(2);
%     %     dyb = -sin(z_t(4))*dz_t(1) + cos(z_t(4))*dz_t(2);
%     %     dzb = dz_t(3);
%     %     dpsi = dz_t(4);
%     % 
%     %     x_t = [z_t(1), z_t(2), z_t(3), dxb, dyb, dzb, z_t(4), dpsi];
%     % 
%     %     ux = (cos(z_t(4))*(ddz_t(1) - kx*dz_t(1)) + sin(z_t(4))*(ddz_t(2) - kx*dz_t(2)))/bx;
%     %     uy = (cos(z_t(4))*(ddz_t(2) - ky*dz_t(2)) + sin(z_t(4))*(-ddz_t(1) + ky*dz_t(1)))/by;
%     %     uz = (ddz_t(3) + g)/bz;
%     %     upsi = (ddz_t(4) - kpsi*dz_t(4))/bpsi;
%     % 
%     %     u_t = [ux, uy, uz, upsi];
%     % 
%     %     % Store the reference flat-output, states and inputs
%     %     Z_ref = [Z_ref; z_t];
%     %     x_ref = [x_ref; x_t];
%     %     u_ref = [u_ref; u_t];   
%     % 
%     %     % Simulate the real system
%     %     if j == 1
%     %         x02 = x_t';
%     %         if i == 1
%     %             x03 = x_t';
%     %         end
%     %     end
%     %     x_sim2 = model.simulate(x02, u_t', Ts);
%     %     x02 = x_sim2;
%     %     x_simul2 = [x_simul2; x_sim2'];
%     %     x_sim3 = model.simulate(x03, u_t', Ts);
%     %     x03 = x_sim3;
%     %     x_simul3 = [x_simul3; x_sim3'];
%     % end

% end













% SECOND APPROACH ──────────────────────────────────────────────────────────────


N_basis = 7;
order   = 2;  

% Basis functions
syms x
basis = sym('x', [1 N_basis]);
for k = 1:N_basis
    basis(k) = x^(k-1);
end

% Small m(t) matrix function
function m = m(t, basis, order)
    syms x
    m = zeros(order+1, length(basis));
    for i = 0:order
        for j = 1:length(basis)
            derivative = diff(basis(j), x, i);
            m(i+1, j) = subs(derivative, x, t);
        end
    end
end

% Trajectory filling
Z_ref = [];
x_ref = [];
u_ref = [];
% x0 = zeros(1,8);
dZ_guide = [zeros(1, nf)];
ddZ_guide = [zeros(1, nf)];
x_simul = [];
x_simul2 = [];
alpha_tensor = zeros(nf, N_basis, length(T_guide) - 1);


for i = 1:10 %length(T_guide) - 1

    disp(i);

    % Time interval and filling points
    t0 = T_guide(i);
    t1 = T_guide(i+1);
    N_filling = ceil((t1 - t0)/Ts) + 1;
    T_filling = linspace(t0, t1, N_filling);

    % Solve linear system without derivative knowledge
    m0 = m(t0, basis, 0);
    m1 = m(t1, basis, 0);
    M0 = kron(eye(nf), m0);
    M1 = kron(eye(nf), m1);
    M = [M0; M1];
    z_bar = [
        Z_guide(i, 1:nf)';
        Z_guide(i+1, 1:nf)';
    ];
    alpha = M\z_bar;
    alpha = reshape(alpha, [], nf)';

    % Reference function z(t) and its derivatives
    z = alpha*basis.';
    dz = diff(z, x);
    ddz = diff(dz, x);

    % Simulate filling points to get derivatives information
    for j = 1:length(T_filling) - 1

        % Evaluate z(t) and derivatives
        z_t = double(subs(z, x, T_filling(j)))';
        z_t(nf) = wrapTo2Pi(z_t(nf)); % wrap the angle state
        dz_t = double(subs(dz, x, T_filling(j)))';
        ddz_t = double(subs(ddz, x, T_filling(j)))';

        % Compute remaining states and inputs
        dxb = cos(z_t(4))*dz_t(1) + sin(z_t(4))*dz_t(2);
        dyb = -sin(z_t(4))*dz_t(1) + cos(z_t(4))*dz_t(2);
        dzb = dz_t(3);
        dpsi = dz_t(4);
        ux = (cos(z_t(4))*(ddz_t(1) - kx*dz_t(1)) + sin(z_t(4))*(ddz_t(2) - kx*dz_t(2)))/bx;
        uy = (cos(z_t(4))*(ddz_t(2) - ky*dz_t(2)) + sin(z_t(4))*(-ddz_t(1) + ky*dz_t(1)))/by;
        uz = (ddz_t(3) + g)/bz;
        upsi = (ddz_t(4) - kpsi*dz_t(4))/bpsi;

        x_t = [z_t(1), z_t(2), z_t(3), dxb, dyb, dzb, z_t(4), dpsi];
        u_t = [ux, uy, uz, upsi];

        % Simulate the real system
        if i == 1 && j == 1
            x0 = x_t';
        end
        x0_tmp = x0;
        x_sim = model.simulate(x0, u_t', Ts);
        x0 = x_sim;
        x_simul = [x_simul; x_sim']; 

    end

    % Obtain derivatives from last simulation value
    dz1 = cos(x_sim(7))*x_sim(4) + sin(x_sim(7))*x_sim(5);
    dz2 = -sin(x_sim(7))*x_sim(4) + cos(x_sim(7))*x_sim(5);
    dz3 = x_sim(6);
    dz4 = x_sim(8);
    u_ = bx*u_t(1);
    v_ = by*u_t(2);
    c_ = cos(x_sim(7));
    s_ = sin(x_sim(7));
    a_ = -kx*dz1;
    b_ = -kx*dz2;
    f_ = -ky*dz2;
    g_ = ky*dz1;
    alpha_ = u_ + c_*a_ + s_*b_;
    beta_ = v_ + c_*f_ - s_*g_;
    A_ = [c_, s_; -s_, c_];
    B_ = [alpha_; beta_];
    ddz_ = A_\B_;
    ddz1 = ddz_(1);
    ddz2 = ddz_(2);
    ddz3 = bz*u_t(3) - g;
    ddz4 = bpsi*u_t(4) + kpsi*x_sim(8);

    % Update next Z_guide with last simulated point and add derivatives
    Z_guide(i+1,:) = [x_sim(1), x_sim(2), x_sim(3), x_sim(7)];
    dZ_guide = [dZ_guide; [dz1, dz2, dz3, dz4]];
    ddZ_guide = [ddZ_guide; [ddz1, ddz2, ddz3, ddz4]];

    % Re-solve the linear system with derivative knowledge
    m0 = m(t0, basis, order);
    m1 = m(t1, basis, order);
    M0 = kron(eye(nf), m0);
    M1 = kron(eye(nf), m1);
    M = [M0; M1];
    z_bar = [
        Z_guide(i, 1:nf)';
        dZ_guide(i, 1:nf)';
        ddZ_guide(i, 1:nf)';
        Z_guide(i+1, 1:nf)';
        dZ_guide(i+1, 1:nf)';
        ddZ_guide(i+1, 1:nf)';
    ];
    alpha = M\z_bar;
    alp = alpha;
    alpha = reshape(alpha, [], nf)';
    alpha_tensor(:, :, i) = alpha;

    % Reference function z(t) and its derivatives
    z = alpha*basis.';
    dz = diff(z, x);
    ddz = diff(dz, x);

    % Generate missing points
    for j = 1:length(T_filling) - 1

        % Evaluate z(t) and derivatives
        z_t = double(subs(z, x, T_filling(j)))';
        z_t(nf) = wrapTo2Pi(z_t(nf)); % wrap the angle state
        dz_t = double(subs(dz, x, T_filling(j)))';
        ddz_t = double(subs(ddz, x, T_filling(j)))';

        % Compute remaining states and inputs
        dxb = cos(z_t(4))*dz_t(1) + sin(z_t(4))*dz_t(2);
        dyb = -sin(z_t(4))*dz_t(1) + cos(z_t(4))*dz_t(2);
        dzb = dz_t(3);
        dpsi = dz_t(4);
        ux = (cos(z_t(4))*(ddz_t(1) - kx*dz_t(1)) + sin(z_t(4))*(ddz_t(2) - kx*dz_t(2))/bx);
        uy = (cos(z_t(4))*(ddz_t(2) - ky*dz_t(2)) + sin(z_t(4))*(-ddz_t(1) + ky*dz_t(1))/by);
        uz = (ddz_t(3) + g)/bz;
        upsi = (ddz_t(4) - kpsi*dz_t(4))/bpsi;

        x_t = [z_t(1), z_t(2), z_t(3), dxb, dyb, dzb, z_t(4), dpsi];
        u_t = [ux, uy, uz, upsi];


        % Simulate the real system
        x0 = x0_tmp;
        x_sim = model.simulate(x0, u_t', Ts);
        x0 = x_sim;
        x_simul2 = [x_simul2; x_sim']; 

        % Store the reference flat-output, states and inputs
        % Z_ref = [Z_ref; z_t];
        % x_ref = [x_ref; x_t];
        Z_ref = [Z_ref; [x_sim(1), x_sim(2), x_sim(3), x_sim(7)]];
        x_ref = [x_ref; x_sim'];
        u_ref = [u_ref; u_t];

    end


end



% Plotting ─────────────────────────────────────────────────────────────────────

% Plot the z(t) fnctions for every interval
figure(2);
for i = 1:length(T_guide) - 1

    % Time interval
    t0 = T_guide(i);
    t1 = T_guide(i+1);

    % Reconstruct z(t) of this interval
    alpha = alpha_tensor(:, :, i);
    z = alpha*basis.';

    % Plot the z(t) function
    ts = linspace(t0, t1, 100);
    zs = double(subs(z, x, ts));
    
    % Subplots for x = z(t)(1), y = z(t)(2), psi = z(t)(4)
    subplot(4, 1, 1);
    plot(ts, zs(1, :), 'blue');
    hold on;
    subplot(4, 1, 2);
    plot(ts, zs(2, :), 'red');
    hold on;
    subplot(4, 1, 3);
    plot(ts, zs(3, :), 'magenta');
    hold on;
    subplot(4, 1, 4);
    plot(ts, zs(4, :), 'green');
    hold on;

end



% % Plot Reference / Guide ───────────────────────────────────────────────────────
figure(1);
arrow_length = 0.01;

% Guide points
guide_points = scatter(Z_guide(:, 1), Z_guide(:, 2), 15, 'filled', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', '#808080');
hold on;

% Reference points
reference_points = scatter(Z_ref(:, 1), Z_ref(:, 2), 5, 'filled', 'MarkerFaceColor', '#808080');
hold on;
for i = 1:length(Z_ref)
    x_arrow = arrow_length * cos(Z_ref(i, 4));
    y_arrow = arrow_length * sin(Z_ref(i, 4));
    quiver(Z_ref(i, 1), Z_ref(i, 2), x_arrow, y_arrow, 'AutoScale', 'off', 'Color', '#808080');
end
hold on;

% Labels
legend([guide_points, reference_points],{'Guide Points', 'Reference trajectory'}, 'Location', 'northwest');
xlabel('x1'); ylabel('x2');
grid on;
axis equal;
hold on;



% % Simulate system to test trajectory ───────────────────────────────────────────
x0 = x_ref(1, :)';
% x = model.simulate(x0, u_ref, Tend);
for i = 1:length(u_ref)
    x_sim = model.simulate(x0, u_ref(i, :)', Ts);
    x0 = x_sim;

    x_line = plot(x_sim(1), x_sim(2), 'green', 'LineWidth', 1);
    x_line.Color(4) = 0.5; % line transparency 50%
    hold on;
    x_points = scatter(x_sim(1), x_sim(2), 5, 'green', 'filled');
    hold on;
    quiver(x_sim(1), x_sim(2), arrow_length * cos(x_sim(7)), arrow_length * sin(x_sim(7)), 'AutoScale', 'off', 'Color', 'green');
    hold on;

    % x_line2 = plot(x_simul2(i, 1), x_simul2(i, 2), 'red', 'LineWidth', 1);
    % x_line2.Color(4) = 0.5; % line transparency 50%
    % hold on;
    % x_points2 = scatter(x_simul2(i, 1), x_simul2(i, 2), 5, 'red', 'filled');
    % hold on;
    % quiver(x_simul2(i, 1), x_simul2(i, 2), arrow_length * cos(x_simul2(i, 7)), arrow_length * sin(x_simul2(i, 7)), 'AutoScale', 'off', 'Color', 'red');
    % hold on;

    % disp('2 ');
    % disp(x_simul2(i, :));

    % x_line3 = plot(x_simul3(i, 1), x_simul3(i, 2), 'blue', 'LineWidth', 1);
    % x_line3.Color(4) = 0.5; % line transparency 50%
    % hold on;
    % x_points3 = scatter(x_simul3(i, 1), x_simul3(i, 2), 5, 'blue', 'filled');
    % hold on;
    % quiver(x_simul3(i, 1), x_simul3(i, 2), arrow_length * cos(x_simul3(i, 7)), arrow_length * sin(x_simul3(i, 7)), 'AutoScale', 'off', 'Color', 'blue');
    % hold on;

    % disp('3');
    % disp(x_simul3(i, :));

    target = scatter(Z_ref(i, 1), Z_ref(i, 2), 20, 'filled', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'red');
    hold on;
    legend([guide_points, reference_points, x_points, x_points2, x_points3, target],{'Guide Points', 'Reference trajectory', 'Real trajectory', 'Real Trajectory 2', 'Real Trajectory 3', 'Target'}, 'Location', 'northwest');
    hold on;

    pause(0.05);
    if i < length(u_ref)
        delete(x_line);
        delete(x_line2);
        delete(x_line3);
        delete(target);
    end
end

