% Trajectory generation Murray method
clear all;
clc;

% Parameters
n = 4; % Number of states
p = 2; % Number of basis functions
q = 0; % Maximum differentiation order

% Initial condition
z0 = [0; 0; 0; 0];

% Times
Tend = 10;
Ts = 0.1;
t = 0:Ts:Tend;
N_guide_points = 24;
guide_points_step = Tend/N_guide_points;
T = 0:guide_points_step:Tend;

% Guide points
Z_guide = zeros(length(T), n);
theta = 0;
delta_theta = 2*pi/N_guide_points;
for i = 1:length(T)
    Z_guide(i, :) = [cos(theta), sin(theta), 0, 0.5 * pi + theta];
    theta = theta + delta_theta;
end

% Matrix-function m(t)
function m = m(t)
    p = 2; % Number of basis functions
    q = 0; % Maximum differentiation order

    % Basis functions
    syms x
    basis = sym('x', [1 p]);
    for k = 1:p
        basis(k) = x^(k-1);
    end

    % Fill the matrix m(t)
    m = zeros(q+1, p);
    for i = 0:q
        for j = 1:p
            derivative = diff(basis(j), x, i);
            m(i+1, j) = subs(derivative, x, t);
        end
    end    
end

% Express z components as linear combination of basis functions
% with coefficients alpha
function z = z(t, alpha)
    p = 2; % Number of basis functions
    n = 4; % Number of states

    % Basis functions
    syms x
    basis = sym('x', [1 p]);
    for k = 1:p
        basis(k) = x^(k-1);
    end
    basis_values = zeros(p, 1);
    for i = 1:p
        basis_values(i) = subs(basis(i), x, t);
    end

    % Reshape alpha
    alpha = reshape(alpha, [], n)';

    % Compute z(t)
    z = alpha * basis_values;

    %// TODO: rewrite z(t) as a symbolic function and evaluate it later...
        % Look at this:
% % Define the symbolic variable
% syms x

% % Define the basis vector
% basis = [1, x];

% % Define the alpha2 matrix
% alpha2 = [
%     0.1822, 0.0818;
%    -6.2117, 0.6212;
%          0,      0;
%     1.5708, 0.6283
% ];

% % Perform the multiplication
% result = alpha2 * basis.';

% % Display the result
% disp(result);

end

% Generate missing trajectory between all guide points
% Z = zeros(length(t), 4);
Z = [];
for i = 1:length(T)-1
    % Matrix M
    mi_0 = m(T(i));
    mi_1 = m(T(i+1));
    Mi0 = kron(eye(n), mi_0);
    Mi1 = kron(eye(n), mi_1);
    M = [Mi0; Mi1];
    
    % Guide points vector
    z_bar = [Z_guide(i, 1:n)'; Z_guide(i+1, 1:n)'];

    % Solve the linear system
    alpha = M\z_bar;

    % Generate missing trajectory between T(i) and T(i+1) with z(t)
    N_filling_points = ceil((T(i+1) - T(i))/Ts);
    ti = linspace(T(i), T(i+1), N_filling_points);
    for j = 1:length(ti) - 1
        % Z((i-1)*length(ti) + j, :) = z(ti(j), alpha)';
        Z = [Z; z(ti(j), alpha)'];
    end
end

% //TODO: recall at the end to wrap all the angles between 0 and 2Pi because it's not done here...
for i = 1:length(Z)
    Z(i, 4) = wrapTo2Pi(Z(i,4));
end

% Model parameters
bx = 2;
by = 2;
bz = 18;
bphi = 110;
kx = -0.5;
ky = -0.5;
kphi = -5;

% Rewrite the other variables w.r.t. the parametrized z(t) and its derivatives

















% % Plot
% figure(1);

% % Reference trajectory
% ref_points = scatter(Z_guide(:, 1), Z_guide(:, 2), 15, 'filled', 'MarkerFaceColor', '#808080');
% hold on;
% arrow_length = 0.03;
% for i = 1:length(T)
%     x_arrow = arrow_length * cos(Z_guide(i, 4));
%     y_arrow = arrow_length * sin(Z_guide(i, 4));
%     quiver(Z_guide(i, 1), Z_guide(i, 2), x_arrow, y_arrow, 'AutoScale', 'off', 'Color', '#808080');
% end
% hold on;

% % Filling the missing trajectory
% ref_points_filling = scatter(Z(:, 1), Z(:, 2), 5, 'filled', 'MarkerFaceColor', '#FF0000');
% hold on;
% for i = 1:length(Z)
%     x_arrow = arrow_length * cos(Z(i, 4));
%     y_arrow = arrow_length * sin(Z(i, 4));
%     quiver(Z(i, 1), Z(i, 2), x_arrow, y_arrow, 'AutoScale', 'off', 'Color', '#FF0000');
% end
% hold on;

% % legend(ref_points,{'Reference trajectory'}, 'Location', 'northwest');
% legend([ref_points, ref_points_filling],{'Reference trajectory', 'Filling trajectory'}, 'Location', 'northwest');

% % Labels
% xlabel('x1'); ylabel('x2');
% grid on;
% axis equal;
% hold on;

