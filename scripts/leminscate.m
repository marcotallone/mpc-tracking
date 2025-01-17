clear all
clc


Ts = 0.1;

% nf = 4;

% N_intervals = 100;
% N_guide = N_intervals + 1;
% Tend = N_intervals * Ts;
% a = 2;
% theta = 0;
% delta = 2*pi/N_intervals;
% m_theta = delta/Ts; % angulat coefficient of theta(t) = m_theta * t
% % radius = 0.5;

% % syms t real;
% % x = (a*sqrt(2)*cos(t))/(sin(t)^2 + 1);
% % y = (a*sqrt(2)*cos(t)*sin(t))/(sin(t)^2 + 1);
% % z = 0;
% % psi = atan2(diff(y, t), diff(x, t));

% % Guide points and derivatives
% T_guide = linspace(0, Tend, N_guide);
% Z_guide = zeros(N_guide, nf);
% for i = 1:N_guide

%     % t_val = wrapTo2Pi(T_guide(i));ù
%     t_val = theta;
%     Z_guide(i, :) = [subs(x, t, t_val), subs(y, t, t_val), z, subs(psi, t, t_val)];

%     % Z_guide(i, :) = [radius*cos(theta), radius*sin(theta), 0, 0.5*pi + theta];
%     theta = theta + delta;
% end 


% Batman parametric equations
function y = periodic_linear(x, y_min, y_max)
    % Ensure y stays within the range [y_min, y_max] periodically
    P = y_max - y_min;  % Period length
    y = y_min + mod(x - y_min, P);

    if y == 0
        y = y + 1e-3;
    end
end

% function x = batman_x(t)
%     ad = @(a1, a2) abs(abs(a1)-a2);
%     x = (abs(t)/t)*( 0.3*ad(t, 0) +0.2*ad(t, 1) +2.2*ad(t, 2) -2.7*ad(t, 3) -3*ad(t, 5) +3.3*ad(t, 7) +5*sin((pi/4)*(ad(t, 3)*ad(t, 4)+1)) +(5/4)*(ad(t, 4)-ad(t, 5)-1)^3 -5.3*cos( ((pi/2) + asin(47/53)) *( (ad(t, 7) - ad(t, 8) - 1)/2) + 2.8));
% end

% function y = batman_y(t)
%     ad = @(a1, a2) abs(abs(a1)-a2);
%     y = (3/2)*ad(t,1) - (3/2)*ad(t, 2) - (29/4)*ad(t,4) +(29/4)*ad(t,5) +(7/16)*( ad(t,2) - ad(t,3) - 1)^4 + 4.5*sin( (pi/4)*( ad(t, 3) - ad(t, 4) -1)) - (3*sqrt(2)/5)*abs( ad(t,5)-abs(t,7))^(5/2) + 6.4*sin( ( (pi/2) + asin(47/53)) * (( ad(t, 7) - ad(t, 8) - 1)/2) + asin(56/64)) + 4.95;
% end

% syms t real;

% % Version 1
% x = (abs(t)/t)*(...
%      0.3*abs(t) ...
%      + 0.2*abs(abs(t)-1) ...
%      + 2.2*abs(abs(t)-2) ...
%      - 2.7*abs(abs(t)-3) ...
%      - 3*abs(abs(t)-5) ...
%      + 3*abs(abs(t)-7) ...
%      + 5*sin((pi/4)*(abs(abs(t)-3)-abs(abs(t)-4)+1)) ...
%      + (5/4)*(abs(abs(t)-4)-abs(abs(t)-5)-1)^3 ...
%      - 5.3*cos(((pi/2) + asin(47/53))*((abs(abs(t)-7)-abs(abs(t)-8)-1)/2) + 2.8)...
%     );

% y = + (3/2)*abs(abs(t)-1) ...
%     - (3/2)*abs(abs(t)-2) ...
%     - (29/4)*abs(abs(t)-4) ...
%     + (29/4)*abs(abs(t)-5) ...
%     + (7/16)*(abs(abs(t)-2) - abs(abs(t)-3) - 1)^4 ...
%     + 4.5*sin((pi/4)*(abs(abs(t)-3) - abs(abs(t)-4) - 1)) ...
%     - (3*sqrt(2)/5)*(abs(abs(abs(t)-5)-abs(abs(t)-7)))^(5/2) ...
%     + 6.4*sin(((pi/2) + asin(47/53))*((abs(abs(t)-7) - abs(abs(t)-8) + 1)/2) + asin(56/64)) ...
%     + 4.95;


% version 2

% formula:
% X(t)=
% ((abs(t))/(t)) 
% (
%  0.3 abs(t)
% +0.2 abs(abs(t)-1)
% +2.2 abs(abs(t)-2)
% -2.7 abs(abs(t)-3)
% -3 abs(abs(t)-5)
% +3 abs(abs(t)-7)
% +5 sin(((π)/(4)) (abs(abs(t)-3)-abs(abs(t)-4)+1))
% +((5)/(4)) (abs(abs(t)-4)-abs(abs(t)-5)-1)^(3)
% -5.3 cos((((π)/(2))+sin^(-1)(((47)/(53))))*((abs(abs(t)-7)-abs(abs(t)-8)-1)/(2)))+2.8)


% formula 
% Y(t)=((3)/(2)) abs(abs(t)-1)-((3)/(2)) abs(abs(t)-2)-((29)/(4)) abs(abs(t)-4)+((29)/(4)) abs(abs(t)-5)+((7)/(16)) (abs(abs(t)-2)-abs(abs(t)-3)-1)^(4)+4.5 sin(((π)/(4)) (abs(abs(t)-3)-abs(abs(t)-4)-1))-3*((sqrt(2))/(5)) (abs(abs(abs(t)-5)-abs(abs(t)-7)))^(((5)/(2)))+6.4 sin((((π)/(2))+sin^(-1)(((47)/(53))))*((abs(abs(t)-7)-abs(abs(t)-8)+1)/(2))+sin^(-1)(((56)/(64))))+4.95

syms t real;

x = (abs(t)/t) * (0.3*abs(t) + 0.2*abs(abs(t)-1) + 2.2*abs(abs(t)-2) ...
    - 2.7*abs(abs(t)-3) - 3*abs(abs(t)-5) + 3*abs(abs(t)-7) ...
    + 5*sin((pi/4)*(abs(abs(t)-3) - abs(abs(t)-4) + 1)) ...
    + (5/4)*(abs(abs(t)-4) - abs(abs(t)-5) - 1)^3 ...
    - 5.3*cos(((pi/2) + asin(47/53)) * ((abs(abs(t)-7) - abs(abs(t)-8) - 1)/2)) ...
    + 2.8);

y = (3/2)*abs(abs(t)-1) - (3/2)*abs(abs(t)-2) - (29/4)*abs(abs(t)-4) ...
    + (29/4)*abs(abs(t)-5) + (7/16)*(abs(abs(t)-2) - abs(abs(t)-3) - 1)^4 ...
    + 4.5*sin((pi/4)*(abs(abs(t)-3) - abs(abs(t)-4) - 1)) ...
    - 3*(sqrt(2)/5) * (abs(abs(abs(t)-5) - abs(abs(t)-7)))^(5/2) ...
    + 6.4*sin(((pi/2) + asin(47/53)) * ((abs(abs(t)-7) - abs(abs(t)-8) + 1)/2) + asin(56/64)) ...
    + 4.95;


z = 0;
psi = atan2(diff(y, t), diff(x, t));

N_guide = 200;
% l = linspace(0, N_guide*Ts, N_guide);
l = linspace(-8, 8, N_guide);
Z_guide = zeros(N_guide, 4);

for i = 1:N_guide
    % t_val = periodic_linear(T_guide(i), -8, 8);
    % t_val = periodic_linear(l(i), -8, 8);
    t_val = l(i);
    % disp(t_val);
    Z_guide(i, :) = [subs(x, t, t_val), subs(y, t, t_val), z, subs(psi, t, t_val)];
end


% Plotting ─────────────────────────────────────────────────────────────────────

% % Plot Reference / Guide ───────────────────────────────────────────────────────
figure(1);
arrow_length = 0.1;

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



% Multiply periodic references for multiple laps
% n_laps = 2;
% Z_guide = repmat(Z_guide(1:end, :), n_laps, 1);
% Tend = Tend*n_laps;
% N_guide = N_intervals * n_laps + 1;



%%%%%%%%%%
% <<<<<<<<
%%%%%%%%%%
% Wait for figure
pause(1);

% Real trajectory
for i = 1:N_guide
    target = scatter(Z_guide(i, 1), Z_guide(i, 2), 20, 'filled', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'red');
    hold on;
    x_arrow = arrow_length * cos(Z_guide(i, 4));
    y_arrow = arrow_length * sin(Z_guide(i, 4));
    arrow = quiver(Z_guide(i, 1), Z_guide(i, 2), x_arrow, y_arrow, 'AutoScale', 'off', 'Color', 'red');
    legend([guide_points, target],{'Guide Points', 'Target'}, 'Location', 'northwest');
    hold on;

    pause(0.05);
    if i < N_guide
        delete(target);
        delete(arrow);
    end
end