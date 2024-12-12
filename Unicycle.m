% Unicycle class
classdef Unicycle
    properties
        r           % wheel radius
        L           % distance between wheels
        Ts          % sampling time
        n = 3;      % number of states
        m = 2;      % number of inputs
        p = 3;      % number of outputs
        max_omega;  % maximum angular velocity
        min_omega;  % minimum angular velocity
        max_x;      % maximum x position
        min_x;      % minimum x position
        max_y;      % maximum y position
        min_y;      % minimum y position
        max_theta;  % maximum heading angle
        min_theta;  % minimum heading angle
        eps_u;      % input constraints matrix
        f_u;        % input constraints vector
        eps_x;      % state constraints matrix
        f_x;        % state constraints vector
    end
    
    methods
        % Constructor to initialize the unicycle parameters and state
        function obj = Unicycle(r, L, Ts, omega_limits, x_limits, y_limits, theta_limits)
            obj.r = r;
            obj.L = L;
            obj.Ts = Ts;
            obj.max_omega = omega_limits(2);
            obj.min_omega = omega_limits(1);
            obj.max_x = x_limits(2);
            obj.min_x = x_limits(1);
            obj.max_y = y_limits(2);
            obj.min_y = y_limits(1);
            obj.max_theta = theta_limits(2);
            obj.min_theta = theta_limits(1);
            obj.eps_u = [
                 1  0;
                -1  0;
                 0  1;
                 0 -1;
            ];
            obj.f_u = [
                +obj.max_omega;
                -obj.min_omega;
                +obj.max_omega;
                -obj.min_omega;
            ];
            obj.eps_x = [
                 1,  0,  0;
                -1,  0,  0;
                 0,  1,  0;
                 0, -1,  0;
                 0,  0,  1;
                 0,  0, -1
            ];
            obj.f_x = [
                +obj.max_x;
                -obj.min_x;
                +obj.max_y;
                -obj.min_y;
                +obj.max_theta;
                -obj.min_theta
            ];
        end

        % Non-linear continuous dynamics / State transition function
        function dxdt = dynamics(obj, t, x, u)
            % x_pos = x(1);                             % x-position
            % y_pos = x(2);                             % y-position
            theta = x(3);                               % orientation angle
            omega1 = u(1);                              % angular velocity of wheel 1
            omega2 = u(2);                              % angular velocity of wheel 2
            v = obj.r / 2 * (omega1 + omega2);          % linear velocity
            omega = obj.r / obj.L * (omega2 - omega1);  % angular velocity
            dx_posdt = v * cos(theta);                  % x velocity
            dy_posdt = v * sin(theta);                  % y velocity
            dthetadt = omega;                           % angular velocity
            dxdt = [dx_posdt; dy_posdt; dthetadt];
        end

        % Simulation function using ode45
        function x_final = simulate(obj, x0, u, T)
            x0(3) = wrapTo2Pi(x0(3));
            [t, x] = ode45(@(t, x) obj.dynamics(t, x, u), [0, T], x0);
            x_final = (x(end, :))';
            x_final(3) = wrapTo2Pi(x_final(3));
        end
        
        % Linearization function using symbolic toolbox
        function [A_lin, B_lin] = linearize(obj, x_bar, u_bar)
            syms x_pos y_pos theta omega1 omega2 real
            syms r L
            x = [x_pos; y_pos; theta];
            u = [omega1; omega2];
            v = r / 2 * (omega1 + omega2);
            omega = r / L * (omega2 - omega1);
            dx_posdt = v * cos(theta);
            dy_posdt = v * sin(theta);
            dthetadt = omega;
            f = [dx_posdt; dy_posdt; dthetadt];
            A = jacobian(f, x);
            B = jacobian(f, u);
            A_lin = double(subs(A, [x; u; r; L], [x_bar; u_bar; obj.r; obj.L]));
            B_lin = double(subs(B, [x; u; r; L], [x_bar; u_bar; obj.r; obj.L]));
        end

        % Discretization function (from linear system)
        function [A_discrete, B_discrete] = discretize(obj, A, B)

            % Exact sampling
            % continuous_sys = ss(A, B, eye(obj.n), zeros(obj.n, obj.m));
            % discrete_sys = c2d(continuous_sys, obj.Ts);
            % A_discrete = discrete_sys.A;
            % B_discrete = discrete_sys.B;

            % Euler discretization
            A_discrete = eye(obj.n) + A * obj.Ts;
            B_discrete = B * obj.Ts;
        end            

    end
end