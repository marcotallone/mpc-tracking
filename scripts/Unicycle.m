% ┌─────────────────────────────────────────────────────────────────────────┐ 
% │                             Unicycle Class                              │ 
% └─────────────────────────────────────────────────────────────────────────┘ 
% Class modelling unicycle non-linear dynamics
%
% Creation
%   Syntax
%     obj = Unicycle(r, L, Ts, x_constraints, u_constraints)
%
%   Input Arguments
%     r - Wheel radius
%         real scalar
%     L - Distance between wheels
%         real scalar
%     Ts - Sampling time
%         real scalar
%     x_constraints - State constraints
%         real matrix
%     u_constraints - Input constraints
%         real vector
%
% Properties
%   n - Number of states
%       real scalar
%   m - Number of inputs
%       real scalar
%   p - Number of outputs
%       real scalar
%   Ts - Sampling time
%       real scalar
%   r - Wheel radius
%       real scalar
%   L - Distance between wheels
%       real scalar
%   min_omega - Minimum angular velocity
%       real scalar
%   max_omega - Maximum angular velocity
%       real scalar
%   min_x - Minimum x position
%       real scalar
%   max_x - Maximum x position
%       real scalar
%   min_y - Minimum y position
%       real scalar
%   max_y - Maximum y position
%       real scalar
%   min_theta - Minimum heading angle
%       real scalar
%   max_theta - Maximum heading angle
%       real scalar
%   eps_x - State constraints matrix
%       real matrix
%   f_x - State constraints vector
%       real vector
%   eps_u - Input constraints matrix
%       real matrix
%   f_u - Input constraints vector
%       real vector
%
% Methods
%   dynamics - State transition function: non-linear continuous dynamics
%   simulate - Simulation function
%   output - Output transformation function
%   linearize - Linearization function
%   fix_angles - Fix reference angles function
%   EKF_step - Extended Kalman Filter (EKF) state estimation step
%
% Examples
%   unicycle = Unicycle(0.1, 0.5, 0.1, [-2 2; -2 2; 0 2*pi], [-30, 30]);
%   x = unicycle.simulate([0; 0; 7*pi/4], [0; 0], 0.1);

classdef Unicycle < DynamicalSystem
    properties

        % Model parameters
        n = 3;      % number of states
        m = 2;      % number of inputs
        p = 2;      % number of outputs
        Ts;         % sampling time
        r;          % wheel radius
        L;          % distance between wheels

        % Constraints
        min_omega;  % minimum angular velocity
        max_omega;  % maximum angular velocity
        min_x;      % minimum x position
        max_x;      % maximum x position
        min_y;      % minimum y position
        max_y;      % maximum y position
        min_theta;  % minimum heading angle
        max_theta;  % maximum heading angle
        eps_x;      % state constraints matrix
        f_x;        % state constraints vector
        eps_u;      % input constraints matrix
        f_u;        % input constraints vector

        % Symbolic properties
        sym_r;      % symbolic wheel radius
        sym_L;      % symbolic distance between wheels
        sym_x;      % symbolic state vector
        sym_u;      % symbolic input vector
        sym_f;      % symbolic state transition function
        sym_g;      % symbolic output measurement
        sym_A;      % symbolic state matrix         (n x n)
        sym_B;      % symbolic input matrix         (n x m)
        sym_C;      % symbolic output matrix        (p x n)

        % EKF-related properties
        P           % State covariance matrix       (n x n)
        Q_tilde     % Process noise covariance      (n x n)
        R_tilde     % Measurement noise covariance  (p x p)
    end
    
    methods
        % Constructor to initialize the unicycle parameters and state
        function obj = Unicycle(r, L, Ts, x_constraints, u_constraints, P0, Q_tilde, R_tilde)
            obj.r = r;
            obj.L = L;
            obj.Ts = Ts;
            obj.min_x = x_constraints(1, 1);
            obj.max_x = x_constraints(1, 2);
            obj.min_y = x_constraints(2, 1);
            obj.max_y = x_constraints(2, 2);
            obj.min_theta = x_constraints(3, 1);
            obj.max_theta = x_constraints(3, 2);
            obj.min_omega = u_constraints(1);
            obj.max_omega = u_constraints(2);
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

            % Initialize symbolic properties
            syms x_pos y_pos theta omega1 omega2 real
            syms sym_r sym_L
            obj.sym_r = sym_r;
            obj.sym_L = sym_L;
            obj.sym_x = [x_pos; y_pos; theta];
            obj.sym_u = [omega1; omega2];
            v = sym_r / 2 * (omega1 + omega2);
            omega = sym_r / sym_L * (omega2 - omega1);
            dx_posdt = v * cos(theta);
            dy_posdt = v * sin(theta);
            dthetadt = omega;
            obj.sym_f = [dx_posdt; dy_posdt; dthetadt];
            obj.sym_g = [x_pos; y_pos];
            obj.sym_A = jacobian(obj.sym_f, obj.sym_x);
            obj.sym_B = jacobian(obj.sym_f, obj.sym_u);
            obj.sym_C = jacobian(obj.sym_g, obj.sym_x);

            % Initialize EKF properties
            obj.P = P0;
            obj.Q_tilde = Q_tilde;
            obj.R_tilde = R_tilde;
        end

        % State transition function: non-linear continuous dynamics
        function dxdt = dynamics(obj, t, x, u)
            % dynamics
            %   State transition function for unicycle model:
            %
            %   dx_posdt = v * cos(theta)
            %   dy_posdt = v * sin(theta)
            %   dthetadt = omega
            %
            % Syntax
            %   dxdt = obj.dynamics(t, x, u)
            %
            % Input Arguments
            %   t - Time
            %       real scalar
            %   x - State vector
            %       real vector
            %   u - Input vector
            %       real vector
            %
            % Output Arguments
            %   dxdt - State transition vector
            %       real vector

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

        % Simulation function
        function x_final = simulate(obj, x0, u, T)
            % simulate
            %   Simulate the unicycle model for a given time period T
            %
            % Syntax
            %   x_final = obj.simulate(x0, u, T)
            %
            % Input Arguments
            %   x0 - Initial state
            %       real vector
            %   u - Input vector
            %       real vector
            %   T - Simulation time
            %       real scalar
            %
            % Output Arguments
            %   x_final - Final state
            %       real vector

            x0(3) = wrapTo2Pi(x0(3));
            [t, x] = ode45(@(t, x) obj.dynamics(t, x, u), [0, T], x0);
            x_final = (x(end, :))';
            x_final(3) = wrapTo2Pi(x_final(3));
        end

        % Output transformation function
        function y = output(obj, x, u)
            % output
            %   Output transformation function for the unicycle model
            %
            % Syntax
            %   y = obj.output(x, u)
            %
            % Input Arguments
            %   x - State vector
            %       real vector
            %   u - Input vector
            %       real vector
            %
            % Output Arguments
            %   y - Output vector
            %       real vector

            y = [x(1); x(2)];
        end
        
        % Linearization function
        function [A_lin, B_lin] = linearize(obj, x_bar, u_bar)
            % linearize
            %   Linearize the unicycle model around the operating point
            %
            % Syntax
            %   [A_lin, B_lin] = obj.linearize(x_bar, u_bar)
            %
            % Input Arguments
            %   x_bar - Operating point of the states
            %       real vector
            %   u_bar - Operating point of the inputs
            %       real vector
            %
            % Output Arguments
            %   A_lin - Linearized state matrix
            %       real matrix
            %   B_lin - Linearized input matrix
            %       real matrix

            A_lin = double(subs(obj.sym_A, [obj.sym_x; obj.sym_u; obj.sym_r; obj.sym_L], [x_bar; u_bar; obj.r; obj.L]));
            B_lin = double(subs(obj.sym_B, [obj.sym_x; obj.sym_u; obj.sym_r; obj.sym_L], [x_bar; u_bar; obj.r; obj.L]));
        end

        % Fix reference angles function
        function x_ref_fixed = fix_angles(obj, x, x_ref)
            % fix_angles
            %   Fix the reference angles w.r.t. the current state to avoid
            %   discontinuities introduced by cuts such as [-pi, pi] or [0, 2*pi]
            %   angles representations
            %
            % Syntax
            %   x_ref_fixed = obj.fix_angles(x, x_ref)
            %
            % Input Arguments
            %   x - Current state
            %       real vector
            %   x_ref - Reference state
            %       real vector
            %
            % Output Arguments
            %   x_ref_fixed - Fixed reference state
            %       real vector

            % Compute the angle between the current and reference states
            delta_theta = atan2(sin(x(3:3:end) - x_ref(3:3:end)), cos(x(3:3:end) - x_ref(3:3:end)));

            % Fix the angular components w.r.t. current/predicted states
            x_ref_fixed = x_ref;
            x_ref_fixed(3:3:end) = x(3:3:end) - delta_theta;
        end

        % Extended Kalman Filter (EKF) step
        function x_hat = EKF_step(obj, x_hat, u, y)
            % EKF_step
            %   Estimates the state of the unicycle model using the Extended Kalman Filter (EKF)
            %   given a past state estimate, input, and output measurements
            %
            % Syntax
            %   x_hat = obj.EKF_step(x_hat, u, y)
            %
            % Input Arguments
            %   x_hat - State estimate
            %       real vector
            %   u - Input vector
            %       real vector
            %   y - Output vector
            %       real vector
            %
            % Output Arguments
            %   x_hat - Updated state estimate
            %       real vector

            % Prediction step

            % State transition matrix
            A = double(subs(obj.sym_A, [obj.sym_x; obj.sym_u; obj.sym_r; obj.sym_L], [x_hat; u; obj.r; obj.L]));

            % Predicted state estimate
            x_hat = obj.simulate(x_hat, u, obj.Ts);

            % Predicted covariance estimate
            P = A * obj.P * A' + obj.Q_tilde;

            % Update step

            % Output transformation matrix
            C = double(subs(obj.sym_C, [obj.sym_x; obj.sym_u; obj.sym_r; obj.sym_L], [x_hat; u; obj.r; obj.L]));

            % Kalman gain
            K = P * C' / (C * P * C' + obj.R_tilde);

            % Updated state estimate
            x_hat = x_hat + K * (y - obj.output(x_hat, u));

            % Updated covariance estimate
            obj.P = (eye(size(P)) - K * C) * P;
            
        end


    end
end