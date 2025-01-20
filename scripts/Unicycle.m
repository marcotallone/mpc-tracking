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

        % Reference trajectory
        x_ref = [];
        u_ref = [];

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

        % Extended Kalman Filter (EKF) state estimation
        function x_hat = EKF_estimate(obj, x_hat, u, y)
            % EKF_estimate
            %   Estimates the state of the unicycle model using the Extended Kalman Filter (EKF)
            %   given a past state estimate, input, and output measurements
            %
            % Syntax
            %   x_hat = obj.EKF_estimate(x_hat, u, y)
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


        % Trajectory generation function
        function [x_ref, u_ref, Tend] = generate_trajectory(obj, N_guide, shape, extra_params)

            % Common generation parameters
            N_intervals = N_guide - 1;
            delta = 2 * pi / N_intervals;
            m_theta = delta / obj.Ts; % angular coefficient of theta(t) = m_theta * t

            % Analytical parametrization
            syms t real;
            sym_theta = m_theta*t;

            % Circular trajectory
            if nargin < 4 && strcmp(shape, 'circle')
                error('Please provide the radius of the circle trajectory.');
            elseif nargin == 4 && strcmp(shape, 'circle')

                % Set radius
                assert(isscalar(extra_params), 'The extra parameter radius must be a scalar value.');
                radius = extra_params;

                % Simulation time and guide time steps
                Tend = N_intervals * obj.Ts;
                T_guide = linspace(0, Tend, N_guide);

                % Analytical definition
                circle_x = radius * cos(sym_theta);
                circle_y = radius * sin(sym_theta);
                z = [circle_x, circle_y, atan2(diff(circle_y, t), diff(circle_x, t))];
                dz = diff(z, t);

                % Analytical definition of remaining states and inputs
                omega1 = (2*sqrt(dz(1)^2 + dz(2)^2) - obj.L*dz(3))/(2*obj.r); 
                omega2 = (2*sqrt(dz(1)^2 + dz(2)^2) + obj.L*dz(3))/(2*obj.r); 
                
                x = [z(1), z(2), z(3)];
                u = [omega1, omega2];

                % Reference trajectory
                obj.x_ref = [];
                obj.u_ref = [];
                for i = 1:N_guide - 1
                    x_t = double(subs(x, t, 1e-6+T_guide(i)));
                    u_t = double(subs(u, t, 1e-6+T_guide(i)));

                    x_t(3) = wrapTo2Pi(x_t(3)); % wrap the angle state

                    obj.x_ref = [obj.x_ref; x_t];
                    obj.u_ref = [obj.u_ref; u_t];
                end
            end

            % Leminscate trajectory
            if nargin < 4 && strcmp(shape, 'leminscate')
                error('Please provide the a parameter of the leminscate trajectory');
            elseif nargin == 4 && strcmp(shape, 'leminscate')

                % Set a parameter
                assert(isscalar(extra_params), 'The extra parameter a must be a scalar value.');
                a = extra_params;

                % Simulation time and guide time steps
                Tend = N_intervals * obj.Ts;
                T_guide = linspace(0, Tend, N_guide);

                % Analytical definition
                lemin_x = (a*sqrt(2)*cos(sym_theta))/(sin(sym_theta)^2 + 1);
                lemin_y = (a*sqrt(2)*cos(sym_theta)*sin(sym_theta))/(sin(sym_theta)^2 + 1);
                % z = [lemin_x, lemin_y, atan2(diff(lemin_y, t)/diff(sym_theta, t), diff(lemin_x, t)/diff(sym_theta, t))]; 
                z = [lemin_x, lemin_y, atan2(diff(lemin_y, t), diff(lemin_x, t))]; 
                dz = diff(z, t);

                % Analytical definition of remaining states and inputs
                omega1 = (2*sqrt(dz(1)^2 + dz(2)^2) - obj.L*dz(3))/(2*obj.r); 
                omega2 = (2*sqrt(dz(1)^2 + dz(2)^2) + obj.L*dz(3))/(2*obj.r); 
                
                x = [z(1), z(2), z(3)];
                u = [omega1, omega2];

                % Reference trajectory
                obj.x_ref = [];
                obj.u_ref = [];
                for i = 1:N_guide - 1
                    x_t = double(subs(x, t, 1e-6+T_guide(i)));
                    % u_t = double(subs(u, t, 1e-6+T_guide(i)));

                    % Compute the average input in [T_guide(i), T_guide(i+1)] between n_samples samples
                    n_samples = 25;
                    u_t_average = zeros(1, obj.m);
                    for j = 1:n_samples
                        u_t = double(subs(u, t, 1e-6+T_guide(i) + j*(T_guide(i+1) - T_guide(i))/n_samples));
                        u_t_average = u_t_average + u_t;
                    end
                    u_t_average = u_t_average / n_samples;

                    u_t = u_t_average;

                    x_t(3) = wrapTo2Pi(x_t(3)); % wrap the angle state

                    obj.x_ref = [obj.x_ref; x_t];
                    obj.u_ref = [obj.u_ref; u_t];
                end
            end

            % Batman trajectory
            if nargin < 4 && strcmp(shape, 'batman')

                % Assert that N_guide is at least 10
                assert(N_guide >= 10, 'The number of guide points must be at least 10 for the Batman trajectory.');

                % If N_guide odd add 1 to make it even
                if mod(N_guide, 2) == 1
                    N_guide = N_guide + 1;
                end

                % Simulation time and guide time steps
                Tend = N_intervals * obj.Ts;
                T_guide = linspace(0, Tend, N_guide);

                % Analytical definition
                batman_x = (abs(t)/t) * (0.3*abs(t) + 0.2*abs(abs(t)-1) + 2.2*abs(abs(t)-2) ...
                    - 2.7*abs(abs(t)-3) - 3*abs(abs(t)-5) + 3*abs(abs(t)-7) ...
                    + 5*sin((pi/4)*(abs(abs(t)-3) - abs(abs(t)-4) + 1)) ...
                    + (5/4)*(abs(abs(t)-4) - abs(abs(t)-5) - 1)^3 ...
                    - 5.3*cos(((pi/2) + asin(47/53)) * ((abs(abs(t)-7) - abs(abs(t)-8) - 1)/2)) ...
                    + 2.8);
                batman_y = (3/2)*abs(abs(t)-1) - (3/2)*abs(abs(t)-2) - (29/4)*abs(abs(t)-4) ...
                    + (29/4)*abs(abs(t)-5) + (7/16)*(abs(abs(t)-2) - abs(abs(t)-3) - 1)^4 ...
                    + 4.5*sin((pi/4)*(abs(abs(t)-3) - abs(abs(t)-4) - 1)) ...
                    - 3*(sqrt(2)/5) * (abs(abs(abs(t)-5) - abs(abs(t)-7)))^(5/2) ...
                    + 6.4*sin(((pi/2) + asin(47/53)) * ((abs(abs(t)-7) - abs(abs(t)-8) + 1)/2) + asin(56/64)) ...
                    + 4.95;
                z = [batman_x, batman_y, atan2(diff(batman_y, t), diff(batman_x, t))];
                dz = diff(z, t);

                % Analytical definition of remaining states and inputs
                omega1 = (2*sqrt(dz(1)^2 + dz(2)^2) - obj.L*dz(3))/(2*obj.r); 
                omega2 = (2*sqrt(dz(1)^2 + dz(2)^2) + obj.L*dz(3))/(2*obj.r); 
                
                x = [z(1), z(2), z(3)];
                u = [omega1, omega2];

                % Distribute reference points more uniformly
                N_head = 6;
                N_left = (N_guide - N_head)/2;
                N_right = N_left;
    
                l_left = linspace(-8, -2, N_left);
                l_head = linspace(-2, 2, N_head + 2);
                l_right = linspace(2, 8, N_right);

                l = unique([l_left, l_head(2:end-1), l_right]);

                % Reference trajectory
                obj.x_ref = [];
                obj.u_ref = [];
                for i = 1:N_guide - 1
                    x_t = double(subs(x, t, 1e-6+l(i)));
                    
                    if i == N_guide-1
                        u_t = double(subs(u, t, 1e-6+l(i)));
                    else

                        % Compute the average input in [T_guide(i), T_guide(i+1)] between n_samples samples
                        n_samples = 10;
                        u_t_average = zeros(1, obj.m);
                        for j = 1:n_samples
                            u_t = double(subs(u, t, 1e-6+l(i) + j*(l(i+1) - l(i))/n_samples));
                            u_t_average = u_t_average + u_t;
                        end
                        u_t_average = u_t_average / n_samples;

                        u_t = u_t_average;
                    end

                    x_t(3) = wrapTo2Pi(x_t(3)); % wrap the angle state

                    obj.x_ref = [obj.x_ref; x_t];
                    obj.u_ref = [obj.u_ref; u_t];
                end
                
            end

            % Arbitrary trajectory with Murray generation method
            if nargin < 4 && strcmp(shape, 'arbitrary')
                error('Please provide a cell array containing {N_points_filling, Z_guide} as extra parameters.');
            elseif nargin == 4 && strcmp(shape, 'arbitrary')

                assert(iscell(extra_params), 'The extra parameters must be a cell array containing {N_points_filling, Z_guide}');

                % Extract the extra parameters
                N_points_filling = extra_params{1};
                Z_guide = extra_params{2};

                assert(isscalar(N_points_filling), 'The number of points to fill must be a scalar value.');
                assert(size(Z_guide, 1) == N_guide, 'The guide points matrix must have N_guide rows: as many as the intervals/guide points.');
                assert(size(Z_guide, 2) == obj.p, 'The guide points matrix must have p columns: as many as the flat outputs.');

                % Generation parameters
                N_intervals = N_guide - 1;
                Tend = N_intervals * N_points_filling * obj.Ts;
                T_guide = linspace(0, Tend, N_guide);
                N_basis = 2;
                order = 0;

                % Basis functions
                syms x
                basis = sym('x', [1 N_basis]);
                for k = 1:N_basis
                    basis(k) = x^(k-1);
                end

                % Trajectory filling
                obj.x_ref = [];
                obj.u_ref = [];
                for i = 1:length(T_guide) - 1

                    % Time interval
                    t0 = T_guide(i);
                    t1 = T_guide(i+1);

                    % Matrix M
                    m0 = obj.m_matrix(t0, basis, order);
                    m1 = obj.m_matrix(t1, basis, 0);
                    M0 = kron(eye(obj.p), m0);
                    M1 = kron(eye(obj.p), m1);
                    M = [M0; M1];

                    % Guide points in the interval
                    z_bar = [
                        Z_guide(i, 1:obj.p)';
                        Z_guide(i+1, 1:obj.p)';
                    ];

                    % Check that M is full column rank
                    if rank(M) < size(M, 2)
                        disp('Matrix M is not full column rank');
                        return;
                    end

                    % Solve the system and reshape
                    alpha = M\z_bar;
                    alpha = reshape(alpha, [], obj.p)';

                    % Reference function z(t) and its derivatives
                    z = alpha*basis.';
                    dz = diff(z, x);

                    % Generate missing points
                    N_filling = ceil((t1 - t0)/obj.Ts) + 1;
                    T_filling = linspace(t0, t1, N_filling);
                    for j = 1:length(T_filling) - 1

                        % Evaluate z(t) and derivatives
                        z_t = double(subs(z, x, T_filling(j)))';
                        z_t(obj.p) = wrapTo2Pi(z_t(obj.p)); % wrap the angle state
                        dz_t = double(subs(dz, x, T_filling(j)))';

                        % Compute remaining states and inputs
                        dtheta = atan2(dz_t(2), dz_t(1));
                        omega1 = (2*sqrt(dz_t(1)^2 + dz_t(2)^2) - obj.L*dtheta)/(2*obj.r); 
                        omega2 = (2*sqrt(dz_t(1)^2 + dz_t(2)^2) + obj.L*dtheta)/(2*obj.r); 
                        
                        x_t = [z_t(1), z_t(2), dtheta];
                        u_t = [omega1, omega2];

                        % Store the reference flat-output, states and inputs
                        obj.x_ref = [obj.x_ref; x_t];
                        obj.u_ref = [obj.u_ref; u_t];
                    end
                end
            end

            % Return reference states and inputs
            x_ref = obj.x_ref;
            u_ref = obj.u_ref;
        end


    end
end