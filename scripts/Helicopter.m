classdef Helicopter < DynamicalSystem
    properties

        % Model parameters
        n = 8;
        m = 4;
        p = 4;
        Ts;
        bx;
        by;
        bz;
        bpsi;
        kx;
        ky;
        kpsi;
        ki;
        g = 9.81;

        % Reference trajectory
        x_ref = [];
        u_ref = [];

        % Constraints
        min_xi;
        max_xi;
        min_yi;
        max_yi;
        min_zi;
        max_zi;
        min_vxb;
        max_vxb;
        min_vyb;
        max_vyb;
        min_vzb;
        max_vzb;
        min_psi;
        max_psi;
        min_vpsi;
        max_vpsi;
        min_ux;
        max_ux;
        min_uy;
        max_uy;
        min_uz;
        max_uz;
        min_upsi;
        max_upsi;
        eps_x;
        f_x;
        eps_u;
        f_u;

        % Symbolic properties
        sym_bx;
        sym_by;
        sym_bz;
        sym_bpsi;
        sym_kx;
        sym_ky;
        sym_kpsi;
        sym_ki;
        sym_x;
        sym_u;
        sym_f;
        sym_g;
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
        function obj = Helicopter(parameters, Ts, x_constraints, u_constraints, P0, Q_tilde, R_tilde)

            % Check that parameters contains exactly 8 elements
            if length(parameters) ~= 8
                error('The parameters vector must be of the form [bx; by; bz; bpsi; kx; ky; kpsi; ki]');
            end

            % Initialize the helicopter parameters
            obj.bx = parameters(1);
            obj.by = parameters(2);
            obj.bz = parameters(3);
            obj.bpsi = parameters(4);
            obj.kx = parameters(5);
            obj.ky = parameters(6);
            obj.kpsi = parameters(7);
            obj.ki = parameters(8);
            obj.Ts = Ts;
            obj.min_xi = x_constraints(1, 1);
            obj.max_xi = x_constraints(1, 2);
            obj.min_yi = x_constraints(2, 1);
            obj.max_yi = x_constraints(2, 2);
            obj.min_zi = x_constraints(3, 1);
            obj.max_zi = x_constraints(3, 2);
            obj.min_vxb = x_constraints(4, 1);
            obj.max_vxb = x_constraints(4, 2);
            obj.min_vyb = x_constraints(5, 1);
            obj.max_vyb = x_constraints(5, 2);
            obj.min_vzb = x_constraints(6, 1);
            obj.max_vzb = x_constraints(6, 2);
            obj.min_psi = x_constraints(7, 1);
            obj.max_psi = x_constraints(7, 2);
            obj.min_vpsi = x_constraints(8, 1);
            obj.max_vpsi = x_constraints(8, 2);
            obj.min_ux = u_constraints(1, 1);
            obj.max_ux = u_constraints(1, 2);
            obj.min_uy = u_constraints(2, 1);
            obj.max_uy = u_constraints(2, 2);
            obj.min_uz = u_constraints(3, 1);
            obj.max_uz = u_constraints(3, 2);
            obj.min_upsi = u_constraints(4, 1);
            obj.max_upsi = u_constraints(4, 2);

            obj.eps_x = kron(eye(obj.n), [1; -1]);
            obj.f_x = [
                +obj.max_xi;
                -obj.min_xi;
                +obj.max_yi;
                -obj.min_yi;
                +obj.max_zi;
                -obj.min_zi;
                +obj.max_vxb;
                -obj.min_vxb;
                +obj.max_vyb;
                -obj.min_vyb;
                +obj.max_vzb;
                -obj.min_vzb;
                +obj.max_psi;
                -obj.min_psi;
                +obj.max_vpsi;
                -obj.min_vpsi;
            ];
            obj.eps_u = kron(eye(obj.m), [1; -1]);
            obj.f_u = [
                +obj.max_ux;
                -obj.min_ux;
                +obj.max_uy;
                -obj.min_uy;
                +obj.max_uz;
                -obj.min_uz;
                +obj.max_upsi;
                -obj.min_upsi;
            ];

            % Initialize symbolic properties
            syms sym_xi sym_yi sym_zi sym_vxb sym_vyb sym_vzb sym_psi sym_vpsi real
            syms sym_ux sym_uy sym_uz sym_upsi
            syms sym_bx sym_by sym_bz sym_bpsi sym_kx sym_ky sym_kpsi sym_ki
            
            obj.sym_bx = sym_bx;
            obj.sym_by = sym_by;
            obj.sym_bz = sym_bz;
            obj.sym_bpsi = sym_bpsi;
            obj.sym_kx = sym_kx;
            obj.sym_ky = sym_ky;
            obj.sym_kpsi = sym_kpsi;
            obj.sym_ki = sym_ki;
            obj.sym_x = [sym_xi; sym_yi; sym_zi; sym_vxb; sym_vyb; sym_vzb; sym_psi; sym_vpsi];
            obj.sym_u = [sym_ux; sym_uy; sym_uz; sym_upsi];

            dxidt = cos(sym_psi) * sym_vxb - sin(sym_psi) * sym_vyb;
            dyidt = sin(sym_psi) * sym_vxb + cos(sym_psi) * sym_vyb;
            dzidt = sym_vzb;
            dvxbdt = sym_bx * sym_ux + sym_kx * sym_vxb + sym_vpsi * sym_vyb;
            dvybdt = sym_by * sym_uy + sym_ky * sym_vyb - sym_vpsi * sym_vxb;
            dvzbd = sym_bz * sym_uz - obj.g;
            dpsidt = sym_vpsi;
            dvpsidt = sym_bpsi * sym_upsi + sym_kpsi * sym_vpsi;

            obj.sym_f = [dxidt; dyidt; dzidt; dvxbdt; dvybdt; dvzbd; dpsidt; dvpsidt];
            obj.sym_g = [sym_xi; sym_yi; sym_zi; sym_psi];
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

            dxidt = cos(x(7)) * x(4) - sin(x(7)) * x(5);
            dyidt = sin(x(7)) * x(4) + cos(x(7)) * x(5);
            dzidt = x(6);
            dvxbdt = obj.bx * u(1) + obj.kx * x(4) + x(8) * x(5);
            dvybdt = obj.by * u(2) + obj.ky * x(5) - x(8) * x(4);
            dvzbd = obj.bz * u(3) - obj.g;
            dpsidt = x(8);
            dvpsidt = obj.bpsi * u(4) + obj.kpsi * x(8);

            dxdt = [dxidt; dyidt; dzidt; dvxbdt; dvybdt; dvzbd; dpsidt; dvpsidt];
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

            % Wrap the angles to [0, 2*pi]
            x0(7) = wrapTo2Pi(x0(7));

            % Simulation
            [t, x] = ode45(@(t, x) obj.dynamics(t, x, u), [0, T], x0);
            x_final = (x(end, :))';

            % Wrap the angles to [0, 2*pi]
            x_final(7) = wrapTo2Pi(x_final(7));
        end

        % Output transformation function
        function y = output(obj, x, u)
            % output
            %   Output transformation function for the Helicopter model
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

            y = [x(1); x(2); x(3); x(7)];
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

            % Combine all symbolic variables and their corresponding values
            sym_vars = [
                obj.sym_x; obj.sym_u; 
                obj.sym_bx; obj.sym_by; obj.sym_bz; obj.sym_bpsi; 
                obj.sym_kx; obj.sym_ky; obj.sym_kpsi; obj.sym_ki;
            ];
            values = [
                x_bar; u_bar; 
                obj.bx; obj.by; obj.bz; obj.bpsi; 
                obj.kx; obj.ky; obj.kpsi; obj.ki;
            ];

            % Ensure the sizes match
            assert(length(sym_vars) == length(values), 'Sizes of symbolic variables and values arrays must match.');

            % Substitute the values into the symbolic expressions
            A_lin = double(subs(obj.sym_A, sym_vars, values));
            B_lin = double(subs(obj.sym_B, sym_vars, values));
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
            %   x - Current state (augmented MPC vector form)
            %       real vector
            %   x_ref - Reference state (augmented MPC vector form)
            %       real vector
            %
            % Output Arguments
            %   x_ref_fixed - Fixed reference state
            %       real vector

            idx = 7;
            step = obj.n;

            % Compute the angle between the current and reference states
            delta_psi = atan2(sin(x(idx:step:end) - x_ref(idx:step:end)), cos(x(idx:step:end) - x_ref(idx:step:end)));

            % Fix the angular components w.r.t. current/predicted states
            x_ref_fixed = x_ref;
            x_ref_fixed(idx:step:end) = x(idx:step:end) - delta_psi;
        end

        % Extended Kalman Filter (EKF) step
        function x_hat = EKF_step(obj, x_hat, u, y)
            % EKF_step
            %   Estimates the state of the Helicopter model using the Extended Kalman Filter (EKF)
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

            % Combine all symbolic variables and their corresponding values
            sym_vars = [
                obj.sym_x; obj.sym_u; 
                obj.sym_bx; obj.sym_by; obj.sym_bz; obj.sym_bpsi; 
                obj.sym_kx; obj.sym_ky; obj.sym_kpsi; obj.sym_ki;
            ];
            values = [
                x_hat; u;
                obj.bx; obj.by; obj.bz; obj.bpsi; 
                obj.kx; obj.ky; obj.kpsi; obj.ki;
            ];

            % Prediction step

            % State transition matrix
            A = double(subs(obj.sym_A, sym_vars, values));

            % Predicted state estimate
            x_hat = obj.simulate(x_hat, u, obj.Ts);

            % Predicted covariance estimate
            P = A * obj.P * A' + obj.Q_tilde;

            % Update step

            % Output transformation matrix
            C = double(subs(obj.sym_C, sym_vars, values));

            % Kalman gain
            K = P * C' / (C * P * C' + obj.R_tilde);

            % Updated state estimate
            x_hat = x_hat + K * (y - obj.output(x_hat, u));

            % Updated covariance estimate
            obj.P = (eye(size(P)) - K * C) * P;
            
        end

        % Trajectory generation function
        function [x_ref, u_ref, Tend] = generate_trajectory(obj, N_guide, shape, extra_params)
            %//TODO: make this model also return Tend which is needed for MPC

            % Common generation parameters
            N_intervals = N_guide - 1;
            delta = 2 * pi / N_intervals;
            m_theta = delta / obj.Ts; % angulat coefficient of theta(t) = m_theta * t

            % Analytical parametrization
            syms t real;
            f_theta = @(t) m_theta * t;
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
                z = [radius * cos(f_theta(t)), radius * sin(f_theta(t)), 0, 0.5 * pi + f_theta(t)];
                dz = diff(z, t);
                ddz = diff(dz, t);

                % Analytical definition of remaining states and inputs
                dxb = cos(z(4)) * dz(1) + sin(z(4)) * dz(2);
                dyb = -sin(z(4)) * dz(1) + cos(z(4)) * dz(2);
                dzb = dz(3);
                dpsi = dz(4);
                ux = (cos(z(4)) * (ddz(1) - obj.kx * dz(1)) + sin(z(4)) * (ddz(2) - obj.kx * dz(2)) / obj.bx);
                uy = (cos(z(4)) * (ddz(2) - obj.ky * dz(2)) + sin(z(4)) * (-ddz(1) + obj.ky * dz(1)) / obj.by);
                uz = (ddz(3) + obj.g) / obj.bz;
                upsi = (ddz(4) - obj.kpsi * dz(4)) / obj.bpsi;

                x = [z(1), z(2), z(3), dxb, dyb, dzb, z(4), dpsi];
                u = [ux, uy, uz, upsi];

                % Reference trajectory
                obj.x_ref = [];
                obj.u_ref = [];
                for i = 1:N_guide - 1
                    x_t = double(subs(x, t, T_guide(i)));
                    u_t = double(subs(u, t, T_guide(i)));

                    x_t(7) = wrapTo2Pi(x_t(7)); % wrap the angle state

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
                z = [lemin_x, lemin_y, 0, atan2(diff(lemin_y, t)/diff(sym_theta, t), diff(lemin_x, t)/diff(sym_theta, t))]; 
                dz = diff(z, t);
                ddz = diff(dz, t);

                % Analytical definition of remaining states and inputs
                dxb = cos(z(4)) * dz(1) + sin(z(4)) * dz(2);
                dyb = -sin(z(4)) * dz(1) + cos(z(4)) * dz(2);
                dzb = dz(3);
                dpsi = dz(4);
                ux = (cos(z(4)) * (ddz(1) - obj.kx * dz(1)) + sin(z(4)) * (ddz(2) - obj.kx * dz(2)) / obj.bx);
                uy = (cos(z(4)) * (ddz(2) - obj.ky * dz(2)) + sin(z(4)) * (-ddz(1) + obj.ky * dz(1)) / obj.by);
                uz = (ddz(3) + obj.g) / obj.bz;
                upsi = (ddz(4) - obj.kpsi * dz(4)) / obj.bpsi;

                x = [z(1), z(2), z(3), dxb, dyb, dzb, z(4), dpsi];
                u = [ux, uy, uz, upsi];

                % Reference trajectory
                obj.x_ref = [];
                obj.u_ref = [];
                for i = 1:N_guide - 1
                    x_t = double(subs(x, t, 1e-6+T_guide(i)));
                    u_t = double(subs(u, t, 1e-6+T_guide(i)));

                    x_t(7) = wrapTo2Pi(x_t(7)); % wrap the angle state

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
                z = [batman_x, batman_y, 0, atan2(diff(batman_y, t), diff(batman_x, t))];
                dz = diff(z, t);
                ddz = diff(dz, t);

                % Analytical definition of remaining states and inputs
                dxb = cos(z(4)) * dz(1) + sin(z(4)) * dz(2);
                dyb = -sin(z(4)) * dz(1) + cos(z(4)) * dz(2);
                dzb = dz(3);
                dpsi = dz(4);
                ux = (cos(z(4)) * (ddz(1) - obj.kx * dz(1)) + sin(z(4)) * (ddz(2) - obj.kx * dz(2)) / obj.bx);
                uy = (cos(z(4)) * (ddz(2) - obj.ky * dz(2)) + sin(z(4)) * (-ddz(1) + obj.ky * dz(1)) / obj.by);
                uz = (ddz(3) + obj.g) / obj.bz;
                upsi = (ddz(4) - obj.kpsi * dz(4)) / obj.bpsi;

                x = [z(1), z(2), z(3), dxb, dyb, dzb, z(4), dpsi];
                u = [ux, uy, uz, upsi];


                % Distribute reference points more uniformly
                N_head = 6; % points for t in [-2,2]
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
                    % x_t = double(subs(x, t, 1e-6+T_guide(i)));
                    % u_t = double(subs(u, t, 1e-6+T_guide(i)));
                    x_t = double(subs(x, t, 1e-6+l(i)));
                    u_t = double(subs(u, t, 1e-6+l(i)));

                    x_t(7) = wrapTo2Pi(x_t(7)); % wrap the angle state

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
                    ddz = diff(dz, x);

                    % Generate missing points
                    N_filling = ceil((t1 - t0)/obj.Ts) + 1;
                    T_filling = linspace(t0, t1, N_filling);
                    for j = 1:length(T_filling) - 1

                        % Evaluate z(t) and derivatives
                        z_t = double(subs(z, x, T_filling(j)))';
                        z_t(obj.p) = wrapTo2Pi(z_t(obj.p)); % wrap the angle state
                        dz_t = double(subs(dz, x, T_filling(j)))';
                        ddz_t = double(subs(ddz, x, T_filling(j)))';

                        % Compute remaining states and inputs
                        dxb = cos(z_t(4))*dz_t(1) + sin(z_t(4))*dz_t(2);
                        dyb = -sin(z_t(4))*dz_t(1) + cos(z_t(4))*dz_t(2);
                        dzb = dz_t(3);
                        dpsi = dz_t(4);

                        x_t = [z_t(1), z_t(2), z_t(3), dxb, dyb, dzb, wrapTo2Pi(z_t(4)), dpsi];

                        ux = (cos(z_t(4))*(ddz_t(1) - obj.kx*dz_t(1)) + sin(z_t(4))*(ddz_t(2) - obj.kx*dz_t(2))/obj.bx);
                        uy = (cos(z_t(4))*(ddz_t(2) - obj.ky*dz_t(2)) + sin(z_t(4))*(-ddz_t(1) + obj.ky*dz_t(1))/obj.by);
                        uz = (ddz_t(3) + obj.g)/obj.bz;
                        upsi = (ddz_t(4) - obj.kpsi*dz_t(4))/obj.bpsi;

                        u_t = [ux, uy, uz, upsi];

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