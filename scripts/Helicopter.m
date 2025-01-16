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
        % x_ref = 0.0;
        % y_ref = 0.0;
        g = 9.81;

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
        % min_xint;
        % max_xint;
        % min_yint;
        % max_yint;
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
        % sym_x_ref;
        % sym_y_ref;
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
            % obj.min_xint = x_constraints(9, 1);
            % obj.max_xint = x_constraints(9, 2);
            % obj.min_yint = x_constraints(10, 1);
            % obj.max_yint = x_constraints(10, 2);
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
                % +obj.max_xint;
                % -obj.min_xint;
                % +obj.max_yint;
                % -obj.min_yint;               
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
            syms sym_xi sym_yi sym_zi sym_vxb sym_vyb sym_vzb sym_psi sym_vpsi real %sym_xint sym_yint real
            syms sym_ux sym_uy sym_uz sym_upsi
            syms sym_bx sym_by sym_bz sym_bpsi sym_kx sym_ky sym_kpsi sym_ki
            % syms sym_x_ref sym_y_ref
            
            obj.sym_bx = sym_bx;
            obj.sym_by = sym_by;
            obj.sym_bz = sym_bz;
            obj.sym_bpsi = sym_bpsi;
            obj.sym_kx = sym_kx;
            obj.sym_ky = sym_ky;
            obj.sym_kpsi = sym_kpsi;
            obj.sym_ki = sym_ki;
            % obj.sym_x_ref = sym_x_ref;
            % obj.sym_y_ref = sym_y_ref;
            obj.sym_x = [sym_xi; sym_yi; sym_zi; sym_vxb; sym_vyb; sym_vzb; sym_psi; sym_vpsi]; % sym_xint; sym_yint];
            obj.sym_u = [sym_ux; sym_uy; sym_uz; sym_upsi];

            dxidt = cos(sym_psi) * sym_vxb - sin(sym_psi) * sym_vyb;
            dyidt = sin(sym_psi) * sym_vxb + cos(sym_psi) * sym_vyb;
            dzidt = sym_vzb;
            dvxbdt = sym_bx * sym_ux + sym_kx * sym_vxb + sym_vpsi * sym_vyb;
            dvybdt = sym_by * sym_uy + sym_ky * sym_vyb - sym_vpsi * sym_vxb;
            dvzbd = sym_bz * sym_uz - obj.g;
            dpsidt = sym_vpsi;
            dvpsidt = sym_bpsi * sym_upsi + sym_kpsi * sym_vpsi;
            % dxintdt = sym_ki * (sym_xi - sym_x_ref);
            % dyintdt = sym_ki * (sym_yi - sym_y_ref);

            obj.sym_f = [dxidt; dyidt; dzidt; dvxbdt; dvybdt; dvzbd; dpsidt; dvpsidt]; %dxintdt; dyintdt];
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
            % dxintdt = obj.ki * (x(1) - obj.x_ref);
            % dyintdt = obj.ki * (x(2) - obj.y_ref);

            dxdt = [dxidt; dyidt; dzidt; dvxbdt; dvybdt; dvzbd; dpsidt; dvpsidt]; %dxintdt; dyintdt];
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
                % obj.sym_x_ref; obj.sym_y_ref
            ];
            values = [
                x_bar; u_bar; 
                obj.bx; obj.by; obj.bz; obj.bpsi; 
                obj.kx; obj.ky; obj.kpsi; obj.ki;
                % obj.x_ref; obj.y_ref
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
            %   x - Current state
            %       real vector
            %   x_ref - Reference state
            %       real vector
            %
            % Output Arguments
            %   x_ref_fixed - Fixed reference state
            %       real vector

            idx = 7;
            step = obj.n;

            % Compute the angle between the current and reference states
            % delta_psi = atan2(sin(x(7:8:end) - x_ref(7:8:end)), cos(x(7:8:end) - x_ref(7:8:end)));
            delta_psi = atan2(sin(x(idx:step:end) - x_ref(idx:step:end)), cos(x(idx:step:end) - x_ref(idx:step:end)));

            % Fix the angular components w.r.t. current/predicted states
            x_ref_fixed = x_ref;
            % x_ref_fixed(7:8:end) = x(7:8:end) - delta_psi;
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

    end
end