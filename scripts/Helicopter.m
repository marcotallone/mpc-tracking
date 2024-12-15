
classdef Helicopter < DynamicalSystem
    properties
        n = 10;
        m = 4;
        p = 10;
        Ts;
        bx;
        by;
        bz;
        bphi;
        kx;
        ky;
        kphi;
        ki;
        g = 9.81;
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
        min_phi;
        max_phi;
        min_vphi;
        max_vphi;
        min_ux;
        max_ux;
        min_uy;
        max_uy;
        min_uz;
        max_uz;
        min_uphi;
        max_uphi;
        eps_x;
        f_x;
        eps_u;
        f_u;
        sym_bx;
        sym_by;
        sym_bz;
        sym_bphi;
        sym_kx;
        sym_ky;
        sym_kphi;
        sym_ki;
        sym_x;
        sym_u;
        sym_f;
        sym_A;
        sym_B;
    end

   methods
        % Constructor to initialize the unicycle parameters and state
        function obj = Helicopter(parameters, Ts, x_constraints, u_constraints)

            % Check that parameters contains exactly 8 elements
            if length(parameters) ~= 8
                error('The parameters vector must be of the form [bx; by; bz; bphi; kx; ky; kphi; ki]');
            end

            % Initialize the helicopter parameters
            obj.bx = parameters(1);
            obj.by = parameters(2);
            obj.bz = parameters(3);
            obj.bphi = parameters(4);
            obj.kx = parameters(5);
            obj.ky = parameters(6);
            obj.kphi = parameters(7);
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
            obj.min_phi = x_constraints(7, 1);
            obj.max_phi = x_constraints(7, 2);
            obj.min_vphi = x_constraints(8, 1);
            obj.max_vphi = x_constraints(8, 2);
            obj.min_ux = u_constraints(1, 1);
            obj.max_ux = u_constraints(1, 2);
            obj.min_uy = u_constraints(2, 1);
            obj.max_uy = u_constraints(2, 2);
            obj.min_uz = u_constraints(3, 1);
            obj.max_uz = u_constraints(3, 2);
            obj.min_uphi = u_constraints(4, 1);
            obj.max_uphi = u_constraints(4, 2);

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
                +obj.max_phi;
                -obj.min_phi;
                +obj.max_vphi;
                -obj.min_vphi;                
            ];
            obj.eps_u = kron(eye(obj.m), [1; -1]);
            obj.f_u = [
                +obj.max_ux;
                -obj.min_ux;
                +obj.max_uy;
                -obj.min_uy;
                +obj.max_uz;
                -obj.min_uz;
                +obj.max_uphi;
                -obj.min_uphi;
            ];

            % Initialize symbolic properties
            syms sym_xi sym_yi sym_zi sym_vxb sym_vyb sym_vzb sym_phi sym_vphi real
            syms sym_ux sym_uy sym_uz sym_uphi real
            syms sym_bx sym_by sym_bz sym_bphi sym_kx sym_ky sym_kphi sym_ki real
            
            obj.sym_bx = sym_bx;
            obj.sym_by = sym_by;
            obj.sym_bz = sym_bz;
            obj.sym_bphi = sym_bphi;
            obj.sym_kx = sym_kx;
            obj.sym_ky = sym_ky;
            obj.sym_kphi = sym_kphi;
            obj.sym_ki = sym_ki;
            obj.sym_x = [sym_xi; sym_yi; sym_zi; sym_vxb; sym_vyb; sym_vzb; sym_phi; sym_vphi];
            obj.sym_u = [sym_ux; sym_uy; sym_uz; sym_uphi];

            dxidt = cos(sym_phi) * sym_vxb - sin(sym_phi) * sym_vyb;
            dyidt = sin(sym_phi) * sym_vxb + cos(sym_phi) * sym_vyb;
            dzidt = sym_vzb;
            dvxbdt = sym_bx * sym_ux + sym_kx * sym_vxb + sym_vphi * sym_vyb;
            dvybdt = sym_by * sym_uy + sym_ky * sym_vyb - sym_vphi * sym_vxb;
            dvzbd = sym_bz * sym_uz - obj.g;
            dphidt = sym_vphi;
            dvphidt = sym_bphi * sym_uphi + sym_kphi * sym_vphi;

            obj.sym_f = [dxidt; dyidt; dzidt; dvxbdt; dvybdt; dvzbd; dphidt; dvphidt];
            obj.sym_A = jacobian(obj.sym_f, obj.sym_x);
            obj.sym_B = jacobian(obj.sym_f, obj.sym_u);
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
            dphidt = x(8);
            dvphidt = obj.bphi * u(4) + obj.kphi * x(8);

            dxdt = [dxidt; dyidt; dzidt; dvxbdt; dvybdt; dvzbd; dphidt; dvphidt];
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
            delta_phi = atan2(sin(x(7:8:end) - x_ref(7:8:end)), cos(x(7:8:end) - x_ref(7:8:end)));

            % Fix the angular components w.r.t. current/predicted states
            x_ref_fixed = x_ref;
            x_ref_fixed(7:8:end) = x(7:8:end) - delta_phi;
        end

    end
end