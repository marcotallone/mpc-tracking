% ┌─────────────────────────────────────────────────────────────────────────┐ 
% │                       Dynamical System Class                            │ 
% └─────────────────────────────────────────────────────────────────────────┘
% Abstract class to define dynamical systems
%
% Usage
%   classdef MySystem < DynamicalSystem
%
% Methods
%   discretize - Discretize a continuous-time linear system 

classdef (Abstract) DynamicalSystem
    properties (Abstract)
        n       % number of states
        m       % number of inputs
        p       % number of outputs
        Ts      % sampling time
        eps_x;  % state constraints matrix
        f_x;    % state constraints vector
        eps_u;  % input constraints matrix
        f_u;    % input constraints vector
    end
    
    methods (Abstract)
        dxdt = dynamics(obj, t, x, u)
        x_final = simulate(obj, x0, u, T)
        [A_lin, B_lin] = linearize(obj, x_bar, u_bar)
        x_ref_fixed = fix_angles(obj, x, x_ref)
    end
   
    methods
        % Discretization function (from linear system)
        function [A_discrete, B_discrete] = discretize(obj, A, B)
            % discretize
            %   Discretize a continuous-time linear system using Euler method
            %
            % Syntax
            %   [A_discrete, B_discrete] = obj.discretize(A, B)
            %
            % Input Arguments
            %   A - State matrix
            %       real matrix
            %   B - Input matrix
            %       real matrix
            %
            % Output Arguments
            %   A_discrete - Discretized state matrix
            %       real matrix
            %   B_discrete - Discretized input matrix
            %       real matrix

            % Euler discretization
            A_discrete = eye(obj.n) + A * obj.Ts;
            B_discrete = B * obj.Ts;

            % #TODO: Check why exact sampling (below) is not working properly
            % continuous_sys = ss(A, B, eye(obj.n), zeros(obj.n, obj.m));
            % discrete_sys = c2d(continuous_sys, obj.Ts);
            % A_discrete = discrete_sys.A;
            % B_discrete = discrete_sys.B;
        end
    end 

end