% ┌─────────────────────────────────────────────────────────────────────────┐ 
% │                  Model Predictive Control (MPC) Class                   │ 
% └─────────────────────────────────────────────────────────────────────────┘ 
% by Marco Tallone, 2024
%
% Class implementing a Model Predictive Control (MPC) algorithm for 
% non-linear systems with either dense or sparse formulation.
%
% Creation
% 	Syntax
% 		obj = MPC(model, x0, Tend, N, Q, R, x_ref, u_ref, preview, formulation, noise, debug)
% 		obj = MPC(model, x0, Tend, N, Q, R, x_ref, u_ref, preview, formulation)
% 		obj = MPC(model, x0, Tend, N, Q, R, x_ref, u_ref, preview)
% 		obj = MPC(model, x0, Tend, N, Q, R, x_ref, u_ref)
% 		obj = MPC(model, x0, Tend, N, Q, R, x_ref)
% 		obj = MPC(model, x0, Tend, N, Q, R)
% 		obj = MPC(model, x0, Tend, N, Q)
% 		obj = MPC(model, x0, Tend, N)
%
% 	Input Arguments
% 		model - System model
% 			DynamicalSystem object
% 		x0 - Initial state
% 			real vector
% 		Tend - Simulation time
% 			real scalar
% 		N - Prediction horizon
% 			integer
% 		Q - State cost
% 			real matrix
% 			optional, default: identity matrix
% 		R - Input cost
% 			real matrix
% 			optional, default: identity matrix
% 		x_ref - Reference states
% 			real matrix
% 			optional, default: zeros matrix (origin)
% 		u_ref - Reference inputs
% 			real matrix
% 			optional, default: zeros matrix
% 		preview - MPC preview flag
% 			integer
% 			0: no preview, 1: preview
% 			optional, default: 1
% 		formulation - MPC formulation flag
% 			integer
% 			0: dense (explicit), 1: sparse (implicit)
% 			optional, default: 0
% 		debug - Debug flag
% 			integer
% 			0: no debug, 1: debug
% 			optional, default: 0
%
% Properties
% 	model - System model
% 		object
% 	Ts - Model sampling time
% 		real scalar
% 	x0 - Initial state
% 		real vector
% 	t - Vector of discrete simulation times with step Ts
% 		real vector
% 	Nsteps - Number of simulations steps
% 		integer
% 	N - Prediction horizon
% 		integer
% 	Q - State cost
% 		real matrix
% 	R - Input cost
% 		real matrix
% 	x_ref - Reference states
% 		real matrix
% 	u_ref - Reference inputs
% 		real matrix
% 	X_REF - Reference states horizon vector
% 		real vector
% 	U_REF - Reference inputs horizon vector
% 		real vector
% 	Z_REF - Reference states and inputs horizon vector
% 		real vector
% 	X_BAR - Linearization states horizon vector
% 		real vector
% 	U_BAR - Linearization inputs horizon vector
% 		real vector
% 	Z_BAR - Linearization states and inputs horizon vector
% 		real vector
% 	x_pred - Prediction states during horizon
% 		real vector
% 	preview - MPC preview flag
% 		0: no preview, 1: preview
% 		integer
% 	formulation - MPC formulation flag
% 		0: dense (explicit), 1: sparse (implicit)
% 		integer
% 	debug - Debug flag
% 		0: no debug, 1: debug
% 		integer
% 	options - Quadprog optimization options
% 		optimoptions object
%
% Methods
% 	preview_reference - Set reference states and inputs with preview
% 	no_preview_reference - Set reference states and inputs without preview
% 	dense_formulation - Set MPC dense formulation matrices
% 	sparse_formulation - Set MPC sparse formulation matrices
% 	dense_constraints - Set MPC constraints matrices for dense formulation
% 	sparse_constraints - Set MPC constraints matrices for sparse formulation
% 	solve_dense - Solve MPC optimization problem with dense formulation
% 	solve_sparse - Solve MPC optimization problem with sparse formulation
% 	optimize - MPC optimization function
%
% Examples
% 	mpc = MPC(model, x0, Tend, N, Q, R, x_ref, u_ref, preview, formulation, debug);
% 	[x, u] = mpc.optimize();

classdef MPC < handle
	properties
		model;      		% system model
		Ts;         		% model sampling time
		x0;         		% initial state
		t;          		% vector of discrete simulation times with step Ts
		Nsteps;     		% number of simulations steps
		N;          		% prediction horizon
		Q;          		% state cost
		R;          		% input cost
		x_ref;      		% reference states
		u_ref;      		% reference inputs
		X_REF;	  			% reference states horizon vector
		U_REF;	  			% reference inputs horizon vector
		Z_REF;	  			% reference states and inputs horizon vector
		X_BAR;	  			% linearization states horizon vector
		U_BAR;	  			% linearization inputs horizon vector
		Z_BAR;	  			% linearization states and inputs horizon vector
		x_pred;     		% prediction states during horizon
		preview;    		% MPC preview flag (0: no preview, 1: preview)
		formulation			% MPC formulation flag (0: dense (explicit), 1: sparse (implicit))
        noise;              % noise flag (0: no noise, 1: noise)
		debug;    			% debug flag (0: no debug, 1: debug)
		set_reference; 		% function handle to set reference states and inputs
		set_formulation;	% function handle to set MPC formulation
		set_constraints; 	% function handle to set MPC constraints
		solve;      		% function handle to solve MPC optimization problem
		options;    		% quadprog optimization options
	end
	
	methods
		% Constructor
		function obj = MPC(model, x0, Tend, N, Q, R, x_ref, u_ref, preview, formulation, noise, debug)

			% Set default values for optional arguments
			if nargin < 12
				debug = 0;
			end
            if nargin < 11
                noise = 0;
            end
			if nargin < 10
				formulation = 0;
			end
			if nargin < 9
				preview = 1;
			end
			if nargin < 8
				u_ref = zeros(length(t), model.m);
			end
			if nargin < 7
				x_ref = zeros(length(t), model.n);
			end
			if nargin < 6
				R = eye(model.m);
			end
			if nargin < 5
				Q = eye(model.n);
			end
			if nargin < 4
				error('Prediction horizon N must be provided');
			end
			if nargin < 3
				error('Simulation time Tend must be provided');
			end
			if nargin < 2
				error('Initial state x0 must be provided');
			end
			if nargin < 1
				error('System model must be provided');
			end

			% Initialize properties
			obj.model = model;
			obj.Ts = obj.model.Ts;
			obj.x0 = x0;
			obj.t = 0:obj.Ts:Tend;
			obj.Nsteps = length(obj.t)-(N+1);
			obj.N = N;
			obj.Q = Q;
			obj.R = R;
			obj.x_ref = reshape(x_ref', [], 1);
			obj.u_ref = reshape(u_ref', [], 1);

			obj.preview = preview;
			if preview
				obj.set_reference = @obj.preview_reference;
			else
				obj.set_reference = @obj.no_preview_reference;
			end

			obj.formulation = formulation;
			if formulation
				obj.set_formulation = @obj.sparse_formulation;
				obj.set_constraints = @obj.sparse_constraints;
				obj.solve = @obj.solve_sparse;
			else
				obj.set_formulation = @obj.dense_formulation;
				obj.set_constraints = @obj.dense_constraints;
				obj.solve = @obj.solve_dense;
			end

			obj.noise = noise;

			obj.debug = debug;
			if obj.debug
				obj.options = optimoptions('quadprog', 'Display', 'iter');  % verbose
			else
				obj.options = optimoptions('quadprog', 'Display', 'none');  % silent
			end
		end

		% Set reference states and inputs with preview
		function preview_reference(obj, k)
			% preview_reference
			% 	Set reference states and inputs with preview option
			%
			% Syntax
			% 	obj.preview_reference(k)
			%
			% Input Arguments
			% 	k - Current time step
			% 		integer

			obj.X_REF = obj.x_ref((k - 1) * obj.model.n + (1:obj.model.n*obj.N));
			obj.U_REF = obj.u_ref((k - 1) * obj.model.m + (1:obj.model.m*obj.N));
		end

		% Set reference states and inputs without preview
		function no_preview_reference(obj, k)
			% no_preview_reference
			% 	Set reference states and inputs without preview option
			%
			% Syntax
			% 	obj.no_preview_reference(k)
			%
			% Input Arguments
			% 	k - Current time step
			% 		integer

			obj.X_REF = kron(ones(obj.N, 1), obj.x_ref((k - 1) * obj.model.n + (1:obj.model.n)));
			obj.U_REF = kron(ones(obj.N, 1), obj.u_ref((k - 1) * obj.model.m + (1:obj.model.m)));
		end

		% Set MPC dense formulation matrices
		function [A, B, D] = dense_formulation(obj, x_last, u_pred)
			% dense_formulation
			% 	Set MPC dense formulation matrices
			%
			% Syntax
			% 	[A, B, D] = obj.dense_formulation(x_last, u_pred)
			%
			% Input Arguments
			% 	x_last - Last state
			% 		real vector
			% 	u_pred - Prediction inputs
			% 		real vector
			%
			% Output Arguments
			% 	A - Big A formulation matrix
			% 		real matrix
			% 	B - Big B formulation matrix
			% 		real matrix
			% 	D - Big D formulation matrix
			% 		real matrix

			% Define number of states, inputs and horizon
			n = obj.model.n;
			m = obj.model.m;
			N = obj.N;

			% Define formulation matrices
			A = zeros(n * N, n);
			B = zeros(n * N, m * N);
			D = zeros(n * N, 1);

			% Initialize temporary variables to fill matrices
			A_prod = eye(n);
			D_sum = zeros(n, 1);
			A(1:n, 1:n) = eye(n);

			% Initialize the prediction states
			obj.x_pred = x_last;

			for i = 1:N - 1

				% Linearize and discretize the system
				x_bar = obj.x_pred(end - n + 1:end);               	 % linearization states 
				u_bar = u_pred((i - 1) * m + (1:m));                 % linearization inputs
				[A_lin, B_lin] = obj.model.linearize(x_bar, u_bar);  % linearization
				[A_k,   B_k]   = obj.model.discretize(A_lin, B_lin); % discretization

				% Update A
				A_prod = A_k * A_prod;
				A(n * i + (1:n), 1:n) = A_prod;

				% Update B
				mapping = [zeros(i,1); ones(N-i, 1)];
				B(:, (i-1)*m + (1:m)) = kron(mapping, B_k);
				if i > 1
					mapping = [zeros(i, i-1); ones(N-i, i-1)];
					for ii = 1:size(mapping, 1)
						for jj = 1:size(mapping, 2)
							if mapping(ii, jj) == 1
								row_start = (ii-1)*n + 1;
								row_end = ii*n;
								col_start = (jj-1)*m + 1;
								col_end = jj*m;
								B_block = B(row_start:row_end, col_start:col_end);
								B(row_start:row_end, col_start:col_end) = A_k * B_block;
							end
						end
					end
				end

				% Update the prediction states
				u_next = u_pred(i * m + (1:m));
				x_next = obj.model.simulate(obj.x_pred(end - n + 1:end), u_next, obj.Ts);
				obj.x_pred = [obj.x_pred; x_next];

				% Update D
				d_k = x_next - A_k * x_bar - B_k * u_bar;
				D_sum = A_k * D_sum + d_k;
				D(n * i + (1:n)) = D_sum;

			end
		end

		% Set MPC sparese formulation matrices
		function [A, B, D] = sparse_formulation(obj, x_last, u_pred)
			% sparse_formulation
			% 	Set MPC sparse formulation matrices
			%
			% Syntax
			% 	[A, B, D] = obj.sparse_formulation(x_last, u_pred)
			%
			% Input Arguments
			% 	x_last - Last state
			% 		real vector
			% 	u_pred - Prediction inputs
			% 		real vector
			%
			% Output Arguments
			% 	A - Big A formulation matrix
			% 		real matrix
			% 	B - Big B formulation matrix
			% 		real matrix
			% 	D - Big D formulation matrix
			% 		real matrix

			% Define number of states, inputs and horizon
			n = obj.model.n;
			m = obj.model.m;
			N = obj.N;

			% Define formulation matrices
			A = zeros(n * N, n * N);
			B = zeros(n * N, m * N);
			D = zeros(n * N, n);
			D(1:n, 1:n) = eye(n);

			% Initialize the prediction states
			obj.x_pred = x_last;

			% Initialize linearization state and input variables
			obj.X_BAR = [];
			obj.U_BAR = [];
			obj.Z_BAR = [];

			for i = 1:N - 1

				% Linearize and discretize the system
				x_bar = obj.x_pred(end - n + 1:end);               	 % linearization states
				u_bar = u_pred((i - 1) * m + (1:m));                 % linearization inputs
				obj.X_BAR = [obj.X_BAR; x_bar];
				obj.U_BAR = [obj.U_BAR; u_bar];
				[A_lin, B_lin] = obj.model.linearize(x_bar, u_bar);  % linearization
				[A_k,   B_k]   = obj.model.discretize(A_lin, B_lin); % discretization

				A(i*n + (1:n), (i-1)*n + (1:n)) = A_k;				% Update A
				B(i*n + (1:n), (i-1)*m + (1:m)) = B_k;				% Update B

				% Update the prediction states
				u_next = u_pred(i*m + (1:m));
				obj.x_pred = [obj.x_pred; obj.model.simulate(obj.x_pred(end - n + 1:end), u_next, obj.Ts)]; % simulate with nonlinear dynamics
				% obj.x_pred = A_k*x_last + B_k*u_next + x_bar - A_k*x_bar - B_k*u_bar; % simulate with linearized dynamics
			end

			% Update linearization state and input variables
			obj.X_BAR = [obj.X_BAR; obj.x_pred(end - n + 1:end)];
			obj.U_BAR = [obj.U_BAR; u_pred((N - 1)*m + (1:m))];
			obj.Z_BAR = [obj.X_BAR; obj.U_BAR];
		end

		% Set MPC constraints matrices for dense formulation
		function [H, f, EPS, F, EPS_eq, F_eq] = dense_constraints(obj, A, B, D, x_last)
			% dense_constraints
			% 	Set MPC constraints matrices for dense formulation
			%
			% Syntax
			% 	[H, f, EPS, F, EPS_eq, F_eq] = obj.dense_constraints(A, B, D, x_last)
			%
			% Input Arguments
			% 	A - Big A formulation matrix
			% 		real matrix
			% 	B - Big B formulation matrix
			% 		real matrix
			% 	D - Big D formulation matrix
			% 		real matrix
			% 	x_last - Last state
			% 		real vector
			%
			% Output Arguments
			% 	H - Quadprog formulation H matrix
			% 		real matrix
			% 	f - Quadprog formulation f vector
			% 		real vector
			% 	EPS - Inequality constraints left-hand side matrix
			% 		real matrix
			% 	F - Inequality constraints right-hand side vector
			% 		real vector
			% 	EPS_eq - Equality constraints left-hand side matrix
			% 		real matrix	
			% 	F_eq - Equality constraints right-hand side vector
			% 		real vector

			% State and input cost matrices
			Q = kron(eye(obj.N), obj.Q);
			R = kron(eye(obj.N), obj.R);

			% Input constraints matrices
			EPS_u = kron(eye(obj.N), obj.model.eps_u);
			F_u   = kron(ones(obj.N, 1), obj.model.f_u);

			% State constraints matrices
			EPS_x = kron(eye(obj.N), obj.model.eps_x);
			F_x   = kron(ones(obj.N, 1), obj.model.f_x);

			% Inequality constraints matrices
			EPS = [
				EPS_u; 
				EPS_x * B
			];
			F = [
				F_u; 
				F_x - EPS_x * A * x_last - EPS_x * D
			];

			% No equality constraints in dense formulation
			EPS_eq = [];
			F_eq = [];

			% Quadprog formulation H matrix and f vector
			H = B' * Q * B + R;
			f = x_last'*A'*Q*B + D'*Q*B - obj.X_REF'*Q*B - obj.U_REF'*R;
		end

		% Set MPC constraints matrices for sparse formulation
		function [H, f, EPS, F, EPS_eq, F_eq] = sparse_constraints(obj, A, B, D, x_last)
			% sparse_constraints
			% 	Set MPC constraints matrices for sparse formulation
			%
			% Syntax
			% 	[H, f, EPS, F, EPS_eq, F_eq] = obj.sparse_constraints(A, B, D, x_last)
			%
			% Input Arguments
			% 	A - Big A formulation matrix
			% 		real matrix
			% 	B - Big B formulation matrix
			% 		real matrix
			% 	D - Big D formulation matrix
			% 		real matrix
			% 	x_last - Last state
			% 		real vector
			%
			% Output Arguments
			% 	H - Quadprog formulation H matrix
			% 		real matrix
			% 	f - Quadprog formulation f vector
			% 		real vector
			% 	EPS - Inequality constraints left-hand side matrix
			% 		real matrix
			% 	F - Inequality constraints right-hand side vector
			% 		real vector
			% 	EPS_eq - Equality constraints left-hand side matrix
			% 		real matrix
			% 	F_eq - Equality constraints right-hand side vector
			% 		real vector

			% State and input cost matrices
			Q = kron(eye(obj.N), obj.Q);
			R = kron(eye(obj.N), obj.R);

			% Input constraints matrices
			EPS_u = kron(eye(obj.N), obj.model.eps_u);
			F_u   = kron(ones(obj.N, 1), obj.model.f_u);

			% State constraints matrices
			EPS_x = kron(eye(obj.N), obj.model.eps_x);
			F_x   = kron(ones(obj.N, 1), obj.model.f_x);

			% Inequality constraints matrices
			EPS = [
				EPS_x zeros(size(EPS_x, 1), size(EPS_u, 2));
				zeros(size(EPS_u, 1), size(EPS_x, 2)) EPS_u
			];
			F = [
				F_x; 
				F_u
			];

			% Equality constraints matrices
			EPS_eq = [eye(obj.model.n*obj.N)-A, -B];
			F_eq = D*(x_last - obj.X_BAR(1:obj.model.n)) + EPS_eq*obj.Z_BAR;

			% Quadprog formulation H matrix and f vector
			H = [
				Q, zeros(obj.model.n*obj.N, obj.model.m*obj.N); 
				zeros(obj.model.m*obj.N, obj.model.n*obj.N), R
			];
			f = -obj.Z_REF'*H;
		end

		% Solve MPC optimization problem with dense formulation
		function UMPC = solve_dense(obj, H, f, EPS, F, EPS_eq, F_eq, k)
			% solve_dense
			% 	Solve MPC optimization problem with dense formulation
			%
			% Syntax
			% 	UMPC = obj.solve_dense(H, f, EPS, F, EPS_eq, F_eq, k)
			%
			% Input Arguments
			% 	H - Quadprog formulation H matrix
			% 		real matrix
			% 	f - Quadprog formulation f vector
			% 		real vector
			% 	EPS - Inequality constraints left-hand side matrix
			% 		real matrix
			% 	F - Inequality constraints right-hand side vector
			% 		real vector
			% 	EPS_eq - Equality constraints left-hand side matrix
			% 		real matrix
			% 	F_eq - Equality constraints right-hand side vector
			% 		real vector
			% 	k - Current time step
			% 		integer

			[UMPC, ~, EXITFLAG] = quadprog((H+H')/2, f, EPS, F, EPS_eq, F_eq, [], [], [], obj.options);
			if EXITFLAG == -2
				disp("Size of H:"); disp(size(H));
				disp("Size of f:"); disp(size(f));
				disp("Size of EPS:"); disp(size(EPS));
				disp("Size of F:"); disp(size(F));
				disp("Size of EPS_eq:"); disp(size(EPS_eq));
				disp("Size of F_eq:"); disp(size(F_eq));
				error('Something wrong with the optimization at iteration ' + string(k))
			end
		end

		% Solve MPC optimization problem with sparse formulation
		function UMPC = solve_sparse(obj, H, f, EPS, F, EPS_eq, F_eq, k)
			% solve_sparse
			% 	Solve MPC optimization problem with sparse formulation
			%
			% Syntax
			% 	UMPC = obj.solve_dense(H, f, EPS, F, EPS_eq, F_eq, k)
			%
			% Input Arguments
			% 	H - Quadprog formulation H matrix
			% 		real matrix
			% 	f - Quadprog formulation f vector
			% 		real vector
			% 	EPS - Inequality constraints left-hand side matrix
			% 		real matrix
			% 	F - Inequality constraints right-hand side vector
			% 		real vector
			% 	EPS_eq - Equality constraints left-hand side matrix
			% 		real matrix
			% 	F_eq - Equality constraints right-hand side vector
			% 		real vector
			% 	k - Current time step
			% 		integer

			[ZMPC, ~, EXITFLAG] = quadprog((H+H')/2, f, EPS, F, EPS_eq, F_eq, [], [], [], obj.options);
			if EXITFLAG == -2
				disp("Iteration: " + string(k));
				disp("Size of H:"); disp(size(H));
				disp("Size of f:"); disp(size(f));
				disp("Size of EPS:"); disp(size(EPS));
				disp("Size of F:"); disp(size(F));
				disp("Size of EPS_eq:"); disp(size(EPS_eq));
				disp("Size of F_eq:"); disp(size(F_eq));
				error('Something wrong with the optimization at iteration ' + string(k))
			end
			UMPC = ZMPC(obj.model.n*obj.N+1:end);
		end

		% MPC optimization function
		function [x, u] = optimize(obj)
			% optimize
			% 	MPC optimization function
			%
			% Syntax
			% 	[x, u] = obj.optimize()
			%
			% Output Arguments
			% 	x - State vector
			% 		real matrix
			% 	u - Input vector
			% 		real matrix		

			% Define number of states, inputs and horizon
			n = obj.model.n;
			m = obj.model.m;

			% Initialize the state and input vectors
			x = obj.x0;
			if obj.noise
				% Noise components generation for Kalman filter
				w = mvnrnd(zeros(obj.model.n, 1), obj.model.Q_tilde, obj.Nsteps);
				v = mvnrnd(zeros(obj.model.p, 1), obj.model.R_tilde, obj.Nsteps);
				x_last = obj.x0 + w(1, :)';
			else
				x_last = obj.x0;
			end
			u = [];
			obj.set_reference(1);
			u_pred = obj.U_REF;

			% Initialize Mean Squared Error (MSE) metrics
			MSE_x = 0;
			MSE_u = 0;
			wMSE_x = 0;
			wMSE_u = 0;

			% Solve optimization problem at each time step
			for k = 1:obj.Nsteps

				% Set reference states and inputs with or without preview
                obj.set_reference(k);

				% Define formulation matrices
				[A, B, D] = obj.set_formulation(x_last, u_pred);

				% Fix angle representation in X_REF
				obj.X_REF = obj.model.fix_angles(obj.x_pred, obj.X_REF);

				% Set augmented reference state Z_REF
				obj.Z_REF = [obj.X_REF; obj.U_REF];

				% Set constraints matrices
				[H, f, EPS, F, EPS_eq, F_eq] = obj.set_constraints(A, B, D, x_last);

				% Solve the optimization problem and get optimal input sequence
				UMPC = obj.solve(H, f, EPS, F, EPS_eq, F_eq);

				% Extract the optimal input and update the input vector
				u_optimal = UMPC(1:m); % only pick the first inputs (of size m)
				u = [u; u_optimal];

				% Update the state vector simulating the system with the optimal input
				if obj.noise
					update_last = obj.model.simulate(x_last, u_optimal, obj.Ts) + w(k,:)';
					measured_output = obj.model.output(update_last, u_optimal) + v(k,:)';
					x_last = obj.model.EKF_estimate(x_last, u_optimal, measured_output);
				else
					x_last = obj.model.simulate(x_last, u_optimal, obj.Ts);
				end
				x = [x; x_last];

				% Update the prediction inputs with the remaining optimization results
				u_remaining = UMPC(m+1:end);
				u_pred = [u_remaining; u_remaining(end-m+1:end)]; % duplicate last inputs

				% Compute distances from reference
				x_distance = norm(x_last - obj.x_ref((k) * n + (1:n)));
				u_distance = norm(u_optimal - obj.u_ref((k) * m + (1:m)));
				weighted_x_distance = (x_last - obj.x_ref((k) * n + (1:n)))' * obj.Q * (x_last - obj.x_ref((k) * n + (1:n)));
				weighted_u_distance = (u_optimal - obj.u_ref((k) * m + (1:m)))' * obj.R * (u_optimal - obj.u_ref((k) * m + (1:m)));

				% Update Mean Squared Error (MSE) metrics
				MSE_x = MSE_x + x_distance^2;
				MSE_u = MSE_u + u_distance^2;
				wMSE_x = wMSE_x + weighted_x_distance;
				wMSE_u = wMSE_u + weighted_u_distance;

				% Display the current iteration results
				fprintf("Iteration: %d/%d, ||x - x_ref|| = %f, ||u - u_ref|| = %f\n", k, obj.Nsteps, x_distance, u_distance);
			end

			% Display the final Mean Squared Error (MSE) metrics
			MSE_x = MSE_x / obj.Nsteps;
			MSE_u = MSE_u / obj.Nsteps;
			wMSE_x = wMSE_x / obj.Nsteps;
			wMSE_u = wMSE_u / obj.Nsteps;
			fprintf("\nFinal Mean Squared Error (MSE) metrics:\n");
			fprintf("         State MSE = %f,          Input MSE = %f\n", MSE_x, MSE_u);
			fprintf("Weighted State MSE = %f, Weighted Input MSE = %f\n", wMSE_x, wMSE_u);

			
			% ------------------------------------------------------------------------------- cut

			% Just for results
            % fprintf("\nhelicopter,lemniscate,%d,%f,%f", obj.N, MSE_x, MSE_u);
            % fprintf("\nunicycle,lemniscate,10,%f,%f", MSE_x, MSE_u);
            % fprintf("\nhelicopter,circle,18,%f,%f",MSE_x, MSE_u);
			% fprintf("\nhelicopter,lemniscate,18,%f,%f",MSE_x,MSE_u);

			% ------------------------------------------------------------------------------- cut

			% Reshape the state and input vectors
			x= (reshape(x, n, length(x) / n))';
			u= (reshape(u, m, length(u) / m))';

		end
	end
end