% Model Predictive Control (MPC) Optimization Class
classdef MPC
	properties
		Q;          % state cost
		R;          % input cost
		N;          % prediction horizon
		model;      % system model
		x_ref;      % reference states
		u_ref;      % reference inputs
		X_REF;	  	% reference states horizon vector
		U_REF;	  	% reference inputs horizon vector
		Z_REF;	  	% reference states and inputs horizon vector
		X_BAR;	  	% linearization states horizon vector
		U_BAR;	  	% linearization inputs horizon vector
		Z_BAR;	  	% linearization states and inputs horizon vector
		x0;         % initial state
		Ts;         % sampling time
		t;          % simulation time vector
		T;          % number of simulations steps
		preview;    % MPC preview flag (0: no preview, 1: preview)
		formulation;% formulation to use for MPC problem (0: dense (explicit), 1: sparse (implicit))
		options;    % quadprog optimization options
	end
	
	methods
		% Constructor to initialize the MPC problem parameters
		function obj = MPC(Q, R, N, model, x_ref, u_ref, x0, Ts, t, T, preview, formulation, options)
			obj.Q = Q;
			obj.R = R;
			obj.N = N;
			obj.model = model;
			obj.x_ref = reshape(x_ref', [], 1);
			obj.u_ref = reshape(u_ref', [], 1);
			obj.x0 = x0;
			obj.Ts = Ts;
			obj.t = t;
			obj.T = length(t)-(N+1);
			obj.preview = preview;
			obj.formulation = formulation;
			obj.options = options;
		end

		% Dense formulation matrices function
		function [A, B, D] = dense_formulation(obj, x_last, u_pred)

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
			x_pred = x_last;

			for i = 1:N - 1
				% Linearize and discretize the system
				x_bar = x_pred(end - n + 1:end);                     % linearization states 
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
				x_next = obj.model.simulate(x_pred(end - n + 1:end), u_next, obj.Ts);
				x_pred = [x_pred; x_next];

				% Update D
				d_k = x_next - A_k * x_bar - B_k * u_bar;
				D_sum = A_k * D_sum + d_k;
				D(n * i + (1:n)) = D_sum;

			end
		end

		% Sparse formulation matrices function
		function [obj, A, B, D] = sparse_formulation(obj, x_last, u_pred)

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
			x_pred = x_last;

			% Initialize linearization state and input variables
			obj.X_BAR = [];
			obj.U_BAR = [];
			obj.Z_BAR = [];

			for i = 1:N - 1
				x_bar = x_pred(end - n + 1:end);                     % linearization states
				u_bar = u_pred((i - 1) * m + (1:m));                 % linearization inputs
				obj.X_BAR = [obj.X_BAR; x_bar];
				obj.U_BAR = [obj.U_BAR; u_bar];

				[A_lin, B_lin] = obj.model.linearize(x_bar, u_bar);  % linearization
				[A_k,   B_k]   = obj.model.discretize(A_lin, B_lin); % discretization

				A(i*n + (1:n), (i-1)*n + (1:n)) = A_k;				% Update A
				B(i*n + (1:n), (i-1)*m + (1:m)) = B_k;				% Update B

				% Update the prediction states
				u_next = u_pred(i*m + (1:m));
				x_pred = [x_pred; obj.model.simulate(x_pred(end - n + 1:end), u_next, obj.Ts)]; % simulate with nonlinear dynamics
				% x_pred = A_k*x_last + B_k*u_next + x_bar - A_k*x_bar - B_k*u_bar; % simulate with linearized dynamics
			end

			% Update linearization state and input variables
			obj.X_BAR = [obj.X_BAR; x_pred(end - n + 1:end)];
			obj.U_BAR = [obj.U_BAR; u_pred((N - 1)*m + (1:m))];
			obj.Z_BAR = [obj.X_BAR; obj.U_BAR];
		end

		% MPC optimization function
		function [x, u] = optimize(obj)

			% Define number of states, inputs and horizon
			n = obj.model.n;
			m = obj.model.m;
			N = obj.N;

			% Initialize the state and input vectors
			x = obj.x0;
			x_last = obj.x0;
			u = [];

			% Solve optimization problem at each time step
			for k = 1:obj.T

				% Handle preview flag
				if obj.preview % Preview enabled: use next N reference states and inputs
					obj.X_REF = obj.x_ref((k - 1) * n + (1:n*N));
					obj.U_REF = obj.u_ref((k - 1) * m + (1:m*N));
					obj.Z_REF = [obj.X_REF; obj.U_REF];
				else % Preview disabled: use same reference state and input N times
					obj.X_REF = kron(ones(N, 1), obj.x_ref((k - 1) * n + (1:n)));
					obj.U_REF = kron(ones(N, 1), obj.u_ref((k - 1) * m + (1:m)));
					obj.Z_REF = [obj.X_REF; obj.U_REF];
				end

				% Select inputs for state prediction and linearization
				if k==1 % First iteration: use reference inputs
					u_pred = obj.U_REF;
				else % Remaining iterations: use last optimization remaining inputs
					u_pred = u_remaining;
				end

				% Define formulation matrices
				if obj.formulation == 0
					[A, B, D] = obj.dense_formulation(x_last, u_pred);
				else
					[obj, A, B, D] = obj.sparse_formulation(x_last, u_pred);
				end

				% State and input cost matrices
				Q = kron(eye(N), obj.Q);
				R = kron(eye(N), obj.R);

				% Cost Function and Constraints

				% Input constraints
				EPS_u = kron(eye(N), obj.model.eps_u);
				F_u   = kron(ones(N, 1), obj.model.f_u);

				% State constraints
				EPS_x = kron(eye(N), obj.model.eps_x);
				F_x   = kron(ones(N, 1), obj.model.f_x);

				% Cost Function and Constraint matrices (formulation dependent)
				if obj.formulation == 0 % Dense (explicit) formulation

					% Inequality constraints matrices
					EPS = [
						EPS_u; 
						EPS_x * B
					];
					F = [
						F_u; 
						F_x - EPS_x * A * x_last - EPS_x * D
					];

					% NOTE: no equality constraints in this case
					EPS_eq = [];
					F_eq = [];

					% Quadprog formulation H matrix and f vector
					H = B' * Q * B + R;
					f = x_last'*A'*Q*B + D'*Q*B - obj.X_REF'*Q*B - obj.U_REF'*R;

				else % Sparse (implicit) formulation

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
					EPS_eq = [eye(n*N)-A, -B];
					F_eq = D*(x_last - obj.X_BAR(1:n)) + EPS_eq*obj.Z_BAR;

					% Quadprog formulation H matrix and f vector
					H = [
						Q, zeros(n*N, m*N); 
						zeros(m*N, n*N), R
					];
					f = -obj.Z_REF'*H;

				end

				% Solve the optimization problem
				[OUT_MPC, FVAL, EXITFLAG] = quadprog((H + H')/2, f, EPS, F, EPS_eq, F_eq, [], [], [], obj.options);
				if EXITFLAG == -2
					disp('Something wrong with the optimization at iteration ' + string(k))
					keyboard
				end
				
				% Extract the optimal input sequence
				if obj.formulation == 0
					UMPC = OUT_MPC;
				else
					UMPC = OUT_MPC(n*N+1:end);
				end

				% Extract the optimal input and update state and input vectors
				u_optimal = UMPC(1:m); % only pick the first inputs (of size m)
				u = [u; u_optimal];
				x_last = obj.model.simulate(x_last, u_optimal, obj.Ts);
				x = [x; x_last];

				% Save the remaining inputs for the next MPC optimization
				u_remaining = UMPC(m+1:end);
				u_remaining = [u_remaining; u_remaining(end-m+1:end)]; % duplicate last inputs

				% Print the current iteration results
				disp('Iteration: ' + string(k) + '/' + string(obj.T) + ' --------------');
				disp('Current state vs reference: ');
				disp(x_last');
				disp(obj.x_ref((k) * n + (1:n))');
				x_distance = norm(x_last - obj.x_ref((k) * n + (1:n)));
				disp('    ||x - x_ref|| = ' + string(x_distance));
				disp(' ');
				disp('Current input vs reference: ');
				disp(u_optimal');
				disp(obj.u_ref((k) * m + (1:m))');
				u_distance = norm(u_optimal - obj.u_ref((k) * m + (1:m)));
				disp('    ||u - u_ref|| =  ' + string(u_distance));
				disp(' ');

			end

			% Reshape the state and input vectors
			x= (reshape(x, n, length(x) / n))';
			u= (reshape(u, m, length(u) / m))';

		end
	end
end