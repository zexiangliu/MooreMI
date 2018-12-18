A = [3, -8; -20, 10];
B = [1;1];
dyn = @(x,u) A*x + B*u;

K = @(x,u) norm(A,2);

% state space [-1,1]x[-1,1], the initial partition has four grids
P1 = [-1 0;  % [lb, ub] of dim 1
      -1 0]; % [lb, ub] of dim 2
P2 = [-1 0;
      0 1];
P3 = [0 1;
      0 1];
P4 = [0 1;
      -1 0];
G = {P1, P2, P3, P4};
label = ["a","b","c","d"];

U = ["a"];

% spec []A ^ <>[]B ^ (AND_i []<> (Ri))
spec.A = [1,2,3,4];
spec.B = [];
spec.R = {};

abs = Abstraction(U, G, spec, label, dyn, K);