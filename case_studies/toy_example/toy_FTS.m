 A = [2, 5; -8, -2];
B = [0;0];
dyn1 = @(x,u) A*x + B*u;

K1 = @(x,u,r) norm(A,'inf');

% state space [-1,1]x[-1,1], the initial partition has four grids
X = [-1 1; -1 1];
P3 = [-1 0;  % [lb, ub] of dim 1
      -1 0]; % [lb, ub] of dim 2
P1 = [-1 0;
      0 1];
P2 = [0 1;
      0 1];
P4 = [0 1;
      -1 0];
G1 = {P1, P2, P3, P4};
label1 = ["1","2","3","4"];

U1 = [0];

% spec []A ^ <>[]B ^ (AND_i []<> (Ri))
spec1.A = [1,2,3,4];
spec1.B = [];
spec1.C = {};

abs_model = Abstraction(X, U1, G1, spec1, label1, dyn1, K1);

abs_model.verifyTransition()