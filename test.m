
n = 3;
m = 2;

A1 = ones(3)-eye(3);
A2 = eye(3);

A = {A1,A2};
Q_name = ["a","b","c"];
U_name = ["1","2"];
Q_label = ["01","02","03"];

G1 = DFA(3,2,A,[],[],[],Q_name,U_name,Q_label);

% test pre
[pre_x,pre_u] = unique(G1.pre(2))
pre_x = G1.pre_xu(2,1)
pre_x = G1.pre_xu("a","1")

% test post
[post_x,post_u] = unique(G1.post(2))
post_x = G1.post_xu(2,2)
post_x = G1.post_xu("a","1")

