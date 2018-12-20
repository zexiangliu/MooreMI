clc;clear all;close all;
X = [0, 100
     0, 80];
U = [0,1];
Target_region = [0.9, 1.8;
                 40, 40.2];
left_parts = Abstraction.surface_diff(X,Target_region);
G = {Target_region};
G(end+1:end+length(left_parts)) = left_parts;
label = {"a","b","b","b","b","b","b","b","b"};
spec.A = [];
spec.B = [];
spec.C = {1};

abs_model = Abstraction(X, U, G, spec, label, @dyn, @Lipschitz,1e-5);
fig = figure;
abs_model.phase_portrait(fig, 20);

abs_model.verifyTransition();
abs_model.to_TransSyst();
[W,CW,~] = abs_model.win_primal();
W