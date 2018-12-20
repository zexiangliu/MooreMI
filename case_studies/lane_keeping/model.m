clc;clear all;close all;
% parameters
param.f0 = 0.1;
param.f1 = 10;%5;
param.f2 = 0.8;%0.25;
param.m = 2000;
param.g = 10;

X = [0, 100
     -50, 50];
F = linspace(-3700,3700,4);
v_y = linspace(-5,5,10);

U = {};
for i = 1:length(F)
    for j = 1:length(v_y)
        U{end+1} = [F(i),v_y(j)];
    end
end

Target_region = [50, 60;
                 -5, 5];
Safe_region = [0,100;
               -15,15];
left_parts1 = Abstraction.surface_diff(X,Safe_region);
left_parts2 = Abstraction.surface_diff(Safe_region,Target_region);
G = {Target_region};

G(end+1:end+length(left_parts2)) = left_parts2;
G(end+1:end+length(left_parts1)) = left_parts1;

label = 1:length(G);
spec.A = [1:length(left_parts2)+1];
spec.B = [1];
spec.C = {};

abs_model = Abstraction(X, U, G, spec, label, @(x,u) dyn(x,u,param)...
    , @(x,u,r) Lipschitz(x,u,r,param),1e-5);

%%
fig = figure(1);
abs_model.phase_portrait(fig, 20);

abs_model.verifyTransition();
abs_model.to_TransSyst();
[W,CW,~] = abs_model.win_primal();

fig2 = figure(2);
abs_model.plot(fig2,0,'k');
fig3 = figure(3);
abs_model.plot(fig3,0,'k');

%% start refinement

W_old = W;
% mapg = {};
% 
% u1 = sub2ind([10,4],10,4);
% u2 = sub2ind([10,4],1,4);
% mapg{1} = {1,2,[u1,u2]};
% 
% u1 = sub2ind([10,4],10,1);
% u2 = sub2ind([10,4],1,1);
% mapg{2} = {1,1,[u1,u2]};

counter = 0;
while(~isempty(CW) || isempty(W))
    grid_vol = zeros(length(CW),1);
    
    for i =1:length(CW)
        PG = abs_model.G{CW(i)};
        grid_vol(i) = prod(PG(:,2)-PG(:,1));
    end
    % pick the one with largest vol to refine
    [~,idx] = max(grid_vol);
    abs_model.refinement(CW(idx(1)));
    W_old = W;
%     abs_model.to_TransSyst(true,mapg);
    abs_model.to_TransSyst();
    [W,CW,cont] = abs_model.win_primal();
    
    abs_model.plot(fig2,CW,'r');
    if ~isempty(W)
        abs_model.plot(fig2,W,'g');
    end
    abs_model.plot(fig3,[],'k');
    drawnow;
    W
    CW
    counter = counter + 1;
    disp("Iteration "+num2str(counter));
    if counter > 20
        break;
    end
end



