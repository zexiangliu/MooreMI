function dx = dyn(x,u,param)
% 2d vehicle dynamics
% v_x ---forward speed
% dv_x --- forward acceleration
% y --- lateral displacement
% v_y --- lateral speed
% states: [v_x;y]

f0 = param.f0;
f1 = param.f1;
f2 = param.f2;
m = param.m;

F = u{1}(1);
v_y = u{1}(2);
v_x = x(1);

% saturation constraints on v_y
if abs(v_x) <= 2
    v_y = sign(v_y)*min(abs(v_y),0.5*v_x);
else
    v_y = sign(v_y)*min(abs(v_y),0.5*2+0.2*v_x);
end

v = v_x^2 + v_y^2;
dv_x = 1/m*(F - f0 - f1*v_x - f2*v_x*sqrt(v));
dx = [dv_x;v_y];
end