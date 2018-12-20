function K = Lipschitz(x,u,range,param)

f1 = param.f1;
f2 = param.f2;
m = param.m;

vx_max = max(abs(range(1,:)));
v_y = u{1}(2);
    
dv_max = sqrt(vx_max^2+v_y^2) + vx_max^2/sqrt(vx_max^2+v_y^2);

if vx_max == 0 && v_y == 0
    dv_max = 0;
end

K = abs(f1)/m + abs(f2)/m*dv_max;
end