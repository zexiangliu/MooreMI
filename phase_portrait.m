A = [3, -8; -20, 10];
f = @(t,Y) A*Y;
y1 = linspace(-1,1,20);
y2 = linspace(-1,1,20);

% creates two matrices one for all the x-values on the grid, and one for
% all the y-values on the grid. Note that x and y are matrices of the same
% size and shape, in this case 20 rows and 20 columns

[x,y] = meshgrid(y1,y2);
size(x);
size(y);

u = zeros(size(x));
v = zeros(size(x));

% we can use a single loop over each element to compute the derivatives at
% each point (y1, y2)
t=0; % we want the derivatives at each point at t=0, i.e. the starting time
for i = 1:numel(x)
    Yprime = f(t,[x(i); y(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
    l = norm([u(i);v(i)],2);
    u(i) = u(i)/l;
    v(i) = v(i)/l;
end
quiver(x,y,u,v,'k'); figure(gcf)
% xlabel('y_1')
% ylabel('y_2')
xticks([-1 -0.5 0.5 1])
yticks([-1 -0.5 0.5 1])
axis tight equal;
% box off;
% axis off
