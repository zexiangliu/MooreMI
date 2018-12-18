% this function is to compute the transition system of the continuous
% system.
% example A=[3, -8; -20, 10]; or nonlinear system with [10*Y(2); -10*Y(1) + 100*Y(2)^2]

global L A
A = [3, -8; -20, 10];
L = norm(A,2);
%L = 1000;
% X space is an example
X1 = [-1, 1];
X2 = [-1, 1];
% G is a cell of vertices of each grid, vertices arrangement is as following
% v1--v2
% |    |
% v3--v4
G = {{[-1;1],[-0.5;1],[-1;0.5],[-0.5;0.5]}, {[-0.5;1], [0.5;1], [-0.5;0.5], [0.5;0.5]},{[0.5;1],[1;1],[0.5;0.5],[1;0.5]},...
    {[-1;0.5],[-0.5;0.5],[-1;-0.5],[-0.5;-0.5]}, {[-0.5;0.5],[0.5;0.5],[-0.5;-0.5],[0.5;-0.5]},{[0.5;0.5],[1;0.5],[0.5;-0.5],[1;-0.5]},...
    {[-1;-0.5],[-0.5;-0.5],[-1;-1],[-0.5;-1]}, {[-0.5;-0.5],[0.5;-0.5],[-0.5;-1],[0.5;-1]}, {[0.5;-0.5],[1;-0.5],[0.5;-1],[1;-1]}};

% representing the connection of vertices, C{i,j}=1 menas grid(i)-> grid(j)
C = eye(length(G));

for i = 1:length(G)
    for j = 1:i-1
        % if grid(i)=grid(j) are neighbours, then do transition function
        [a, P] = if_neighbor(G{i},G{j});
        if a == 2
            trans = transition(P{1}, P{2});
            pin = prod(mean(cell2mat(G{i}),2)>=mean(cell2mat(G{j}),2));
            % pin is a flag for the statement G{i}>G{j}
            i_big = i*pin + j*(1-pin);
            i_sml = i*(1-pin) + j*pin;
            C(i_big,i_sml) = ismember(-1, trans);
            C(i_sml,i_big) = ismember(1, trans);
        elseif a==1
            % do vert_transition function, which takes two grid, and
            % returns either trans=0 meaning no transition between the two,
            % or trans=1 meaning flow from higher grid to lower grid, or
            % trans=-1 meaning flow from lower grid to higher grid
            trans = vert_transition(P, G{i}, G{j});
            mid1 = mean(cell2mat(G{i}),2); mid2 = mean(cell2mat(G{j}),2);
            pin = mid1(2)>mid2(2); % pin is flag for G{i} is a higher grid
            i_hgh = i*pin + j*(1-pin);
            i_low = i*(1-pin) + j*pin;
            C(i_low,i_hgh) = (trans==1);
            C(i_hgh,i_low) = (trans==-1);
        end
    end
end


% this function compute if the two grid are neighbours
function [a, P] = if_neighbor(G1, G2)
    % a=2 means two cell have common edge, a=1 means two cell have common
    % vertices, a=0 means two cell are not neighbors
    G1 = cell2mat(G1)';
    G2 = cell2mat(G2)';
    P = num2cell(intersect(G1, G2,'rows')', 1);
    a = length(P);
end


% this function aims to compute the transition between the two grid, and
% takes two common vertices as input, no order requirement.
function trans = transition(p1, p2)
% output: trans=1 means flow to the positive direction
%         trans=-1 means flow to the negative direction
%         trans=[1,-1] means flow is bi-directional
global list T
T = [];
% change p1 and p2 position
while prod(p1>=p2)
    p_hold = p2;
    p2 = p1;
    p1 = p_hold;
end
list = {{p1, p2}};
while length(list)&&(~(ismember(1, T)&&ismember(-1,T)))
    [t, c1, c2, finish] = edge(list{1}{1}, list{1}{2});
    T = [T, t];
end
trans = unique(T);
end


% this function aims to compute the cover size, flow direction of that
% cover, and returns the uncovered region
% this function is a subfunction of transition function
function [t, c1, c2, finish] = edge(p1, p2)
    global L list
    % p1 = [x1; x2]; p2 = [y1; y2];
    % c1 = co-ordinate of one side of the cover; c2 defined similarly
    % example: p1-----c1---mid---c2-----p2
    list = list(2:end);
    mid = (p1+p2)/2;
    n = double([p1(1)==p2(1); p1(2)==p2(2)]); % normal vector
    v = dynamic(mid);
    r = abs(n'*v)/(norm(n, 2)*L); % radius
    finish = norm(mid - p1, 2) <= r; % flag, finish=1 means covers the whole interval
    if ~finish
        c1 = mid - [p1(1)~= p2(1); p1(2)~=p2(2)]*r; % compute coordinate of c1
        c2 = mid + [p1(1)~= p2(1); p1(2)~=p2(2)]*r; % compute coordinate of c2
        % add (p1,c1), (c2,p2) to list
        list{end+1} = {p1,c1};
        list{end+1} = {c2,p2};
    else
        c1 = p1;
        c2 = p2;
    end
    t = sign(n'*v);
end

% this function aims to compute the transition when grids have one common
% vertice
% trans=0 means no flow between the two grid;
% trans=1 means flow from lower grid to higher grid
% trans=-1 means flow from higher grid to lower grid
function trans = vert_transition(p, G1, G2)
    p = cell2mat(p);
    mid1 = mean(cell2mat(G1),2);
    mid2 = mean(cell2mat(G2),2);
    pos = prod(sign(mid1-mid2)); 
    % pos=1 means 0..1| pos=-1 means 1..0
    %             1..0|              0..1
    n1 = [1;0]; n2 = [0;1];
    v = dynamic(p);
    flow = sign(v'*n1)*sign(v'*n2); % flow=0 means at p, flow is horizontal or vertical
    % flow=1 means flow is in "/" direction, flow=-1 means flow is in "\" direction.
    trans = 0;
    if flow==pos
        trans = sign(v'*n2);
    end
end


% this is the query step
function v = dynamic(Y)
    global A
    %v = [10*Y(2); -10*Y(1) + 100*Y(2)^2];
    v = A*Y;
end