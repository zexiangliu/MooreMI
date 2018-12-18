classdef Abstraction < handle
    properties (SetAccess = protected)
        n; % number of grids
        dim; % dimension of state space
        m; % number of actions
        G; % partitions, cell of grids
        X; % the state space, e.g. [-1,1; -1,1]
        U; % array of actions
        A; % adjacent matrix of the grids
        spec; % the verification/synthesis specification
        label; % labeling of grids
        ts; % transition system corresponding to the partitions
        Win; % winning set
        dyn; % function handle, dynamics of the system, take state ...
             % and action as inputs
        K; % Lipschitz constant, a function of state and action
        mini_width; % the minimal width of grid
    end
    
    methods
        function obj = Abstraction(X, U, G, spec, label, dyn, K, mini_width)
            if nargin == 7
                mini_width = 1e-6;
            end
            obj.X = X;
            obj.U = U;
            obj.G = G;
            obj.dim = size(G{1},1);
            obj.n = length(G);
            obj.m = length(U);
            obj.spec = spec;
            obj.label = label;
            obj.dyn = dyn;
            obj.K = K;
            obj.initialize_A();
            obj.ts = [];
            obj.Win = [];
            obj.mini_width = mini_width;
        end
        
        function obj = verifyTransition(obj, mini_width)
            if nargin == 1
                mini_width = obj.mini_width;
            end
            for idx_u = 1:obj.m
                for G1 = 1:obj.n
                    for G2 = 1:G1-1
                        % verify if there is any transition between i and j
                        if obj.A{idx_u}(G1,G2) == 1 || obj.A{idx_u}(G2,G1) == 1
                            [~, ~, cmn_surface] = has_cmn_surface(...
                                                            obj, G1, G2);
                            direct = is_transit(obj, G1, G2, ...
                                cmn_surface,idx_u,mini_width);

                            switch direct
                                case 0
                                    obj.A{idx_u}(G1,G2) = 1;
                                    obj.A{idx_u}(G2,G1) = 1;
                                case 1
                                    obj.A{idx_u}(G1,G2) = 1;
                                    obj.A{idx_u}(G2,G1) = 0;
                                case 2
                                    obj.A{idx_u}(G1,G2) = 0;
                                    obj.A{idx_u}(G2,G1) = 1;
                                case 3
                                    obj.A{idx_u}(G1,G2) = 0;
                                    obj.A{idx_u}(G2,G1) = 0;
                            end
                        end
                    end
                end
            end
        end
        
        % verify the direction of the flow on the common surface.
        function [direct, list_cover, cover_direct] =...
                is_transit(obj, G1, G2, cmn_surface, idx_u, mini_width)
        % input: cmn_surface --- output from has_cmn_surface(obj,G1,G2)
        %        u --- action
        % output: direct --- 0: indeterminate or i <---> j
        %                    1: i ---> j
        %                    2: j ---> i
        %                    3: no transition i <--x--> j
        %         list_cover --- cell of covers of the surface
        %         cover_direct --- direction of covers in the list_cover,
        %                        the interpretation is the same as direct.
            if nargin == 5
                mini_width = obj.mini_width;
            end
            
            list_cover = {};
            cover_direct = {};
            flag_cover = nargout > 1;
            u = obj.U(idx_u);
            % compute the normal vector of the common surface from i to j
            normal = double(cmn_surface(:,2)==cmn_surface(:,1));
            center_G1 = mean(obj.G{G1},2);
            center_G2 = mean(obj.G{G2},2);
            normal = sign(center_G2-center_G1).*normal;
            
            queue = {cmn_surface};
            flag_12 = false;
            flag_21 = false;
            while ~isempty(queue)
                % queue pop
                surface = queue{1};
                queue(1) = [];
                % compute the cover and decide flow direction
                x = mean(surface,2);
                [flow, cover] = obj.compute_cover(x,u,normal);
                norm_flow = flow.*normal;
                norm_flow = norm_flow(normal~=0);
               
                % if surface and cover are too small, skip
                if sum(normal~=0)~= obj.dim... % see if a point or not
                        && max(surface(:,2) - surface(:,1)) <= mini_width ...
                        && any(surface(:,1) < cover(:,1) & ...
                        surface(:,2) > cover(:,2)) 
                    list_cover{end+1} = surface;
                    cover_direct{end+1} = 0; % indeterminate
                    flag_12 = true;
                    flag_21 = true;
                    if flag_cover
                        continue;
                    else
                        break;
                    end
                end
                
                % if need list_covers as output
                if flag_cover
                    list_cover{end+1} = cover;
                    if all(norm_flow > 0)
                        cover_direct{end+1} = 1;
                    elseif all(norm_flow < 0)
                        cover_direct{end+1} = 2;
                    else
                        cover_direct{end+1} = 3;
                    end
                end
                
                if ~flag_12 && all(norm_flow > 0)
                    flag_12 = true;
                elseif ~flag_21 && all(norm_flow < 0)
                    flag_21 = true;
                end
                
                residual = obj.surface_diff(surface,cover);
                % queue push
                queue(end+1:end+length(residual)) = residual;
                
                % if do not need list_covers, can be stopped
                if ~flag_cover && flag_12 && flag_21
                    break;
                end
            end
            
            if flag_12 && flag_21
                direct = 0;
            elseif flag_12
                direct = 1;
            elseif flag_21
                direct = 2;
            else
                direct = 3;
            end
        end
        
        % compute the radius of the cover
        function [flow, cover] = compute_cover(obj,x,u,normal)
        % inputs: x --- state you want to query
        %         u --- action applied
        %         normal --- normal vector of the common surface
            flow = obj.dyn(x,u);
            r = min(abs(flow(normal~=0)))/obj.K(x,u);
            cover = [x-ones(obj.dim,1)*r,x+ones(obj.dim,1)*r];
        end
        
        % compute residuals after one iteration  
        function residual = surface_diff(obj,surface, cover)
        % compute the difference (surface - cover) and represent the 
        % residual using rectangles
            list_interval = cell(obj.dim,1);
            
            % the idx of cover in the new partition, needs to be removed
            idx_cover = cell(obj.dim,1);
            % dim of list_interval
            list_dim = zeros(obj.dim,1);
            
            for i = 1:obj.dim
                surf_i = surface(i,:);
                cover_i = cover(i,:);
                
                if cover_i(1)>surf_i(1) && cover_i(2)<surf_i(2)
                    % ---xxx---
                    int1 = [surf_i(1),cover_i(1)];
                    int2 = [cover_i(1),cover_i(2)];
                    int3 = [cover_i(2),surf_i(2)];
                    list_interval{i} = {int1,int2,int3};
                    idx_cover{i} = 2;
                elseif surf_i(1)>=cover_i(1) && ...
                        surf_i(2)<=cover_i(2)
                    list_interval{i} = {surf_i};
                    idx_cover{i} = 1;
                elseif surf_i(1)>=cover_i(1) && ...
                        surf_i(1)<cover_i(2)
                    % xxx---- 
                    int1 = [surf_i(1),cover_i(2)];
                    int2 = [cover_i(2),surf_i(2)];
                    list_interval{i} = {int1,int2};
                    idx_cover{i} = 1;
                elseif surf_i(2)>cover_i(1) && ...
                        surf_i(2)<=cover_i(2)
                    % ----xxx
                    int1 = [surf_i(1),cover_i(1)];
                    int2 = [cover_i(1),surf_i(2)];
                    list_interval{i} = {int1,int2};
                    idx_cover{i} = 2;
                else
                    error("surface and cover does not intersect!");
                end
                list_dim(i) = length(list_interval{i});
            end
            
            total = prod(list_dim);
            residual = cell(total-1,1);
            for i = 1:total
                sub = ind2sub2(list_dim,i);
                surf = zeros(obj.dim,2);
                for j = 1:obj.dim
                    surf(j,:) = list_interval{j}{sub(j)};
                end
                residual{i} = surf;
            end
            idx_cover = sub2ind2(list_dim,idx_cover);
            residual(idx_cover) = [];
        end
        % check if the two grids have common surface
        function [bool, cmn_dim, cmn_surface] = has_cmn_surface(obj, G1, G2)
        % output: bool --- true: has common surface
        %                  false: has no common surface
        %         cmn_dim --- num of dimension of the space
        %                  the surface lies in; cmn_dim = -1 <---> no
        %                  common surface.
        %         cmn_surface --- the common surface represented by dimx2
        %                         matrix
            PG1 = obj.G{G1};
            PG2 = obj.G{G2};
            bool = false;
            cmn_dim = -1;
            cmn_surface = [max([PG1(:,1) PG2(:,1)],[],2)...
                            min([PG1(:,2) PG2(:,2)],[],2)];
            bnd_diff = cmn_surface(:,2) - cmn_surface(:,1);
            if all(bnd_diff>=0)
                cmn_dim = sum(bnd_diff>0);
                bool = cmn_dim < obj.dim;
            end
        end
        
        % abstraction mapping: input a state, return the idx of the grid.
        function idx_G = abstract_mapping(obj,state)
        % input: state
        % output: idx_G --- = []: state is not in any grid
        %                   = n : state is in grid n, 
        %                   = [n1,n2,...]: if state is in
        %                   mulitple grids, return an array of indices.
            idx_G = zeros(obj.n,1,'logical');
            for i = 1:obj.n
                if all(state>=obj.G{i}(:,1) & state<=obj.G{i}(:,2))
                    idx_G(i) = true;
                end
            end
            idx_G = find(idx_G);
        end
        
        % plot phase portrait
        function phase_portrait(obj,fig, dense)
        % input: fig --- a figure handle, created by `fig = figure();`
            if nargin == 2
                dense = 20;
            end
            if obj.dim == 2
                figure(fig); hold on;
                for idx_u = 1:obj.m
                    y1 = linspace(obj.X(1,1),obj.X(1,2),dense);
                    y2 = linspace(obj.X(2,1),obj.X(2,2),dense);
                    [x,y] = meshgrid(y1,y2);
                    size(x);
                    size(y);
                    u = zeros(size(x));
                    v = zeros(size(x));
                    for i = 1:numel(x)
                        Yprime = obj.dyn([x(i); y(i)],obj.U(idx_u));
                        u(i) = Yprime(1);
                        v(i) = Yprime(2);
                        u(i) = u(i);
                        v(i) = v(i);
                    end
                    quiver(x,y,u,v,'color',rand(1,3));
                    xticks(linspace(obj.X(1,1),obj.X(1,2),5));
                    yticks(linspace(obj.X(2,1),obj.X(2,2),5));
                    axis tight equal;
                end
            else
                warning("Phase portrait is not available for dim > 2!");
            end
        end
        
        function ts = to_TransSyst(obj)
            
        end
        plot(obj);
    end
    methods(Access = protected)
        function initialize_A(obj)
            obj.A = cell(obj.m,1);
            Adj = eye(obj.n);
            for i = 1:obj.n
                for j = 1:i-1
                    if has_cmn_surface(obj,i,j)
                        Adj(i,j) = 1;
                        Adj(j,i) = 1;
                    end
                end
            end
            for i = 1:obj.m
                obj.A{i} = Adj;
            end
        end
    end
end