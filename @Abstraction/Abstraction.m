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
        dyn; % function handle, dynamics of the system, take state ...
             % and action as inputs
        K; % Lipschitz constant, a function of state and action
        mini_width; % the minimal width of grid
        
        ts; % transition system corresponding to the partitions
        Win; % winning set
        CWin; % candidate winning set
        cont; % controller
        ts_enabled; % if ts is computed
        transit_cerf; % the trasit direction for each grid
                      % format: {action}{direction}{+,-}{s1,s2,s3...}
        is_conservative; % whether re-verify transit state (incomplete)
        specical_action_group; % save the multi-action progress group
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
            obj.ts_enabled = false;
            obj.transit_cerf = cell(obj.m,1);
            for i = 1:obj.m
                obj.transit_cerf{i} = cell(obj.dim,1);
                for j = 1:obj.dim
                    obj.transit_cerf{i}{j} = {[],[]};
                end
            end
            obj.is_conservative = true;
            obj.specical_action_group = {};
        end
        
        % verify the transitions in A
        function obj = verifyTransition(obj, mini_width, G2_idx)
            % input: G2_idx --- optional, you can select some idx to update
            %                   partially.
           
            if nargin == 1 || isempty(mini_width)
                mini_width = obj.mini_width;
            end
            if nargin <= 2
                flag_part = false;
                G2_idx = [];
            else
                flag_part = true;
                G2_idx = sort(G2_idx);
            end
            for idx_u = 1:obj.m
                for G1 = 1:obj.n
                    if flag_part % if verify part of the states
                        list_G2 = setdiff(G2_idx,G1);
                    else
                        list_G2 = 1:G1-1;
                    end
                    for G2 = list_G2
                        % verify if there is any transition between i and j
                        if obj.A{idx_u}(G1,G2) == 1 || obj.A{idx_u}(G2,G1) == 1
                            [~, ~, cmn_surface] = has_cmn_surface(...
                                                            obj, G1, G2);
                            direct = has_transition(obj, G1, G2, ...
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
                    
                    if ~flag_part || flag_part && ismember(G1,G2_idx)
                        if obj.A{idx_u}(G1,end) == 1 || obj.A{idx_u}(end,G1)
                            direct = obj.is_bnd_inv(G1,idx_u,mini_width);
                            switch direct
                                case 0 % both
                                    obj.A{idx_u}(G1,end) = 1;
                                    obj.A{idx_u}(end,G1) = 1;
                                case 1 % in
                                    obj.A{idx_u}(G1,end) = 0;
                                    obj.A{idx_u}(end,G1) = 1;
                                case 2 % out
                                    obj.A{idx_u}(G1,end) = 1;
                                    obj.A{idx_u}(end,G1) = 0;
                                otherwise
                                    error("direct is not correct!");
                            end
                        end
                        
                        if obj.A{idx_u}(G1,G1) == 1 
                            [bool, flag_transit] = ...
                                obj.is_transit(G1,idx_u);
                            if bool
                                obj.A{idx_u}(G1,G1) = 0;
                                % add G1 to transit_certificate
                                idx_transit = find(flag_transit~=0);
                                for t = 1:length(idx_transit)
                                    idx = idx_transit(t);
                                    direction = 1+(flag_transit(idx)+1)/2;
                                    obj.transit_cerf{idx_u}...
                                        {idx}{direction}(end+1) = G1;
                                end
                            end
                        end
                    end
                end
            end
        end
        % verify the self loop of a state
        function [bool,flag_transit] = is_transit(obj,G1,idx_u,mini_width)
            if nargin <=3
                mini_width = obj.mini_width;
            end
            
            PG1 = obj.G{G1};
            u = obj.U(idx_u);
            x = mean(PG1,2);
            flag_transit = sign(obj.dyn(x,u));
            [~, cover] = obj.compute_cover(x,u,flag_transit,1, PG1);
            residual = obj.surface_diff(PG1,cover);
            
            queue = residual;
            
            bool = false;
            while ~isempty(queue)
                % queue pop
                grid = queue{1};
                
                queue(1) = [];
                % compute the cover and decide flow direction
                x = mean(grid,2);
                [flow, cover] = obj.compute_cover(x,u,flag_transit,1,grid);
                
                % if the grid is too small or cover is too small, stop
                if max(grid(:,2)-grid(:,1)) <= mini_width...
                         && any(grid(:,1) < cover(:,1) & ...
                        grid(:,2) > cover(:,2)) 
                    flow(flag_transit==0) = inf;
                    flow(flow == 0) = inf;
                    % discard the dimension that has minimal flow
                    [~,discard_idx] = min(abs(flow));
                    flag_transit(discard_idx) = 0;
                    queue{end+1} = grid; % re-push the grid into the queue
                    continue;
                end
                
                idx = sign(flow) ~= flag_transit;
                flag_transit(idx) = 0;
                residual = obj.surface_diff(grid,cover);
                % queue push
                queue(end+1:end+length(residual)) = residual;
                
                if all(flag_transit == 0) || all(flow == 0)
                    return;
                end
            end 
            bool = true;
        end
        
        % verify if the grid on the boundary of the state space is inv,
        % means that the flow points into the grid/state space
        function direct = is_bnd_inv(obj, ...
                G1, idx_u, mini_width)
            % output: direct --- =0 indeterminate
            %                    =1 in
            %                    =2 out
            if nargin == 3
                mini_width = obj.mini_width;
            end
            PG = obj.G{G1};
            idx_lb = find(PG(:,1)==obj.X(:,1));
            idx_ub = find(PG(:,2)==obj.X(:,2));
            u = obj.U(idx_u);
            direct = []; % collect directions of common surfaces on bnd.
            % lower or left bnd case
            for i = 1:length(idx_lb)
                idx = idx_lb(i);
                normal = zeros(obj.dim,1);
                normal(idx) = 1;
                cmn_surface = PG;
                cmn_surface(idx,2) = cmn_surface(idx,1);
                direct(end+1) =  verify_surface(...
                obj, u, cmn_surface,normal,mini_width);
                if direct(end) == 0 || ...
                        ismember(1,direct) && ismember(2,direct)
                    direct = 0;
                    return;
                end
            end
            
            % upper or right bnd case
            for i = 1:length(idx_ub)
                idx = idx_ub(i);
                normal = zeros(obj.dim,1);
                normal(idx) = -1;
                cmn_surface = PG;
                cmn_surface(idx,1) = cmn_surface(idx,2);
                direct(end+1) =  verify_surface(...
                    obj, u, cmn_surface,normal,mini_width);
                if direct(end) == 0 || ...
                        ismember(1,direct) && ismember(2,direct)
                    direct = 0;
                    return;
                end
            end 
            
            if all(direct == 1)
                direct = 1;
            elseif all(direct == 2)
                direct = 2;
            else
                direct = 0;
            end
        end
        
        % verify the direction of the flow on the common surface.
        function direct =...
                has_transition(obj, G1, G2, cmn_surface, idx_u, mini_width)
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
            u = obj.U(idx_u);
            % compute the normal vector of the common surface from i to j
            normal = double(cmn_surface(:,2)==cmn_surface(:,1));
            center_G1 = mean(obj.G{G1},2);
            center_G2 = mean(obj.G{G2},2);
            normal = sign(center_G2-center_G1).*normal;
            direct =  verify_surface(obj, u, cmn_surface,...
                    normal,mini_width);
        end
        
        % helper function of has_transition
        function [direct, list_cover, cover_direct] =  verify_surface(...
                obj, u, cmn_surface, normal, mini_width)
        % output: direct --- = 0 indeterminate
        %                    = 1 along the normal
        %                    = 2 along the negative normal
            flag_cover = nargout > 1;
            list_cover = {};
            cover_direct = {};
            queue = {cmn_surface};
            flag_12 = false;
            flag_21 = false;
            while ~isempty(queue)
                % queue pop
                surface = queue{1};
                queue(1) = [];
                % compute the cover and decide flow direction
                x = mean(surface,2);
                [flow, cover] = obj.compute_cover(x,u,normal,0,surface);
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
        function [flow, cover] = compute_cover(obj,x,u,normal, mode, range)
        % inputs: x --- state you want to query
        %         u --- action applied
        %         normal --- normal vector of the common surface
        %         mode --- = 0 find cover for surface
        %                  = 1 find cover for grid
        %         range --- the range you want to compute Lipschitz
        %                   contant
            if nargin <= 4 || isempty(mode)
                mode = 0;
            end
            flow = obj.dyn(x,u);
            if mode == 0
                r = min(abs(flow(normal~=0)))/obj.K(x,u,range);
            elseif mode == 1
                normal(flow == 0) = 0;
                if all(normal==0)
                    r=0;
                else
                    % here we make the r a little bit smaller such that the
                    % flow has strictly positive component by dividing
                    % (1+1e-6).
                    r = min(abs(flow(normal~=0)))/obj.K(x,u,range)/(1+1e-6);
                end
            end
            if isnan(r)
                keyboard();
            end
            cover = [x-ones(obj.dim,1)*r,x+ones(obj.dim,1)*r];
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
            if G1 == G2
                bool = true;
                cmn_dim = 0;
                cmn_surface =[];
                return;
            end
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
            if nargin <= 1
                fig = figure;
            end
            if nargin <= 2
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
        
        % convert the abstraction to TransSyst format
        function obj = to_TransSyst(obj, magp_enabled, mapg_pairs)
        % mapg_pairs --- {dim,direct,{u_list}}
            if nargin == 1
                magp_enabled = false;
                mapg_pairs = [];
            end
            if magp_enabled == true && isempty(mapg_pairs)
                mapg_pairs = ones(obj.dim,2);
            end
            n_s = obj.n+1;
            n_a = obj.m;
            
            s1 = [];
            s2 = [];
            a = [];
            for idx_u = 1:obj.m
                for i = 1:obj.n+1
                    for j = 1:obj.n+1
                        if obj.A{idx_u}(i,j)
                            s1(end+1) = i;
                            s2(end+1) = j;
                            a(end+1) = idx_u;
                        end
                    end
                end
            end
            obj.ts = TransSyst(n_s, n_a);
            obj.ts.add_transition(s1,s2,a);
            
            % add special multi-action progess group
            if magp_enabled
%                 obj.specical_action_group = {};
%                 for i = 1:obj.dim
%                     for j = 1:2
%                         if magp_direction(i,j)
%                             mapg_x = [];
%                             xu_mapping = zeros(obj.n,obj.m,'logical');
%                             for idx_u = 1:obj.m
%                                 tmp_x = obj.transit_cerf{idx_u}{i}{j};
%                                 if ~isempty(tmp_x)
%                                     mapg_x = union(mapg_x,tmp_x);
%                                     xu_mapping(tmp_x,idx_u) = 1;
%                                 end
%                             end
%                             mapg_x = unique(mapg_x);
%                             xu_mapping = xu_mapping(mapg_x,:);
%                             % lookup table x---> available u
%                             xu_list = cell(1,length(mapg_x));
%                             for idx_x = 1:length(mapg_x)
%                                 xu_list{idx_x} = ...
%                                     find(xu_mapping(idx_x,:)==1);
%                             end
%                             SUB = ndgrid2(xu_list);
%                             mapg_u = [];
%                             new_u_mapping = zeros(obj.n,1);
%                             for k = 1:length(SUB{1})
%                                 obj.ts.add_action();
%                                 new_u = obj.ts.n_a;
%                                 mapg_u = [mapg_u;new_u];
%                                 % add transition for new input
%                                 for l = 1:length(mapg_x)
%                                     x = mapg_x(l);
%                                     u = SUB{l}(k);
%                                     new_u_mapping(x) = u;
%                                     s2 = find(obj.A{u}(x,:));
%                                     a = s2*0 + new_u;
%                                     s1 = s2*0 + x;
%                                     obj.ts.add_transition(s1,s2,a);
%                                 end
%                                 obj.specical_action_group{end+1}=new_u_mapping;
%                             end
%                             obj.ts.add_progress_group(mapg_u,mapg_x);
%                         end
%                     end
%                 end
                % simpler implementation
                for i = 1:length(mapg_pairs)
                    d = mapg_pairs{i}{1};
                    dir = mapg_pairs{i}{2};
                    u_list = mapg_pairs{i}{3};
                    mapg_x = ones(obj.n,1,'logical');
                    for j = 1:length(u_list)
                        u = u_list(j);
                        tmp_x = zeros(obj.n,1,'logical');
                        tmp_x(obj.transit_cerf{u}{d}{dir})=1;
                        mapg_x =  mapg_x & tmp_x;
                    end
                    obj.ts.add_progress_group(u_list,find(mapg_x));
                end
            end
            
            obj.ts.trans_array_enable();
            obj.ts_enabled = true;
        end
        
        % set spec property
        function obj = set_spec(obj, A, B, C_list, quant1, quant2)
        % spec supported is []A & <>[]B & (AND_i []<>C_i)
        % inputs: A --- a subset of 1:n
        %         B --- a subset of 1:n
        %         C_list --- a cell of subsets of 1:n
            obj.spec.A = A;
            obj.spec.B = B;
            obj.spec.C = C_list;
            if nargin == 4
                obj.spec.quant1 = 'exists';
                obj.spec.quant2 = 'forall';
            else
                obj.spec.quant1 = quant1;
                obj.spec.quant2 = quant2;
            end
        end
        
        % given labeling, return the states with that label.
        function states = label_inverse(obj, label)
            states = [];
            for i = 1:obj.n
                if ismember(label, obj.label{i})
                    states(end+1) = i;
                end
            end
        end
        
        % synthsis
        function [V, Cv, control, obj] = win_primal(obj, A, B, C_list, ...
                quant1, quant2, V)
            if obj.ts_enabled == false
                obj.to_TransSyst();
            end
            if nargin == 1
                A = obj.spec.A;
                B = obj.spec.B;
                C_list = obj.spec.C;
            end
            if nargin <= 4
                try
                    quant1 = obj.spec.quant1;
                    quant2 = obj.spec.quant2;
                catch
                    quant1 = 'exists';
                    quant2 = 'forall';
                end
            end
            if nargin <= 6
                V = [];
            end
            [V, Cv, control] = obj.ts.win_primal(...
                                A, B, C_list, quant1, quant2, V);
            Cv = setdiff(Cv,obj.n+1); % remove dummy state
            obj.Win = V;
            obj.CWin = Cv;
            obj.cont = control;
        end
        
        % refine grid G1 according to the covers on the boundary
        function obj = refinement(obj, G1)
            % disable synthesis
            obj.ts_enabled = false;
            new_grids = obj.grid_refine_naive(G1);
            num_new = length(new_grids);
            new_idx = [G1,obj.n+1:obj.n+num_new-1];
            % update partition
            obj.G{G1} = new_grids{1};
            obj.G(end+1:end+num_new-1)=new_grids(2:end);
            % update labeling
            obj.label(end+1:end+num_new-1) = obj.label(G1);
            % update n
            obj.n = length(obj.G);
            % update A
            obj.update_A(new_idx);
            
            % verify new A
            obj.verifyTransition([],new_idx);
            
            % update transit_cerf
            for i = 1:obj.m
                for j = 1:obj.dim
                    for k = 1:2
                        if obj.is_conservative
                            if ismember(G1,obj.transit_cerf{i}{j}{k})
                                obj.transit_cerf{i}{j}{k}...
                                    (end+1:end+num_new-1) = new_idx(2:end);
                            end
                        else
                            if ismember(G1,obj.transit_cerf{i}{j}{k})
                                obj.transit_cerf{i}{j}{k} = setdiff(...
                                    obj.transit_cerf{i}{j}{k},G1);
                            end
                        end
                    end
                end
            end
                    
            
            % update spec
            try
                if ismember(G1,obj.spec.A)
                    obj.spec.A(end+1:end+num_new-1) = new_idx(2:end);
                end
                if ismember(G1,obj.spec.B)
                    obj.spec.B(end+1:end+num_new-1) = new_idx(2:end);
                end
                
                for i = 1:length(obj.spec.C)
                    if ismember(G1,obj.spec.C{i})
                        obj.spec.C{i}(end+1:end+num_new-1) = new_idx(2:end);
                    end
                end
            catch
                disp("No spec yet.");
            end
            
        end
        
        % naive refinement strategy
        function new_grids = grid_refine_naive(obj,grid)
            PG = obj.G{grid};
            new_grids = cell(2^obj.dim,1);
            
            % divide each dimension into two intervals
            list_intervals = cell(obj.dim,2);
            DIM = cell(obj.dim,1);
            for i = 1:obj.dim
                cnt_pt = mean(PG(i,:));
                list_intervals{i}{1} = [PG(i,1) cnt_pt];
                list_intervals{i}{2} = [cnt_pt PG(i,2)];
                DIM{i} = 1:2;
            end
            % take the combinations as new grids
            SUB = ndgrid2(DIM);
            for i = 1:length(SUB{1})
                new_grid = zeros(obj.dim,2);
                for j = 1:obj.dim
                    idx_int = SUB{j}(i);
                    new_grid(j,:) = list_intervals{j}{idx_int};
                end
                new_grids{i} = new_grid;
            end
        end
        
        % check if a grid is on the boundary.
        function bool = is_boundary(obj,s)
            PG = obj.G{s};
            bool = false;
            if any(PG(:,1)==obj.X(:,1)) || any(PG(:,2)==obj.X(:,2))
                bool = true;
            end
        end
        
        % plot the partition of the abstraction
        function plot(obj,fig,list_grid,color)
        % input: fig --- figure handle
        %        list_grid --- you can specify some grids to plot
            if nargin <= 1
                figure;
            else
                figure(fig);
            end
            hold on;
            if nargin <= 2 || isempty(list_grid)
                list_grid = 1:obj.n;
            end
            if nargin <= 3
                color = rand(1,3);
            end

            
            for i = 1:length(list_grid)
                idx = list_grid(i);
                if idx == 0
                    grid = obj.X;
                else
                    grid = obj.G{idx};
                end
                v = [grid(1,1) grid(2,1);
                     grid(1,1) grid(2,2);
                     grid(1,2) grid(2,2);
                     grid(1,2) grid(2,1)];
                f = [1 2 3 4];
                patch('Faces',f,'Vertices',v,...
                'EdgeColor',color,'FaceColor','none','LineWidth',2);
            end
        end
    end
    methods(Access = protected)
        function initialize_A(obj)
            obj.A = cell(obj.m,1);
            Adj = eye(obj.n+1,'logical');% we implicitly add one dummy state in A to
                                 %  indicate the region beyond X
            for i = 1:obj.n
                for j = 1:i-1
                    if has_cmn_surface(obj,i,j)
                        Adj(i,j) = 1;
                        Adj(j,i) = 1;
                    end
                end
                if is_boundary(obj,i)
                    Adj(i,end) = 1;
                    Adj(end,i) = 1;
                end
            end
            
            for i = 1:obj.m
                obj.A{i} = Adj;
            end
        end
        
            
        function obj = update_A(obj, new_idx)
            A_pad = zeros(obj.n+1,'logical');
            for G1 = new_idx
                for G2 = 1:obj.n
                    if obj.has_cmn_surface(G1,G2)
                        A_pad(G1,G2) = 1;
                        A_pad(G2,G1) = 1;
                    end
                end
                if is_boundary(obj,G1)
                    A_pad(G1,end) = 1;
                    A_pad(end,G1) = 1;
                end
            end
            
            siz_A = size(obj.A{1},1);
            for i = 1:obj.m
                A_i = obj.A{i};
                obj.A{i} = zeros(obj.n+1,'logical');
                % copy main part of old A
                obj.A{i}(1:siz_A-1,1:siz_A-1)= A_i(1:end-1,1:end-1);
                % reset G1
                obj.A{i}(new_idx(1),:) = 0;
                obj.A{i}(:,new_idx(1)) = 0;
                % copy old dummy transitions
                obj.A{i}(end,1:siz_A-1) = A_i(end,1:siz_A-1);
                obj.A{i}(1:siz_A-1,end) = A_i(1:siz_A-1,end);
                % reset dummy transition of G1
                obj.A{i}(end,new_idx(1)) = 0;
                obj.A{i}(new_idx(1),end) = 0;
                % add self-loop to dummy state
                obj.A{i}(end,end) = 1;
                % add A_pad to obj.A{i}
                obj.A{i} = obj.A{i} | A_pad;
                % self loop of new states: if original grid has no self
                % loop, then the new states have no self loop
                if A_i(new_idx(1),new_idx(1))==0 
                    if obj.is_conservative
                        for j = 1:length(new_idx)
                            obj.A{i}(new_idx(j),new_idx(j)) = 0;
                        end
                    else
                        obj.A{i}(new_idx(1),new_idx(1))=1;
                    end
                end
            end
        end
    end
    
    methods(Static)
        % compute residuals after one iteration  
        function residual = surface_diff(surface, cover)
        % compute the difference (surface - cover) and represent the 
        % residual using rectangles
            dim = size(surface,1);
            list_interval = cell(dim,1);
            
            % the idx of cover in the new partition, needs to be removed
            idx_cover = cell(dim,1);
            % dim of list_interval
            list_dim = zeros(dim,1);
            
            for i = 1:dim
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
                surf = zeros(dim,2);
                for j = 1:dim
                    surf(j,:) = list_interval{j}{sub(j)};
                end
                residual{i} = surf;
            end
            idx_cover = sub2ind2(list_dim,idx_cover);
            residual(idx_cover) = [];
        end 
    end
end