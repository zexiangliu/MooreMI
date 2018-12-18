classdef Abstraction < handle
    properties (SetAccess = protected)
        n; % number of grids
        dim; % dimension of state space
        m; % number of actions
        G; % partitions, cell of grids
        U; % array of actions
        A; % adjacent matrix of the grids
        spec; % the verification/synthesis specification
        label; % labeling of grids
        ts; % transition system corresponding to the partitions
        dyn; % function handle, dynamics of the system, take state ...
             % and action as inputs
        K; % Lipschitz constant, a function of state and action
    end
    
    methods
        function obj = Abstraction(U, G, spec, label, dyn, K)
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
        end
        
        function obj = verifyTransition(obj)
            for i = 1:length(obj.n)
                for j = 1:i-1
                    % if grid(i)=grid(j) are neighbours, then do transition function
                    [a, P] = if_neighbor(obj.G{i},obj.G{j});
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
        end
        
        % check if two grids are neighbors
        function bool = is_neighbor(obj, G1, G2)
            PG1 = obj.G{G1};
            PG2 = obj.G{G2};
            bool = false;
            if any(PG1(:,1)==PG2(:,2))||any(PG1(:,2)==PG2(:,1))
                bool = true;
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
        
        ts = to_TransSyst(obj);
        plot(obj);
    end
    methods(Access = protected)
        function initialize_A(obj)
            obj.A = eye(obj.n);
            for i = 1:obj.n
                for j = i+1:obj.n
                    if is_neighbor(obj,i,j)
                        obj.A(i,j) = 1;
                        obj.A(j,i) = 1;
                    end
                end
            end
        end
    end
end