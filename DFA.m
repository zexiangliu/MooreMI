classdef DFA
    properties (SetAccess = protected)
        Q; % states [1,2,3,4,...]
        Q_name; % name of states, string array ["ep","ab","abc"]
        Q_label;
        U; % actions {'a','b','c','d',...}
        U_name;
        n; % number of states
        m; % number of actions
        % kx1 vectors to save transitions s1--a-->s2
        state1;
        action;
        state2;
        % mx1 cell to save adjacency matrices for every action
        A;
    end
    methods(Access = public)
        % constructor
        function obj = DFA(n,m,A,s1,a,s2,Q_name,U_name,Q_label)
            if(nargin==3)
                s1 = [];                
                a = [];
                s2 = [];
            end
            if(nargin <= 6)
                Q_name = [];
                U_name = [];
            end
            if(nargin <= 8)
                Q_label = [];
            end
            
            obj.n = n;
            obj.Q = 1:n;
            obj.m = m;
            obj.U = 1:m;
            obj.Q_name = Q_name;
            obj.U_name = U_name;
            obj.Q_label= Q_label;
            
            if ~isempty(A)
                obj.A = A;
                [s1,a,s2] = A_to_sas(obj);
                obj.state1 = s1;
                obj.action = a;
                obj.state2 = s2;
            elseif ~isempty(s1)
                obj.state1 = s1;
                obj.action = a;
                obj.state2 = s2;
                obj.A = sas_to_A(obj);
            end
        end

        % santiy check consistency
        function check(obj)
            assert(lenght(obj.Q) == obj.n);
            assert(lenght(obj.U) == obj.m);
            if ~isempty(obj.Q_name)
                assert(length(obj.Q)==length(obj.Q_name));
            end
            if ~isempty(obj.Q_label)
                assert(length(obj.Q)==length(obj.Q_label));
            end
            if ~isempty(obj.U_name)
                assert(length(obj.U) == length(obj.U_name));
            end
            assert(length(obj.U) == length(obj.A));
            assert(length(obj.state1) == length(obj.state2) && ...
                length(obj.state2) == length(obj.action));
        end
        
        % get state name
        function [x_name] = get_x_name(obj,x)
            x_name = obj.Q_name(x);
        end
        
        % get input name
        function [u_name] = get_u_name(obj,u)
            u_name = obj.U_name(u);
        end
        
        % get state label
        function [x_label] = get_x_label(obj,x)
            x_label = obj.Q_label(x);
        end
        
        % get state index
        function [x_idx] = get_x_idx(obj,x)
            x_idx = zeros(length(x),1);
            for i = 1:length(x_idx)
                x_idx(i) = find(obj.Q_name == x(i));
            end
        end
        
        % get input index
        function [u_idx] = get_u_idx(obj,u)
            u_idx = zeros(length(u),1);
            for i = 1:length(u_idx)
                u_idx(i) = find(obj.U_name == u(i));
            end
        end
        
         % predecessors of state x
        function [pre_x, pre_u] = pre(obj,x)
            % input x could be either index or string name
            if isa(x,"string")
                x = get_x_idx(obj,x);
            end
            
            idx = obj.state2==x;
            pre_x = obj.state1(idx);
            pre_u = obj.action(idx);
        end
        
        % get the predecessor of x under action u
        function [pre_x] = pre_xu(obj,x,u)
            % input x could be either index or string name
            if isa(x,"string")
                x = get_x_idx(obj,x);            
            end
            if isa(u,"string")
                u = get_u_idx(obj,u);            
            end
            
            idx = obj.state2==x & obj.action==u;
            pre_x = obj.state1(idx);
        end
        
        % successor of state x
        function [post_x,post_u] = post(obj,x)
            if isa(x,"string")
                x = get_x_idx(obj,x);
            end
            
            idx = obj.state1==x;
            post_x = obj.state2(idx);
            post_u = obj.action(idx);
        end
        
        % successor of state x
        function post_x = post_xu(obj,x,u)
            if isa(x,"string")
                x = get_x_idx(obj,x);            
            end
            if isa(u,"string")
                u = get_u_idx(obj,u);            
            end
            
            idx = obj.state1==x & obj.action==u;
            post_x = obj.state2(idx);
        end
        
        % controlled predecessors of a set of state X
        function preX = CPre(obj,X)
            if isa(X(1),"string")
                X = get_x_idx(obj.X);
            end
            preX = zeros(obj.n,1,'logical');
            for i = 1:obj.m
                idx = sum(obj.A{i}(:,X),2) == sum(obj.A{i},2);
                preX = preX|idx;
            end
            preX = find(idx);
        end
    end
    
    
    methods(Access = protected)
        function A = sas_to_A(obj)
            s1 = obj.state1;
            a =  obj.action;
            s2 = obj.state2;
            A = cell(obj.m,1);
            for i = 1:length(s1)
                if isempty(A{a(i)})
                    A{a(i)} = zeros(obj.n,obj.n);
                end
                A{a(i)}(s1(i),s2(i)) = 1;
            end
        end
        
        function [s1,a,s2] = A_to_sas(obj)
            s1 = [];
            a = [];
            s2 = [];
            
            for i = 1:obj.m
                [s1_i,s2_i] = ind2sub([obj.n,obj.n],find(obj.A{i}~=0));
                a_i = i*ones(size(s1_i));
                s1 = [s1;s1_i];
                a = [a;a_i];
                s2 = [s2;s2_i];
            end
        end
    end
end