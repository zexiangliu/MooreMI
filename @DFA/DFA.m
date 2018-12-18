classdef DFA<handle
    properties (SetAccess = protected)
        Q; % states [1,2,3,4,...]
        Q0;
        Q_name; % name of states, string array ["ep","ab","abc"]
        Q_label;
        Q_final; % final states or accepting states
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
        function obj = DFA(n,m,A,s1,a,s2,Q_name,U_name,Q0,Q_final,Q_label)
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
                Q0 = 1;
            end
            if(nargin <= 9)
                Q_final = [];
            end
            if(nargin <= 10)
                Q_label = [];
            end
            
            obj.n = n;
            obj.Q = 1:n;
            obj.m = m;
            obj.U = 1:m;
            obj.Q_name = Q_name;
            obj.U_name = U_name;
            if isa(Q0,"string")
                Q0 = obj.get_x_idx(Q0);
            end
            obj.Q0 = Q0;
            if ~isempty(Q_final)
                if isa(Q_final(1),'string')
                    obj.Q_final = get_x_idx(obj,Q_final);
                else
                    obj.Q_final = Q_final;
                end
            end
            
            obj.Q_label= Q_label;
            
            if ~isempty(A)
                obj.A = A;
                obj.update_sas();
            elseif ~isempty(s1)
                if isa(s1,"string")
                    s1 = get_x_idx(obj,s1);
                end
                if isa(s2,"string")
                    s2 = get_x_idx(obj,s2);
                end
                if isa(a,"string")
                    a = get_u_idx(obj,a);
                end
                len_s1 = length(s1);
                obj.state1 = zeros(len_s1,1);
                obj.action = zeros(len_s1,1);
                obj.state2 = zeros(len_s1,1);
                obj.state1(1:end) = s1;
                obj.action(1:end) = a;
                obj.state2(1:end) = s2;
                obj.A = sas_to_A(obj);
            end
            obj.check();
        end
        
        function G = copy(obj)
            G = DFA(obj.n,obj.m,obj.A,obj.state1,...
                obj.action,obj.state2,obj.Q_name,obj.U_name,...
                obj.Q0,obj.Q_final,obj.Q_label);
        end
        
        % santiy check consistency
        function check(obj)
            assert(length(obj.Q) == obj.n);
            assert(length(obj.U) == obj.m);
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
            assert(length(obj.Q0)==1);
        end
        % add final states
        function obj = add_final(obj, Q_final)
            if isa(Q_final, "string")
                Q_final = obj.get_x_idx(Q_final);
            end
            obj.Q_final(end+1:end+length(Q_final)) = Q_final;
        end
        
        % set final states
        function obj = set_final(obj, Q_final)
            if isa(Q_final, "string")
                Q_final = obj.get_x_idx(Q_final);
            end
            obj.Q_final = Q_final;
        end
        
        % set state label
        function obj = set_label(obj, Q_label)
            if length(Q_label) ~= obj.n
                error("The number of Q_label and states mismatch.");
            end
            obj.Q_label = Q_label;
        end
        
        % add new states
        function obj = add_x(obj, x_name, x_label, isFinal)
            if ismember(x_name, obj.Q_name)
                error(x_name + " has existed!");
            end
            obj.n = obj.n+1;
            obj.Q(end+1) = obj.n;
            obj.Q_name(end+1) = x_name;
            if ~isempty(obj.Q_label)
                obj.Q_label(end+1) = x_label;
            end
            if isFinal
                obj.Q_final(end+1) = obj.n;
            end
            
            for i = 1:obj.m
                obj.A{i} = padarray(obj.A{i},[1,1],'post');
            end
        end
        
        % remove state x
        function obj = remove_x(obj,x)
            if isa(x,"string")
                x = get_x_idx(obj,x);
            end
            
            obj.n = obj.n-1;
            obj.Q(end) = [];
            obj.Q_name(x) = [];
            if ~isempty(obj.Q_label)
                obj.Q_label(x) = [];
            end
            idx = obj.Q_final==x;
            obj.Q_final(idx) = [];
            % need to re-index the states after x
            idx = obj.Q_final>x;
            obj.Q_final(idx) = obj.Q_final(idx)-1;
            % update the transition matrix
            for i = 1:obj.m
                obj.A{i} = obj.A{i}([1:x-1,x+1:end],[1:x-1,x+1:end]);
            end
            obj.update_sas();
        end
        
        % add new action
        function obj = add_u(obj,u_name)
           obj.m = obj.m+1;
           obj.U(end+1) = obj.m;
           obj.U_name(end+1) = u_name;
        end
        
        % remove action u
        function obj = remove_u(obj, u)
            if isa( u, "string")
                u = get_u_idx(obj,u);
            end
            
            obj.m = obj.m-1;
            obj.U(end) = [];
            obj.U_name(u) = [];
            idx = obj.action==u;
            obj.state1(idx) = [];
            obj.state2(idx) = [];
            obj.action(idx) = [];
            obj.A(u)= [];
            % need to re-index actions after u
            idx = obj.action > u;
            obj.action(idx) = obj.action(idx) - 1;
        end
        
        % add new transitions
        function obj = add_trans(obj, s1, a, s2)
            if isa(s1,"string")
                s1 = obj.get_x_idx(s1);
            end
            if isa(s2,"string")
                s2 = obj.get_x_idx(s2);
            end
            if isa(a, "string")
                a = obj.get_u_idx(a);
            end
            
            obj.state1(end+1:end+length(s1)) = s1;
            obj.state2(end+1:end+length(s2)) = s2;
            obj.action(end+1:end+length(a)) = a;
            obj.A = sas_to_A(obj);
        end
        
        % remove transitions
        function obj = remove_trans(obj, s1, a, s2)
            if isa(s1,"string")
                s1 = get_x_idx(obj,s1);
            end
            if isa(s2,"string")
                s2 = get_x_idx(obj,s2);
            end
            if isa(a,"string")
                a = get_u_idx(obj,a);
            end
            
            for i = 1:length(s1)
                obj.A{a(i)}(s1(i),s2(i)) = 0;
                idx = obj.state1==s1(i) & ...
                    obj.state2==s2(i) & obj.action==a(i);
                obj.state1(idx) = [];
                obj.state2(idx) = [];
                obj.action(idx) = [];
            end
        end
        
        % get state name
        function [x_name] = get_x_name(obj,x)
            if isa(x,"string")
                x_name = x;
                return;
            end
            x_name = obj.Q_name(x);
        end
        
        % get input name
        function [u_name] = get_u_name(obj,u)
            if isa(u,"string")
                u_name = u;
                return;
            end
            u_name = obj.U_name(u);
        end
        
        % get state label
        function [x_label] = get_x_label(obj,x)
            if isa(x,"string")
                x = get_x_idx(x);
            end
            
            x_label = obj.Q_label(x);
        end
        
        % get state index
        function [x_idx] = get_x_idx(obj,x)
            if isa(x,"double")
                x_idx = x;
                return;
            end
            
            x_idx = zeros(length(x),1);
            for i = 1:length(x_idx)
                x_idx(i) = find(obj.Q_name == x(i));
            end
        end
        
        % get input index
        function [u_idx] = get_u_idx(obj,u)
            if isa(u,"double")
                u_idx = u;
                return;
            end
    
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
            post_x = [];
            post_u = [];
            for i =1:length(x)
                idx = obj.state1==x(i);
                post_x = [post_x;obj.state2(idx)];
                post_u = [post_u;obj.action(idx)];
            end
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
        
        % check if a state is accepting
        function bool = isaccepting(obj,x)
            if isa(x,"string")
                x = get_x_idx(obj,x);
            end
            bool = ismember(x,obj.Q_final);
        end
        
        % check if a run is feasible or reaching an accepting state
        function [status,q,idx_u] = run(obj,inputs)
        % output: status = -1 infeasible
        %                =  0 feasible, not accepting
        %                =  1 feasible and accepting
        %         q --- the last feasible state, if status = -1
        %         idx_u --- idx of the first infeasible u if status = -1
            if isa(inputs,"string")
                inputs = get_u_idx(obj,inputs);
            end
            q = obj.Q0;
            idx_u = [];
            for i = 1:length(inputs)
                u = inputs(i);
                q_hist = q;
                q = find(obj.A{u}(q,:));
                if(isempty(q))
                    q = q_hist;
                    idx_u = i;
                    status = -1;
                    return;
                end
                % in case non-deterministic, always pick the first one
                q = q(1);
            end
            if ismember(q,obj.Q_final)
                status = 1;
            else
                status = 0;
            end
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
        
        function obj = merge(obj,s1,s2)
            if isa(s1,"string")
                s1 = get_x_idx(obj,s1);
            end
            if isa(s2,"string")
                s2 = get_x_idx(obj,s2);
            end
            
            [q_u,a_u] = obj.pre(s2);
            if length(q_u)~=1
                error("The parents of s2 are not unique!");
            end
            
            obj.A{a_u}(q_u,s2) = 0;
            obj.A{a_u}(q_u,s1) = 1;
            
            merge_stack_s1 = s1;
            merge_stack_s2 = s2;
            
            while(~isempty(merge_stack_s1))
                q1 = merge_stack_s1(end);
                q2 = merge_stack_s2(end);
                merge_stack_s1(end) = [];
                merge_stack_s2(end) = [];
                
                if q1==q2
                    continue;
                end
                if q1~=s1 && q2~=s2 && obj.isLess(q2,q1)
                    tmp_q1 = q1;
                    q1 = q2;
                    q2 = tmp_q1;
                end
                if ismember(q2,obj.Q_final) && ~ismember(q1,obj.Q_final)
                    obj.Q_final(end+1) = q1;
                end
                
                for u = 1:obj.m
                    q1_new = find(obj.A{u}(q1,:));
                    q2_new = find(obj.A{u}(q2,:));
                    
                    if ~isempty(q2_new)
                        if  ~isempty(q1_new)
                            merge_stack_s1(end+1) = q1_new;
                            merge_stack_s2(end+1) = q2_new;
                        else
                            obj.A{u}(q1,:) = obj.A{u}(q2,:);
                        end
                    end
                end
                % remove q2 and obj.s1, a, s2 are updated in remove_x
                obj.remove_x(q2); 
                % update indices in merge_stack
                idx = merge_stack_s1>q2;
                merge_stack_s1(idx) = merge_stack_s1(idx)-1;
                idx = merge_stack_s2>q2;
                merge_stack_s2(idx) = merge_stack_s2(idx)-1;
                % update s1 and s2
                if s1>q2
                    s1 = s1-1;
                end
                if s2>q2
                    s2 = s2-1;
                end
            end
        end

        function obj = complete(obj)
            for i = 1:obj.m
                sub = find(sum(obj.A{i},2)==0);
                idx = sub2ind([obj.n,obj.n],sub,sub);
                obj.A{i}(idx) = 1;
            end
            obj.update_sas();
        end
    end
    
    
    methods(Access = protected)
        function A = sas_to_A(obj)
            A = cell(obj.m,1);
            for i = 1:obj.m
                A{i} = zeros(obj.n,obj.n,'logical');
            end
            s1 = obj.state1;
            a =  obj.action;
            s2 = obj.state2;
            for i = 1:length(s1)
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
        
        function obj = update_sas(obj)
            [s1,a,s2] = obj.A_to_sas();
            obj.state1 = s1;
            obj.state2 = s2;
            obj.action = a;
        end
        
        % decide whether q1<q2
        function bool = isLess(obj,q1,q2)
            
            bool = false;
            
            if isa(q1,"double")
                q1 = obj.get_x_name(q1);
            end
            if isa(q2,"double")
                q2 = obj.get_x_name(q2);
            end
            eps0 = obj.get_x_name(obj.Q0);
            if q2 == eps0
                return;
            end
            
            if q1 == eps0 || strlength(q1) < strlength(q2) ||...
                    strlength(q2)==strlength(q1) && q1<q2
                bool = true;
            end
        end

    end
end