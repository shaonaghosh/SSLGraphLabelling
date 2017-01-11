 function[f] = quadraticprog(partlab1, partlab2, nodelabi, nodelabj, Adj, D, L, num, f_j)
% Author: Shaona Ghosh
% Date: 11/01/2014
% 
% 

n = num;

% A.  Creating the constraint A matrix in the inequailty constraint Ax <= b :
%Dealing with Coefficient matrix for ui-uj = label difference
% m = nnz(Adj);       %No of edges
% 
% ids = find(Adj);
% [r,c] = ind2sub(size(Adj),ids);
% f_j;
% r = r';
% c = c';
% % Create column index vector
% iv = 1:1:length(r);
% 
% % Handle calculation of v1 and v2 separately, then add them
% M1 = zeros(n,length(r));
% M2 = zeros(n,length(c));
% 
% % This will generate vectors A and B, which are the 
% % matrices we require (just in vector form)
% A = zeros(1,n*length(r));
% B = zeros(1,n*length(c));
% 
% B(sub2ind(size(M2), c, iv)) = -1;
% A(sub2ind(size(M1), r, iv)) = 1;
% 
% % Turn them into matrices
% %mat1 = vec2mat(A,length(r))
% %mat2 = vec2mat(B,length(c))
% 
% mat1 = reshape(A,n,length(r));
% mat2 = reshape(B,n,length(c));
% 
% % Make final matrix
% Amat = mat1 + mat2;
% 
% %Creating the b matrix in the inequality constraint Ax <= b
% b1 = zeros(length(r),1);
% b2 = zeros(length(c),1);
% 
% iv = ones(length(r),1);
% %length(iv)
% %length(r)
% 
% %b1 = f_j(sub2ind(size(b1),r,iv)); %TO FIX ISSUE ON INDEX OUT OF RANGE
% b1 = f_j(1,r);
% b2 = f_j(1,c);
% b = b1-b2;



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           METHOD 1 : using quadprog solver

% %Other parameters for the solver
H = full(L);             %The Q matrix here is the Laplacian
g = zeros(n,1);          %The linear part of the objective
A = zeros(n,n);          %When considering no constraint as in A above
b = zeros(n,1);

% Equality constraints u_s = 0, u_t = -1
f_j;
labids = find(f_j);
colids = labids;
Aeq = zeros(n,n);
Aeq(sub2ind(size(Aeq),labids,colids)) = 1;
beq = f_j;       

% Upper and lower bounds on the labels
ub = ones(n,1);       
lb = -1*ones(n,1);
% 
ff_j = quadprog(H,g,A,b,Aeq,beq,lb,ub,[]); %[],'Algorithm', 'active-set','Display', 'on');
ff_j;
f = ff_j;
% exitflag


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            FAST  METHOD 2 : using MOSEK solver

% i = ones (3 ,1);
% j = 1:3;
% vij = ones (3 ,1);
% bux = sparse ( i , j , vij )
%
% i = ones (3 ,1);
% j = 1:3;
% vij = -1* ones (3 ,1);
% blx = sparse ( i , j , vij )

%%%%%%%% METHOD 1: Calling the optimizer through wrapper - Mosekqpopt

% f = zeros(n,1);         %The linear part of the objective
% Aeq = diag(ones(n,1));  % Equality constraints u_s = 0, u_t = -1
% beq = f_j;
% 
% blc = [];                % Upper and lower bounds on the constraints
% buc = beq;
% 
% bux = ones(n,1);         % Upper and lower bounds on the labels
% blx = -1*ones(n,1);
% 
% [res] = mskqpopt(full(L),f,Aeq,blc,buc,bux,blx);
% sol = res.sol;
% f = sol.itr.xx';

%%%%%% METHOD 2: Calling the optimizer directly - Mosekopt

 end