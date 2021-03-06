function[f] = quadraticprogflowcons(partlab1, partlab2, nodelabi, nodelabj, Adj, D, L, num, f_j, mincut)
%Date:12.01.2014

%Calculating max flow outside
%Let there be a lambda_ij for each edge

%No of edges
[rw, col] = find(Adj); %to keep reference of the nodes
ids = find(Adj);
length(ids);
m = nnz(Adj);

%Weight vector
w = zeros(m,1);
w(ids) = 1;

n = num;
%ui-uj <= lambdaij  m x (m+n)
diagval = -1*ones(m,1);
Amatlambda  = diag(diagval);
%Amatu = zeros(m,n);

[r,c] = ind2sub(size(Adj),ids);
f_j;
r = r';
c = c';
% Create column index vector
iv = 1:1:length(r);

% Handle calculation of v1 and v2 separately, then add them
M1 = zeros(n,length(r));
M2 = zeros(n,length(c));

% This will generate vectors A and B, which are the 
% matrices we require (just in vector form)
A = zeros(1,n*length(r));
B = zeros(1,n*length(c));

B(sub2ind(size(M2), c, iv)) = -1;
A(sub2ind(size(M1), r, iv)) = 1;

% Turn them into matrices
%mat1 = vec2mat(A,length(r))
%mat2 = vec2mat(B,length(c))

mat1 = reshape(A,n,length(r));
mat2 = reshape(B,n,length(c));

% Make final matrix
Amatu = mat1 + mat2;
Amatu = Amatu';
size(Amatu);
size(Amatlambda);
Amat = horzcat(Amatu,Amatlambda);
size(Amat);
b = zeros(m,1);
size(b);

%Min cut constraint  TEST CODE
mincutvect = zeros(1,n+m);
mincutvect(1,n+1:end) = ones(1,m);
bformincutvect = mincut;



Amat = vertcat(Amat,mincutvect)
b = vertcat(b,bformincutvect)
%Min cut constraint  TEST CODE


ub = zeros(n+m,1);
ub(1:n) = 1*ones(n,1);
ub(n+1:end) = inf*ones(m,1)

lb = zeros(n+m,1);
lb(1:n) = -1*ones(n,1);
lb(n+1:end) = zeros(m,1)


%%%%%%%%%%%%%%%%%%%%%%%%%Equalities
labids = find(f_j);
colids = labids;
Aeq = zeros(n,n);
Aeq(sub2ind(size(Aeq),labids,colids)) = 1;
beq = f_j'  
Aeqpart = zeros(n,m);
Aeq = horzcat(Aeq,Aeqpart)

% B(sub2ind(size(M2), c, iv)) = 1;
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
% Aeq = mat1 + mat2
% beq = f_j'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Linear objective coefficient vector for summation of lambda in objective
%h = ones(1,m+n);
%h(1:n) = zeros(1,n);

%TEST CODE
h = zeros(m+n,1);

Laug = zeros(n+m,n+m);
Laug(1:n,1:n) = L

%Linear program formulation
%TEST CODE
%flowbasedlabels = quadprog(Laug,h,Amat,b,Aeq,beq,lb,ub,[],[]);
flowbasedlabels = quadprog(Laug,h,Amat,b,Aeq,beq,lb,ub,[],[]);
f = flowbasedlabels(1:n);

%size(Amat)
%size(b)
%size(Aeq)
%size(beq)

%res = msklpopt(w,A,blc,buc,blx,bux);

%flow = res.sol

end