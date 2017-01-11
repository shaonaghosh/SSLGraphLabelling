function[f] = harmonicinverse(partlab1, partlab2, nodelabi, nodelabj, A, D, L, num, f_j)

%Function to calculate the minimum harmonic energy minimization : Ref -Semi
%Supervised Learning using Gaussian Fields and Harmonic Functions 


%Author:Shaona Ghosh
%Date:10.01.2014

%Check for number of arguments
%error(nargchk(1, 9, nargin'));


%Copy labelled before changing order
fcpy = zeros(length(f_j),1);
id1s = find(f_j);
len1 = length(nodelabi);
len2 = length(nodelabj);


fcpy(1:len1) = partlab1;
fcpy(len1+1:len1+len2) = partlab2; %the rest does not matter as all 0s


%Find the permutation order 
ids = find(not(f_j));

%Check if col vector
if size(ids,1) ~= 1
    ids = ids';
end

p  = [nodelabi nodelabj ids];

%no of labelled and unlabelled points
u = numel(ids);
l = numel(find(fcpy));
f_l = fcpy(1:l);

%Permute the weight matrix
W_j = sparse(num,num);
%W_j(p,:) = A(:,p);
W_j = A(p,p); %Faster than the above method
full(W_j);

%Permute the degree matrix
D_j = sparse(num,num);
D_j = D(p,p);
%D_j(p,:) = D(:,p);

Duu = D_j(l+1:end,l+1:end);
Wuu = W_j(l+1:end,l+1:end);

Wul = W_j(l+1:end,1:l);

%Test code
% Wul = full(Wul)
% f_l 
% Duu = full(Duu)
% Wuu = full(Wuu)
% 
% interm2 = Wul*f_l
% interm2 = full(interm2)
% 
% finterm = inv(Duu-Wuu)
% f_u = finterm*(Wul*f_l)

P = D\A;
P_j = P(p,p);

Puu = P_j(l+1:end,l+1:end);
Pul = P_j(l+1:end,1:l);

interm = (speye(size(Puu)) - Puu);
f_u = (interm)\(Pul*f_l);


ff = vertcat(f_l,f_u);
%Reorder the vertices in the original order
f = zeros(num,1);
f(p) = ff;



end