% %Computing max flow min cut
% 
% % directed graph
% Adir = Asp;
% 
% [i,j,V] = find(Asp)
% I = horzcat(i,j)
% Icpy = I;
% it = 1;
% cnt = 1;
% while  (size(Icpy,1))
%    ij = wrev(Icpy(it,:));
% 
%    if (ismember(ij, I, 'rows'))
%       [~,id] = ismember(ij, Icpy, 'rows');
%       rw = Icpy(id,1);
%       cl = Icpy(id,2);
%       Adir(rw,cl) = 0; %To connected graph now
%    end
%    Icpy(it,:) = [];
%    Icpy(id-1,:) = [];
%  
%    cnt = cnt + 1;
% end
% 
% 
% Adir
%  
% %no of edges
% m = nnz(Adir)
% 
% % Find min cut
% capacity = rand(m,1)
% 
% [M, F, C] = graphmaxflow(Adir,1,7, 'Capacity', capacity, 'Method', 'Edmonds')
% 
%Through linear programming
%weights matrix
[~,~,vals] = find(D(1,:))
edfromsrc = length(vals)
f = zeros(m,1);
f(1) = -1;
if(edfromsrc > 1)
    f(2:edfromsrc) = -1;
end
b = rand(1,m)
lb = zeros(m,1);
ub = ones(m,1);
x = linprog(f,W,b,[],[],lb,ub)
cut = f'*x

