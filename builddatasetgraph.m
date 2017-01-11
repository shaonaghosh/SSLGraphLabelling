function [A,D,L] = builddatasetgraph(n,data, distancemat1, K)

%Function to compute the graph corresponding to the dataset
%Author: Shaona Ghosh
%Date: 10.01.2014

%Check for number of arguments
%error(nargchk(1, 3, nargin'));

if nargin < 4
    K = 5;
end

%Test Case
% data = [1,2,3;4,5,6;7,8,9]
% n = 3;
% K = 2


%Adjacency matrix 
A = sparse(n,n);

% Nearest neighbour method 1 :Faster than 2 but bug in construcing A Todo
% for i = 1:n
%     %for every node, connect to its k most nearest neighbours
%     nodei = data(i,:);
%     [weight, idx] = sort(distancemat1(i,:),2,'ascend');
%     %[weight, idx] = sort(distancemat1(i,:),'ascend')%TEST CODE
%     %Ignore the first index
%     idx(1) = [];
%     idx(K+1:end) = [];
%     idx
%     for j = 1:K  %Short loop 
%          A(i,idx(j)) = weight(j);
%          A(idx(j),i) = weight(j);
%     end
% end

%Build graph with k nearest neighbours using optimized toolbox- Method 2
tic
for i = 1:n
   nodei = data(i,3:end);
   [idx] = nearestneighbour(nodei',data(:,3:end)','NumberOfNeighbours',K );
   %logindices = sub2ind(data',idx,idx)
   % Ignore the first index as that will be nearest to itself
   idx(1) = [];
   for j = 1:length(idx) %Shorter loop
       A(i,idx(j)) = distancemat1(i,idx(j));  %Todo not efficient need to sparse access
       A(idx(j),i) = distancemat1(i,idx(j));
   end
end
toc



%Save the weighted graph to file
%filesavepth = fullfile('C:\Users\','graphnearest.mat');
csvwrite(fullfile(pwd,'graphnearest.mat'),full(A));

%fprintf( 'Nearest Neighbours Graph saved\n');




%Making unweighted - ToDo check requirement
A = spones(A);
% A = bsxfun(@times,A,-1);
% A = exp(-A);


%Save adjacency as CSV for visualization
%filesavepth = fullfile('C:\Users','dataadj.csv');
%csvwrite(filesavepth, full(A), 1, 1); %Starts writing at column offset 1
%end


%Create the Degree matrix
degrees = sum(A,2);

%Create degree sparse matrix
D = sparse(1:size(A,1), 1:size(A,2), degrees);
%D  = full(D);  %Redundant as had made sparse before ToDo Check to see exploiting sparsity

%Calculate laplacian
 L = D - A;

%Check to see all the eigenvalues nonnegative for positive definiteness
%eigvals = eigs((L+L')/2);
eigvals = eigs(L);
ids = find(not(eigvals));
if ~isempty(ids)
    warning(1,'L not PSD');
end

end