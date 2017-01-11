function [distancemat1] = distances( data, n, d )

%function to compute the euclidean distances
%Author: Shaona Ghosh
% Date: 10.01.2014

%Check for number of arguments
error(nargchk(1, 3, nargin'));

%Test Code
% data = [1,2,3,6,2;4,5,6,2,8;7,8,9,1,3;2,0,5,4,1]
% n = size(data,1)
% K = 2

if (nargin < 3)
    d = 256; %Default dimensionality 
    %d = 617;
elseif (3 == nargin)
    %reserve first two columns for labels
    data = data(:,3:end);
    
end


%Distance matrix
distancemat1 = zeros(n,n);
%distancemat2 = zeros(n,n);
 
% sigmad  = 380;%Reference :Semi supervised learning using Gaussian fields and Harmonic functions
% sigmad = 380^2;
sigmad = 1; %We need unweighted graph


%Find squared euclidean distance between each pair - Faster than faster
%alternative below
tic
for i = 1:n
     nodei = data(i,:);
     for j = 1:n  %TODO remove for loop use bsxfun, at least remove one loop
        nodej = data(j,:);
        interm = (nodei - nodej).^2;
        interm = interm/sigmad;
        %distancemat1(i, j) = sum((nodei - nodej).^2); 
        distancemat1(i, j) = sum(interm);
     end
end
toc


%Commented as not faster than method 1 even if should be - need to check
%data = data';
%tic
%Alternative - calculate euclidean distances fast using norm and dot
%product
%D = bsxfun(@plus,dot(data,data,1)',dot(data,data,1))-2*(data'*data)
%toc

%Test code
distancemat1;


%Alternative method3 for large dataset, can be reduced time complexity if
%not space, need to test for large dataset though
% 
% tic
% datacpy = zeros(3,3);
% sums = zeros(n,1);
% for i = 1:n
%    nodei  = data(i,:);
%    datacpy = bsxfun(@minus,data,nodei);
%  
%    datacpy = bsxfun(@power,datacpy,2);
     %datacpy = bsxfun(@,datacpy,2)
     
%    sums = sqrt(sum(datacpy,2));
     %sums = sum(datacpy,2);
%    distancemat2(i,:) = sums';
% end
% 
% distancemat2
% toc
%end

%Save the diatances to file:
%filesavepth = fullfile('C:\Users\','distance.mat');
save(fullfile(pwd,'distance.mat'), 'distancemat1');

fprintf( 'Distances saved\n')

end
