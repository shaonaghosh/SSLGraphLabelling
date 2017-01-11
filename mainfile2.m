function[] = mainfile2()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function: Main control flow, does the avergaing runs over the experiments
%on different sizes of datasets with harmonic energy minimization using
%different methods and max flow-min cut respecting labelling methods
%
%Author: Shaona Ghosh
% Date: 11-01-2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%clear all;
close all;

%Seed the random number generator 
rng(555);

%parse the file 
%[data, labels, N, d] = parsedata('dtrain123.dat');
[data, labels, N, d] = parsedata('ziptrain.dat');

%Set parameters
%No of nodes to sample
nofsets = 800 %[200,400,800,1600];

%No of labelled nodes
lset = [10,20,40,80,160,320]; %For different labeled sets
%l = 6%50;


%No of neighbours
%K = [3,4,5,6,7,8,9];
k = 3;
labtype1 = 3;
labtype2 = 6;
%Type of optimization methods
M = 5;


%Experiments
navgruns = 10;


c = 1;
%while c <= length(nofsets)
while c <= length(lset)
%while c <= length(K)
    
    %for diff n sets
    %n = nofsets(c);
    
    %for diff l sets
    n = nofsets(1);
    l = lset(c);
    
    %for diff K's
%     n = nofsets(1);
%     k = K(c);
    
    %prediction matrix
    f_hat_mat = zeros(n,M);%will be larger than actual n
    
    for run = 1: navgruns
        
        %divide into random sets
        [num,fullset] = sampledata( labels, data, n, labtype1, labtype2 );   %Arguments:labels,data,n, labtype1, labtype2
        
        %distance computation
        [distancemat1] = distances( fullset, num, 258 );  %Third parameter default to 258
        
        %build Graph
        [A,D,L] = builddatasetgraph(num, fullset, distancemat1, k);
        
        %sample two classes randomly
        [partlab1, partlab2, nodelabi, nodelabj, f_j] = samplelabels(num,fullset(:,1), l);
        
        %seminorm interpolation
        [f_hat_mat(1:num,1)]  = optimizer('harmonicinv', partlab1, partlab2, nodelabi, nodelabj, A, D, L, num, f_j);
        
        [f_hat_mat(1:num,2)]  = optimizer('harmonicquad', partlab1, partlab2, nodelabi, nodelabj, A, D, L, num, f_j);
        
        [f_hat_mat(1:num,3)]  = optimizer('maxflowlinear', partlab1, partlab2, nodelabi, nodelabj, A, D, L, num, f_j);
        
        [f_hat_mat(1:num,4)]  = optimizer('harplusmaxflowcons', partlab1, partlab2, nodelabi, nodelabj, A, D, L, num, f_j);
        
        [f_hat_mat(1:num,5)]  = optimizer('harplusmaxflowobj', partlab1, partlab2, nodelabi, nodelabj, A, D, L, num, f_j, 2*(num^2));
        
        %Save results
        saveresults(run, num, l, f_hat_mat, fullset(:,1) );
        
        run
    end
    f_j'
    f_hat_mat
    %Increment counter
    c = c+1;
end


end
%Plot results 
%plotresults(noavgruns,noofsetscnt,modes); 






