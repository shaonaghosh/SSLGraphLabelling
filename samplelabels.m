function[partlab1, partlab2, nodesforlab1s, nodesforlab2s, f_j] = samplelables(num,fullset,l)

%function: to sample two labels from the two sets randomly
%Author:Shaona Ghosh
%Date:10.01.2014


%Check for number of arguments
error(nargchk(1, 4, nargin'));


%%For label 1 vs label 2
%%Sample partial labels randomly with replacement

%Without replacement
% nodesforlab1s = randperm(floor(num/2),floor(l/2));  
% 
% %With replacement 
% nodesforlab2s = randi([floor(num/2)+1,num],floor(l/2));
% %nodesforlab2s = nodesforlab2s(1,:)
% nodesforlab2s = unique(nodesforlab2s)';
% while(length(nodesforlab2s) < floor(l/2))
%     nodesforlab2s = randi([floor(num/2)+1,num],floor(l/2));
%     %nodesforlab2s = nodesforlab2s(1,:)
%     nodesforlab2s = unique(nodesforlab2s)';
% end
% nodesforlab2s = nodesforlab2s(1:floor(l/2));



%For even vs odd labels
%sample balanced labels
% nodesforlab1s = 1:1:60; %there will be 5 even
% nodesforlab2s = floor(num/2)+1:1:floor(num/2)+1+59; %will be 5 odd 




%Superbalanced trimmed
nsuperbal = floor(num/10);
numelebal = 1;
nodesforlab1s = zeros(1,l/2);
nodesforlab2s = zeros(1,l/2);
%ii = 1:1:13;
quot = (l/2)/5;
for i = 1:quot:(l/2)
    %LAZY CODE BELOW TODO REMOVE
    for j = 1:quot
        nodesforlab1s(i+(j-1)) = numelebal+(j-1);
        nodesforlab2s(i+(j-1)) = numelebal+nsuperbal+(j-1);
    end
    numelebal = numelebal + 2*nsuperbal;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Superbalanced trimmed  for Isolet
% nsuperbal = floor(num/26);
% numelebal = 1;
% nodesforlab1s = zeros(1,l/2);%Usually for superbalanced no of labels should be 26 in total or at least multiple of 26
% nodesforlab2s = zeros(1,l/2);
% %ii = 1:1:13;
% %quot = (l/2)/5;
% quot = floor((l/2)/13);%Try 13 assuming 5 was for 5 odd 5 even in total no of digits in orig so here 13 odd and 13 even
% for i = 1:quot:(l/2)
%     for j = 1:quot
%         nodesforlab1s(i+(j-1)) = numelebal+(j-1);%one for odd
%         nodesforlab2s(i+(j-1)) = numelebal+nsuperbal+(j-1);%one for even
%     end
%     numelebal = numelebal + 2*nsuperbal;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



partlab1 = fullset(nodesforlab1s(1),1)
partlab2 = fullset(nodesforlab2s(1),1)


%partially labelled vector
f_j = zeros(num,1);
f_j(nodesforlab1s,1) =  fullset(nodesforlab1s,1);
f_j(nodesforlab2s,1) =  fullset(nodesforlab2s,1);

f_j = sign(f_j);

% sprintf('Labels selected \n');


end