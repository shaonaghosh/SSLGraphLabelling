function [num,fullset] = sampledata(labels,data,n,labtype1,labtype2,nfeat)

%Finds equal sets of labels 1 and 2, then divides the sets into smaller blocks for experiments
%Author:Shaona Ghosh
%Date:10/01/2014

%Check for number of arguments
error(nargchk(1, 6, nargin'));

%For label 1 vs label 2
% id1s = find(labels==labtype1);
% id2s = find(labels==labtype2);



%For even vs odds - this section just gurantees that there are enough
%balanced labels to sample from even if instances will not be balanced
% id1s = find(not(mod(labels,2))); %For evens versus odds
% id2s = find(mod(labels,2));
% %Need to check if at least one of each label there
% digiteven = [0,2,4,6,8];
% digitodd = [1,3,5,7,9];
% 
% nodeforlab1s = zeros(1,60);
% nodeforlab2s = zeros(1,60);
% %%% For even label vs odd labels, extra code at the end
% numel = 1;
% for i = 1:length(digiteven)
% idievens = find(labels == digiteven(i)); %For evens versus odds
% idiodds = find(labels == digitodd(i));
% randideven = randi([1 length(idievens)],12);
% randideven = randideven(1,:);
% randidodd = randi([1 length(idiodds)],12);
% randidodd = randidodd(1,:)
% nodeforlab1s(numel:numel+12-1) = idievens(randideven);
% idievens(randideven);
% labels(idievens(randideven));
% nodeforlab2s(numel:numel+12-1) = idiodds(randidodd);
% idiodds(randidodd);
% labels(idiodds(randidodd));
% numel = numel+12;
% end


% %For super balanced even vs odd
nsuperbal = floor(n/10);
digits= [0,1,2,3,4,5,6,7,8,9];
minlen  = nsuperbal;
fullsetbal = zeros(n,258);
numelembal = 1;
superbalancedsamp()

    function superbalancedsamp
      
        for ii = 1:length(digits)
            ids = find(labels == digits(ii));
            len = length(ids);
            if len <= minlen
                minlen = len;
            else
                p = randperm(nsuperbal)';
                randids = ids(p);
                fullsetbal(numelembal:numelembal+nsuperbal-1,3:end) = data(randids,2:end);
               
                if not(mod(digits(ii),2))%even
                    fullsetbal(numelembal:numelembal+nsuperbal-1,1) = 1*data(randids,1);
                else
                    fullsetbal(numelembal:numelembal+nsuperbal-1,1) = -1*data(randids,1);
                end
                numelembal = numelembal + nsuperbal;
            end
        end
        fullsetbal(1:numelembal-1,1);
        fullsetbal(1:numelembal-1,1) = sign(fullsetbal(1:numelembal-1,1));
        % %For even zeros
        idzeros = find(fullsetbal(1:numelembal-1,1) == 0);
        fullsetbal(idzeros) = 1;
        fullsetbal(1:numelembal-1,1);
    end
if minlen < nsuperbal
    %Repeat the balanced sampling
    superbalancedsamp()
end
fullset = fullsetbal;

%Have balanced labels at beginning of sets
fullset(:,1) = sign(fullset(:,1));
num = n;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Checking superbalanced for isolet
% % For super balanced even vs odd
% nsuperbal = floor(n/26);
% % digits= [0,1,2,3,4,5,6,7,8,9];
% digits = 1:1:26;
% minlen  = nsuperbal;
% % fullsetbal = zeros(n,619); %process isolet added first col of labels here adding one more col for preds
% fullsetbal = zeros(nsuperbal*26,619);  % will not be size n as depends on size on n and no of letter
% numelembal = 1;
% superbalancedsamp()
% 
%     function superbalancedsamp
%       
%         for ii = 1:length(digits)
%             ids = find(labels == digits(ii));
%             len = length(ids);
%             if len <= minlen
%                 minlen = len;
%             else
%                 randidforsel = randperm(length(ids))';
%                 randid = ids(randidforsel);
%                 %select nsuperbal no of such randomly sorted ids for letter
%                 %ii
%                 randids = randid((1:nsuperbal)');
%                 fullsetbal(numelembal:numelembal+nsuperbal-1,3:end) = data(randids,2:end);
%                
%                 if not(mod(digits(ii),2))%even
%                     fullsetbal(numelembal:numelembal+nsuperbal-1,1) = 1*data(randids,1);
%                 else
%                     fullsetbal(numelembal:numelembal+nsuperbal-1,1) = -1*data(randids,1);
%                 end
%                 numelembal = numelembal + nsuperbal;
%             end
%         end
%         fullsetbal(1:numelembal-1,1);
%         fullsetbal(1:numelembal-1,1) = sign(fullsetbal(1:numelembal-1,1));
%         % %For even zeros
%         idzeros = find(fullsetbal(1:numelembal-1,1) == 0);
%         fullsetbal(idzeros) = 1;
%         fullsetbal(1:numelembal-1,1);
%     end
% if minlen < nsuperbal
%     %Repeat the balanced sampling
%     superbalancedsamp()
% end
% fullset = fullsetbal;
% 
% %Have balanced labels at beginning of sets
% fullset(:,1) = sign(fullset(:,1));
% %num = n; % not n anymore
% [num,d] = size(fullset);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% if the num of ids are less than n/2
% numel1 = numel(id1s);
% numel2 = numel(id2s);
% 
% 
% if (numel1 < (n/2)) || (numel2 < (n/2))
%     warning('SampleData', 'Cannot find labelled sets of given sizes');
%     %Uniformly sample each labels proprotional to as many available
%     
%     %Choose the minimum of the two sizes of unbalanced labeled sets
%     %available
%     numelmin = min(numel1,numel2);
%     numelmin = floor(numelmin/2);
%     index = randperm(numelmin)';
%     randid1s = id1s(index,1);
%     index = randperm(numelmin)';
%     randid2s = id2s(index,1);
%     numel1 = numelmin;
%     numel2 = numelmin;
% else
%     %Randomly sample n/2 1's and 2's with replacement
%     index = randperm(n/2)';
%     randid1s = id1s(index,1);
%     index = randperm(n/2)';
%     randid2s = id2s(index,1);
%     numel1 = n/2;
%     numel2 = n/2;
% end
% 
% 
% set1 = zeros(numel1,nfeat); %First column for true labels, second for predicted labels --->
% set2 = zeros(numel1,nfeat);
% 
% %set1(1:1+numel1-1,3:end) = data(randid1s,:);%always commented
% %set2(1:1+numel2-1,3:end) = data(randid2s,:);
% 
% 
% 
% %For even and odd - this ensures that the labels at least will be balanced
% %even if instances are not when sampled in samplelabels -Even odd enable
% % randid1s(1:length(nodeforlab1s)) = nodeforlab1s';
% % randid2s(1:length(nodeforlab2s)) = nodeforlab2s';
% 
% 
% 
% set1(1:numel1,3:end) = data(randid1s,2:end);
% set2(1:numel1,3:end) = data(randid2s,2:end);
% 
% set1(1:numel1,1) = 1*data(randid1s,1); %First type of label in first column of dataset
% set2(1:numel1,1) = -1*data(randid2s,1); %Second type of labels
% 
% %For even zeros -Enable for even vs odd
% % idzeros = find(set1(1:numel1,1) == 0);%Evens will be in set1
% % set1(idzeros) = 1;
% 
% 
% set1(1:10,1:10)
% set2(1:10,1:10)
% 
% 
% %TEST CODE
% %set2 = set2(1:numel1,:);
% 
% %Also return full set
% fullset = zeros(numel1+numel2,nfeat);
% fullset = vertcat(set1,set2);
% fullset(:,1) = sign(fullset(:,1));
% 
% %fullset(100:115,1:10);
% %Permute the set :Todo
% num = numel1+numel2;

end