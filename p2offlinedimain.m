function[] = p2online()

varK = [8,9]; 
varL = [8,16,32,64,128]; 
varN = [200];
lab1 = [1];%[1,2,3];%,4,6];
lab2 = [2];%[2,3,8];%,7,9];
% varK = [7,8,9,10,11,12];
% varL = [16,32,64,128,256];
% varN = [600,800,1000,1200,1400];
%modVarN = [988];



exppath = fullfile(pwd,'july1408');%uncomment when not newsgroup
%exppath = fullfile(pwd,'20news-bydate','matlab','class1Vs2','newsgroupn400');

for class = 1:length(lab1)

for niter = 1:length(varN)
    
for kiter = 1:length(varK)

for liter = 1:length(varL)
    
% % %labtype1 = 1;  
% % %labtype2 = 2;
    
navgruns = 10;

%parameters
k = varK(kiter);
l = varL(liter);
n = varN(niter);
M = 1;

labtype1 = lab1(class);
labtype2 = lab2(class);
    
f_hat_mat = zeros(n,M);  %----------------------REMOVE FOR ISOLET COMMENTED


%parameter settings directory
classstr = sprintf('class%dVs%d', labtype1, labtype2);
settingstr = sprintf('n%dl%dk%d',n,l,k);
settingpath = fullfile(exppath,classstr,settingstr); %uncomment for not
%newsgroup
%settingpath = fullfile(exppath,settingstr);



%error logging information
cumumistakesovertrials = zeros(n,navgruns);
cumumistakeoverruns = zeros(navgruns,1);
errormat = zeros(navgruns,M,1);

for run = 1: navgruns
    
   %cumumistakesperrun = 0
   %read the graph and labels
   datasetfilename = sprintf('%s%d%d%s', 'datasetgraph', n, run, '.csv'); 
   
   
   %for isolet
   %datasetfilename = sprintf('%s%d%d%s', 'datasetgraph', modVarN(niter), run, '.csv'); 
   
   
   
   
   
   
   datasetfilepath = fullfile(settingpath,datasetfilename); 
   try
   fid = fopen(datasetfilepath,'r+');
   %C = textscan(fid, '%f', 'delimiter', '\t');
   %C = reshape(cell2mat(C),[n n]);
   %C = csvread(datasetfilepath);
   data = textscan(fid, repmat('%d',1,3),  'delimiter', '\t', 'CollectOutput', 1);
   catch exception
       datasetfilepath
       continue;
   end
   
   data = data{1};
  
   %drop the weights
   E = data(:,1:2);
   E = double(E);
   A = sparse(E(:,1),E(:,2),1);
   degrees = sum(A,2);

    %Create degree sparse matrix
    D = sparse(1:size(A,1), 1:size(A,2), degrees);
    %D  = full(D);  %Redundant as had made sparse before ToDo Check to see exploiting sparsity

    %Calculate laplacian
    L = D - A;
    eigvals = eigs(L);
    ids = find(not(eigvals));
    if ~isempty(ids)
        warning(1,'L not PSD');
    end
    
    %read all the labels and already labelled nodes
    labsetfilename = sprintf('%s%d%d%d%s', 'labset', n, run, 0, '.csv'); 
    
       
    
    
    %for isolet
    
    %labsetfilename = sprintf('%s%d%d%d%s', 'labset', modVarN(niter), run, 0, '.csv'); 
    
    
    
    
    
    labsetfilepath = fullfile(settingpath,labsetfilename);
    labels = csvread(labsetfilepath);
    labels = labels(:,2);
    labels(labels==0) = -1;
    labels(labels==1) = 1;
    
    availabfilename = sprintf('%s%d%d%s', 'trainingset', n, run,'.csv'); 
    
    
    
    
    %for isolet
    %availabfilename = sprintf('%s%d%d%s', 'trainingset', modVarN(niter), run,'.csv'); 
    
    
    
    
    
    
    availlabfilepath = fullfile(settingpath, availabfilename);
    availlab = csvread(availlabfilepath);
    %availlab = availlab(:,1);
    nlab = length(availlab);
    
    nlabhalf = nlab/2;
    availlabels = availlab(:,2);
    availabnodes = availlab(:,1);
    
    partlab1 = availlabels(1:nlabhalf)'
    partlab2 = availlabels(nlabhalf+1:nlab)'
    nodelab1 = availabnodes(1:nlabhalf)'
    nodelab2 = availabnodes(nlabhalf+1:nlab)'
    
    f_j = zeros(1,length(labels));
    f_j(nodelab1) = partlab1;
    f_j(nodelab2) = partlab2;
    
    
    %online prediction
    %for predrnd = 1: length(f_j)
       %find an unlabelled node
%        if f_j(predrnd) ~= 0
%            continue;
%        end
       %call harmonic semi norm to predict
       
       
       
       
       [f_hat_mat(1:n,1)]  = optimizer('harmonicinv', partlab1, partlab2, nodelab1, nodelab2, A, D, L, n, f_j);
       
      



   pred = f_hat_mat(1:n,1);
   %ISOLET
   %pred = f_hat_mat(1:modVarN(niter),1);
   
   
    ytrue = labels;

%Calculate accuracy
pos = (sign(ytrue) == +1);     
        
        p = nnz(pos);
neg = (sign(ytrue) == -1);
noofneg = nnz(neg);

tpos = (sign(pred) == +1);

%test code
%testpos = (pred == 1); %want true positive not all positive
%testposcnt = nnz(testpos);

length(tpos);
length(pos);
addpos = bsxfun(@plus,tpos,pos);
trposcnt = (addpos == 2);
fposcnt = (addpos == 1);
fposcnt2 = nnz(fposcnt)
trposcnt = nnz(trposcnt)


tneg = (sign(pred) == -1);

addneg = bsxfun(@plus,tneg,neg);
trnegcnt = (addneg == 2);
fnegcnt = (addneg == 1)
fnegcnt2 = nnz(fnegcnt);
trnegcnt = nnz(trnegcnt)


error = (trposcnt + trnegcnt)/(p+noofneg);


errormat(run,M,1) = error;
fclose(fid);  
end
    resultsfilename = sprintf('results-n%dl%dk%d.csv',n,l,k);
    resultfilepath1 = fullfile(exppath,classstr,'combinedp2edit');
if ~exist(resultfilepath1, 'dir')
    mkdir(resultfilepath1);
end
resultsfilepath = fullfile(resultfilepath1,resultsfilename);
errormat

    kerrormat = errormat(:,:,1);
    meanmat(1,:) = mean(kerrormat,1);
    stdmat(1,:) = std(kerrormat,1);
    


resultmat = [meanmat,stdmat];

csvwrite(resultsfilepath, resultmat);




end
end
end
cd ('C:\Users\Shaona\Documents\MATLAB\');
end
end