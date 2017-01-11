function[] = p2online()

%varK = [8,9,10,11]; 
%varL = [26,52,78,104,130]; 
%varN = [600,800,1000];
varK = [4];
varL = [40,80,160];
varN = [600];
%modVarN = [988];

exppath = fullfile(pwd,'p2rerunexp','EvenVsOddHighLab');

for niter = 1:length(varN)
    
for kiter = 1:length(varK)

for liter = 1:length(varL)
    
%labtype1 = 1;
%labtype2 = 2;
    
navgruns = 10;

%parameters
k = varK(kiter);
l = varL(liter);
n = varN(niter);
M = 1;
    
f_hat_mat = zeros(n,M);  %----------------------REMOVE FOR ISOLET COMMENTED


%parameter settings directory
settingstr = sprintf('n%dl%dk%d',n,l,k);
settingpath = fullfile(exppath,settingstr);

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
   
   fid = fopen(datasetfilepath,'r+');
   %C = textscan(fid, '%f', 'delimiter', '\t');
   %C = reshape(cell2mat(C),[n n]);
   %C = csvread(datasetfilepath);
   data = textscan(fid, repmat('%d',1,3),  'delimiter', '\t', 'CollectOutput', 1);
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
       
       %ISOLET ONLY
       %f_hat_mat = zeros(modVarN(niter),M);
       %ISOLET ONLY
       %[f_hat_mat(1:modVarN(niter),1)]  = optimizer('harmonicinv', partlab1, partlab2, nodelab1, nodelab2, A, D, L, modVarN(niter), f_j);
       
       
       
       
       
       
%        prediction = f_hat_mat(1:n,1);
%        predbyp2 = prediction(predrnd);
%        
%        %receive true label
%        f_j(predrnd) = labels(predrnd);
%        %update
%        if 1 == labels(predrnd)
%            partlab1 = horzcat(partlab1,labels(predrnd));
%            nodelab1 = horzcat(nodelab1,predrnd);
%            
%        else
%            partlab2 = horzcat(partlab2,labels(predrnd));
%            nodelab2 = horzcat(nodelab2,predrnd);
%        end
%        %mistakes
%        if predbyp2 ~= labels(predrnd)
%            cumumistakesperrun = cumumistakesperrun+1;
%        end
%        cumumistakesovertrials(predrnd,run) = cumumistakesperrun; 
       
    %end
    %find error of prediction
%     error = cumumistakesperrun/n;
%     cumumistakeoverruns(run) = error;
%     
    
      



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
    
end
    resultsfilename = sprintf('results-n%dl%dk%d.csv',n,l,k);
resultfilepath1 = fullfile(exppath,'combinedp2edit');
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


% 
%     %save to file
%     p2onlineresdirname = 'p2onlineresults';
%     %p2onlineresstr = sprintf('n%dl%dk%d',n,l,k);
%     p2onlinerespath = fullfile(pwd,p2onlineresdirname);
%     if ~exist(p2onlinerespath, 'dir')
%         mkdir(p2onlinerespath);
%     end
%     filenamesave1 = sprintf('resultsonlinematlab-n%dl%dk%d.csv', n, l, k);
%     filepathsave1 = (fullfile(p2onlinerespath, filenamesave1));
%     meanerror = mean(cumumistakeoverruns);
%     stderror = std(cumumistakeoverruns);
%     materror = [meanerror stderror];
%     csvwrite(filepathsave1,materror);
%     
%     filenamesave2 = sprintf('resultsonlinematlabmistakestrials-n%dl%dk%d.csv', n, l, k);
%     filepathsave2 = (fullfile(p2onlinerespath, filenamesave2));
%     %take mean of mistakes over trials
%     meanmistakestrials = max(cumumistakesovertrials,2);
%     meanmistakestrials = ceil(meanmistakestrials);
%     stddevmistakestrials = std(cumumistakesovertrials, 1, 2);
%     %mistakesovertrials = [meanmistakestrials,stddevmistakestrials];
%     csvwrite(filepathsave2,meanmistakestrials);

end
end
end