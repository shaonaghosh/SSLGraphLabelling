function[] = saveresults(runno, n, l, ypred, ytrue)
%Author:Shaona Ghosh

%Check for number of arguments
error(nargchk(1, 5, nargin'));

%Crop prediction matrix
if size(ypred,1) ~= n
    ypred(size(ytrue,1)+1:end,:) = [];  %ypred was allocated outside loop in main
end

field1 = 'runno';  value1 = runno;
field2 = 'num';  value2 = n; 
field3 = 'nopartlab';  value3 = l;
field4 = 'ypred';  value4 = ypred; 
field5 = 'ytrue'; value5 = ytrue;

results = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5);


probsizestr  = sprintf('nodes%d',n); %nodes and labels can interchange
%probsizestr  = sprintf('nodes%d',l);

%cd(fullfile('C:\Users\Shaona\Documents\MATLAB'));  %need to remove hardcoding
%mkdir('datasetsize',probsizestr);
filename = sprintf('resultfile-%d-%d.mat',l,runno);
filepath = (fullfile(pwd, 'datasetsize', probsizestr, filename)); %ToDo find way to create directory
if ~exist(fullfile(pwd, 'datasetsize', probsizestr), 'dir')
  mkdir(fullfile(pwd, 'datasetsize', probsizestr));
end
save(filepath, 'results');






end