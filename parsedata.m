function [data, labels, N, d] = parseData(filename)

%Function to load the dataset, create the distance matrix
%Author: Shaona Ghosh
%Date: 10.01.2014


%No of records and values per record
novals = 257;

%Load the file :TODO change filepath to relative path
filepath = fullfile(pwd,filename);

try
    fid = fopen(filepath,'r+')
    data = textscan(fid, repmat('%f',1,novals), 'CollectOutput', 1);
    
    %Retrieve the data from the cell array in which the file was loaded
    data = data{1};
    
    %Have the true labels
    labels = data(:,1);
    
    %data(:,1) = [];
    
    [N, d] = size(data)
    
    fclose(fid);
catch
   if -1 == fid, error('Cannot open file for reading.'); end


end