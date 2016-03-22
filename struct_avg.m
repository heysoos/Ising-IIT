function [ avgMat ] = struct_avg( struct, fSrch1, srchStr1, fSrch2, srchStr2, fAvg )
%   struct_avg
% 
%   Input:
%   struct - struct to look at
% 
%   fSearch - the field whose values will be looked at to decide which
%   struct indices will be averaged.
% 
%   fAvg - the field whose indices will be averages
% 
%   srchStr - string which the field fSrch will be searched for

%   Output:
%   avgMat = 2D matrix of averages over the fields fAvg through indices
%   given b fSrch
%   

%   Function that averages the fields (2D matrices) of a struct over the
%   indices of another field. The indices are picked depending on the value
%   of the values in the 2nd field.

avgMat = [];
fSrchVal1 =  cat(1,struct.(fSrch1));
for i = 1:numel(struct)
    fSrchVal2{i,:} =  struct(i).(fSrch2);
end
avgInd = find((ismember(fSrchVal1,srchStr1,'rows') + ismember(fSrchVal2,srchStr2))==2);

numDim = ndims(struct(1).(fAvg));
% Context: timeseries data
% Averages are done only over the matrix with the smallest length
minDim = min(cell2mat((arrayfun(@(x) size(x.(fAvg)), struct,'UniformOutput',false))'));

% avgMat = cat(numDim + 1, struct(avgInd).(fAvg)(1:minDim(1),1:minDim(2)));

for i = 1:length(avgInd)
%     avgInd(i)
    avgMat = cat(numDim + 1, avgMat, struct(avgInd(i)).(fAvg)(1:minDim(1),1:minDim(2)));
%     size(avgMat)
end

avgMat = mean(avgMat,numDim + 1);

end

