function [numSelect, valueSelect] = removeConsective(num, value)

% input:
%         *num: an index serial, e.g. [1 2 3 6 7 9],
%         *value: the value w.r.t num, e.g. [0.1, 0.2, 0.3, 4, 1, pi];
% This function removes the consective index 
% and keep the one with smallest value:
%        *num_rm: result num, e.g. [1 7 9];
%        *value_rm = result value, e.g. [0.1, 1, pi];

% sort the vectors first
[numSort, numSortIdx] = sort(num);
valueSort = value(numSortIdx);

indDiffJump = find(diff(numSort)~=1);
indJump = [1, indDiffJump+1];

% construct cell
indEnd = length(numSort);
% sections = cell(length(indJump));
numSelect = [];
valueSelect = [];
for i = 1:length(indJump)
    if i == length(indJump)
        sections_i = indJump(i):indEnd;
    else
        sections_i = indJump(i):(indJump(i+1)-1);
    end
    % find the smallest 
    [minI, idxMinI] = min(valueSort(sections_i));
    idxSelectI = sections_i(idxMinI);
    numSelect = [numSelect, numSort(idxSelectI)];
    valueSelect = [valueSelect, valueSort(idxSelectI)];
end

