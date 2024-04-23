function [index] = findindex(v,val)
%FINDINDEX Finds index of v where value is closest to val
%   Detailed explanation goes here

index = find(abs(v-val) == min(abs(v-val)),1);
end

