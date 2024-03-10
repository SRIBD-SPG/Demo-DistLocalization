function [v] = rowVec(v)
%ROWVEC 此处显示有关此函数的摘要
%   此处显示详细说明

if iscolumn(v)
    v = v.';
end

end

