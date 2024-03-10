function [dis] = distWGS(p1, p2)
%DISTWGS 此处显示有关此函数的摘要
%   此处显示详细说明

pt1 = locationToCartesian(p1(1), p1(2), p1(3));
pt2 = locationToCartesian(p2(1), p2(2), p2(3));

dis = norm(pt1 - pt2);
end

