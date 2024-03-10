function [st] = detcTimeCorr(x1, x2)
% INPUT
% x1 -- signal samples of node 1
% x2 -- signal samples of node 2
% nfft -- number of fft point

% OUTPUT
% st -- test statistic value
xc = xcorr(x1,x2,'biased');

st = max(abs(xc));

end