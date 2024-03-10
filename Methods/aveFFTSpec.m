function [sfft] = aveFFTSpec(sig, nfft)
%   This is a brief description of the function
%   计算时域信号FFT频谱能量均值
%   [OUTPUT1, OUTPUT2, ...] = FUNCTIONNAME(INPUT1, INPUT2, VARARGIN)
%
%   Detailed description of the function goes here, if needed.
%
%   Input Arguments:
%   x -- signal samples in time domain
%   nfft -- number of fft point
%   varargin - optional input arguments
%
%   Output Arguments:
%   st -- test statistic value
%   ...
%
%   Example:
%   [out1, out2] = functionName(in1, in2, 'option1', value1, 'option2', value2);
%
%   See also: relatedFunctionName
%
%   Author: Wenqiang PU
%   Email: wenqiangpu@cuhk.edu.cn

% Overwrite default values with user-specified values
n = floor(length(sig)/nfft);
S = zeros(n, nfft);

for i = 1:n
    S(i,:) = fft(sig(1+(i-1)*nfft:i*nfft), nfft);
end

sfft = mean(S);

end