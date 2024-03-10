function [st] = detcEnergy(sig, nfft, varargin)
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

sig = rowVec(sig);
sig = sig - mean(sig);
% 解析可选参数
p = inputParser;
addParameter(p, 'CalType', 'AVE');

parse(p, varargin{:});
% 提取可选参数
calstr = p.Results.CalType;

switch upper(calstr)

    case 'MAXLEN'
        nfft = 2^nextpow2(length(sig));
        nf = nfft/2+1;
        sfft = (fft(sig, nfft));
        st = mean(abs(sfft(1:nf)).^2);

    case 'AVE'
        nf = nfft/2+1;
        sfft = aveFFTSpec(sig, nfft);
        st = mean(abs(sfft(1:nf)).^2);

    otherwise
        nf = nfft/2+1;
        sfft = (fft(sig, nfft));
        st = mean(abs(sfft(1:nf)).^2);

end





end

