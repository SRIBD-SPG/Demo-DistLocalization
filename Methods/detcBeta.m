function [st] = detcBeta(sig, nfft, varargin)
%   This is a brief description of the function
%   计算基于 beta 散度的假设检验统计量
%   [OUTPUT1, OUTPUT2, ...] = FUNCTIONNAME(INPUT1, INPUT2, VARARGIN)
%
%   Detailed description of the function goes here, if needed.
%
%   Input Arguments:
%   x -- signal samples in time domain
%   nfft -- number of fft point
%   VARARGIN - optional input arguments
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
% 解析可选参数
p = inputParser;
addParameter(p,'betaVal', 7);
addParameter(p,'etaVal', 1);

parse(p, varargin{:});

% 提取可选参数
beta = p.Results.betaVal;
eta = p.Results.etaVal;


%%  计算功率谱熵
sig = rowVec(sig);
sig = sig - mean(sig);
window = hamming(nfft);
noverlap = int32(0.5*nfft); %数据重叠
np = nfft/2 + 1;

nend = floor(length(sig)/noverlap) * noverlap - noverlap;

[Px,~] = pwelch(sig(1:nend), window, noverlap, nfft);
Px = (Px(1:np));

Px = medfilt1(Px,13);
Pxn = length(Px)*Px./sum(Px)*eta;


tmp = ones(size(Pxn))*eta;
st = betadiv(Pxn,tmp, beta)/np;% /Para_Stru.M  /length(Ps_scale)



end

