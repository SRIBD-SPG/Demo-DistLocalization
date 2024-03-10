function d = DPD(S, delay, fs, varargin)
%   This is a brief description of the function
%   DPD方法定位，对搜索空间划分网格，根据搜索位置作时延补偿，搜索空间出现峰值的位置即为辐射源的位置
%
%
%   Input Arguments:
%   S - complex vector 各个感知节点接收到的数据, n-by-m
%   fs - number, sample rate\Hz(120KHz)
%   delay - vector  搜索位置到各个感知节点的时延\s,向量长度等于感知节点的个数
%   Nstart - number  截断位置对数据进行处理
%
%   Output Arguments:
%   sfft -- averged fft spectrum
%   ...
%
%   Example:
%   [out1, out2] = functionName(in1, in2, 'option1', value1, 'option2', value2);
%
%   See also: relatedFunctionName
%
%   Author: Wenqiang PU
%   Email: wenqiangpu@cuhk.edu.cn

% 解析可选参数
p = inputParser;
addParameter(p,'Nstart',1);
parse(p,varargin{:});

% 提取可选参数
NN = p.Results.Nstart; 


delay = rowVec(delay);
Ns = size(S, 1); % number of samples
NFFT = 2^(nextpow2(Ns)-1); % number of DFT points

idx = (0:NFFT-1)';
wk = 2*pi*fs/NFFT;
Tp = exp(-1j*wk*idx*delay);
Y = fft(S(NN:NN+NFFT-1, :), NFFT, 1);

L = Y.* Tp;
D = L'*L;

[~, v] = eig(D);
d = max(diag(v));

end