% Generalized Likelihood Ratio (GLT) Detector
function  statistics = detcGlrMul(X, varargin)
% GLR检测方法计算检验统计量
%
% STATISTICS = GLRMUL(X, R, NOISE_POWER_RANGE) 
% Reference: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=5755211
% 
% 输入参数：
% X - complex vector 各个感知节点接收到的数据 （60000*5）
% R - number 观测数据的协方差矩阵可以分解为一个秩为1的矩阵加一个对角阵 （1）
% NOISE_POWER_RANGE - vector 用于描述各个感知节点的接收噪声水平，向量长度为感知节点个数
%
% 输出参数：
% STATISTICS - number GLR检测方法，计算出来的检验统计量
%
% Example:
% n = 60000  p = 4
% X = complex(randn(n, p), randn(n, p))
% r = 1
% noise_power_range = randn(1,p) 
% statistics = GLRMUL(X, r, noise_power_range)

%解析可选参数
p = inputParser;
  
addParameter(p,'r',1);
addParameter(p,'noise_power_range',[-0.95,-0.95]);
p.parse(varargin{:});
      
% 提取可选参数
r = p.Results.r;
noise_power_range = p.Results.noise_power_range;

n = size(X, 1);
S = 1/n * X' * X;
S = S/trace(S);


max_iter = 1e2;
ptol = 1e-6;
ftol = 1e-6;
[H, Sgm, ~] = FAD(S, r, noise_power_range, max_iter, ptol, ftol);
% disp(objs);

% compute statistics
H1_M = H*H'+Sgm;
H0_N = diag(diag(S));
% aa=log(det(H0_N))+trace(H0_N\S);
% bb=log(det(H1_M))+trace(H1_M\S);
statistics = real((trace(H0_N\S)-trace(H1_M\S)+log(det(H0_N))-log(det(H1_M))));
end
