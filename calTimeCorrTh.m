function th = calTimeCorrTh(nlevel, pf, n, varargin) 
% INPUT
% nlevel -- noise level in dBm
% pf -- false alarm prob
% n -- number of samples

% OUTPUT
% th -- threshold

% 解析可选参数
p = inputParser;
addParameter(p,'MonteCarlo',2000);

parse(p,varargin{:});

% 提取可选参数
MonteCarlo = p.Results.MonteCarlo;

stCorr = zeros(MonteCarlo,1);
for i = 1:MonteCarlo
    noise = sqrt(db2pow(nlevel)/2)*(randn(2,n) + 1j*randn(2,n));
    stCorr(i) = detcTimeCorr(noise(1,:),noise(2,:));
end
SortstCorr = sort(stCorr,'descend');
th = SortstCorr(round(MonteCarlo*pf));
end