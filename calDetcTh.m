function th = calDetcTh(nlevel, pf, n, varargin)
% INPUT
% nlevel -- noise level in dBm
% pf -- false alarm prob
% n -- number of sample point
% method -- detection method, i.e., psd energy, entropy, beta detect

% OUTPUT
% th -- threshold
% 解析可选参数
p = inputParser;
addParameter(p, 'CalType', 'AVE');
addParameter(p,'MonteCarlo',2000);
addParameter(p,'etaVal',1);
addParameter(p,'betaVal',7);
addParameter(p,'NumFFT',1024);
addParameter(p, 'method', 'ED');

parse(p,varargin{:});

% 提取可选参数
eta = p.Results.etaVal;
beta = p.Results.betaVal;
MonteCarlo = p.Results.MonteCarlo;
method = p.Results.method;
calstr = p.Results.CalType;
nfft = p.Results.NumFFT;


tmp = zeros(MonteCarlo,1);

switch upper(method)
    case 'ED'
        for i = 1:MonteCarlo
            noise = sqrt(db2pow(nlevel)/2)*(randn(1,n) + 1j*randn(1,n));
            tmp(i) = detcEnergy(noise, nfft, 'CalType', calstr); %calDetcStatisc(noise, n, 'method',method,'CalType',calstr );
        end

    case 'ENTROPY'
        for i = 1:MonteCarlo
            noise = sqrt(db2pow(nlevel)/2)*(randn(1,n) + 1j*randn(1,n));
            tmp(i) = detcEntropy(noise, nfft); %calDetcStatisc(noise, nfft,'method',method);
        end

    case 'BETA'
        for i = 1:MonteCarlo
            noise = sqrt(db2pow(nlevel)/2)*(randn(1,n) + 1j*randn(1,n));
            tmp(i) = detcBeta(noise, nfft, 'etaVal', eta, 'betaVal', beta); %calDetcStatisc(noise, nfft,'method',method ,'eta',eta,'beta',beta);
        end

end

sort_st = sort(tmp,'descend');
th = sort_st(round(MonteCarlo*pf));

end
