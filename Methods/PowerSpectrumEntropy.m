function [Px, Pxn, Pfft, Efft, Hfft, Ex, Hx] = PowerSpectrumEntropy(Smat,Para, nFFT)

Np = Para.nfft/2+1;%Para_Stru.nfft;
Nf = nFFT/2+1;%Para_Stru.nfft;

% Np = Para.nfft;%Para_Stru.nfft;
% Nf = nFFT;%Para_Stru.nfft;

Smat = Smat - mean(Smat); % normalize 

SmatM = reshape(Smat, [nFFT/Para.nfft, Para.nfft]);
%%  计算功率谱熵
window = hamming(Para.nfft);
noverlap = 0.5*Para.nfft; %数据重叠
[Px,~] = pwelch(Smat,window,noverlap,Para.nfft,Para.Fs);
Px = (Px(1:Np)); 
Px = medfilt1(Px,13);
Pxn = length(Px)*Px./sum(Px);


%%  计算频谱 FFT
Sfft = (fft(Smat, nFFT));
Pfft = abs(Sfft(1:Nf)).^2;

% SfftM = (fft(SmatM, Para.nfft,2));
SfftM = fft(Smat, Para.nfft);
Ex = sum(abs(SfftM(1:Np)).^2);

%% metric
% 基于频谱 FFT 
Efft = sum(Pfft); % energy
tmp = Pfft/sum(Pfft);
Hfft = sum(tmp.*(log2(tmp))); % entropy

% 基于功率谱
% Ex = sum(Px); % energy

tmp = Px/sum(Px);
Hx = sum(tmp.*(log2(tmp))); % entropy

end
