function  [data,noise,H] = channelTrans(fc, fd, samp_rate, modu_signal, snr, tau, sens,FDmax)
% 
%   Modulated Signal Transmission over Channel
%
%   [DATA,Noise,H] = channelTrans(FC, FD, SAMPLE_RATE, MODU_SIGNAL, SNR, TAU, SENS)
%   Input Arguments:
%   FC - number, center frequency\Hz(1.42e9Hz)
%   FD - number, Received Doppler frequency\Hz, FD=FC_HOP-FC+fD_TEMP
%   假设辐射源信号载载频FC(1.42e9Hz),频波动区间FC_RANGE(15e6Hz),波动区间个数FC_HOP_NUM(10)
%   跳频频点FC_HOP = FC-FC_RANGE/2+randperm(FC_HOP_NUM,1)*FC_RANGE/FC_HOP_NUM
%   FD_TEMP=-FC_HOP/c*(Vs-Vr)'*(Ps-Pr)/norm(Ps-Pr) c=physconst('lightspeed')
%   SAMPLE_RATE - number, sample rate\Hz(120KHz)
%   MODU_SIGNAL - complex Vector modulation signal
%   SNR - number, received signal-to-noise ratio at sensing node\dB
%   SENS - number, sensitivity\dBm (-100dBm) 
%
%   Output Arguments:
%   DATA - complex Vector, The signal received by the sense node(60000*1)
%   Noise - complex Vector, The noise received by the sense node(60000*1)
%   H - complex Vector, Channel parameters,Single radiation source.
%       Each element is a vector of identical complex numbers.(60000*1)
%   
%   Example:
% 
%
%   NOTE: 
%
%   Author: Zehui Zhang, Wenqiang PU 
%   Email: wenqiangpu@cuhk.edu.cn
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2023-07-09: Initial Version 0.1 
% 2023-07-27: Version 0.11, add Doppler effect 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nlevel = db2pow(sens); %noise level in miliwatt
SigGenClass = helperModClassTestChannel(...
         'SampleRate', samp_rate, ...
         'SNR', snr, ...
         'ALPHA', nlevel,...
         'PathDelays', tau ,...
         'KFactor', 4, ...
         'AveragePathGains', 0 ,...
         'DirectPathDopplerShift',fd,...
         'DopplerShift',0,...
         'MaximumDopplerShift', samp_rate/10, ...
         'MaximumClockOffset', 0, ...
         'PathGainsOutputPort',1,...
         'CenterFrequency',fc);
[data,noise,H] = SigGenClass(modu_signal); 
end
%samp_rate/10
