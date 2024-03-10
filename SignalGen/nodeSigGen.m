%%  各感知节点接收数据生成
% 输入参数
%modu_type string  调制方式选择 modulationTypes = categorical(["BPSK", "QPSK", "8PSK", ...
%  "16QAM", "64QAM", "PAM4", "GFSK", "CPFSK", ...
%  "B-FM", "DSB-AM", "SSB-AM"]);  
%baund_rate number (20k)  波特率\baund
%samp_rate  number (120k) 采样率\Hz
%observe_time number (0.5s)  观测时间\s
%tau  矢量（长度为感知节点个数） 传输延迟\s
%snr  矢量（长度为感知节点个数） 感知节点处接收信噪比\dB
%alpha  矢量(长度为感知节点个数) 代表每个节点的噪声水平 
% fc   载频
% 输出参数
% allNodeData（struct）  感知节点处接收的数据 包含背景噪声
function sigGen = nodeSigGen(fc, samp_rate, modu_signal, tau, fd, snr, sens)
       
       [sig,noise,h] = channelTrans(fc, fd, samp_rate, modu_signal, snr, tau, sens);
        
       sigGen.signal = sig;
       sigGen.noise = noise;
       sigGen.ch = h;
end