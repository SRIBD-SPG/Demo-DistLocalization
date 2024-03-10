%% 根据采样率以及辐射源位置到感知节点时延进行时延补偿，按照四舍五入计算需移动的采样点个数，实现各感知节点接收数据大致同步
% 输入参数
% data  各感知节点接收到的数据
%noise  各感知节点所处的背景噪声
% sample_rate  采样率\Hz
%delay  辐射源位置到各感知节点的传输延迟
% NN  信道延迟会对数据进行补零，因此截取数据的一部分进行同步处理
%输出参数
%data_syn   noise_syn  经过同步处理后的各感知节点的数据
function allNodeDataSyn = roundDelayAndFdSyn(allNodeData,sample_rate,delay,fd,NN)
% 根据时延以及多普勒频率对接收数据作补偿
%
% ALLNODEDATASYN = ROUNDDELAYANDFDSYN(ALLNODEDATA,SAMPLE_RATE,DELAY,FD,NN)
%
% 输入参数：
% ALLNODEDATA - struct（complex vector） 各感知节点接收到的数据包括噪声、信号+噪声、信道 allNodeData.signal、allNodeData.noise、allNodeData.ch (60000*1)
% SAMPLE_RATE - number  采样率\Hz  (120e3)
% DELAY - vector  辐射源到各个感知节点的时延\s,向量长度等于感知节点的个数
% FD - vector 各个感知节点的多普勒频率，向量长度与感知节点的个数有关\Hz
% NN - number  截断位置对异步数据进行处理 （1000）
%
% 输出参数：
% ALLNODEDATASYN - struct （complex vector）经过时间补偿与多普勒补偿的数据 allNodeDataSyn.signal（60000*1），allNodeDataSyn.noise（60000*1）
%
% Example
% n = 60000  p = 4
% allNodeData.signal = complex(randn(n, p), randn(n, p))
% allNodeData.noise = complex(randn(n, p), randn(n, p))
% allNodeData.ch = complex(randn(1,p),randn(1,p))*ones(n,p)
% sample_rate = 120e3
% delay = randn(1,p)
% fd = randn(1,p)
% NN = 1000
% allNodeDataSyn = roundTimeSyn(allNodeData,sample_rate,delay,NN)

data = allNodeData.signal;
noise = allNodeData.noise;
% move = round((delay - delay(1))*sample_rate);
Len = length(data(:,1)) - 2*NN;
t=(0:Len-1)./sample_rate;
data_syn_delay = zeros(Len,size(data,2));
noise_syn_delay = zeros(Len,size(data,2));
data_syn_delay_and_fd = zeros(Len,size(data,2));
noise_syn_delay_and_fd = zeros(Len,size(data,2));
for i=1:size(data,2)
    idx = NN + round(delay(i)*sample_rate);
    data_syn_delay(:,i) = data( idx: idx + Len-1, i);
    data_syn_delay_and_fd(:,i) =data_syn_delay(:,i).*exp(-1j*2*pi*fd(i)*t.');
    noise_syn_delay(:,i) = noise( idx: idx + Len-1, i);
    noise_syn_delay_and_fd(:,i)=noise_syn_delay(:,i).*exp(-1j*2*pi*fd(i)*t.');
end
allNodeDataSyn.signal = data_syn_delay_and_fd;
allNodeDataSyn.noise = noise_syn_delay_and_fd;

end