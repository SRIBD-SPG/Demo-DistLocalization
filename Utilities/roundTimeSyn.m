%% 根据采样率以及辐射源位置到感知节点时延进行时延补偿，按照四舍五入计算需移动的采样点个数，实现各感知节点接收数据大致同步
% 输入参数
% allNodeData  struct 各感知节点接收到的数据 allNodeData.signal、allNodeData.noise、allNodeData.ch    
%noise  各感知节点所处的背景噪声
%sample_rate number  (120KHz)  采样率\Hz 
%delay  辐射源位置到各感知节点的传输延迟
% NN  信道延迟会对数据进行补零，因此截取数据的一部分进行同步处理
%输出参数
%data_syn   noise_syn  经过同步处理后的各感知节点的数据
%  move  各感知节点数据需移动的采样点数以粗略达到时间同步
function data_syn = roundTimeSyn(allNodeData,sample_rate,delay,NN)
% 根据采样率以及辐射源位置到感知节点时延进行时延补偿，按照四舍五入计算需移动的采样点个数，实现各感知节点接收数据大致同步
%
% ALLNODEDATASYN = ROUNDTIMESYN(ALLNODEDATASYN,SAMPLE_RATE,DELAY,NN)
%
% 输入参数
% ALLNODEDATA - struct（complex vector） 各感知节点接收到的数据包括噪声、信号+噪声、信道  allNodeData.signal、allNodeData.noise、allNodeData.ch
% SAMPLE_RATE - number  采样率\Hz  (120e3)
% DELAY - vector  辐射源到各个感知节点的时延\s,向量长度等于感知节点的个数
% NN - number  截断位置对异步数据进行处理 （1000）
%
% 输出参数
% ALLNODEDATASYN - struct （complex vector）经过四射五入同步处理后的数据 allNodeDataSyn.signal（60000*1），allNodeDataSyn.noise（60000*1）
% 
% Example:
% n = 60000  p = 3
% allNodeData.signal = complex(randn(n, p), randn(n, p))
% allNodeData.noise = complex(randn(n, p), randn(n, p))
% allNodeData.ch = complex(randn(1,p),randn(1,p))*ones(n,p)
% sample_rate = 120e3
% delay = randn(1,p)
% NN = 1000
% allNodeDataSyn = roundTimeSyn(allNodeData,sample_rate,delay,NN)
data = allNodeData;
move = round((delay - delay(1))*sample_rate);
Len = length(data(:,1)) - 2*NN;

data_syn = zeros(Len,size(data,2));

for t=1:size(data,2)
    idx = NN + move(t);
    data_syn(:,t) = data( idx: idx + Len-1, t);
end


end