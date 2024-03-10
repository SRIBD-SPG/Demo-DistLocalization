function [ss, snr, delay, fd] = chFreeSpace(TxNode, RxNode)
%CHANNELFREESPACE 此处显示有关此函数的摘要
%   此处显示详细说明
% 
%  function of relaizing freespace channel model
%  
%
%   Example:
% 
%
%   NOTE: 
%
%   Author: Wenqiang PU 
%   Email: wenqiangpu@cuhk.edu.cn
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2023-07-16: Initial Version 0.1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\

PM = propagationModel("freespace");
ss = sigstrength(RxNode.Rx, TxNode.Tx, PM,'Map','none'); % dBm
snr = ss - RxNode.Rx.ReceiverSensitivity; % dB

p1 = [RxNode.Rx.Latitude,RxNode.Rx.Longitude,RxNode.Rx.AntennaHeight];
p2 = [TxNode.Tx.Latitude,TxNode.Tx.Longitude,TxNode.Tx.AntennaHeight];
delay = distWGS(p1,p2)/physconst("LightSpeed");
% delay = distance(RxNode.Rx, TxNode.Tx)/physconst("LightSpeed"); % second

% Doppler shift, Hz
fd = -TxNode.Tx.TransmitterFrequency/physconst("LightSpeed")...
    *(TxNode.currentPosition-RxNode.currentPosition)'...
    *(TxNode.currentSpeed-RxNode.currentSpeed)...
    /norm(TxNode.currentPosition-RxNode.currentPosition);

end

