function [TxNode] = singleAntTx(name, txlat, txlon, alt, pow, txfreq)
%SINGLEANTTX 此处显示有关此函数的摘要
%   此处显示详细说明
% 
%  function of relaizing tx of signle antenna
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
% 2023-07-09: Initial Version 0.1 
% 2023-09-09: Version 0.2, revise parameter initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


geo_locat = [txlat; txlon;alt]; % geo location
cart_locat =  locationToCartesian(txlat, txlon, alt); %Cartesian location

% Single Antenna Tx
Tx = txsite( 'Name', name, 'Latitude', txlat, 'Longitude', txlon, 'TransmitterFrequency', txfreq,...
    'AntennaHeight', alt, 'TransmitterPower', pow, 'Antenna', 'isotropic');


TxNode = txNode('Name', name, 'initialGeoLocation', geo_locat,...
    'initialPosition', cart_locat, 'Tx', Tx);
end

