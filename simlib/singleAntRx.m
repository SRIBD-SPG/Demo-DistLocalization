function [RxNode] = singleAntRx(name, txlat, txlon, alt, sensitivity)
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
% 2023-07-16: Initial Version 0.1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

paraRx = struct();
paraRx.name = name;
paraRx.initialLocation = [txlat; txlon];
paraRx.initialAltitude = alt;

paraRx.initialPosition =  locationToCartesian(txlat, txlon, alt);
paraRx.initialSpeed = 0 * randn(3,1);
paraRx.motionDensity = 0 * rand();

paraRx.Rx = rxsite('Name', name,  'Latitude', txlat, 'Longitude', txlon, ...
    'AntennaHeight', alt, 'Antenna', 'isotropic', 'ReceiverSensitivity', sensitivity);
RxNode = rxNode(paraRx);
end

