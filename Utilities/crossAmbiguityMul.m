function [fd_diff,delay_diff] = crossAmbiguityMul(all_node_data, samp_rate, fd_max, ngrid)
%Cross ambiguity function, using the correlation of data from each
%receiving node to calculate time difference and frequency difference
%The positions corresponding to the peaks of the mutual ambiguity function are the time and frequency 
%differences of the measurements
%
% 
%   [FD_DIFF,DELAY_DIFF] = CROSSAMBIGUITY(ALL_NODE_DATA,SAMP_RATE,fd_range)
%   Input Arguments:
%   ALL_NODE_DATA - The data transmitted by each perception node to the fusion center. a T × N matrix,
%   T represents the length of data received by the sense node, which is related to the sampling rate 
%   and observation time.(obs_time/samp_rate).N represents the number of sensing nodes.
%   SAMP_RATE - number, sample rate\Hz(120KHz) 
%   FD_RANGE - Frequency difference range of search element.(Maximum search range：-B:B,B indicates the bandwidth of the signal)
%
%   Output Arguments
%   FD_DIFF - number,Frequency difference calculated using the mutual ambiguity function \Hz
%   DELAY_DIFF - number,Time difference calculated using the mutual ambiguity function  \s

% Initialize time difference frequency difference
delay_diff = zeros(1,size(all_node_data,2));
fd_diff = zeros(1,size(all_node_data,2));

% the length of data received by the sense node
N = size(all_node_data,1);

% obs_time
t=(0:N-1)/samp_rate;
fd_range = linspace(-fd_max, fd_max, ngrid);
% Time difference search range
tau_range=(-(N-1):(N-1))/samp_rate;
% fd_range=-100:2:100;
% fd_range=-BW:BW/10:BW-BW/10;
% fd_range=samp_rate/N*(-(N-1):(N-1));
for m = 2:size(all_node_data,2)
    % Value Initialization for Mutual ambiguity Function
    X = zeros(length(fd_range),length(tau_range));
    for b = 1:length(fd_range)
        x1 = all_node_data(:,1).*exp(1j*2*pi*fd_range(b)*t.');
        x2 = all_node_data(:,m);
        X(b,:) = xcorr(x2,x1,'biased');
    end
[indx,indy] = find(abs(X) == max(max(abs(X))));
delay_diff(m) = tau_range(indy);
fd_diff(m) = fd_range(indx);

% figure
% mesh(abs(X));
% colorbar
% hold on
% plot(indy,indx,'k.','markersize',20)   %标记一个黑色的圆点 

end
end
