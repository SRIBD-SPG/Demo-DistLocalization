function [indicSig] = alignTxTime(t_start, t_end, TxNode)
%ALIGNTXTIME 此处显示有关此函数的摘要
%   此处显示详细说明

% t_start = t_start*1e3; % ms
% t_end = t_end*1e3; 
% t_resol = t_resol*1e3; 
Tdur = t_end - t_start;
t_resol = TxNode.traPatPara.resolTime;
switch TxNode.transPattern

    case "Cyclical"
        n_resol = int32(Tdur/t_resol);
        % number of periods over the  time segment
%         n_period = floor(Tdur/TxNode.transPara.periodTime);
        
        % total number of segment over one time period
        n_dur = int32(TxNode.traPatPara.periodTime/t_resol);
        m_dur = int32(TxNode.traPatPara.durTime/t_resol);
        
        indicSig = false(n_resol, 1);
        
        for i = 1:n_dur:n_resol
            start_idx = randi(n_dur-m_dur);
            indicSig(i+start_idx:i+start_idx+m_dur-1) = true(m_dur, 1);
        end



    otherwise

end

end

