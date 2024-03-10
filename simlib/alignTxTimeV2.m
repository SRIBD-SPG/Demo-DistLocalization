function [indSig] = alignTxTimeV2(TxNode)
%ALIGNTXTIME 此处显示有关此函数的摘要
%   此处显示详细说明

T_perid = TxNode.traPatPara.periodTime;
t_pice  = TxNode.traPatPara.resolTime;
t_dur   = TxNode.traPatPara.durTime;

% number of samples per pice
NSampFr = int32(t_pice * TxNode.sampRate);
switch TxNode.transPattern

    case "Cyclical"

        % total number of pice over one period
        n_resol                               = int32(T_perid/t_pice);
        n_dur                                 = int32(t_dur/t_pice);

        indicSig                              = false(n_resol, 1);
        if n_dur < n_resol
            start_idx                             = randi(n_resol-n_dur);
        else
            start_idx                             = 1;
        end

        indicSig(start_idx:start_idx+n_dur-1) = true(n_dur, 1);

        indSig                                = repmat(indicSig, [1, NSampFr]);
        indSig                                = reshape(indSig', [numel(indSig), 1]);


    otherwise

end

end

