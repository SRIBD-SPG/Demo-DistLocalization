function  [] = calTh(SENS,PF,LEN_SIG,NFFT,nNode,BANDWIDTH)
DetcResult = struct();
DetcResult.Th.ED = calDetcTh(SENS, PF, LEN_SIG, 'NumFFT', NFFT, 'method', 'ED', 'CalType', 'AVE');
% DetcResult.Th.Beta = calDetcTh(SENS, PF, LEN_SIG, 'NumFFT', NFFT, 'method', 'Beta');
DetcResult.MulTh.Xcor = calTimeCorrTh(SENS, PF,LEN_SIG);
DetcResult.MulTh.GLR = calMulDetch(SENS,PF,nNode,LEN_SIG);
DetcResult.MulTh.Energy = calMulDetch(SENS,PF,nNode,LEN_SIG,'method','ED');

save(['DetcResult' num2str(round(BANDWIDTH/1e3)) 'K.mat'],"DetcResult")
end