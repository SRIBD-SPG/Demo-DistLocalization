function [nsc,nov,nff] = myPwelch1(x,fs,titleStr)
%MYPWELCH1 根据信号x长度计算一组[nsc,nov,nff],用于画信号功率谱。
%   pwelch画信号功率谱，并返回pwelch使用的参数[nsc,nov,nff]
    Nx = length(x);
    nsc = floor(Nx/4.5);
    nov = floor(nsc/2);
    nff = max(1024,2^nextpow2(nsc));
    pwelch(x, hamming(nsc), nov, nff, fs,'centered');
    title(titleStr);
end

