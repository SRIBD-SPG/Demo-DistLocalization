function  TimeFrequencySpectrum(x,RxNode,ax)
 
 Nx = length(x);
 nsc = floor(Nx/8);
 nov = floor(nsc/2);
 nff = max(1024,2^nextpow2(nsc));
 fs = RxNode.sampRate;
 [s,f,t,p]=spectrogram(x,hamming(nsc),nov,nff,fs,'yaxis','centered');
[tt,ff] = meshgrid(t,f);
surf(ax,tt,ff,10*log10(abs(p)),'EdgeColor','none');
colorbar(ax);
view(ax,2);
end