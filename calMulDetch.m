function  th = calMulDetch(nlevel,pf,m,len,varargin)
% INPUT
% nlevel -- noise level in dBm
% pf -- false alarm prob
% len -- number of samples
% m -- number of nodes
% method -- detection method, i.e., xxx

% OUTPUT
% th -- threshold
% 解析可选参数
  p = inputParser;
  addParameter(p,'MonteCarlo',2000);
  addParameter(p, 'method', 'GLR'); 
  addParameter(p,'r',1);
  addParameter(p,'noise_power_range',[-0.95,-0.95]);
  parse(p,varargin{:});
  
% 提取可选参数
  MonteCarlo = p.Results.MonteCarlo;
  method = p.Results.method;
  r = p.Results.r; 
  noise_power_range = p.Results.noise_power_range;
  

stAGM = zeros(MonteCarlo,1);
stED = zeros(MonteCarlo,1);
stEMR = zeros(MonteCarlo,1);
stGLR = zeros(MonteCarlo,1);
stMME = zeros(MonteCarlo,1);
stSLE = zeros(MonteCarlo,1);


switch upper(method)
    case 'AGM'
       for i = 1:MonteCarlo
           noise = sqrt(db2pow(nlevel)/2)*(randn(len,m) + 1j*randn(len,m));
           stAGM(i) = detcAgmMul(noise);
       end 
       SortstAGM = sort(stAGM,'descend');
       th= SortstAGM(round(MonteCarlo*pf));
        
    case 'ED'
        for i = 1:MonteCarlo
            noise = sqrt(db2pow(nlevel)/2)*(randn(len,m) + 1j*randn(len,m));
            stED(i) = detcEnergyMul(noise);
%             stED(i) = GLR(noise,1,[-0.95,-0.95]);
        end
       SortstED = sort(stED,'descend');
       th= SortstED(round(MonteCarlo*pf));
        
    case 'EMR'
        for i = 1:MonteCarlo
            noise = sqrt(db2pow(nlevel)/2)*(randn(len,m) + 1j*randn(len,m));
            stEMR(i) = detcEmrMul(noise);
        end
       SortstEMR = sort(stEMR,'descend');
       th= SortstEMR(round(MonteCarlo*pf));
       
    case 'GLR'
        for i = 1:MonteCarlo
            noise = sqrt(db2pow(nlevel)/2)*(randn(len,m) + 1j*randn(len,m));
            stGLR(i) = detcGlrMul(noise,'r',r,'noise_power_range',noise_power_range);
        end
       SortstGLR = sort(stGLR,'descend');
       th= SortstGLR(round(MonteCarlo*pf));   
       
    case 'MME'
        for i = 1:MonteCarlo
            noise = sqrt(db2pow(nlevel)/2)*(randn(len,m) + 1j*randn(len,m));
            stMME(i) = detcMmeMul(noise);
        end
       SortstMME = sort(stMME,'descend');
       th= SortstMME(round(MonteCarlo*pf));

    case 'SLE'
        for i = 1:MonteCarlo
            noise = sqrt(db2pow(nlevel)/2)*(randn(len,m) + 1j*randn(len,m));
            stSLE(i) = detcSleMul(noise);
        end
       SortstSLE = sort(stSLE,'descend');
       th= SortstSLE(round(MonteCarlo*pf));
       
        
end 

end