function Max_chara_value = DPDMUL(data,varargin)

% DPD方法定位，对搜索空间划分网格，根据搜索位置作时延补偿，搜索空间出现峰值的位置即为辐射源的位置 
%
% MAX_CHARA_VALUE = DPDMUL(DATA, SAMPLE_RATE, TAU, NN)
%
% 输入参数:         
% DATA - complex vector 各个感知节点接收到的数据 （60000*5）
% SAMPLE_RATE - number, sample rate\Hz(120KHz)
% TAU - vector  搜索位置到各个感知节点的时延\s,向量长度等于感知节点的个数
% NN - number  截断位置对异步数据进行处理 （1000）
%
% 输出参数：
% MAX_CHARA_VALUE - complex number 
%
% Example:
% n = 60000  p = 4
% data = complex(randn(n, p), randn(n, p))
% sample_rate = 120e3
% tau = randn(1,p) NN = 1000
% Max_chara_value = DPDMUL(data,sample_rate,tau,NN)

% 解析可选参数
  p = inputParser;
  addParameter(p,'sampRate',[]);
  addParameter(p,'delay',[]);
  addParameter(p,'Nstart',500);
  parse(p,varargin{:});
  
% 提取可选参数
  sample_rate = p.Results.sampRate;
  tau = p.Results.delay;
  NN = p.Results.Nstart;
  

      node_number=size(data,2);
      D=zeros(node_number,node_number);
      for i=1:node_number
          for j=1:node_number
               choise=[i,j];
               D(i,j)=receive_data_process(tau,choise,data,sample_rate,NN);
          end
      end
      [~, v] = eig(D);
      Max_chara_value=max(diag(v));      
end

function S = receive_data_process(tau,choise,data,sample_rate,NN)
% 对不同节点的数据作时延补偿，再作内积并取平均
%
% S = RECEIVE_DATA_PROCESS(TAU,CHOISE,DATA,SAMPLE_RATE,NN)
%
% 输入参数：
% TAU - vector  搜索位置到各个感知节点的时延\s,向量长度等于感知节点的个数
% CHOISE - vector [i,j] 表示对第 i 个节点和第j个节点的数据进行处理
% DATA - complex vector 各个感知节点接收到的数据 （60000*5）
% SAMPLE_RATE - number, sample rate\Hz(120KHz)
% NN - number  截断位置对异步数据进行处理 （1000）
%
% 输出参数：
% S - complex number 对时延补偿后的两个节点的数据作内积并取平均
%
% Example:
% n = 60000  p = 4
% data = complex(randn(n, p), randn(n, p))
% sample_rate = 120e3
% tau = randn(1,p) NN = 1000
% choise = rand(1,2)
% S = RECEIVE_DATA_PROCESS(tau, choise, data, sample_rate, NN)

  move1 = round(tau(choise(1))*sample_rate);
  move2 = round(tau(choise(2))*sample_rate);
  Len = length(data(:,1)) - 2*NN;
  idx1 = NN + move1;
  idx2 = NN + move2;
  R1 = data( idx1: idx1 + Len-1, choise(1)).';
  R2 = data( idx2: idx2 + Len-1, choise(2)).';
  S = mean(R1*R2');
end