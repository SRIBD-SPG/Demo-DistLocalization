function st = Position(data,RxNode,varargin)

%解析可选参数
      p = inputParser;
      
      addParameter(p,'latlim',[22.758955 22.817597]);
      addParameter(p,'lonlim',[114.186553 114.259567]);
      addParameter(p,'method','GLR');
      addParameter(p,'ngrid',100);
      addParameter(p,'P_TRA',100);
      addParameter(p,'CAR_F',1.42e9);
      addParameter(p,'Height',1000);
      addParameter(p,'sampRate',120e3);
      addParameter(p,'NStart',500);
      addParameter(p,'r',1);
      addParameter(p,'noise_power_range',[-0.95,-0.95]);
      p.parse(varargin{:});
      
      % 提取可选参数
      latlim = p.Results.latlim;
      lonlim = p.Results.lonlim;
      method = p.Results.method;
      ngrid = p.Results.ngrid;
      P_TRA =  p.Results.P_TRA ;
      CAR_F =  p.Results.CAR_F;
      Height = p.Results.Height;
      sampRate = p.Results.sampRate;
      NStart = p.Results. NStart;
      r      = p.Results.r;
      noise_power_range = p.Results.noise_power_range;
      
      tgtlatv = linspace(latlim(1),latlim(2),ngrid);
      tgtlonv = linspace(lonlim(1),lonlim(2),ngrid);
      [tgtlatv1,tgtlonv1] = meshgrid(tgtlatv,tgtlonv);
      st = zeros(length(tgtlatv1),length(tgtlonv1));
      delay =zeros(1,length(RxNode));
      for m=1:length(tgtlatv1)
          for n=1:length(tgtlonv1)
                 TxNode = singleAntTx("TxNode",tgtlatv1(1,m), tgtlonv1(n,1), Height, P_TRA, CAR_F);
                 for i = 1:length(RxNode)
                     delay(i) = distance(RxNode{i}.Rx, TxNode.Tx)/physconst("LightSpeed");
                 end
                         
                 % 根据时延补偿各感知节点信号
                 allNodeDataSyn= roundTimeSyn(data,sampRate,delay,NStart);
                 switch upper(method)
                     case 'GLR'               
                         st(m,n) = detcGlrMul(allNodeDataSyn,'r',r,'noise_power_range',noise_power_range);
                     case 'DPD'
                         st(m,n) = DPDMUL(allNodeDataSyn,'delay',delay,'sampRate',sampRate,'NStart',NStart);
                 end
          end
      end
st_max=max(max(st));  
[id_xmax,id_ymax]=find(st==st_max);
tgtlatv1max=tgtlatv1(1,id_xmax).';
tgtlonv1max=tgtlonv1(id_ymax,1);
position = [mean(tgtlatv1max),mean(tgtlonv1max)].'; 
figure
imagesc(tgtlatv1(1,:),tgtlonv1(:,1),abs(st.'))
colorbar
hold on
plot(position(1),position(2),'k.','markersize',20)   %标记一个黑色的圆点
axis([latlim(1) latlim(2) lonlim(1) lonlim(2) ]);
xlabel('x (m)')
ylabel('y (m)')
% for mm=1:length(tgtlatv1max)
%     text(tgtlatv1max(mm),tgtlonv1max(mm),['x=',num2str(tgtlatv1max(mm)),newline,'y=',num2str(tgtlonv1max(mm))]);
% end

end