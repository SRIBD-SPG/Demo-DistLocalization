function [position,grid_point] = locatDirectPos(data, RxNode,varargin)

%解析可选参数
      p = inputParser;
      addParameter(p,'TruePosi',[0,0]);
%       addParameter(p,'latlim',[25.560634 26.742902]);
%       addParameter(p,'lonlim',[118.186553 120.852217]);
      addParameter(p,'method','GLRT');
      addParameter(p,'ngrid',100);
      addParameter(p,'Error',0.025);
      addParameter(p,'Height',1000);
      addParameter(p,'sampRate',120e3);
      addParameter(p,'NStart',500);
      addParameter(p,'r',1);
      addParameter(p,'noise_power_range',[-0.95,-0.95]);
      addParameter(p,'grid',[])
      p.parse(varargin{:});
      
      % 提取可选参数
%       latlim = p.Results.latlim;
%       lonlim = p.Results.lonlim;
      TruePosi = p.Results.TruePosi;
      latlim = [TruePosi(1)-p.Results.Error,TruePosi(1)+p.Results.Error];
      lonlim = [TruePosi(2)-p.Results.Error,TruePosi(2)+p.Results.Error];
      method = p.Results.method;
      ngrid = p.Results.ngrid;
      Height = p.Results.Height;
      sampRate = p.Results.sampRate;
      NStart = p.Results. NStart;
      r      = p.Results.r;
      noise_power_range = p.Results.noise_power_range;
      grid = p.Results.grid;
      nNode = size(RxNode,1);

      tgtlatv = linspace(latlim(1),latlim(2),ngrid);
      tgtlonv = linspace(lonlim(1),lonlim(2),ngrid);
      [tgtlatv1,tgtlonv1] = meshgrid(tgtlatv,tgtlonv);
      if double(isempty(grid)) == 0
         grid_point = geoshow(grid,tgtlatv1(:),  tgtlonv1(:),  'DisplayType', 'point', 'Marker', '.', 'Color', 'r', 'MarkerSize',8); 
      else
         grid_point = [];
       
      end

      st = zeros(length(tgtlatv1),length(tgtlonv1));
      delay = zeros(1,nNode);
      for m=1:length(tgtlatv)
          for n=1:length(tgtlonv)
              p1 = [tgtlatv(m), tgtlonv1(n), Height];
              for i = 1:nNode
                  p2 = [RxNode(i,1), RxNode(i,2), RxNode(i,3)];
                  delay(i) = distWGS(p1,p2)/physconst("LightSpeed");
              end
%                  TxNode = singleAntTx("TxNode",tgtlatv1(1,m), tgtlonv1(n,1), Height, P_TRA, CAR_F);
%                  for i = 1:length(RxNode)
%                      delay(i) = distance(RxNode{i}.Rx, TxNode.Tx)/physconst("LightSpeed");
%                  end
%                     delay = distance(recsite(1,1),recsite(2,1),25.560634,118.186553)/physconst("LightSpeed");     
                 % 根据时延补偿各感知节点信号
                 allNodeDataSyn= roundTimeSyn(data,sampRate,delay,NStart);
                 switch upper(method)
                     case 'GLRT'               
                         st(m,n) = detcGlrMul(allNodeDataSyn,'r',r,'noise_power_range',noise_power_range);
                     case 'DPD'
                         st(m,n) = DPDMUL(data,'delay',delay,'sampRate',sampRate,'NStart',NStart);
                 end
          end
      end
st_max=max(max(st));  
[id_xmax,id_ymax]=find(0.9*st_max<st);
tgtlatv1max=tgtlatv1(1,id_xmax).';
tgtlonv1max=tgtlonv1(id_ymax,1);
position = [mean(tgtlatv1max),mean(tgtlonv1max),Height]; 
% figure
% imagesc(tgtlatv1(1,:),tgtlonv1(:,1),abs(st.'))
% colorbar
% hold on
% plot(position(1),position(2),'k.','markersize',20)   %标记一个黑色的圆点
% axis([latlim(1) latlim(2) lonlim(1) lonlim(2) ]);
% xlabel('x (m)')
% ylabel('y (m)')
% for mm=1:length(tgtlatv1max)
%     text(tgtlatv1max(mm),tgtlonv1max(mm),['x=',num2str(tgtlatv1max(mm)),newline,'y=',num2str(tgtlonv1max(mm))]);
% end

end