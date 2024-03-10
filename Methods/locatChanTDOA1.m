function  position = locatChanTDOA1(varargin)
% chan 两步法时差定位 (二维空间)
%
% POSITION_TDOA = CHANTDOAMUL(AMDELAY, ALL_NODE_POSITION)
% 
% 输入参数：
% AMDELAY - vector 由互模糊函数估计各个感知节点相对于第一个感知节点的时延差
% ALL_NODE_POSITION - matric 感知节点的位置，行数与感知节点的个数有关，列数为2
%
% 输出参数：
% POSITION_TDOA - two-dimensional vector  由Chan算法估计出来的辐射源的位置
%
% Example:
% n = 5 
% amdelay = randn(1,n)
% all_node_position = randn(n,2)
% Position_TDOA=ChanTDOAMUL(amDelay,all_node_position)

% 解析可选参数
  p = inputParser;
  addParameter(p,'amDelay',[]);
  addParameter(p,'all_node_position',[]);
  addParameter(p,'Height',1000);
  parse(p,varargin{:});
  
% 提取可选参数
  amDelay = p.Results.amDelay;
  Height = p.Results.Height;
  all_node_position_temp = p.Results.all_node_position; % 经纬度坐标
%   temp = all_node_position_temp(:,3)+linspace(-20,20,size(all_node_position_temp,1)).';
%   all_node_position_temp = [all_node_position_temp(:,1:2),temp];
  lla0 =[all_node_position_temp(1,1:2),0];
  all_node_position1 = lla2enu(all_node_position_temp, lla0, "flat");

%   all_node_position =locationToCartesian(all_node_position_temp(:,1),all_node_position_temp(:,2),all_node_position_temp(:,3));
  all_node_position = all_node_position1(:,1:2);

  
       N = size(all_node_position,1);
       c = physconst('lightspeed');
%      delay_theory_ori = zeros(node_number,1);
%      data_length=size(all_node_data.signal,1);
%      tau=(-(data_length-1):(data_length-1))/sample_rate;
%         for refinx = 2:node_number
%             x1 = all_node_data.signal(:,1);
%             x2 = all_node_data.signal(:,refinx);
%             xcorrMat = xcorr(x2,x1,'biased');
%             [~,ind] = max(abs(xcorrMat));
%             delay_theory_ori(refinx) =tau(ind);
% %         delay_theory_ori(refinx) = (length(all_node_data.signal(:,1))-ind )/sample_rate;
% %        (  tau_1(refinx) - tau_1(1)  )*600e3
%         end
        delay_theory_ori = amDelay;
        BS = all_node_position - repmat(all_node_position(1,:),N,1);
        R_theory_ori = zeros(1,N);
        Kj_ori       = zeros(1,N);
        
        for i = 2:N
            R_theory_ori(i) = delay_theory_ori(i)*c;
%             Kj_ori(i)       = BS(i,1)^2 + BS(i,2)^2 +BS(i,3)^2; 
            Kj_ori(i)       = sum(BS(i,:).^2,2);
        end
        
        H_tdoa_ori = BS(2:N,:);
        C_tdoa_ori = -1*R_theory_ori(2:N)';
        D_tdoa_ori = 0.5*(Kj_ori(2:N)-R_theory_ori(2:N).^2)';
        K_tdoa_ori = inv(H_tdoa_ori'*H_tdoa_ori) *(H_tdoa_ori')*C_tdoa_ori;
        Q_tdoa_ori = inv(H_tdoa_ori'*H_tdoa_ori) *(H_tdoa_ori')*D_tdoa_ori;
        
        %计算得到abc三个参数
        as_ori     = sum(K_tdoa_ori.^2) -1;
        bs_ori     = sum(K_tdoa_ori.*Q_tdoa_ori);
        cs_ori     = sum(Q_tdoa_ori.^2);
%         Position_TDOA = zeros(size(all_node_position_temp,2),2);
%         position =Position_TDOA;
        if (bs_ori^2-as_ori*cs_ori)>=0
           R1_ori     = (-bs_ori + sqrt(bs_ori^2-as_ori*cs_ori))/as_ori;
           xyChan_ori = K_tdoa_ori.*R1_ori + Q_tdoa_ori;
           Position_TDOA = xyChan_ori+ all_node_position(1,:).';
%          R1_ori     = (-bs_ori - sqrt(bs_ori^2-as_ori*cs_ori))/as_ori;
%          xyChan_ori = K_tdoa_ori.*R1_ori + Q_tdoa_ori;
%          Position_TDOA(:,2) = xyChan_ori+ all_node_position(1,:).';
%          xyChan_ori + all_node_position(1,:)'
          %计算估计得到XY坐标点
                   
%            X_theory_ori = xyChan_ori(1) + all_node_position(1,1);
%            Y_theory_ori = xyChan_ori(2) + all_node_position(1,2);
%            Z_theory_ori = xyChan_ori(3) + all_node_position(1,3);
%            Position_TDOA=[X_theory_ori,Y_theory_ori,Z_theory_ori]; 
        else
           Position_TDOA= 1e9*ones(size(all_node_position,2),1);
        end
        position_1 = [Position_TDOA(1),Position_TDOA(2), all_node_position1(1,3)];
        position = enu2lla(position_1,lla0,'flat');
%         position = cartesianToLocation(Position_TDOA(1),Position_TDOA(2),Position_TDOA(3));
%         position(3) = Height;
%         position(:,2) =cartesianToGeoLocation(Position_TDOA(1,2),Position_TDOA(2,2),Position_TDOA(3,2));
end