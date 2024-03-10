function  [] = scDynamicTx(ax,ED,xcor,glr,posi,txt_re,txt_re_1,button,Front,Back)


load(['map30Km.mat']);
txt = [];
Contain = struct();
%% 演示节奏控制
% num_step         = 100;  % 仿真段落数

Tsim_start = 0; % start time
Tsim_step        = 1;  % 仿真的每段持续时间
% Tsim_end         = num_step * Tsim_step;  % 仿真总时间

pauseTime        = 0.01;  % 播放的时间间隔（与仿真效果无关）

% fprintf('仿真起始时间: %d 秒, 每帧时间长度 %f 秒 \n', Tsim_start, Tsim_step);
txt1 = sprintf('仿真起始时间: %d 秒, 每帧时间长度 %f 秒 \n', Tsim_start, Tsim_step);
txt = [{txt1};txt];

%% scenario related basic settings
CAR_F = 1.42e9; % carraier frequency, Hz
BANDWIDTH = 500e3; %Hz
SAMP_RATE = BANDWIDTH/2*8; % sampling rate, Hz
FD_MAX = 400/physconst('LightSpeed')*CAR_F; % maximum doppler shift, 100 m/s
NFFT = 1024; % number of fft point
SENS = -90; % dBm
nNode = 6;

%%%%%%%  高信噪比场景下的仿真 P_TRA = 5w  ,低信噪比场景下的仿真 P_TRA = 0.8*1e-2
P_TRA = 0.8*1e-2; % power in watt
PF = 1e-2; % false alarm probability
LEN_SIG = Tsim_step/10 * SAMP_RATE;
load(['DetcResult' num2str(round(BANDWIDTH/1e3)) 'K.mat'])

% fprintf('辐射信号载频 %f GHz, 信号采样率 %d kHz \n', CAR_F/1e9, SAMP_RATE/1e3);
txt1 = sprintf('辐射信号载频 %f GHz, 信号采样率 %d kHz \n', CAR_F/1e9, SAMP_RATE/1e3);
txt = [{txt1};txt];
%% 定义路径点，包括经纬度和高度
%%%%%%低信噪比辐射源运动轨迹%%%%%%%%%
% if P_TRA<1e-2
   path_enemy = [29.94897,122.64,1000;
       29.915698, 122.539497,1000;
       29.891439,122.507425, 1000;
       29.87896, 122.507981,1000;
       29.845698, 122.569497,1000;
      29.825458, 122.636762,1000;    
    ];
% else
% %%%%%% 高信噪比辐射源运动轨迹%%%%%%%%%
% path_enemy = [29.94897,122.64,1000;
%      29.917698,122.61597,1000;
% %      29.915698, 122.539497,1000;
%      29.925698,122.539497,1000;
%      29.9962521,122.507425, 1000;
%      30.0062521,122.567425, 1000;
%      30.042482, 122.626762,1000;    
% ];
% end


% 移动速度
speed_enemy = 300;  % 1000 m/s

% 生成轨迹
track_enemy = uniformMotion(path_enemy, -1, speed_enemy, Tsim_step);
% track_enemy = circshift(track_enemy,-49);

% 更新 step 数量
num_step = size(track_enemy, 1)-1;

% 检查 "设定的仿真时间下" 走完 "轨迹" 的平均速度
ave_speed = calculateTotalDistance(path_enemy) / (num_step * Tsim_step);
% fprintf("信号源的大概平均速度为：%2f m/s\n", ave_speed);
% fprintf("仿真模拟时间总计：%f second\n", Tsim_step * num_step);

txt1 = sprintf('信号源的大概平均速度为：%2f m/s', ave_speed);
txt = [{txt1};txt];
txt1 = sprintf('仿真模拟时间总计：%f second', Tsim_step * num_step);
txt = [{txt1};txt];

txt_re.Value = txt ;

% latlim = [25.560634 26.742902];
% lonlim = [118.186553 120.852217];
% latlim = [22.758955 22.817597];
% lonlim = [114.186553 114.259567];
latlim = [29.785795 30.088752];
lonlim = [122.322216 122.664003];


%% 生成演示地图
% % 创建一个地图容器: figure
% imageHeight = 800;
% imageWidth = 800;
% fig = figure(1);
% pos = [0 0 imageHeight imageWidth];
% set(fig, 'Position', pos);




%% 设置无人机节点构型
center = [22.783614 114.207968 1000];
drone_positions =[
29.860016,122.393303,1000;
29.90309,122.401622,1000;
29.93189,122.397188,1000;
29.947428,122.359096,1000;
29.89047,122.465459,1000;
29.92461,122.464597,1000
];
% drone_positions =[29.860016,122.393303,0;29.947428,122.359096,0;29.879872,122.341273,0;29.90309,122.401622,0;29.936224,122.429192,0;29.960258,122.416831,0];

% center = [25.783614 119.207968 0];
% drone_positions = [25.95,119.65,0;26.24,119.65,0;26.5,119.25,0;26.3,119.35,0;25.9,119.35,0;25.7,119.25,0];

% drone_positions = [25.95,119.65,0;26.24,119.65,0;26.5,119.25,0;26.3,119.35,0;26.1,119.4,0;25.9,119.35,0;25.7,119.25,0];
% drone_positions = droneFormation(center, nNode, 1/80, 'Trapezoid', 50);
% drone_positions =[22.77,114.25 0;22.78,114.2,0;22.79,114.25,0;22.81,114.21,0;22.8,114.19,0];
% drone_positions = [22.7669,114.208,0;22.7725,114.19,0;22.7781,114.2,0;22.7836,114.22,0;...
%     22.7892,114.2,0;22.7947,114.209,0;22.8003,114.208,0];

% figure;hold on;
% grid on; set(gca,'GridLineStyle','--');
% plot(track_enemy(:,1),track_enemy(:,2),'h')
% hold on
% plot(drone_positions(:,1),drone_positions(:,2),'o')
% hold on
% plot(25,119.9,'*')
% % axis([22.758955 22.817597 114.186553 114.259567 ])
% % axis equal
% xlabel('x (m)');
% ylabel('y (m)');
% legend('Emitter Position','Sense Node Position ','Location', 'northeast', 'Interpreter', 'latex')

% 创建一个地图

geoshow(ax,A, R, 'DisplayType', 'texturemap');
alley_points = geoshow(ax,drone_positions(:, 1),  drone_positions(:, 2),  'DisplayType', 'point', 'Marker', '.', 'Color', 'r', 'MarkerSize', 10);
enemy_points = geoshow(ax,track_enemy(1, 1),   track_enemy(1, 2), 'DisplayType', 'point', 'Marker', '.', 'MarkerEdgeColor', 'blue', 'MarkerSize', 10);

pause(pauseTime);
%% setup transmiter
TxNode = singleAntTx("TxNode", 22.806537, 114.250512, 1000, P_TRA, CAR_F);

% baseband signal setting at Tx
% baseband signal parameters
TxNode.moduType = "BPSK"; % modulation type
TxNode.baudRate = BANDWIDTH/2; % baud rate in Hz
TxNode.sampRate = SAMP_RATE; % sample rate in Hz

TxNode.transPattern = "Cyclical";
TxNode.traPatPara.periodTime = 0.1; % time period of one circle
TxNode.traPatPara.durTime = 1*TxNode.traPatPara.periodTime; % duration time of the signal
TxNode.traPatPara.resolTime = 0.25*TxNode.traPatPara.periodTime; % reselu time of the signal

Contain.static.moduType = TxNode.moduType;
Contain.static.baudRate = TxNode.baudRate;
Contain.static.fc = CAR_F;
%% modulator class -- modulator with symbol as input
SFS = round(TxNode.sampRate/TxNode.baudRate); % samples per symbol
ModuSig = helperModClassGetModulator(TxNode.moduType, SFS, TxNode.sampRate);

% source class -- symbol sequence
NST = round(TxNode.traPatPara.periodTime*TxNode.sampRate); % samples per frame
SigSrc = helperModClassGetSource(TxNode.moduType, SFS, NST, TxNode.sampRate);
%% setup receiver

rxlat = drone_positions(:,1);
rxlon = drone_positions(:,2);
rxalt = 1000*ones(nNode, 1);
recsite = [rxlat'; rxlon'; rxalt'];


RxNode = cell(nNode,1);
RxPosi = cell(nNode,1);
for i = 1:nNode
    name = "RxNode " + (i);
    txlat = recsite(1, i);
    txlon = recsite(2, i);
    alt = recsite(3, i);
    sensitivity = SENS; % dBm

    RxNode{i} = singleAntRx(name, txlat, txlon, alt, sensitivity);
    RxNode{i}.sampRate = SAMP_RATE;
    Contain.static.position{i,1} =  [txlat txlon alt];
    Contain.static.sampleRate{i,1} = SAMP_RATE;    
end



%% forloop for time

% Tsim_start = 0; % start time
Tsim_end = Tsim_step*num_step; % end time
% Tsim_step = 0.05; % time per simulation frame

curseg_time = 0; % indicate the end time of current node signal
len_seg = int32(TxNode.traPatPara.periodTime*TxNode.sampRate); % length of signal segment
len_step = int32(Tsim_step/10*TxNode.sampRate); % length of simul s

% temp signal data matrix, need to be saved
sigtmp_rx = zeros(nNode, len_step); % temp signal matrix of tx
sigtmp_tx = zeros(1, len_seg); % temp signal matrix of rx





% DetcResult = struct(); % save detection result
LocaResult = struct(); % save localization result
LocaResult.Chan = cell(1,num_step);
LocaResult.GLR = cell(1,num_step);
% LocaResult.TdoaTrack = cell(1,num_step);
% LocaResult.TdoaTrack1 = cell(1,num_step);
LocaResult.GLRTrack = cell(1,num_step);
% LocaResult.GLRTrack1 = cell(1,num_step);
% LocaResult.ChanTrack = cell(1,num_step);
% LocaResult.ChanTrack1 = cell(1,num_step);




Contain.Dynamic.SNR = cell(nNode,num_step);
Contain.Dynamic.Fd =  cell(nNode,num_step);
Contain.Dynamic.TxPosi = cell(1,num_step);
Contain.Dynamic.TxVeloci = cell(1,num_step);
% Contain.Data = cell(1,num_step);
Contain.Dynamic.delay_real_diff = cell(1,num_step);
Contain.Dynamic.fd_real_diff = cell(1,num_step);
Contain.Dynamic.amDelay = cell(1,num_step);
Contain.Dynamic.amFd = cell(1,num_step);
% % 计算检测门限
% DetcResult.Th.ED = calDetcTh(SENS, PF, LEN_SIG, 'NumFFT', NFFT, 'method', 'ED', 'CalType', 'AVE');
% % DetcResult.Th.Beta = calDetcTh(SENS, PF, LEN_SIG, 'NumFFT', NFFT, 'method', 'Beta');
% DetcResult.MulTh.Xcor = calTimeCorrTh(SENS, PF,LEN_SIG);
% DetcResult.MulTh.GLR = calMulDetch(SENS,PF,nNode,LEN_SIG);
% DetcResult.MulTh.Energy = calMulDetch(SENS,PF,nNode,LEN_SIG,'method','ED');


%计算检验统计量
% Chan_points = [];
grid_point = [];
for t = Tsim_start:Tsim_step:Tsim_end - Tsim_step

%     disp(['当前仿真时间区间: ', num2str(t),' to ', num2str(t+Tsim_step), ' second']);
    txt1 = ['当前仿真时间区间: ', num2str(t),' to ', num2str(t+Tsim_step), ' second'];
    txt = [{txt1};txt];
    %% update location and velocity of Tx node
    countT = int32(t/Tsim_step) + 1;
    while button.Value == 0
        pause(0.1);
        if button.Value ==1
           break;
        end

    end
    
    % 计算当前点的速度（笛卡尔坐标系）
    alley_speed = calculateVelocity(drone_positions, drone_positions, Tsim_step);
    enemy_speed = calculateVelocity(track_enemy(countT, :), track_enemy(countT+1, :), Tsim_step);
%     fprintf("信号源的当前速度为：%f m/s\n", norm(enemy_speed));
    txt1 = sprintf('信号源的当前速度为：%f m/s\n', norm(enemy_speed));
    txt = [{txt1};txt];
    txt_re.Value = txt ;
    % update Tx node geo location and velocity
    TxNode.setGeoLocation(track_enemy(countT, 1),   track_enemy(countT, 2), track_enemy(countT, 3));
    TxNode.setSpeed(enemy_speed)
    Contain.Dynamic.TxPosi{1,countT} = [track_enemy(countT, 1), track_enemy(countT, 2), track_enemy(countT, 3)];
    Contain.Dynamic.TxVeloci{1,countT} = enemy_speed;
 
    %% generate baseband signal data
    % check simul time frame within signal segment
    if t + Tsim_step <= curseg_time-1e-8 % within signal segment

        % find signal index
        idxtmp = int32(TxNode.sampRate * (t + Tsim_step - curseg_time + ...
            TxNode.traPatPara.periodTime) + 1); % start index
        % add signal to data class
        for i = 1:nNode
            sigtmp_rx(i,:) = rowVec(RecSigTmp(i).sig(idxtmp:idxtmp+len_step-1));
 %             %% ----------- 单节点检测统计量 -----------
            [st] = detcEnergy(sigtmp_rx(i,:), NFFT, 'CalType', 'AVE');
            DetcResult.RxNode(i).st_energy(countT) = st;
%             [st] = detcBeta(sigtmp_rx(i,:), NFFT);
%             DetcResult.RxNode(i).st_beta(countT) = st;
            Contain.Dynamic.sigstr{i,countT} = Contain.Dynamic.sigstr{i,countT-1};
            Contain.Dynamic.SNR{i,countT} = Contain.Dynamic.SNR{i,countT-1};
            Contain.Dynamic.delay{i,countT} = Contain.Dynamic.delay{i,countT-1};
            Contain.Dynamic.Fd{i,countT} = Contain.Dynamic.Fd{i,countT-1};
        end
%         Contain.Data{1,countT} = sigtmp_rx;

    else % beyound signal segment

        curseg_time = curseg_time + TxNode.traPatPara.periodTime;


        % generate new baseband signal received at Rx
        sigsrc = SigSrc(); % generate new symbol signal
        modu_signal = ModuSig(sigsrc); % generate modulated signal

        % generate signal with desired tx pattern
        indSig = alignTxTimeV2(TxNode);
        sigtmp_tx = modu_signal .* indSig;

        RecSigTmp = struct();
        fd_real = zeros(1,nNode);
        delay_real = zeros(1,nNode);

      
        for i = 1:nNode

            % calculate current channel parameters
            [sigstr, snr, delay, fd] = chFreeSpace(TxNode, RxNode{i});
%             delay = 0;
%             snr = SNR(countT);
%             fd = 0;
%             fprintf("接收节点 %d：信号强度 %.5f dBm, 信噪比 %.5f dB, 时延 %.5f ms， 多普勒 %.5f Hz \n", ...
%                 i, sigstr, snr, 1e3*delay, fd);
            txt1 = sprintf("接收节点 %d：信号强度 %.5f dBm, 信噪比 %.5f dB, 时延 %.5f ms， 多普勒 %.5f Hz \n", ...
                i, sigstr, snr, 1e3*delay, fd);
            txt = [txt1;txt];
            Contain.Dynamic.sigstr{i,countT} = sigstr;
            Contain.Dynamic.SNR{i,countT} = snr;
            Contain.Dynamic.delay{i,countT} = delay;
            Contain.Dynamic.Fd{i,countT} = fd;
            fd_real(i)= fd;
            delay_real(i)=delay;
            % generate signal according to parameters
            sig_data = nodeSigGen(TxNode.Tx.TransmitterFrequency, RxNode{i}.sampRate, sigtmp_tx, ...
                delay, fd, snr, RxNode{i}.Rx.ReceiverSensitivity);

            RecSigTmp(i).sig = rowVec(sig_data.signal);
            sigtmp_rx(i,:) = RecSigTmp(i).sig(1:len_step);
            
            %% ----------- 单节点检测统计量 -----------
            [st] = detcEnergy(sigtmp_rx(i,:), NFFT, 'CalType', 'AVE');
            DetcResult.RxNode(i).st_energy(countT) = st;
%             [st] = detcBeta(sigtmp_rx(i,:), NFFT);
%             DetcResult.RxNode(i).st_beta(countT) = st;

        end
        txt_re.Value = txt ;
%         Contain.Data{1,countT} = sigtmp_rx; 
    end
        %% ----------- 时差频差测试 -------------
        fd_diff_real = fd_real-fd_real(1);
        delay_diff_real = delay_real-delay_real(1);
        [amFd, amDelay] = crossAmbiguityMul(sigtmp_rx(:,1:20000).', SAMP_RATE, FD_MAX, 201);
        Contain.Dynamic.delay_real_diff{countT} = delay_diff_real;
        Contain.Dynamic.fd_real_diff{countT} = fd_diff_real;
        Contain.Dynamic.amDelay{countT} = amDelay;
        Contain.Dynamic.amFd{countT} = amFd ;       



        %% ----------- 对接收信号进行多普勒补偿
        DataDopplerCompensate = DopplerCompensate(sigtmp_rx.',SAMP_RATE, fd_diff_real);  
        allNodeDataSyn= roundTimeSyn(DataDopplerCompensate,SAMP_RATE,delay_diff_real,500);

        %% ----------- 互相关检测统计量 -----------
         [st] = detcTimeCorr(allNodeDataSyn(:,5),allNodeDataSyn(:,6));
         DetcResult.st_xcor(countT) = st;

        %% ----------- 多节点检测统计量 -----------
        [st] = detcGlrMul(allNodeDataSyn);
        DetcResult.simulMul.st_GLR(countT) = st;
        [st] = detcEnergyMul(sigtmp_rx.');
        DetcResult.simulMul.st_energy(countT) = st;
%         [st] = detcAgmMul(sigtmp_rx.') ;
%         DetcResult.simulMul.st_Agm(countT) = st;
%         [st] = detcMmeMul(sigtmp_rx.') ;
%         DetcResult.simulMul.st_Mme(countT) = st;
%         [st] = detcEmrMul(sigtmp_rx.');
%         DetcResult.simulMul.st_Emr(countT) = st;
%         [st] = detcSleMul(sigtmp_rx.');
%         DetcResult.simulMul.st_Sle(countT) = st; 
        

        %% ----------- chan 定位，GLRT 定位 -----------
     % chan定位      
       posiChan = locatChanTDOA('amDelay',amDelay,'all_node_position',recsite.');
       LocaResult.Chan{1,countT}= posiChan;
     % GLRT定位
       delete(grid_point);
       [posiGLR,grid_point] = locatDirectPos(DataDopplerCompensate,recsite.','ngrid',10,'sampRate',SAMP_RATE,'TruePosi',track_enemy(countT,:),'grid',[]);
       LocaResult.GLR{1,countT} = posiGLR;
        %% ----------- 跟踪算法 -----------
%        drone_positions_xyzecef = lla2ecef(drone_positions);
       track_enemy_xyz_ecef = lla2ecef(track_enemy);
       priori_ecef_Z = track_enemy_xyz_ecef(:,3).';
%        tdoa_signal = cell(1,countT);
       LocaResult_GLR_xyz_ecef = cell(1,countT);
%        LocaResult_Chan_xyz_ecef = cell(1,countT);
       for i = 1:countT
%            tdoa_signal{i} = Contain.Dynamic.amDelay{i}(2:end).';
%            LocaResult_Chan_xyz_ecef{i} = lla2ecef(LocaResult.Chan{1,countT}).';
           LocaResult_GLR_xyz_ecef{i} = lla2ecef(LocaResult.GLR{1,countT}).';
       end
       %%% 将GLR网格搜索直接定位结果（从1时刻到当前时刻的定位结果）传入PHD_Rand_Linmeas跟踪算法作为量测（地心地固坐标系，Z作为先验信息）
        GLRTrack = PHD_Rand_Linmeas(LocaResult_GLR_xyz_ecef,priori_ecef_Z(1:countT));
        LocaResult.GLRTrack{1,countT} = [GLRTrack(:,end);1000];
%             LocaResult.GLRTrack1{1,countT} = GLRTrack;

%             ChanTrack = PHD_Rand_Linmeas(LocaResult_Chan_xyz_ecef,priori_ecef_Z(1:countT));
%             LocaResult.ChanTrack{1,countT} = [ChanTrack(:,end);1000];
%             LocaResult.ChanTrack1{1,countT} = ChanTrack;


        


%%   plot
    delete(alley_points);
%     delete(enemy_points);
    alley_points = geoshow(ax,drone_positions(:, 1),  drone_positions(:, 2),  'DisplayType', 'point', 'Marker', '.', 'Color', 'r', 'MarkerSize', 15);
    enemy_points = geoshow(ax,track_enemy(countT, 1),   track_enemy(countT, 2),  'DisplayType', 'point', 'Marker', '.', 'MarkerEdgeColor', 'blue', 'MarkerSize', 15);
   % 节点1  节点5 时频图
    TimeFrequencySpectrum(sigtmp_rx(1,:),RxNode{1},Back);
    TimeFrequencySpectrum(sigtmp_rx(5,:),RxNode{5},Front);

   if countT>=5
%      delete(Chan_points);
        
    if DetcResult.st_xcor(countT-3:countT)>DetcResult.MulTh.Xcor
           txt4 = sprintf('GLR融合检测在工作\n');
       if DetcResult.RxNode(5).st_energy(countT-3:countT)>DetcResult.Th.ED
           txt4 = sprintf('单节点能量检测在工作\n');
       end
    else
           txt4 = sprintf('互相关检测在工作\n');
    end
    txt_re_1.Value = txt4 ;


%        if DetcResult.st_xcor(countT-3:countT)>DetcResult.MulTh.Xcor 
%          
%          
%            if calDistance(LocaResult.Chan{1,countT}(1:2),track_enemy(countT,1:2))>3000
%               Chan_points = geoshow(ax,track_enemy(countT,1)-0.05,  track_enemy(countT,2)-0.05,  'DisplayType', 'point', 'Marker', '.', 'MarkerEdgeColor', 'black', 'MarkerSize', 1);   
%            else 
%               Chan_points = geoshow(ax,LocaResult.Chan{1,countT}(1),  LocaResult.Chan{1,countT}(2),  'DisplayType', 'point', 'Marker', '.', 'MarkerEdgeColor', 'black', 'MarkerSize', 30);  
%            end
%            geoshow(ax,LocaResult.TdoaTrack{1,countT}(1),  LocaResult.TdoaTrack{1,countT}(2),  'DisplayType', 'point', 'Marker', '.', 'MarkerEdgeColor', 'yellow', 'MarkerSize', 30);  
%  
%        end
    if P_TRA >=1e-2         %%%高信噪比下才在地图显示Chan定位结果;低信噪比下不显示定位结果，因为低信噪比时差估计误差大，定位误差太大
       if DetcResult.st_xcor(countT-3:countT)>DetcResult.MulTh.Xcor 

           if calDistance(LocaResult.Chan{1,countT}(1:2),track_enemy(countT,1:2))<= 6000
              geoshow(ax,LocaResult.Chan{1,countT}(1),  LocaResult.Chan{1,countT}(2),  'DisplayType', 'point', 'Marker', '.', 'MarkerEdgeColor', 'black', 'MarkerSize', 30);  
           end   
       end
    end
    

     if DetcResult.st_xcor(countT-3:countT)>DetcResult.MulTh.Xcor
        geoshow(ax,LocaResult.GLR{1,countT}(1),  LocaResult.GLR{1,countT}(2),  'DisplayType', 'point', 'Marker', '.', 'MarkerEdgeColor', 'green', 'MarkerSize', 20);
        geoshow(ax,LocaResult.GLRTrack{1,countT}(1),  LocaResult.GLRTrack{1,countT}(2),  'DisplayType', 'point', 'Marker', '.', 'MarkerEdgeColor', 'yellow', 'MarkerSize', 30);  

     end
     


    txt1 = sprintf('能量检测法/检验统计量: %d, 能量检测法/检测阈值: %d, 能量检测法/判决结果: %d\n',DetcResult.RxNode(5).st_energy(countT),DetcResult.Th.ED,double(DetcResult.RxNode(5).st_energy(countT)>=DetcResult.Th.ED));
    txt2 = sprintf('互相关检测法/检验统计量: %d, 互相关检测法/检测阈值: %d, 互相关检测法/判决结果: %d\n',DetcResult.st_xcor(countT),DetcResult.MulTh.Xcor,double(DetcResult.st_xcor(countT)>=DetcResult.MulTh.Xcor));
%     txt3 = sprintf('MultiED/st: %d, MultiED/th: %d, MultiED/judge: %d\n',DetcResult.simulMul.st_energy(countT),DetcResult.MulTh.Energy,double(DetcResult.simulMul.st_energy(countT)>=DetcResult.MulTh.Energy));
    
    if DetcResult.st_xcor(countT-3:countT)>DetcResult.MulTh.Xcor
       txt4 = sprintf('GLR融合检测工作\n');
       if DetcResult.RxNode(5).st_energy(countT-3:countT)>DetcResult.Th.ED
           txt4 = sprintf('单节点能量检测在工作\n');
       end
       txt5 = sprintf('GLR融合检测/检验统计量: %d, GLR融合检测/检测阈值: %d, GLR融合检测/判决结果: %d\n',DetcResult.simulMul.st_GLR(countT),DetcResult.MulTh.GLR,double(DetcResult.simulMul.st_GLR(countT)>=DetcResult.MulTh.GLR));
    else 
       txt4 = sprintf('互相关检测工作\n');
       txt5 = sprintf('<<<<<<<Waiting for>>>>>>\n');
    end
    txt_temp = {['Information' num2str(countT) ':'];txt1;txt2;txt4;txt5} ;
    txt = [txt_temp;txt];
    txt_re.Value = txt ;
    


    %%%%%单节点能量检测法检验统计量随时间的变化曲线
    h1 = plot(ED,countT-1:countT,DetcResult.Th.ED*ones(1,2)/DetcResult.Th.ED,'b-.');
    legend(ED,'检测门限（1% 虚警）')
    if DetcResult.RxNode(5).st_energy(countT)>DetcResult.Th.ED
       hold(ED,'on')
       h2 = plot(ED,countT-1:countT,DetcResult.RxNode(5).st_energy(countT-1:countT)/DetcResult.Th.ED,'r-.');
       legend(ED,[h1,h2],'检测门限（1% 虚警）','检测统计量')
    else 
       hold(ED,'on')
       h2 = plot(ED,countT-1:countT,DetcResult.RxNode(5).st_energy(countT-1:countT)/DetcResult.Th.ED,'black-.');    
       legend(ED,[h1,h2],'检测门限（1% 虚警）','检测统计量')
    end
    
    %%%%plot  互相关检验统计量随时间的变化曲线
    h1 = plot(xcor,countT-1:countT,DetcResult.MulTh.Xcor*ones(1,2)/DetcResult.MulTh.Xcor,'b-.');
    legend(xcor,'检测门限（1% 虚警）')
    if DetcResult.st_xcor(countT)>DetcResult.MulTh.Xcor
       hold(xcor,'on')  
       h2 = plot(xcor,countT-1:countT,DetcResult.st_xcor(countT-1:countT)/DetcResult.MulTh.Xcor,'r-.');         
       legend(xcor,[h1,h2],'检测门限（1% 虚警）','检测统计量')
    else
       hold(xcor,'on')  
       h2 = plot(xcor,countT-1:countT,DetcResult.st_xcor(countT-1:countT)/DetcResult.MulTh.Xcor,'black-.');         
       legend(xcor,[h1,h2],'检测门限（1% 虚警）','检测统计量')
    end


    %%%%plot glr检验统计量随时间的变化曲线  互相关检测连续4个时刻的检验统计量超过检测门限，GLR检测开启
    h1 = plot(glr,countT-1:countT,DetcResult.MulTh.GLR*ones(1,2)/DetcResult.MulTh.GLR,'b-.');
    legend(glr,'检测门限（1% 虚警）')
    if DetcResult.st_xcor(countT-3:countT)>DetcResult.MulTh.Xcor

        if DetcResult.simulMul.st_GLR(countT)>DetcResult.MulTh.GLR
           hold(glr,'on')
           h2 = plot(glr,countT-1:countT,DetcResult.simulMul.st_GLR(countT-1:countT)/DetcResult.MulTh.GLR,'r-.');
           
           legend(glr,[h1,h2],'检测门限（1% 虚警）','检测统计量')
        else
           hold(glr,'on')
           h2 = plot(glr,countT-1:countT,DetcResult.simulMul.st_GLR(countT-1:countT)/DetcResult.MulTh.GLR,'black-.');
           legend(glr,[h1,h2],'检测门限（1% 虚警）','检测统计量')
        end 
    
    end 

    
    

    %%%plot定位误差曲线  互相关检测连续4个时刻的检验统计量超过检测门限，定位算法开启
    if DetcResult.st_xcor(countT-3:countT)>DetcResult.MulTh.Xcor
       Rms = zeros(3,countT);
       for i = 1:countT
           Rms(1,i) = distWGS(LocaResult.Chan{1,i},track_enemy(i,:));
           Rms(2,i) = distWGS(LocaResult.GLR{1,i},track_enemy(i,:));
           Rms(3,i) = distWGS(LocaResult.GLRTrack{1,i},track_enemy(i,:));
       end
       Rms(find(abs(Rms)>=3000))=3000;
       plot(posi,countT-1:countT,Rms(1,countT-1:countT),'r-.')
       hold(posi,'on')
       plot(posi,countT-1:countT,Rms(2,countT-1:countT),'blue-.')
       hold(posi,'on')
       plot(posi,countT-1:countT,Rms(3,countT-1:countT),'black-.')
       legend(posi,'两步法Chan时差定位','GLR网格搜索定位','GLR直接定位结果用于PHD跟踪算法')
    end

   end
   % 刷新图片并等待
    drawnow
    pause(pauseTime)
end
folderName = datestr(now,'yyyymmdd');
if ~exist(folderName,'dir')
   mkdir(folderName)
end
save([folderName ,'/data_hight_snr.mat'],"DetcResult","LocaResult","Contain");
end

  