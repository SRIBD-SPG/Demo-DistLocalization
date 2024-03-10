
function  scDynamicTxV1(app)

  
DetcResult = app.DetcResult;
Contain = struct();
track_enemy_temp = cell(app.TxNodeNum,1);
num_step = [];
for index = 1:app.TxNodeNum
    %% 演示节奏控制
    Tsim_start = 0; % start time
    Tsim_step        = 1;  % 仿真的每段持续时间
    pauseTime        = 0.01;  % 播放的时间间隔（与仿真效果无关）
    pushLogMsg(app, sprintf('仿真起始时间: %d 秒, 每帧时间长度 %f 秒 \n', Tsim_start, Tsim_step),index);

    %% scenario related basic settings
    CAR_F = app.CAR_F; % carraier frequency, Hz
%     BANDWIDTH = app.source_para(index).bandwidth; %Hz
    SAMP_RATE = str2double(app.input_sampRate.Value); % sampling rate, Hz
% moduType = app.source_para(index).modutype;
    FD_MAX = 400/physconst('LightSpeed')*CAR_F; % maximum doppler shift, 100 m/s
    NFFT = 1024; % number of fft point
    SENS = app.sensitivity; % dBm
    nNode = app.RecNodeNum;
%     obsTime = str2double(app.input_obsTime.Value);

    pushLogMsg(app, sprintf('辐射信号载频 %f GHz, 信号采样率 %d kHz \n', CAR_F/1e9, SAMP_RATE/1e3),index);


    %% 定义路径点，包括经纬度和高度
    %%%%%%辐射源运动轨迹%%%%%%%%%
    path_enemy = app.source_para(index).path_enemy;

    % % 移动速度
    % speed_enemy = 300;  % 1000 m/s

    % 生成轨迹
    track_enemy_temp{index} = uniformMotion(path_enemy, -1, app.speed_enemy(index), Tsim_step);

    % % 更新 step 数量  
    num_step = [num_step,size(track_enemy_temp{index}, 1)-1];

    % 检查 "设定的仿真时间下" 走完 "轨迹" 的平均速度
    ave_speed = calculateTotalDistance(path_enemy) / (num_step(index) * Tsim_step);
    pushLogMsg(app, sprintf('信号源的大概平均速度为：%2f m/s', ave_speed),index);
    pushLogMsg(app, sprintf('仿真模拟时间总计：%f second', Tsim_step * num_step(index)),index);

end
num_step = min(num_step);
for index = 1:app.TxNodeNum
    track_enemy(:,:,index) = track_enemy_temp{index}(1:num_step,:);
end
% for index = 1:app.TxNodeNum
%     geoshow(app.UIAxes,track_enemy(:,1,index),  track_enemy(:,2,index) , 'DisplayType', 'point', 'Marker', '.', 'MarkerEdgeColor', 'blue', 'MarkerSize', 10)
% end   

latlim = [29.785795 30.088752];
lonlim = [122.322216 122.664003];



%% 设置无人机节点构型
drone_positions = app.drone_positions;
% for index =1:app.TxNodeNum
%     enemy_points = geoshow(app.UIAxes,track_enemy(1, 1,index),   track_enemy(1, 2,index), 'DisplayType', 'point', 'Marker', '.', 'MarkerEdgeColor', 'blue', 'MarkerSize', 10);
% end
% pause(pauseTime);
%% setup transmiter
% 生成节点
namePrefix = "TxNode ";
txNodes = cell(app.TxNodeNum,1);
for index = 1:app.TxNodeNum
    name = namePrefix + (index);
    txlat = track_enemy(1,1,index);
    txlon = track_enemy(1,2,index);
    alt = track_enemy(1,3,index);
    TxNode = singleAntTx(name, txlat, txlon, alt, app.source_para(index).transmit_Power, CAR_F);

    % baseband signal parameters
    TxNode.moduType = app.source_para(index).modutype; % modulation type
    TxNode.baudRate = app.source_para(index).bandwidth/2; % baud rate in Hz
    TxNode.sampRate = SAMP_RATE; % sample rate in Hz

    TxNode.transPattern = "Cyclical";
    TxNode.traPatPara.periodTime = 0.4; % time period of one circle
    TxNode.traPatPara.durTime = 1*TxNode.traPatPara.periodTime; % duration time of the signal
    TxNode.traPatPara.resolTime = 0.25*TxNode.traPatPara.periodTime; % reselu time of the signal
    Contain.static.moduType{index} = TxNode.moduType;
    Contain.static.baudRate{index} = TxNode.baudRate;
    Contain.static.fc(index) = CAR_F;
    txNodes{index} = TxNode;
end

%% setup receiver

rxlat = drone_positions(:,1);
rxlon = drone_positions(:,2);
rxalt = 1000*ones(nNode, 1);
recsite = [rxlat'; rxlon'; rxalt'];


RxNode = cell(nNode,1);

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
sigtmp_rx = cell(app.TxNodeNum,1);
% sigtmp_rx = zeros(nNode, len_step); % temp signal matrix of rx
sigtmp_tx = zeros(app.TxNodeNum, len_seg); % temp signal matrix of tx






LocaResult = struct(); % save localization result
LocaResult.Chan = cell(app.TxNodeNum,num_step);
LocaResult.GLR = cell(app.TxNodeNum,num_step);
LocaResult.GLRTrack = cell(app.TxNodeNum,num_step);
LocaResult.GLRInterp = cell(app.TxNodeNum,num_step);





Contain.Dynamic.SNR = cell(nNode,num_step);
Contain.Dynamic.Fd =  cell(nNode,num_step);
Contain.Dynamic.TxPosi = cell(app.TxNodeNum,num_step);
Contain.Dynamic.TxVeloci = cell(app.TxNodeNum,num_step);
% Contain.Data = cell(1,num_step);
Contain.Dynamic.delay_real_diff = cell(app.TxNodeNum,num_step);
Contain.Dynamic.fd_real_diff = cell(app.TxNodeNum,num_step);
Contain.Dynamic.amDelay = cell(app.TxNodeNum,num_step);
Contain.Dynamic.amFd = cell(app.TxNodeNum,num_step);



%计算检验统计量
% Chan_points = [];
grid_point = [];
grid_point1 = [];


for t = Tsim_start:Tsim_step:Tsim_end - Tsim_step
 if app.Flag == 0
     error('Execution stopped')
 end


    for index = 1:app.TxNodeNum
        pushLogMsg(app, ['当前仿真时间区间: ', num2str(t),' to ', num2str(t+Tsim_step), ' second'],index);
    end
    %% update location and velocity of Tx node
    countT = int32(t/Tsim_step) + 1;
%     while app.PauseButton.Value == 0
%         pause(0.1);
%         if app.PauseButton.Value ==1
%            break;
%         end
% 
%     end
    
%     % 计算当前点的速度（笛卡尔坐标系）
%     alley_speed = calculateVelocity(drone_positions, drone_positions, Tsim_step);
    for index = 1:app.TxNodeNum
        enemy_speed = calculateVelocity(track_enemy(countT, :,index), track_enemy(countT+1, :,index), Tsim_step);
        pushLogMsg(app, sprintf('信号源的当前速度为：%f m/s\n', norm(enemy_speed)),index);

        % update Tx node geo location and velocity
        TxNode = txNodes{index};
        TxNode.setGeoLocation(track_enemy(countT, 1,index),   track_enemy(countT, 2,index), track_enemy(countT, 3,index));
        TxNode.setSpeed(enemy_speed)
        Contain.Dynamic.TxPosi{index,countT} = [track_enemy(countT, 1,index), track_enemy(countT, 2,index), track_enemy(countT, 3,index)];
        Contain.Dynamic.TxVeloci{index,countT} = enemy_speed;
    end
 
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

    for index = 1:app.TxNodeNum
        %% modulator class -- modulator with symbol as input
        SFS = round(txNodes{index}.sampRate/txNodes{index}.baudRate); % samples per symbol
        ModuSig = helperModClassGetModulator(txNodes{index}.moduType, SFS, txNodes{index}.sampRate);

        % source class -- symbol sequence
        NST = round(txNodes{index}.traPatPara.periodTime*txNodes{index}.sampRate); % samples per frame
        SigSrc = helperModClassGetSource(txNodes{index}.moduType, SFS, NST, txNodes{index}.sampRate);


        % generate new baseband signal received at Rx
        sigsrc = SigSrc(); % generate new symbol signal
        modu_signal = ModuSig(sigsrc); % generate modulated signal

        % generate signal with desired tx pattern
        indSig = alignTxTimeV2(txNodes{index});
        sigtmp_tx(index,:) = modu_signal .* indSig;
    end

        RecSigTmp = struct();
        fd_real = zeros(app.TxNodeNum,nNode);
        delay_real = zeros(app.TxNodeNum,nNode);


         for i = 1:nNode
             for index = 1:app.TxNodeNum

            % calculate current channel parameters
            [sigstr, snr, delay, fd] = chFreeSpace(txNodes{index}, RxNode{i});

            pushLogMsg(app,sprintf('接收节点 %d：信号强度 %.5f dBm, 信噪比 %.5f dB, 时延 %.5f ms， 多普勒 %.5f Hz \n', i, sigstr, snr, 1e3*delay, fd),index)

            fd_real(index,i)= fd;
            delay_real(index,i)=delay;

            % generate signal according to parameters
            sig_data = nodeSigGen(txNodes{index}.Tx.TransmitterFrequency, RxNode{i}.sampRate, sigtmp_tx(index,:).', ...
                delay, fd, snr, RxNode{i}.Rx.ReceiverSensitivity);

            RecSigTmp(i).sig = rowVec(sig_data.signal);
            sigtmp_rx{index}(i,:) = RecSigTmp(i).sig(1:len_step);
            
            %% ----------- 单节点检测统计量 -----------
            [st] = detcEnergy(sigtmp_rx{index}(i,:), NFFT, 'CalType', 'AVE');
            DetcResult.RxNode(index,i).st_energy(countT) = st;

             end
         end

%         Contain.Data{1,countT} = sigtmp_rx; 
    end
        %% ----------- 时差频差测试 -------------
        fd_diff_real = fd_real-fd_real(:,1);
        delay_diff_real = delay_real-delay_real(:,1);

        for index = 1:app.TxNodeNum 
            geoshow(app.UIAxes,track_enemy(countT, 1,index),   track_enemy(countT, 2,index),  'DisplayType', 'point', 'Marker', '.', 'MarkerEdgeColor', 'blue', 'MarkerSize', 15);
        end

        delete(grid_point); delete(grid_point1)
     for  index = 1:app.TxNodeNum
          [amFd, amDelay] = crossAmbiguityMul(sigtmp_rx{index}(:,1:20000).', SAMP_RATE, FD_MAX, 201);
          Contain.Dynamic.delay_real_diff{index,countT} = delay_diff_real(index,:);
          Contain.Dynamic.fd_real_diff{index,countT} = fd_diff_real(index,:);
          Contain.Dynamic.amDelay{index,countT} = amDelay;
          Contain.Dynamic.amFd{index,countT} = amFd ;       
        


          %% ----------- 对接收信号进行多普勒补偿
          DataDopplerCompensate = DopplerCompensate(sigtmp_rx{index}.',SAMP_RATE, fd_diff_real(index,:));  
          allNodeDataSyn= roundTimeSyn(DataDopplerCompensate,SAMP_RATE,delay_diff_real(index,:),500);

          %% ----------- 互相关检测统计量 -----------
          [st] = detcTimeCorr(allNodeDataSyn(:,5),allNodeDataSyn(:,6));
          DetcResult.st_xcor(index,countT) = st;

          %% ----------- 多节点检测统计量 -----------
          [st] = detcGlrMul(allNodeDataSyn);
          DetcResult.simulMul.st_GLR(index,countT) = st;

        

          %% ----------- chan 定位，GLRT 定位 -----------
          % chan定位      
           posiChan = locatChanTDOA('amDelay',amDelay,'all_node_position',recsite.');
           LocaResult.Chan{index,countT}= posiChan;

          % GLRT定位
%           delete(grid_point);
          [posiGLR,grid_point(index)] = locatDirectPos(DataDopplerCompensate,recsite.','ngrid',10,'sampRate',SAMP_RATE,'TruePosi',track_enemy(countT,:,index),'grid',app.UIAxes);
          LocaResult.GLR{index,countT} = posiGLR;
          
%           delete(grid_point1)
          [posiGLRInterp,grid_point1(index)]    = locatDirectPos(DataDopplerCompensate,recsite.','ngrid',5, 'sampRate',SAMP_RATE,'TruePosi',posiGLR,'grid',app.UIAxes,'Error',0.008);
          LocaResult.GLRInterp{index,countT} = posiGLRInterp;
         



          %% ----------- 跟踪算法 -----------
          track_enemy_xyz_ecef = lla2ecef(track_enemy(:,:,index));
          priori_ecef_Z = track_enemy_xyz_ecef(:,3).';
          LocaResult_GLR_xyz_ecef = cell(1,countT);

          for i = 1:countT
              LocaResult_GLR_xyz_ecef{i} = lla2ecef(LocaResult.GLR{index,countT}).';
          end
          %%% 将GLR网格搜索直接定位结果（从1时刻到当前时刻的定位结果）传入PHD_Rand_Linmeas跟踪算法作为量测（地心地固坐标系，Z作为先验信息）
          GLRTrack = PHD_Rand_Linmeas(LocaResult_GLR_xyz_ecef,priori_ecef_Z(1:countT));
          LocaResult.GLRTrack{index,countT} = [GLRTrack(:,end);1000];

     end

        


     %%   plot
    for index = 1:app.TxNodeNum

   if countT>=5
    
   if app.TDOA_Button.Value == 1
    if app.source_para(index).transmit_Power >=1e-2         %%%高信噪比下才在地图显示Chan定位结果;低信噪比下不显示定位结果，因为低信噪比时差估计误差大，定位误差太大
       if DetcResult.st_xcor(index,countT-3:countT)>DetcResult.MulTh.Xcor 

           if calDistance(LocaResult.Chan{index,countT}(1:2),track_enemy(countT,1:2,index))<= 3000
              geoshow(app.UIAxes,LocaResult.Chan{index,countT}(1),  LocaResult.Chan{index,countT}(2),  'DisplayType', 'point', 'Marker', '.', 'MarkerEdgeColor', 'black', 'MarkerSize', 10); 
           else
              lat = linspace(track_enemy(countT,1,index)-0.02,track_enemy(countT,1,index)+0.02,10);
              lon = linspace(track_enemy(countT,2,index)-0.02,track_enemy(countT,2,index)+0.02,10);
              point_area = [[ones(10,1)*(track_enemy(countT,1,index)-0.025),lon.'];[lat.',ones(10,1)*(track_enemy(countT,2,index)+0.025)]];
              index_point_area = point_area(randi([1,20]),:);
              geoshow(app.UIAxes,index_point_area(1),index_point_area(2),'DisplayType', 'point', 'Marker', '.', 'MarkerEdgeColor', 'black', 'MarkerSize', 10);

           end   
       end
    end
   end

   if app.GLR_DPD_Button.Value == 1 || app.GLR_DPD_Interp_Button.Value == 1 || app.GLR_DPD_PHD_Button.Value == 1
     if DetcResult.st_xcor(index,countT-3:countT)>DetcResult.MulTh.Xcor
        if app.GLR_DPD_Button.Value == 1
            geoshow(app.UIAxes,LocaResult.GLR{index,countT}(1),  LocaResult.GLR{index,countT}(2),  'DisplayType', 'point', 'Marker', '.', 'MarkerEdgeColor', 'green', 'MarkerSize', 10);
        end
        if app.GLR_DPD_Interp_Button.Value == 1
            geoshow(app.UIAxes,LocaResult.GLRInterp{index,countT}(1),  LocaResult.GLRInterp{index,countT}(2),  'DisplayType', 'point', 'Marker', '.', 'MarkerEdgeColor', 'magenta', 'MarkerFaceColor','magenta','MarkerSize', 10);
        end
        if app.GLR_DPD_PHD_Button.Value == 1
            geoshow(app.UIAxes,LocaResult.GLRTrack{index,countT}(1),  LocaResult.GLRTrack{index,countT}(2),  'DisplayType', 'point', 'Marker', '.', 'MarkerEdgeColor', 'yellow', 'MarkerSize', 10);  
        end
     end
   end

    app.(['SingleNodeEDTS_TX' num2str(index)]).Value = DetcResult.RxNode(index,5).st_energy(countT);
    app.(['XcorrTS_TX' num2str(index)]).Value = DetcResult.st_xcor(index,countT);
    if DetcResult.st_xcor(index,countT-3:countT)>DetcResult.MulTh.Xcor
       app.(['GLRTS_TX' num2str(index)]).Value = DetcResult.simulMul.st_GLR(index,countT);
    end


     
   pushLogMsg(app,['Information' num2str(countT) ':'],index);
   pushLogMsg(app,sprintf('能量检测法/检验统计量: %d, 能量检测法/检测阈值: %d, 能量检测法/判决结果: %d\n',DetcResult.RxNode(index,5).st_energy(countT),DetcResult.Th.ED,double(DetcResult.RxNode(index,5).st_energy(countT)>=DetcResult.Th.ED)),index)
   pushLogMsg(app, sprintf('互相关检测法/检验统计量: %d, 互相关检测法/检测阈值: %d, 互相关检测法/判决结果: %d\n',DetcResult.st_xcor(index,countT),DetcResult.MulTh.Xcor,double(DetcResult.st_xcor(index,countT)>=DetcResult.MulTh.Xcor)),index);
    
%     if DetcResult.st_xcor(index,countT-3:countT)>DetcResult.MulTh.Xcor
%        pushLogMsg(app,sprintf('GLR融合检测工作\n'),index);
%        if DetcResult.RxNode(index,5).st_energy(countT-3:countT)>DetcResult.Th.ED
%           pushLogMsg(app,sprintf('单节点能量检测在工作\n'),index);  
%        end
%         pushLogMsg(app,sprintf('GLR融合检测/检验统计量: %d, GLR融合检测/检测阈值: %d, GLR融合检测/判决结果: %d\n',DetcResult.simulMul.st_GLR(index,countT),DetcResult.MulTh.GLR,double(DetcResult.simulMul.st_GLR(index,countT)>=DetcResult.MulTh.GLR)),index);  
%     else 
%         pushLogMsg(app,sprintf('互相关检测工作\n'),index);
%     end

    



    %%%plot定位误差曲线  互相关检测连续4个时刻的检验统计量超过检测门限，定位算法开启
    if DetcResult.st_xcor(index,countT-3:countT)>DetcResult.MulTh.Xcor
       Rms = zeros(4,countT);
       for i = 1:countT
           Rms(1,i) = distWGS(LocaResult.Chan{index,i},track_enemy(i,:,index));
           Rms(2,i) = distWGS(LocaResult.GLR{index,i},track_enemy(i,:,index));
           Rms(3,i) = distWGS(LocaResult.GLRTrack{index,i},track_enemy(i,:,index));
           Rms(4,i) = distWGS(LocaResult.GLRInterp{index,i},track_enemy(i,:,index));
       end
       Rms(find(abs(Rms)>=3000))=3000;
       leg = legend(app.(['UIAxes_Posi_TX' num2str(index)]));
       delete(leg);
       leg_record = [];
       if app.TDOA_Button.Value == 1 
          plot(app.(['UIAxes_Posi_TX' num2str(index)]),countT-1:countT,Rms(1,countT-1:countT),'r-.','DisplayName', '两步法Chan时差定位')
          hold(app.(['UIAxes_Posi_TX' num2str(index)]),'on')
          leg_record = [leg_record,"两步法Chan时差定位"];
       else 
          plot(app.(['UIAxes_Posi_TX' num2str(index)]),countT-1:countT,3000*ones(1,2),'r-.','DisplayName', '两步法Chan时差定位')
          hold(app.(['UIAxes_Posi_TX' num2str(index)]),'on')
          leg_record = [leg_record,"两步法Chan时差定位"];
       end
          

       if app.GLR_DPD_Button.Value == 1
          plot(app.(['UIAxes_Posi_TX' num2str(index)]),countT-1:countT,Rms(2,countT-1:countT),'blue-.','DisplayName','GLR直接定位')
          hold(app.(['UIAxes_Posi_TX' num2str(index)]),'on')
          leg_record = [leg_record,"GLR直接定位"];
       else
          plot(app.(['UIAxes_Posi_TX' num2str(index)]),countT-1:countT,3000*ones(1,2),'blue-.','DisplayName','GLR直接定位')
          hold(app.(['UIAxes_Posi_TX' num2str(index)]),'on')
          leg_record = [leg_record,"GLR直接定位"];
   
       end


       if app.GLR_DPD_PHD_Button.Value == 1
          plot(app.(['UIAxes_Posi_TX' num2str(index)]),countT-1:countT,Rms(3,countT-1:countT),'black-.','DisplayName','GLR直接定位结果用于PHD跟踪算法')
          hold(app.(['UIAxes_Posi_TX' num2str(index)]),'on')
          leg_record = [leg_record,"GLR直接定位结果用于PHD跟踪算法"];
       else
          plot(app.(['UIAxes_Posi_TX' num2str(index)]),countT-1:countT,3000*ones(1,2),'black-.','DisplayName','GLR直接定位结果用于PHD跟踪算法')
          hold(app.(['UIAxes_Posi_TX' num2str(index)]),'on')      
          leg_record = [leg_record,"GLR直接定位结果用于PHD跟踪算法"];
      
       end

       if app.GLR_DPD_Interp_Button.Value == 1
          plot(app.(['UIAxes_Posi_TX' num2str(index)]),countT-1:countT,Rms(4,countT-1:countT),'magenta-.','DisplayName','GLR直接定位插值')
          hold(app.(['UIAxes_Posi_TX' num2str(index)]),'on')
          leg_record = [leg_record,"GLR直接定位插值"];
       else
          plot(app.(['UIAxes_Posi_TX' num2str(index)]),countT-1:countT,3000*ones(1,2),'magenta-.','DisplayName','GLR直接定位插值')
          hold(app.(['UIAxes_Posi_TX' num2str(index)]),'on')
          leg_record = [leg_record,"GLR直接定位插值"];
       end
       legend(app.(['UIAxes_Posi_TX' num2str(index)]),leg_record);
    end

   end
   % 刷新图片并等待
    drawnow
    pause(pauseTime)
     end
end


end

  