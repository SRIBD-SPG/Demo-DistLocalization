function [track_X_lla_Linmeas] = PHD_Rand_Linmeas(pos_meas,z_truth)
%基于TDOA量测的多目标跟踪，用GM-PHD实现
% clc; clear; close all;
%% PHD滤波器
%将TDOA信息用于EKF滤波器
% formatMatlabCode

%% 基本参数（可修改）
monteCarlo = 1; %蒙特卡洛次数
Ttotal = size(pos_meas, 2); T = 0.1; %总时间和采样间隔
x_dim = 6; %运动状态维度
c_dim = 1; %类信息维度
z_dim = size(pos_meas{1, 1},1); %观测维度
lambda_c = 20e-2000; %杂波平均个数
pdf_c = 1e-6;
F1 = [1 T 0 0 0 0; 
      0 1 0 0 0 0; 
      0 0 1 T 0 0; 
      0 0 0 1 0 0;
      0 0 0 0 1 T;
      0 0 0 0 0 1]; %cv模型

F_ex = [F1 zeros(x_dim, c_dim); 
    zeros(c_dim, x_dim) ones(c_dim)]; %扩展转移矩阵（同时转移运动状态和类信息）
Q = diag([30 20 30 20 1e-8 1e-8]'/1) .^ 2; %过程噪声
sigma = 4; %Z轴先验已知
R = diag(sigma .^ 2); %观测噪声
P_S = 0.95; %存活概率
P_D = 0.95; %检测概率
C = [0.82 0.03 0.03 0.03 0.03 0.03 0.03]; %混淆矩阵
l_C = [max(C, [], 2) sum(C, 2) - max(C, [], 2)]; %类似然

% 新生参数
L_bir = 1; %高斯分量个数
L_C = 1; 
w_bir = zeros(L_bir, 1); 
m_bir = zeros(x_dim + c_dim, L_bir); 
B_bir = zeros(x_dim, x_dim, L_bir); 
P_bir = zeros(x_dim, x_dim, L_bir); 
C_bir = zeros(L_C, 1); 

%% 产生真实状态
nbirths = 1; 
t_sur = zeros(2, nbirths); %记录各目标的存活时间起始

%% 滤波器主体
est_X = cell(Ttotal, monteCarlo); 
est_P = cell(Ttotal, monteCarlo); 
est_N = zeros(Ttotal, monteCarlo); 
est_L = cell(Ttotal, monteCarlo); 
Z = cell(Ttotal, monteCarlo); 

for mc = 1 : monteCarlo
    w_bir(1, 1) = 0.1; 
    m_bir(:, 1) = [-2982588.36647978; 0; 4658069.44143427; 0; 3165384.62762881; 0; 1];
    B_bir(:, :, 1) = diag([ 30; 20; 30; 20; 1e-8; 1e-8]/1)*1e0;%Z轴先验已知
    P_bir(:, :, 1) = B_bir(:, :, 1) * B_bir(:, :, 1)'; 
    C_bir(1, 1) = 1;
    %一些初始值
    w_up = []; 
    m_up = []; 
    P_up = []; 
    l_up = [];
    Traj = cell(1, 0); %提取的航迹
    Tst = zeros(1, 0); %航迹起始时刻
    IDt = 0; 
    m_dis = cell(1, IDt); 
    assot = cell(1, IDt); 
    x_alltruth2 = cell(1, Ttotal); %存储每个时刻的所有状态(量测到达时间不定)

    %确定各个目标对应量测的到达时刻
    t_arr = cell(1, nbirths);
    t_max = 5;%最长量测到达时间间隔
    for i = 1 : nbirths
        t_arr{1, i} = t_sur(1, i);
        while t_arr{1, i}(1, end) <= t_sur(2, i)
              t_diff = randi(t_max);
              t_arr{1, i}(1, end + 1) = t_arr{1, i}(1, end) + t_diff;
        end
        t_end = find(t_arr{1, i} <= t_sur(2, i), 1, 'last');
        t_arr{1, i} = t_arr{1, i}(1, 1 : t_end);
    end
    
    for t = 1 : Ttotal
        for i = 1 : nbirths
            if ismember(t, t_arr{1, i})
               x_alltruth2{1, t} = [x_alltruth2{1, t} x_truth{1, i}(:, t - t_sur(1, i) + 1)]; 
            end
        end
    end


    asso = zeros(0, Ttotal); %用于表示每一个接收到量测的时刻，量测信息的关联情况
    num = 0;
    dis_threshold = 500000;
    w_temp = cell(1, 0);
    m_temp = cell(1, 0);
    P_temp = cell(1, 0);
    %递归滤波
    for t = 1 : Ttotal
       fprintf('Monte Carlo: %d, Time: %d\n', mc, t); 
       if ~isempty(pos_meas{1, t})
           Z{t, mc}{1, 1}(1:3, :) = pos_meas{1, t};
           Z{t, mc}{1, 1}(4, :) = 1;
       else
           Z{t, mc} = []; 
       end
           

        %% 预测
        %存活        
        if isempty(m_up)
            w_pre = []; m_pre = []; P_pre = []; 
        else
            w_pre = P_S * w_up; 
            m_pre = F_ex * m_up; 
            P_pre = zeros(size(P_up)); 
            for j = 1 : size(P_up, 3)
                P_pre(:, :, j) = F1 * P_up(:, :, j) * F1' + Q; 
            end
            m_pre(5,:) = z_truth(1,t);
        end

        if t == 1
            w_pre = w_bir; 
            m_pre = m_bir; 
            P_pre = P_bir; 
        end

        %% 更新
        w_up = w_pre; 
        m_up = m_pre; 
        P_up = P_pre; 

        if t == 1
            l_up = zeros(2, size(m_up, 2));
        end
           

        for s = 1 : size(Z{t, mc}, 1)
            w_up1 = zeros(0, 1); 
            m_up1 = zeros(x_dim + c_dim, 0); 
            P_up1 = zeros(x_dim, x_dim, 0); 

            %观测更新
            if ~isempty(Z{t, mc}{s, 1})
                w_up2 = zeros(size(m_up, 2), size(Z{t, mc}{s, 1}, 2)); 
                m_up2 = cell(size(m_up, 2), size(Z{t, mc}{s, 1}, 2)); 
                P_up2 = cell(size(m_up, 2), size(Z{t, mc}{s, 1}, 2)); 
                M_dis = zeros(size(m_up, 2), size(Z{t, mc}{s, 1}, 2));
                like = zeros(size(m_up, 2), size(Z{t, mc}{s, 1}, 2)); %似然值
                sum_like = 0;
                H = zeros(z_dim, x_dim, size(m_up, 2));
                for i = 1 : size(m_up, 2)
                    H(:, :, i) = zeros(z_dim, x_dim); %观测矩阵
                    H(1, 1, i) = 1; 
                    H(2, 3, i) = 1; 
                    H(3, 5, i) = 1;
    
                    Z_pre(:, i) = [m_up(1,1);m_up(3,1);m_up(5,1)]; 
                    S(:, :, i) = H(:, :, i) * P_up(:, :, i) * H(:, :, i)' + R; 
                    K(:, :, i) = P_up(:, :, i) * H(:, :, i)' * pinv(S(:, :, i)); 
    
                    for j = 1 : size(Z{t, mc}{s, 1}, 2)
                        M_dis(i, j) = (Z{t, mc}{s, 1}(1:3, j) - Z_pre(:, i))' * (Z{t, mc}{s, 1}(1:3, j) - Z_pre(:, i)); 
                        like(i, j) = 1 / sqrt(det(2 * pi * S(:, :, i))) * 1 / ( M_dis(i, j) ); 
                        is_corr = m_up(end, i) ~= Z{t, mc}{s, 1}(4, j); 
                        is_corr = is_corr + 1; 
                        w_up2(i, j) = w_up(i) * (P_D / sqrt(det(2 * pi * S(:, :, i))) * ...
                                      exp(-0.5 * (Z{t, mc}{s, 1}(1:3, j) - Z_pre(:, i))' * pinv(S(:, :, i)) * (Z{t, mc}{s, 1}(1:3, j) - Z_pre(:, i))) + 1e-299)...
                                      * l_C(m_up(end, i), is_corr); 
                    end
                end
                COSTmat = log(M_dis / dis_threshold^2);
                [uasses, ~] = mbestwrap_updt_custom(COSTmat, 1); 
    
                for i = 1 : size(m_up, 2)
                    %如果该目标和某个量测之间的距离小于门限，则用该量测更新目标的信息
                     if uasses(i) ~= 0
                         w_up(i) = w_up2(i, uasses(i)) /(lambda_c * pdf_c + sum(w_up2(:, uasses(i))));
                         m_up(1 : end - 1, i) = m_up(1 : end - 1, i) + K(:, :, i) * (Z{t, mc}{s, 1}(1:3, uasses(i)) - Z_pre(:, i)); 
                         P_up(:, :, i) = (eye(size(P_up(:, :, i))) - K(:, :, i) * H(:, :, i)) * P_up(:, :, i);
                 
                         if l_up(1, i) == 0 && l_up(2, i) == 0
                             %第一次和观测关联上
                             num = num + 1;
                             l_up(:, i) = [t; num];
                             asso(num, t) = 1;
                         else
                             %第n次和观测关联上
                             asso(l_up(2, i), t) = 1; 
                         end 
                     end
                end
           end
        end

        for l = 1 : size(l_up, 2)
            if  t > l_up(1, l) && l_up(2, l) ~= 0 
                w_temp{1, l_up(2, l)} = cat(2, w_temp{1, l_up(2, l)}, w_up(l_up(2, l)));
                m_temp{1, l_up(2, l)} = cat(2, m_temp{1, l_up(2, l)}, m_up(:, l_up(2, l)));
                P_temp{1, l_up(2, l)} = cat(3, P_temp{1, l_up(2, l)}, P_up(:, :, l_up(2, l)));
            elseif t == l_up(1, l) && l_up(2, l) ~= 0 
                w_temp{1, l_up(2, l)} = w_up(l_up(2, l));
                m_temp{1, l_up(2, l)} = m_up(:, l_up(2, l));
                P_temp{1, l_up(2, l)} = P_up(:, :, l_up(2, l));
            end
        end

        for l = 1 : size(asso, 1)
            t_asso = find(asso(l, 1 : t) == 1, 1, 'last');
            len = length(w_temp{1, l_up(2, l)});
            if ~isempty(t_asso) && t == t_asso && len <= 10
               idx = find(w_temp{1, l_up(2, l)} > 0.5); 
                if w_temp{1, l_up(2, l)}(end) > 0.5
                    for j = 1 : length(w_temp{1, l_up(2, l)})
                        repeat_num_targets = round(w_temp{1, l_up(2, l)}(end)); 
                        est_X{t - len + j, mc} = cat(2, est_X{t - len + j, mc}, repmat(m_temp{1, l_up(2, l)}(:, j), [1, repeat_num_targets]) ); 
                        est_P{t - len + j, mc} = cat(3, est_P{t - len + j, mc}, repmat(P_temp{1, l_up(2, l)}(:, :, j), [1, 1, repeat_num_targets]) ); 
                        est_N(t - len + j, mc) = est_N(t - len + j, mc) + repeat_num_targets; 
                        est_L{t - len + j, mc} = cat(2, est_L{t - len + j, mc}, l_up(:, l)); 
                    end
                end
                w_temp{1, l_up(2, l)} = [];
                m_temp{1, l_up(2, l)} = [];
                P_temp{1, l_up(2, l)} = [];
            end
        end
    end
end

%% 输出
%输出跟踪轨迹
est_X_sum = cell(Ttotal,1);
for i = 1 : Ttotal
    est_X_sum{i,1} = zeros(7,1);
end
for i = 1 : Ttotal
    for j = 1 : mc
        est_X_sum{i,1} = est_X_sum{i,1} + est_X{i,mc};
    end
    est_X_sum{i,1} = est_X_sum{i,1}/mc;
end

for i = 1 : Ttotal
    track_X(1,i) = est_X_sum{i,1}(1,:);
    track_X(2,i) = est_X_sum{i,1}(3,:);
    track_X(3,i) = est_X_sum{i,1}(5,:);
end
track_lla = ecef2lla(track_X')';
track_X_lla_Linmeas = track_lla(1:2,:);
end