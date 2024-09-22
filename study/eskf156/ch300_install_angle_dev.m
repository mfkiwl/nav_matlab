close all;
clear;
clc;

format long g;
format compact;

N = 17;             % ESKF维度
R2D = 180/pi;       % Rad to Deg
D2R = pi/180;       % Deg to Rad
GRAVITY = 9.8;      % 重力加速度

%%
% Cb2n: b系到n系的旋转矩阵 Vn = Cb2n * Vb

% 切换到当前工作目录
scriptPath = mfilename('fullpath');
scriptFolder = fileparts(scriptPath);
cd(scriptFolder);

%% 数据载入
load('dataset\240922\240922B3.mat');

%单位国际化
data.imu.acc =  data.imu.acc*GRAVITY;
data.imu.gyr =  data.imu.gyr*D2R;
data.gnss.lat = data.gnss.lat*D2R;
data.gnss.lon = data.gnss.lon*D2R;
data.dev.lat = data.dev.ins_lat*D2R;
data.dev.lon = data.dev.ins_lon*D2R;

%加入仿真噪声
% data.imu.acc(:,2) = data.imu.acc(:,2) + 0.1*GRAVITY;
% data.imu.gyr(:,3) = data.imu.gyr(:,3) + 0.5*D2R;

att = [0 0 0]*D2R; %初始安装角
Cb2v_simulate = att2Cnb(att); % matlab仿真设置的安装角误差
Cb2v = eye(3);

%定义变量
ESKF156_FB_A = bitshift(1,0); %反馈失准角
ESKF156_FB_V = bitshift(1,1); %反馈速度
ESKF156_FB_P = bitshift(1,2); %反馈位置
ESKF156_FB_W = bitshift(1,3); %反馈陀螺零篇
ESKF156_FB_G = bitshift(1,4); %反馈加计零篇
ESKF156_FB_CBV = bitshift(1,5); %反馈安装角误差
gnss_vel_R = 0;
gnss_pos_R = 0;
gnss_lost_elapsed = 0;
gnss_last_valid_time = 0;
chi_lambda_vel = 999;
chi_lambda_pos = 999;
zupt_detect_time = 0;
is_zupt = 0;
%% 说明
% KF 状态量: 失准角(3) 速度误差(3) 位置误差(3) 陀螺零偏(3) 加计零偏(3)

%% 相关选项及参数设置
opt.alignment_time = 1;          % 初始对准时间(s)
opt.gnss_outage = 0;             % 模拟GNSS丢失
opt.outage_start = 600;         % 丢失开始时间(s)
opt.outage_stop = 700;          % 丢失结束时间(s)
opt.nhc_enable = 1;              % 车辆运动学约束
opt.nhc_lever_arm = 0*[0.35,0.35,-1.35]; %nhc杆臂长度 b系下（右-前-上）240816测试杆臂,仅测试天向，其他两轴为目测值
opt.nhc_R = 3;                % 车载非完整性约束噪声
opt.gnss_min_interval = 0;    % 设定时间间隔，例如0.5秒
opt.gnss_delay = 0;              % GNSS量测延迟 sec
opt.gnss_lever_arm = 0*[0.45;0.45;-1.34]; %GNSS杆臂长度 b系下（右-前-上）240816测试杆臂
opt.has_install_esti = 1;       %% can close or open ;1:esti insatllangle ; 0:no esti

opt.zupt_enable = 1;              % 启用ZUPT
opt.zupt_vel_threshold = 2.5;     % ZUPT速度阈值(m/s)
opt.zupt_gyr_threshold = 0.5*D2R; % ZUPT角速度阈值(rad/s)
opt.zupt_time_threshold = 1.0;    % ZUPT持续时间阈值(s)
opt.zupt_R = diag([0.1, 0.1, 0.1].^2); % ZUPT测量噪声协方差

% 初始状态方差:    姿态       ENU速度  水平位置      陀螺零偏                加速度计零偏        安装俯仰角 安装航向角
opt.P0 = diag([[2 2 10]*D2R, [1 1 1], [5 5 5], 0.1*D2R*ones(1,3), 1e-2*GRAVITY*ones(1,3), 2*D2R*ones(1,2) ])^2;
N = length(opt.P0);
% 系统误差:         角度随机游走          速度随机游走                     角速度随机游走            加速度随机游走
opt.Q = diag([(0.1*D2R)*ones(1,3), (0.01)*ones(1,3), 0*ones(1,3), 1/3600*D2R*ones(1,3), 0*GRAVITY*ones(1,3), [0 0] ])^2;

imu_len = length(data.imu.tow);
dev_len = length(data.dev.tow);

imu_dt = mean(diff(data.imu.tow));
gnss_dt = mean(diff(data.gnss.tow));

% 开启静止时间
indices = find(data.imu.tow <= opt.alignment_time);  % 找到在x_seconds时间范围内的所有索引
acc_align0 = data.imu.acc(indices, :);
gyr_align0 = data.imu.gyr(indices, :);
gyr_bias0 = mean(gyr_align0);

% 结束静止时间
start_tow = data.imu.tow(end) - opt.alignment_time;  % 将x秒转换为毫秒，并从结束时间戳中减去
indices = find(data.imu.tow >= start_tow & data.imu.tow <= data.imu.tow(end));
gyr_bias_end = mean(data.imu.gyr(indices, :));

fprintf("gyro起始时刻bias估计:%7.3f,%7.3f,%7.3f deg/s\n", gyr_bias0(1)*R2D, gyr_bias0(2)*R2D, gyr_bias0(3)*R2D);
fprintf("gyro结束时刻bias估计:%7.3f,%7.3f,%7.3f deg/s\n", gyr_bias_end(1)*R2D, gyr_bias_end(2)*R2D, gyr_bias_end(3)*R2D);
fprintf("IMU帧平均间隔:%.3fs\n", imu_dt);
fprintf("GNSS帧平均间隔:%.3fs\n", gnss_dt);

%% 经纬度转换为当地东北天坐标系
lat0 = data.gnss.lat(1);
lon0 = data.gnss.lon(1);
h0 = data.gnss.msl(1);

time_sum = 0;
distance_sum = 0;
gnss_enu = zeros(length(data.gnss.tow), 3);
log.vel_norm = zeros(length(data.gnss.tow), 1);

% 初始化上一次融合的时间
last_gnss_fusion_time = -inf;

inital_gnss_idx = 1;
% 根据速度 获得初始航向角
for i=1:length(data.gnss.tow)
    if norm(data.gnss.vel_enu(i,:)) > 3
        opt.inital_yaw = atan2(data.gnss.vel_enu(i,1),data.gnss.vel_enu(i,2));
        if(opt.inital_yaw < 0)
            opt.inital_yaw =  opt.inital_yaw + 360*D2R;
        end

        inital_gnss_idx = i;
        diff_values = abs(data.imu.tow - data.gnss.tow(inital_gnss_idx));
        [~, inital_imu_idx] = min(diff_values);

        fprintf("初始航向角:%.2f°, 从GNSS数据:%d开始, IMU数据:%d\r\n",  opt.inital_yaw*R2D, inital_gnss_idx, inital_imu_idx);
        break;
    end
end
if i == length(data.gnss.tow)
    opt.inital_yaw = 0;
    fprintf("无法通过速度矢量找到初始航向角，设置为:%.2f°\r\n",  opt.inital_yaw*R2D);
    inital_imu_idx = 1;
end


for i=inital_gnss_idx : length(data.gnss.tow)
    [gnss_enu(i,1), gnss_enu(i,2), gnss_enu(i,3)] =  ch_LLA2ENU(data.gnss.lat(i), data.gnss.lon(i), data.gnss.msl(i), lat0, lon0, h0);
    log.vel_norm(i) = norm(data.gnss.vel_enu(i, :));
    distance_sum = distance_sum + norm(data.gnss.vel_enu(i, :))*gnss_dt;
    time_sum = time_sum + gnss_dt;
end

%% MCU结果转换为当地东北天坐标系
dev_pos_enu = zeros(dev_len, 3);
for i=1:dev_len
    [dev_pos_enu(i,1), dev_pos_enu(i,2), dev_pos_enu(i,3)] =  ch_LLA2ENU(data.dev.lat(i), data.dev.lon(i),  data.dev.ins_msl(i), lat0, lon0, h0);
end

%% 初始参数设置
% 粗对准
g_b = - mean(acc_align0)';
g_b = g_b/norm(g_b);
pitch0 = asin(-g_b(2));
roll0 = atan2( g_b(1), -g_b(3));
yaw0 = opt.inital_yaw;
pitch_sins = pitch0;
roll_sins = roll0;
yaw_sins = yaw0;
Qb2n = angle2quat(-yaw0, pitch0, roll0, 'ZXY');
Qb2n_sins = angle2quat(-yaw0, pitch0, roll0, 'ZXY');
vel = [0 0 0]';
pos = [0 0 0]';

X = zeros(N,1);
X_temp = X;
gyro_bias = X(10:12);
acc_bias = X(13:15);

P = opt.P0;

log.pitch = zeros(imu_len, 1);
log.roll = zeros(imu_len, 1);
log.yaw = zeros(imu_len, 1);
log.vel = zeros(imu_len, 3);
log.pos = zeros(imu_len, 3);
log.P = zeros(imu_len, N);
log.X = zeros(imu_len, N);
log.gyro_bias = zeros(imu_len, 3);
log.acc_bias = zeros(imu_len, 3);
log.sins_att = zeros(imu_len, 3);
log.vb = zeros(imu_len, 3);
log.installangle = zeros(imu_len, 3);
tic;
last_time = toc;
gnss_idx = inital_gnss_idx;
imucnt= 0;j=0;
imu_after_cal = nan(imu_len, 7);

for i=inital_imu_idx:imu_len
    imucnt=imucnt+1;
    FB_BIT = 0; %反馈标志
    curr_time = toc;
    if curr_time - last_time >= 1 % 如果自上次更新已经过去了至少1秒
        fprintf('已完成 %.2f%%\n', (i / imu_len) * 100);
        last_time = curr_time; % 更新上次的时间
    end

    %% 捷联更新
    % 单子样等效旋转矢量法
    w_b = data.imu.gyr(i,:)';
    w_b = Cb2v_simulate*w_b;
    w_b = w_b - gyro_bias;

    f_b = data.imu.acc(i,:)';
    f_b = Cb2v_simulate*f_b;
    f_b = f_b - acc_bias;

    j = j+1;
    imu_after_cal(j,:)=[data.imu.tow(i),w_b',f_b'];

    % 捷联更新
    [Qb2n, pos, vel, ~] = ins(w_b, f_b, Qb2n, pos, vel, GRAVITY, imu_dt);

    % 纯捷联姿态
    [Qb2n_sins, ~, ~, ~] = ins(w_b, f_b, Qb2n_sins, pos, vel, GRAVITY, imu_dt);

    f_n = ch_qmulv(Qb2n, f_b);
    a_n = f_n + [0; 0; -GRAVITY];
    Cb2n = ch_q2m(Qb2n); %更新Cb2n阵
    Cn2b = Cb2n'; %更新Cn2b阵

    log.vb(i, :) = (Cn2b * vel)';

    %% 卡尔曼滤波
    F = zeros(N);
    F(4,2) = -f_n(3); %f_u天向比力
    F(4,3) = f_n(2); %f_n北向比力
    F(5,1) = f_n(3); %f_u天向比力
    F(5,3) = -f_n(1); %f_e东向比力
    F(6,1) = -f_n(2); %f_n北向比力
    F(6,2) = f_n(1); %f_e东向比力
    F(7:9, 4:6) = eye(3);
    F(1:3, 10:12) = -Cb2n;
    F(4:6, 13:15) =  Cb2n;
   % F(16:17, 16:17) = -1/3600;%% 相关时间
    % 状态转移矩阵F离散化
    F = eye(N) + F*imu_dt;

    % 卡尔曼时间更新
    X = F*X;
    P = F*P*F' + opt.Q*imu_dt;

    current_time = data.imu.tow(i);
    gnss_lost_elapsed = current_time - gnss_last_valid_time;

    %% GNSS量测更新
    if gnss_idx <= length(data.gnss.tow) && abs(data.imu.tow(i) - data.gnss.tow(gnss_idx)) < 0.01 % threshold 是允许的最大差异

        if data.gnss.solq_pos(gnss_idx) > 0

            gnss_vel_R = diag(ones(3,1) * data.gnss.gnss_vel_std_n(gnss_idx))^2 * 1;
            gnss_pos_R = diag(ones(3,1) * data.gnss.gnss_pos_std_n(gnss_idx))^2 * 1;

            % gnss_vel_R = diag([0.1 0.1 0.1])^2;
            % gnss_pos_R = diag([1.2 1.2 1.2])^2;

            Hvel = zeros(3,N);
            Hvel(1:3, 4:6) = eye(3);
            Zvel = vel - data.gnss.vel_enu(gnss_idx,:)';
            Zvel = Zvel - a_n*opt.gnss_delay; % GNSS量测延迟补偿
            Zvel = Zvel + (-Cb2n*v3_skew(w_b))*opt.gnss_lever_arm; % GNSS天线杆壁效应补偿

            Hpos = zeros(3,N);
            Hpos(1:3, 7:9) = eye(3);
            Zpos = pos - gnss_enu(gnss_idx,:)';

            Zpos = Zpos - vel*opt.gnss_delay; % GNSS量测延迟补偿
            Zpos = Zpos + (-Cb2n)*opt.gnss_lever_arm; % GNSS天线杆壁效应补偿

            if (opt.gnss_outage == 0 || (opt.gnss_outage == 1 && (current_time < opt.outage_start || current_time > opt.outage_stop))) ...
                    && (current_time - last_gnss_fusion_time >= opt.gnss_min_interval)

                % 速度位置更新：  Reference:  一种基于软卡方检测的自适应Ｋａｌｍａｎ滤波方法
                A = Hvel * P * Hvel' + gnss_vel_R;
                Kvel = P * Hvel' / A;
                innov_vel = Zvel - Hvel * X;

                chi_lambda_vel = innov_vel'*A^(-1)*innov_vel;
                log.lambda_vel(gnss_idx, :) = chi_lambda_vel;

                % 位置更新
                A = Hpos * P * Hpos' + gnss_pos_R;
                Kpos = P * Hpos' / (A);
                innov_pos = Zpos - Hpos * X;

                chi_lambda_pos = innov_pos'*A^(-1)*innov_pos;
                log.lambda_pos(gnss_idx, :) = chi_lambda_pos;
                
                if((chi_lambda_vel < 1.5 && chi_lambda_pos < 0.15) || (data.gnss.gnss_vel_std_n(gnss_idx) < 0.5 && data.gnss.gnss_pos_std_n(gnss_idx) < 5) || gnss_lost_elapsed > 3)
                    X = X + Kvel * innov_vel;
                    P = (eye(N) - Kvel * Hvel) * P;

                    X = X + Kpos * innov_pos;
                    P = (eye(N) - Kpos * Hpos) * P;

                    gnss_last_valid_time = data.gnss.tow(gnss_idx);
                end

                FB_BIT = bitor(FB_BIT, ESKF156_FB_A);
                FB_BIT = bitor(FB_BIT, ESKF156_FB_V);
                FB_BIT = bitor(FB_BIT, ESKF156_FB_P);
                FB_BIT = bitor(FB_BIT, ESKF156_FB_G);
                FB_BIT = bitor(FB_BIT, ESKF156_FB_W);
                
                % 更新上一次融合的时间
                last_gnss_fusion_time = current_time;
            end
        end
        gnss_idx = gnss_idx + 1;
    end

    %% NHC 约束
    if(opt.nhc_enable)
        norm_vel = norm(vel);

        if norm_vel > 0.2 && norm(w_b) < 20*D2R
            if gnss_lost_elapsed > 0.5
                H = zeros(2,N);
                A = [1 0 0; 0 0 1];
                H(:,4:6) = A*Cb2v*Cn2b;
                Z = 0 + (A*Cb2v*Cn2b)*vel;
                R = diag(ones(1, size(H, 1))*opt.nhc_R)^2;

                % 卡尔曼量测更新
                K = P * H' / (H * P * H' + R);
                X = X + K * (Z - H * X);
                P = (eye(N) - K * H) * P;
                FB_BIT = bitor(FB_BIT, ESKF156_FB_A);
                FB_BIT = bitor(FB_BIT, ESKF156_FB_V);
            else
                %% v frame 速度误差作为安装角估计的量测方程 [Vvx_ins,Vvy_ins,Vvz_ins]=dVv = Vv'-Vv = -Cb2v*Cn2b*v3_skew(Vn)*dphi + v3_skew(Cb2v*Cn2b*Vn)*dinstallangle + Cb2v*Cn2b*dVn
                %% paper ref:车载MIMUGNSS组合高精度无缝导航关键技术研究--sunzhenqian
                if (~opt.has_install_esti)
                    Cb2v = eye(3);
                end
                M1 = - Cb2v * Cn2b * v3_skew(vel);% v frame vel
                % M = Cn2b * v3_skew(data.gnss.vel_enu(gnss_idx,:));
                M2 = Cb2v * Cn2b;%M2 = Cb2v*Cn2b
                M3 = v3_skew(Cb2v * log.vb(i,:)');% M3=v3_skew(Cb2v*Cn2b*Vn)
                H = zeros(2, N);

                H(1, 1:3) = M1(1,:);
                H(1, 4:6) = M2(1,:);
                H(1, 17)  = M3(1,3);
                H(2, 1:3) = M1(3,:);
                H(2, 4:6) = M2(3,:);
                H(2, 16)  = M3(3,1);

                Vn2v_sins = Cb2v * log.vb(i, 1:3)';

                Z = Vn2v_sins([1,3]);

                R = diag(ones(1, size(H, 1))*opt.nhc_R)^2;

                % 卡尔曼量测更新
                K = P * H' / (H * P * H' + R);
                X = X + K * (Z - H * X);
                P = (eye(N) - K * H) * P;

                FB_BIT = bitor(FB_BIT, ESKF156_FB_A);
                FB_BIT = bitor(FB_BIT, ESKF156_FB_V);
                if (opt.has_install_esti)
                    FB_BIT = bitor(FB_BIT, ESKF156_FB_CBV);
                end
            end
        end
    end

    %% ZUPT检测和更新
    if opt.zupt_enable
        vel_norm = norm(vel);
        gyr_norm = norm(w_b);

        if vel_norm < opt.zupt_vel_threshold && gyr_norm < opt.zupt_gyr_threshold
            zupt_detect_time = zupt_detect_time + imu_dt;
            if zupt_detect_time >= opt.zupt_time_threshold
                is_zupt = true;
            end
        else
            zupt_detect_time = 0;
            is_zupt = false;
        end

        if is_zupt
            % ZUPT更新
            is_zupt = false;
            zupt_detect_time = 0;

            H_zupt = zeros(3, N);
            H_zupt(:, 4:6) = eye(3);  % 速度误差状态

            z_zupt = vel;  % 观测值：当前速度（应该接近零）

            % 卡尔曼滤波更新
            K_zupt = P * H_zupt' / (H_zupt * P * H_zupt' + opt.zupt_R);
            X = X + K_zupt * (z_zupt - H_zupt * X);
            P = (eye(N) - K_zupt * H_zupt) * P;
            
            FB_BIT = bitor(FB_BIT, ESKF156_FB_V);

        end
    end



    % 状态暂存
    X_temp = X;

    if bitand(FB_BIT, ESKF156_FB_A)
        rv = X(1:3);
        rv_norm = norm(rv);
        if rv_norm > 0
            qe = [cos(rv_norm/2); sin(rv_norm/2)*rv/rv_norm]';
            Qb2n = ch_qmul(qe, Qb2n);
            Qb2n = ch_qnormlz(Qb2n); %单位化四元数
            Cb2n = ch_q2m(Qb2n); %更新Cb2n阵
            Cn2b = Cb2n'; %更新Cn2b阵
            X(1:3) = 0;
        end

    end

    if bitand(FB_BIT, ESKF156_FB_V)
        vel = vel - X(4:6);
        X(4:6) = 0;
    end

    if bitand(FB_BIT, ESKF156_FB_P)
        pos = pos - X(7:9);
        X(7:9) = 0;
    end

    % 零偏反馈
    if bitand(FB_BIT, ESKF156_FB_W)
        gyro_bias = gyro_bias + X(10:12);
        X(10:12) = 0;
    end

    if bitand(FB_BIT, ESKF156_FB_G)
        acc_bias = acc_bias + X(13:15);
        X(13:15) = 0;
    end

    if bitand(FB_BIT, ESKF156_FB_CBV)
        cvv = att2Cnb([X(16),0,X(17)]);
        Cb2v = cvv*Cb2v;
        X(16:17) = 0;
        att = m2att(Cb2v);
    end



    %% 信息存储
    [pitch, roll, yaw] = q2att(Qb2n);
    log.pitch(i,:) = pitch;
    log.roll(i,:) = roll;
    log.yaw(i,:) = yaw;
    log.installangle(i,:)=att'*180/pi;
    log.vel(i,:) = vel';
    log.pos(i,:) = pos';
    log.X(i, :) = X_temp';
    log.P(i, :) = sqrt(diag(P))';
    log.gyro_bias(i, :) = gyro_bias;
    log.acc_bias(i, :) = acc_bias;
    log.tow(i,:) = current_time;

    % 纯惯性信息存储
    [pitch_sins, roll_sins, yaw_sins] = q2att(Qb2n_sins);
    log.sins_att(i,:) = [pitch_sins roll_sins yaw_sins];
end
imu_after_cal(j+1:end,:) = [];

fprintf('数据处理完毕，用时%.3f秒\n', toc);

set(groot, 'defaultAxesXGrid', 'on');
set(groot, 'defaultAxesYGrid', 'on');
set(groot, 'defaultAxesZGrid', 'on');
%% 姿态与航向估计曲线
% plot_att(data.imu.tow,log.att, mcu_time,mcu.att, data.imu.tow,log.sins_att, data.imu.tow,[bl_pitch bl_yaw]);
%
% figure('name', "MGNSS组合导航航向与双天线航向");
% subplot(2,1,1);
% plot(data.imu.tow, mcu.att(:,3),  data.imu.tow, bl_yaw, '.-'); grid on;
% xlim([data.imu.tow(1) data.imu.tow(end)]);
% ylim([-10 370]);
% yticks(0:45:360);
% legend("MCU YAW", "GNSS DUAL YAW");
% subplot(2,1,2);
% plot(data.imu.tow, atand(tand(mcu.att(:,3) - bl_yaw)));  grid on;
% xlim([data.imu.tow(1) data.imu.tow(end)]);
% set(gcf, 'Units', 'normalized', 'Position', [0.025, 0.05, 0.95, 0.85]);

navigation_plots.axis_comparison('XLabel', '时间(s)', ...
                                      'YLabel', {{'East (m)', 'North (m)', 'Up (m)'}, {'Std East (m)', 'Std North (m)', 'Std Up (m)'}}, ...
                                      'Title', '位置对比', ...
                                      'Time', {data.imu.tow, data.imu.tow, data.imu.tow, data.imu.tow}, ...
                                      'Values', {log.pos, data.dev.pos_enu, log.P(:,7:9), data.dev.kf_p_pos}, ...
                                      'Labels', {'Matlab', 'HI30', 'Matlab_P', 'HI30_P'}, ...
                                      'Subplots', [1, 1, 2, 2]);


navigation_plots.axis_comparison('XLabel', '时间(s)', ...
                                      'YLabel', {{'East (m/s)', 'North (m/s)', 'Up (m/s)'}, {'Std East (m/s)', 'Std North (m/s)', 'Std Up (m/s)'}}, ...
                                      'Title', '速度对比', ...
                                      'Time', {data.imu.tow, data.imu.tow, data.imu.tow, data.imu.tow}, ...
                                      'Values', {log.vel, data.dev.vel_enu, log.P(:,4:6), data.dev.kf_p_vel}, ...
                                      'Labels', {'Matlab', 'HI30', 'Matlab_P', 'HI30_P'}, ...
                                      'Subplots', [1, 1, 2, 2]);

navigation_plots.axis_comparison('XLabel', '时间(s)', ...
                                      'YLabel', {{'Roll (deg)', 'Pitch (deg)', 'Yaw (deg)'}, {'Std Roll (deg)', 'Std Pitch (deg)', 'Std Yaw (deg)'}}, ...
                                      'Title', '姿态对比', ...
                                      'Time', {data.imu.tow, data.imu.tow, data.imu.tow, data.imu.tow}, ...
                                      'Values', {[log.roll log.pitch log.yaw], [data.dev.roll data.dev.pitch data.dev.yaw], log.P(:,1:3)*R2D, data.dev.kf_p_att}, ...
                                      'Labels', {'Matlab', 'HI30', 'Matlab_P', 'HI30_P'}, ...
                                      'Subplots', [1, 1, 2, 2]);


navigation_plots.axis_comparison('XLabel', '时间(s)', ...
                                      'YLabel', {{'X (deg)', 'Y (deg)', 'Z (deg)'}, {'Std X (deg)', 'Std Y (deg)', 'Std Z (deg)'}}, ...
                                      'Title', '陀螺零偏对比', ...
                                      'Time', {data.imu.tow, data.imu.tow, data.imu.tow, data.imu.tow, data.imu.tow}, ...
                                      'Values', {log.gyro_bias*R2D, repmat(gyr_bias0*R2D, length(data.imu.tow), 1), data.dev.kf_wb*R2D, log.P(:, 10:12)*R2D, data.dev.kf_p_wb*R2D}, ...
                                      'Labels', {'Matlab', 'bias0', 'HI30', 'Matlab_P', 'HI30_P'}, ...
                                      'Subplots', [1, 1, 1, 2, 2], ...
                                      'LineStyles', {'-', '-.', ':', '-'});

navigation_plots.axis_comparison('XLabel', '时间(s)', ...
                                      'YLabel', {{'X (mG)', 'Y (mG)', 'Z (mG)'}, {'Std X (mG)', 'Std Y (mG)', 'Std Z (mG)'}}, ...
                                      'Title', '加计零偏对比', ...
                                      'Time', {data.imu.tow, data.imu.tow, data.imu.tow, data.imu.tow, data.imu.tow}, ...
                                      'Values', {log.acc_bias*1000/GRAVITY, data.dev.kf_gb*1000/GRAVITY, log.P(:, 13:15)*1000/GRAVITY, data.dev.kf_p_gb*1000/GRAVITY}, ...
                                      'Labels', {'Matlab', 'HI30', 'Matlab_P', 'HI30_P'}, ...
                                      'Subplots', [1, 1, 2, 2]);



%% 状态量曲线
figure('name','状态量曲线');
subplot(2,2,1);
plot(data.imu.tow, log.X(:, 1) * R2D, 'c', 'linewidth', 1.5); hold on; grid on;
plot(data.imu.tow, log.X(:, 2) * R2D, 'm', 'linewidth', 1.5);
plot(data.imu.tow, log.P(:, 1) * R2D * 1, 'r-.', 'linewidth', 1);
plot(data.imu.tow, log.P(:, 2) * R2D * 1, 'g-.', 'linewidth', 1);
plot(data.imu.tow, log.P(:, 1) * R2D * -1, 'r-.', 'linewidth', 1);
plot(data.imu.tow, log.P(:, 2) * R2D * -1, 'g-.', 'linewidth', 1);
xlim tight;
ylim([-0.5 0.5]);
xlabel('时间(s)'); ylabel('平台失准角(°)'); legend('Pitch', 'Roll', 'Orientation','horizontal');

subplot(2,2,3);
plot(data.imu.tow, log.X(:, 3) * R2D, 'c', 'linewidth', 1.5); hold on; grid on;
plot(data.imu.tow, log.P(:, 3) * R2D * 1, 'b-.', 'linewidth', 1);
plot(data.imu.tow, log.P(:, 3) * R2D * -1, 'b-.', 'linewidth', 1);
xlim tight;
ylim([-5 5]);
xlabel('时间(s)'); ylabel('平台失准角(°)'); legend('Yaw', 'Orientation','horizontal');

subplot(2,2,2);
plot(data.imu.tow, log.X(:, 4:6), 'linewidth', 1.5); hold on; grid on;
plot(data.imu.tow, log.P(:, 4:6)  * 1, '-.', 'linewidth', 1);
plot(data.imu.tow, log.P(:, 4:6)  * -1, '-.', 'linewidth', 1);
xlim tight;
ylim([-10 10]);
xlabel('时间(s)'); ylabel('速度误差(m/s)'); legend('E', 'N', 'U', 'Orientation','horizontal');

subplot(2,2,4);
plot(data.imu.tow, log.X(:, 7:9), 'linewidth', 1.5); hold on; grid on;
plot(data.imu.tow, log.P(:, 7:9) * 1, '-.', 'linewidth', 1);
plot(data.imu.tow, log.P(:, 7:9) * -1, '-.', 'linewidth', 1);
xlim tight;
ylim([-100 100]);
xlabel('时间(s)'); ylabel('位置误差(m)'); legend('E', 'N', 'U', 'Orientation','horizontal');

set(gcf, 'Units', 'normalized', 'Position', [0.025, 0.05, 0.95, 0.85]);

figure('name', 'GNSS接收机给出的信息');
subplot(2,2,1);
plot(data.gnss.tow, data.gnss.gnss_vel_std_n); hold on; grid on;
plot(data.gnss.tow, log.lambda_vel);
xlabel('时间(s)'); legend('gnss_vel_std_n', 'log.lambda_vel', 'Orientation','horizontal');

subplot(2,2,2);
plot(data.gnss.tow, data.gnss.gnss_pos_std_n); hold on; grid on;
plot(data.gnss.tow, log.lambda_pos);
xlabel('时间(s)'); legend('gnss_pos_std_n', 'log.lambda_pos', 'Orientation','horizontal');

subplot(2,2,3);
plot(data.gnss.tow, data.gnss.hdop); grid on;
xlabel('时间(s)'); legend('hdop', 'Orientation','horizontal');

subplot(2,2,4);
plot(data.gnss.tow, data.gnss.nv); grid on;
xlabel('时间(s)'); legend('卫星数', 'Orientation','horizontal');

%% 安装误差角在线估计结果
figure('name', '安装误差角在线估计结果');
subplot(2,2,1);
plot(data.imu.tow, log.installangle(:,1), 'LineWidth', 1.5); grid on;
xlabel('时间(s)'); ylabel('俯仰安装角(°)'); xlim tight;
subplot(2,2,3);
plot(data.imu.tow, log.P(:, 16)*R2D, 'LineWidth', 1.5); grid on;
xlabel('时间(s)'); ylabel('俯仰安装角方差(°)'); xlim tight;

subplot(2,2,2);
plot(data.imu.tow, log.installangle(:,3), 'LineWidth', 1.5); grid on;
xlabel('时间(s)'); ylabel('航向安装角(°)'); xlim tight;
subplot(2,2,4);
plot(data.imu.tow, log.P(:, 17)*R2D, 'LineWidth', 1.5); grid on;
xlabel('时间(s)'); ylabel('航向安装角方差(°)'); xlim tight;
%
% subplot(111);
% plot(data.imu.tow, log.installangle(:,2), 'LineWidth', 1.5); grid on;
% xlabel('时间(s)'); ylabel('滚动安装角(°)'); xlim tight;

set(gcf, 'Units', 'normalized', 'Position', [0.025, 0.05, 0.95, 0.85]);

navigation_plots.trajectory_2d_plot({gnss_enu, log.pos, dev_pos_enu}, {'GNSS', 'MATLAB', 'DEV嵌入式设备轨迹'});

%% 数据统计
time_duration = seconds(time_sum);
time_duration.Format = 'hh:mm:ss';

fprintf('行驶时间: %s\n', char(time_duration));
fprintf('行驶距离: %.3fkm\n', distance_sum/1000);
fprintf('最高时速: %.3fkm/h\n', max(log.vel_norm)*3.6);

