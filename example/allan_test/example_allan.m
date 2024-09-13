clear;
clc;
close all;

%% : read matlabe internal dataset 
% load('LoggedSingleAxisGyroscope', 'omega', 'Fs')
% gyroReading = omega;

%% ADIS16488 dataset
% load('adis16488_gyr.mat');
% gyroReading = omega(:,1);
% gyroReading = deg2rad(gyroReading);

%%  ch00   deg/s    G

% load('./allan_plot_ch100/ch100.mat');
% gyroReading = imu.gyr(:,2);
% gyroReading = deg2rad(gyroReading);
% accelReading = imu.acc(:,2)*9.80665;
% Fs = 400;

% load('./allan_plot_ch100/hi226_hi229.mat');
% gyroReading = imu.gyr(:,1);
% gyroReading = deg2rad(gyroReading);
% accelReading = imu.acc(:,1)*9.80665;
% Fs = 400;

% 
% load('./stim300.mat');
% gyroReading = imu.gyr;
% accelReading = imu.acc;
% 
% 
% gyroReading = gyroReading(:,3);
% gyroReading = deg2rad(gyroReading);
% accelReading = accelReading(:,3)*9.80665;
% Fs = 100;

% glvs
% 
% ts = 1 / Fs;
% 
% dd = load('stim300.txt');
% 
% imu = [[imu.gyr*glv.dps,  imu.acc*glv.g0]*ts, (1:length(dd))'*ts];
% 
% avars(imu(:,1:3)/ts/glv.dph, ts);



%% https://github.com/Aceinna/gnss-ins-sim
% VRW ��λ:                    m/s/sqrt(hr)           'accel_vrw': np.array([0.3119, 0.6009, 0.9779]) * 1.0,
% ���ٶ���ƫ���ȶ���:     m/s^(2)                 'accel_b_stability': np.array([1e-3, 3e-3, 5e-3]) * 1.0,
% ARW ��λ:                    deg/sqrt(hr)           'gyro_arw': np.array([0.25, 0.25, 0.25]) * 1.0,
% ���ٶ���ƫ�ȶ���:        deg/h                     'gyro_b_stability': np.array([3.5, 3.5, 3.5]) * 1.0,  

%���ݵ�λ: rad/s, m/s^(2)
M = csvread('gyro-0.csv',1 ,0);
gyroReading = M(:,1);
gyroReading = deg2rad(gyroReading); 

M = csvread('accel-0.csv',1 ,0);
Fs = 100;
accelReading = M(:,1);


 

            
%% ���� allan  ��λ����ת��Ϊ ����: rad/s,  ���ٶ�:m/s^(2s)
[avar1, tau1 , N, K, B] = ch_allan(gyroReading , Fs);
fprintf('GYR: ��ƫ���ȶ���                                                             %frad/s                    ��   %fdeg/h \n', B, rad2deg(B)*3600);
fprintf('GYR: �����ܶ�(�Ƕ��������, ARW, Noise density)              %f(rad/s)/sqrt(Hz)     ��  %f deg/sqrt(h)\n', N, rad2deg(N)*3600^(0.5));
fprintf('GYR: �������������, bias variations ("random walks")       %f(rad/s)sqrt(Hz)       ��  %f deg/h/sqrt(h)\n', K, rad2deg(K) * (3600^(1.5)));



%% ���ٶȼ� allan
[avar1, tau1 , N, K, B] = ch_allan(accelReading, Fs);

fprintf('ACC: ��ƫ���ȶ���                                                                                       %fm/s^(2)                       ��   %fmg  ��  %fug\n', B, B / 9.80665 *1000,  B / 9.80665 *1000*1000);
fprintf('ACC: �����ܶ�(�����������,VRW, Noise Density, Rate Noise Density)          %f(m/s^(2))/sqrt(Hz)        ��   %f m/s/sqrt(hr)\n', N, N * 3600^(0.5));
fprintf('ACC: ���ٶ�������ߣ�bias variations ("random walks")                               %f(m/s^(2)sqrt(Hz))\n',  K);

