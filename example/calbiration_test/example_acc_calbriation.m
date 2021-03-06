%% Example: 最小二乘法校准加速度计
clc;
clear;
close all;
format short;
%% 校准数据
input =[
    1.0556    0.0194   -0.0651
    0.0195    0.9952   -0.0430
    0.0714   -0.0175    0.9652
    -0.9272   -0.0181   -0.0296
    0.0921   -0.9974   -0.0491
    0.0503    0.0159   -1.0561
    ];

[C, B] = acc_calibration(input);

fprintf('校准矩阵:');
C
fprintf('零偏:');
B

%% 计算校准后的数据
output(1,:) = C*(input(1,:)') - B;
output(2,:) = C*(input(2,:)') - B;
output(3,:) = C*(input(3,:)') - B;

output(4,:) = C*(input(4,:)') - B;
output(5,:) = C*(input(5,:)') - B;
output(6,:) = C*(input(6,:)') - B;



%% 计算校准前后误差
X =  input - [1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1];
error_input =  sum(sum(abs(X).^2, 2).^(1/2));

X =  output - [1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1];
error_ouput =  sum(sum(abs(X).^2, 2).^(1/2));
fprintf('校准前误差: %f    校准后误差: %f\n', error_input, error_ouput);

%% plot
grid on;
plot3(input(:,1), input(:,2), input(:,3), 'or');
hold on;
plot3(output(:,1), output(:,2), output(:,3), '*b');
axis equal

legend('输入', '校准后');

