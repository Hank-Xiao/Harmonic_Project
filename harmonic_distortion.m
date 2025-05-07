clear;

Fs=512;%采样频率
Ts=1/Fs;%采样周期
N=512;%N—DFT
% L=20;%原信号序列长度
t=(0:N-1).*Ts;%时域自变量
f = [50.1 100 150 200 250 300 350 400 450 500 550]; % Frequency of signal components (Hz)
A = [220 1.4 6.3 1.0 5.5 0.6 3.2 0.2 1.5 0.1 0.9]; % Amplitude of signal components
ph = deg2rad([30 39.3 44.5 122.3 60 75.5 20 143.2 88.3 122.3 90]); % Phase of signal components (rad)
x = zeros(1,length(t));
for k = 1:length(f)
      x = x + A(k)*cos(2*pi*f(k)*t+ph(k));
end

%% 4项3阶Nuttall窗函数参数
a0 = 0.355768;
a1 = 0.487396;
a2 = 0.144232;
a3 = 0.012604;
N = 512;                % 窗长
n = 0:N-1;              % 时间索引

%% 计算Nuttall窗
w1 = a0 - a1*cos(2*pi*n/(N-1)) + a2*cos(4*pi*n/(N-1)) - a3*cos(6*pi*n/(N-1)); %4项3阶Nuttall窗函数
y=x.*w1;
% x= 3.2*cos(25*2*pi*t+pi/9)+380*cos(50*2*pi*t+pi/3)+20*cos(150*2*pi*t+pi/6)+2.8*cos(180*2*pi*t+25/180*pi)+16.4*cos(250*2*pi*t+100/180*pi)+...
            % +1.8*cos(320*2*pi*t+2*pi/3)+12*cos(350*2*pi*t+pi*5/6)+1.4*cos(380*2*pi*t+pi)+10*cos(450*2*pi*t+7*pi/6); 
 % plot(t,y,'k','LineWidth',1.2);
y=fft(y);
plot(t,y,'k');
% xlabel("Time");
% ylabel('Amplitude');
% % legend('参考信号','估计信号');
% noise = 0.01*randn(size(x));
% x = x + noise;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %加窗效果有待验证%上和%下的代码都可以跑通
% win = nuttall(N, 'numTerms', 4, 'Order', 3);
% win = win/sum(win); % Ensures that the sum of the window is 1
% % Convolution of windows (double windowing)
% win_conv = conv(win,win);
% win_conv = win_conv/sum(win_conv); % Ensures that the sum of the 
% % window is 1
% % Single windowing of last N samples for FFT
% x_win = x(end-N+1:end).*win;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% legend("估计信号")
% Y=fft(y)/N*2.8;%fft变换
% P=angle(Y);
% Y=abs(Y);%实际幅值变换
% f=(0:N-1)*Fs./N;%实际频率变换
% % subplot(211);
% Y(2:end) = 2*Y(2:end);
% plot(f(1:N/2),Y(1:N./2),'r',LineWidth=1.1);
% xlabel("f/Hz")
% grid on
% f=f-Fs./2;%移位
% subplot(212);
% stem(f,fftshift(Y));%移位
% title("N-DFT变换幅频响应双边")
% xlabel("f/Hz")
% grid on

% figure
% subplot(211)
% plot(t,x)
% title("原信号")
% xlabel("t/s")
% grid on
% Y=fft(x);
% xx=ifft(Y);
% subplot(212)
% plot(t,xx)
% title("ifft原信号")
% xlabel("t/s")
% grid on

% % 
% t=0:0.2/1023:0.2;
% x=100*cos(50*2*pi*t)+30*cos(100*2*pi*t)+80*cos(150*2*pi*t)+15*cos(200*2*pi*t)+20*cos(250*2*pi*t)...
%     +8*cos(300*2*pi*t)+10*cos(350*2*pi*t)+2*cos(400*2*pi*t)+5*cos(450*2*pi*t)+cos(500*2*pi*t); %参考项,10次谐波
% y=100*cos(50*1.2*2*pi*t)+30*cos(100*0.9*2*pi*t)+80*cos(150*1.15*2*pi*t)+15*cos(200*0.98*2*pi*t)+20*cos(1.02*250*2*pi*t)...
%     +8*cos(300*1.04*2*pi*t)+10*cos(350*2*pi*t)+2*cos(400*1.1*2*pi*t)+5*cos(450*2*pi*t)+cos(0.99*500*2*pi*t); %参考项,10次谐波
% plot(t,x,'b');
% grid on
% hold on
% plot(t,y,'ro-');
% 
% % Parameters
% f0 = 50;             % Fundamental frequency in Hz
% A0 = 100;            % Amplitude of the fundamental wave
% harmonics = 10;      % Number of harmonics to include
% fs = 5000;           % Sampling frequency in Hz
% N = 1024;            % Number of sampling points
% t = (0:N-1) / fs;    % Time vector based on sampling frequency
% 
% % Initialize waveform with the fundamental frequency
% y = A0 * sin(2 * pi * f0 * t);
% 
% % Add harmonics
% for k = 2:harmonics
%     Ak = A0 / k;           % Amplitude of the k-th harmonic
%     fk = k * f0;           % Frequency of the k-th harmonic
%     y = y + Ak * sin(2 * pi * fk * t);
% end
% 
% % Perform DFT using FFT
% Y = fft(y, N);
% 
% % Perform inverse DFT using IFFT
% y_inv = ifft(Y, N);
% 
% % Plot the original and inverse transformed waveforms
% figure;
% subplot(2,1,1);
% plot(t, y);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Original Waveform');
% grid on;
% 
% subplot(2,1,2);
% plot(t, real(y_inv));
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Waveform after Inverse DFT');
% grid on;

% % 参数设置
% f0 = 50;             % 基波频率 (Hz)
% A0 = 100;            % 基波幅度
% harmonics = 10;      % 谐波数量
% fs = 5000;           % 采样频率 (Hz)
% N = 1024;            % 采样点数
% t = (0:N-1) / fs;    % 时间向量
% 
% % 构造原始谐波波形
% y_orig = A0 * sin(2 * pi * f0 * t);  
% for k = 2:harmonics
%     Ak = A0 / k;                    % k次谐波幅度
%     fk = k * f0;                    % k次谐波频率
%     y_orig = y_orig + Ak * sin(2 * pi * fk * t);
% end
% 
% % 初始猜测：间谐波幅度与频率
% initial_params = [5 * ones(1, harmonics), 0.1 * f0 * ones(1, harmonics)];
% 
% % 全局变量用于记录目标函数值及谐波参数
% global objective_values amplitude_history frequency_history;
% objective_values = [];
% amplitude_history = []; % 存储每次迭代的谐波幅度
% frequency_history = []; % 存储每次迭代的谐波频率
% 
% % 优化目标函数定义
% objective_func = @(params) compute_error(params, y_orig, f0, t, harmonics);
% 
% % 使用fmincon进行优化
% options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
% lb = zeros(1, 2 * harmonics);  % 参数下界（幅度和频率均大于0）
% ub = [A0 * ones(1, harmonics), 2 * f0 * ones(1, harmonics)];  % 参数上界
% 
% optimal_params = fmincon(objective_func, initial_params, [], [], [], [], lb, ub, [], options);
% 
% % 计算最优参数下的间谐波波形
% [y_inter_opt] = generate_waveform(optimal_params, f0, t, harmonics);
% 
% % 绘制结果
% figure;
% subplot(3,1,1);
% plot(t, y_orig);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Original Harmonic Waveform');
% grid on;
% 
% subplot(3,1,2);
% plot(t, y_inter_opt);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Optimized Waveform with Interharmonics');
% grid on;
% 
% subplot(3,1,3);
% plot(objective_values);
% xlabel('Iteration');
% ylabel('Objective Function Value');
% title('Objective Function Value per Iteration');
% grid on;

% % 绘制谐波参数变化
% figure;
% for k = 1:harmonics
%     subplot(harmonics,2,2*k-1);
%     plot(amplitude_history(:, k));
%     xlabel('Iteration');
%     ylabel(['Amplitude of Harmonic ', num2str(k)]);
%     title(['Amplitude of Harmonic ', num2str(k)]);
%     grid on;
% 
%     % subplot(harmonics,2,2*k);
%     % plot(frequency_history(:, k));
%     % xlabel('Iteration');
%     % ylabel(['Frequency of Harmonic ', num2str(k)]);
%     % title(['Frequency of Harmonic ', num2str(k)]);
%     % grid on;
% end
% 
% %% 目标函数计算：原始波形与间谐波波形之间的误差平方和
% function error = compute_error(params, y_orig, f0, t, harmonics)
%     global objective_values amplitude_history frequency_history;
%     y_inter = generate_waveform(params, f0, t, harmonics);
%     error = sum((y_orig - y_inter).^2);  % 计算误差的平方和
%     objective_values = [objective_values; error]; % 记录每次迭代的目标函数值
% 
%     % 分解参数，记录每个谐波的幅度和频率
%     amplitude_history = [amplitude_history; params(1:harmonics)];
%     frequency_history = [frequency_history; params(harmonics+1:end)];
% end
% 
% %% 根据参数生成间谐波波形
% function y_inter = generate_waveform(params, f0, t, harmonics)
%     y_inter = zeros(size(t));
%     for k = 1:harmonics
%         Ak = params(k);                    % 第k个间谐波的幅度
%         fk = k * f0 + params(harmonics + k);  % 第k个间谐波的频率
%         y_inter = y_inter + Ak * sin(2 * pi * fk * t);
%     end
% end
% 
% 
% 
