clc; clear; close all;

%% 仿真参数
SNR_dB = 0:2:30;                % Es/N0 范围 (dB)
numSymbols = 2e5;               % 每个SNR点的符号数（200000）
mod_list = {'QPSK','8PSK','16PSK','16QAM','32QAM'};
M_list   = [4, 8, 16, 16, 32];
type_list= {'psk','psk','psk','qam','qam'};

% 格雷码映射表
grayQPSK  = [0 1 3 2];
gray8PSK  = [0 1 3 2 7 6 4 5];
gray16QAM = [0 1 3 2 7 6 4 5 15 14 12 13 8 9 11 10];
gray32QAM = [0 1 3 2 7 6 4 5 15 14 12 13 8 9 11 10 31 30 28 29 24 25 27 26 16 17 19 18 23 22 20 21];
% 16PSK格雷码
M16 = 16;
gray16PSK = zeros(1, M16);
for i = 0:M16-1
    gray16PSK(i+1) = bitxor(i, bitshift(i, -1));
end

% 预分配BER矩阵
ber_sim = zeros(length(mod_list), length(SNR_dB));

%% 主仿真循环
for m = 1:length(mod_list)
    M = M_list(m);
    mod_type = type_list{m};
    
    % 选择对应的格雷码表
    switch mod_list{m}
        case 'QPSK'
            gray_map = grayQPSK;
            phase_offset = pi/4;
            bits_per_sym = 2;
        case '8PSK'
            gray_map = gray8PSK;
            phase_offset = pi/8;
            bits_per_sym = 3;
        case '16PSK'
            gray_map = gray16PSK;
            phase_offset = pi/16;
            bits_per_sym = 4;
        case '16QAM'
            gray_map = gray16QAM;
            bits_per_sym = 4;
        case '32QAM'
            gray_map = gray32QAM;
            bits_per_sym = 5;
    end
    
    for idx = 1:length(SNR_dB)
        EsN0_lin = 10^(SNR_dB(idx)/10);
        N0 = 1 / EsN0_lin;               % 噪声功率谱密度（信号功率=1）
        
        % --- 生成符号索引（列向量）---
        dataIn = randi([0 M-1], numSymbols, 1);
        
        % --- 格雷映射 ---
        dataGray = gray_map(dataIn + 1);   % 仍然是列向量
        
        % --- 调制 ---
        if strcmp(mod_type, 'psk')
            tx = pskmod(dataGray, M, phase_offset);
        else
            tx = qammod(dataGray, M, 'UnitAveragePower', true);
        end
        % 确保 tx 是列向量
        tx = tx(:);
        
        % --- 生成平坦瑞利衰落系数 h ~ CN(0,1) ---
        h = (randn(numSymbols, 1) + 1j*randn(numSymbols, 1)) / sqrt(2);
        h = h(:);   % 确保列向量
        
        % --- 生成AWGN噪声（列向量）---
        noise = sqrt(N0/2) * (randn(numSymbols, 1) + 1j*randn(numSymbols, 1));
        noise = noise(:);
        
        % --- 关键检查：防止意外矩阵 ---
        if ~isvector(tx) || ~isvector(h) || ~isvector(noise)
            error('tx, h, 或 noise 不是向量！');
        end
        if length(tx) ~= numSymbols || length(h) ~= numSymbols || length(noise) ~= numSymbols
            error('变量长度不等于 numSymbols！');
        end
        
        % --- 接收信号（点乘，不产生矩阵）---
        rx = h .* tx + noise;
        
        % --- 已知信道均衡 ---
        rx_eq = rx ./ h;
        
        % --- 解调 ---
        if strcmp(mod_type, 'psk')
            rx_demod = pskdemod(rx_eq, M, phase_offset);
        else
            rx_demod = qamdemod(rx_eq, M, 'UnitAveragePower', true);
        end
        
        % --- 格雷逆映射 ---
        [~, dataOut] = ismember(rx_demod, gray_map);
        dataOut = dataOut - 1;
        
        % --- 计算误比特率 ---
        bits_tx = de2bi(dataIn, bits_per_sym);
        bits_rx = de2bi(dataOut, bits_per_sym);
        [~, ber_sim(m, idx)] = biterr(bits_tx, bits_rx);
        
        % 进度显示
        fprintf('已完成: %s, SNR=%.1fdB, BER=%.2e\n', mod_list{m}, SNR_dB(idx), ber_sim(m,idx));
    end
end

%% 理论BER计算（瑞利衰落，1阶分集）
EbN0_dB_QPSK   = SNR_dB - 10*log10(2);
EbN0_dB_8PSK   = SNR_dB - 10*log10(3);
EbN0_dB_16PSK  = SNR_dB - 10*log10(4);
EbN0_dB_16QAM  = SNR_dB - 10*log10(4);
EbN0_dB_32QAM  = SNR_dB - 10*log10(5);

ber_QPSK_theory   = berfading(EbN0_dB_QPSK,  'psk', 4,  1);
ber_8PSK_theory   = berfading(EbN0_dB_8PSK,  'psk', 8,  1);
ber_16PSK_theory  = berfading(EbN0_dB_16PSK, 'psk', 16, 1);
ber_16QAM_theory  = berfading(EbN0_dB_16QAM, 'qam', 16, 1);
ber_32QAM_theory  = berfading(EbN0_dB_32QAM, 'qam', 32, 1);

%% 绘图
figure;
semilogy(SNR_dB, ber_sim(1,:), 'bo-'); hold on;
semilogy(SNR_dB, ber_sim(2,:), 'r*-');
semilogy(SNR_dB, ber_sim(3,:), 'ks-');
semilogy(SNR_dB, ber_sim(4,:), 'gs-');
semilogy(SNR_dB, ber_sim(5,:), 'md-');

semilogy(SNR_dB, ber_QPSK_theory,   'b--');
semilogy(SNR_dB, ber_8PSK_theory,   'r--');
semilogy(SNR_dB, ber_16PSK_theory,  'k--');
semilogy(SNR_dB, ber_16QAM_theory,  'g--');
semilogy(SNR_dB, ber_32QAM_theory,  'm--');

xlabel('Es/N0 (dB)');
ylabel('BER');
legend('QPSK仿真', '8PSK仿真', '16PSK仿真', '16QAM仿真', '32QAM仿真', ...
       'QPSK理论', '8PSK理论', '16PSK理论', '16QAM理论', '32QAM理论', ...
       'Location', 'southwest');
title('平坦瑞利衰落信道下不同调制方式的BER性能对比');
grid on;
hold off;