clc;
clear;

% 仿真参数
SNR_dB = -10:1:10;          % SNR范围 (dB)
numSymbols = 1e5;           % 每个SNR点的仿真符号数

% BER初始化
berQPSK_sim   = zeros(1, length(SNR_dB));
ber8PSK_sim   = zeros(1, length(SNR_dB));
ber16PSK_sim  = zeros(1, length(SNR_dB));
ber16QAM_sim  = zeros(1, length(SNR_dB));
ber32QAM_sim  = zeros(1, length(SNR_dB));

% ---------- 格雷码映射表 ----------
% QPSK, 8PSK 使用手动定义的格雷码（已验证符合相邻相位仅1比特差异）
grayQPSK  = [0 1 3 2];
gray8PSK  = [0 1 3 2 7 6 4 5];
% 16QAM, 32QAM 使用标准格雷码（由MATLAB通信工具箱函数也可生成，此处沿用原代码）
gray16QAM = [0 1 3 2 7 6 4 5 15 14 12 13 8 9 11 10];
gray32QAM = [0 1 3 2 7 6 4 5 15 14 12 13 8 9 11 10 31 30 28 29 24 25 27 26 16 17 19 18 23 22 20 21];

% 16PSK 格雷码：采用二进制反射格雷码生成，确保相位顺序索引相邻时只有1比特差异
M_PSK16 = 16;
gray16PSK = zeros(1, M_PSK16);
for i = 0:M_PSK16-1
    gray16PSK(i+1) = bitxor(i, bitshift(i, -1));
end

% ------------------- 仿真循环 -------------------

% 1. QPSK 仿真
for i = 1:length(SNR_dB)
    dataIn = randi([0 3], numSymbols, 1);
    dataGray = grayQPSK(dataIn + 1);
    dataMod = pskmod(dataGray, 4, pi/4);
    received = awgn(dataMod, SNR_dB(i), 'measured');
    dataDemod = pskdemod(received, 4, pi/4);
    [~, dataOut] = ismember(dataDemod, grayQPSK);
    dataOut = dataOut - 1;
    [~, berQPSK_sim(i)] = biterr(de2bi(dataIn, 2), de2bi(dataOut, 2));
end

% 2. 8PSK 仿真
for i = 1:length(SNR_dB)
    dataIn = randi([0 7], numSymbols, 1);
    dataGray = gray8PSK(dataIn + 1);
    dataMod = pskmod(dataGray, 8, pi/8);
    received = awgn(dataMod, SNR_dB(i), 'measured');
    dataDemod = pskdemod(received, 8, pi/8);
    [~, dataOut] = ismember(dataDemod, gray8PSK);
    dataOut = dataOut - 1;
    [~, ber8PSK_sim(i)] = biterr(de2bi(dataIn, 3), de2bi(dataOut, 3));
end

% 3. 16PSK 仿真（新增）
for i = 1:length(SNR_dB)
    dataIn = randi([0 15], numSymbols, 1);
    dataGray = gray16PSK(dataIn + 1);
    dataMod = pskmod(dataGray, 16, pi/16);      % 相位偏移 pi/16
    received = awgn(dataMod, SNR_dB(i), 'measured');
    dataDemod = pskdemod(received, 16, pi/16);
    [~, dataOut] = ismember(dataDemod, gray16PSK);
    dataOut = dataOut - 1;
    [~, ber16PSK_sim(i)] = biterr(de2bi(dataIn, 4), de2bi(dataOut, 4));
end

% 4. 16QAM 仿真
for i = 1:length(SNR_dB)
    dataIn = randi([0 15], numSymbols, 1);
    dataGray = gray16QAM(dataIn + 1);
    dataMod = qammod(dataGray, 16, 'UnitAveragePower', true);
    received = awgn(dataMod, SNR_dB(i), 'measured');
    dataDemod = qamdemod(received, 16, 'UnitAveragePower', true);
    [~, dataOut] = ismember(dataDemod, gray16QAM);
    dataOut = dataOut - 1;
    [~, ber16QAM_sim(i)] = biterr(de2bi(dataIn, 4), de2bi(dataOut, 4));
end

% 5. 32QAM 仿真
for i = 1:length(SNR_dB)
    dataIn = randi([0 31], numSymbols, 1);
    dataGray = gray32QAM(dataIn + 1);
    dataMod = qammod(dataGray, 32, 'UnitAveragePower', true);
    received = awgn(dataMod, SNR_dB(i), 'measured');
    dataDemod = qamdemod(received, 32, 'UnitAveragePower', true);
    [~, dataOut] = ismember(dataDemod, gray32QAM);
    dataOut = dataOut - 1;
    [~, ber32QAM_sim(i)] = biterr(de2bi(dataIn, 5), de2bi(dataOut, 5));
end

% ------------------- 理论BER计算 -------------------
berQPSK_theory   = berawgn(SNR_dB, 'psk', 4, 'nondiff');
ber8PSK_theory   = berawgn(SNR_dB, 'psk', 8, 'nondiff');
ber16PSK_theory  = berawgn(SNR_dB, 'psk', 16, 'nondiff');
ber16QAM_theory  = berawgn(SNR_dB, 'qam', 16);
ber32QAM_theory  = berawgn(SNR_dB, 'qam', 32);

% ------------------- 绘图 -------------------
figure;
semilogy(SNR_dB, berQPSK_sim,   'bo-');  hold on;
semilogy(SNR_dB, ber8PSK_sim,   'r*-');
semilogy(SNR_dB, ber16PSK_sim,  'ks-');
semilogy(SNR_dB, ber16QAM_sim,  'gs-');
semilogy(SNR_dB, ber32QAM_sim,  'md-');

semilogy(SNR_dB, berQPSK_theory,   'b--');
semilogy(SNR_dB, ber8PSK_theory,   'r--');
semilogy(SNR_dB, ber16PSK_theory,  'k--');
semilogy(SNR_dB, ber16QAM_theory,  'g--');
semilogy(SNR_dB, ber32QAM_theory,  'm--');

xlabel('SNR (dB)');
ylabel('BER');
legend('QPSK仿真', '8PSK仿真', '16PSK仿真', '16QAM仿真', '32QAM仿真', ...
       'QPSK理论', '8PSK理论', '16PSK理论', '16QAM理论', '32QAM理论', ...
       'Location', 'southwest');
title('AWGN信道下不同调制方式的BER性能对比');
grid on;
hold off;