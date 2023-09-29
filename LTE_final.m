%% LTE全系统仿真
%系统基本组成：信源——>CRC（循环冗余校验）——>信道编码——>交织——>符号调制——>插入导频——>成帧——>上变频——>过信道（EVA信道，瑞利信道）
% ——>下变频——>时间同步——>频率同步——>数据分离——>信道估计——>符号检测——>解交织——>信道译码——>解CRC——>输出

%系统条件：信号带宽 2Mhz 信号载频 2Ghz 信号频偏 1Khz 信道条件 EVA信道(多径衰落)，移动速度60km/h

%仿真要求：自行设计帧结构，调制、编码方式，进行总信源数据量不少于1e6的仿真
%分析和测试系统的时间同步捕获概率，频率估计精度，系统误码率，信息速率和频谱效率

hError = comm.ErrorRate;
%% 基本参数设置和初始化
%信源速率 
Rs = 2e6;
%传输一帧 
rand('state',123456);
msg_source = randi([0 1],25572,1);%一个帧25572个bit
M = 2;%2 bit ——> 1 symbol
Rb = Rs*M;
%% CRC
%msg_afer_crc = lte_crc(msg_source);
crcgenerator = comm.CRCGenerator('Polynomial',[1 1 zeros(1, 16) 1 1 0 0 0 1 1]);     %采用的是CRC-24B
msg_after_crc = crcgenerator(msg_source);%在后面添加24位CRC——>变成一个帧25596个bit
%% 信道编码
%msg_afer_channelcoding = lte_channelcoding(msg_afer_crc);
frmLen = 25596;%必须和msg_after_crc相同
rng default;
intrlvrIndices = randperm(frmLen);
turboenc = comm.TurboEncoder('TrellisStructure',poly2trellis(4,[13 15],13),'InterleaverIndices',intrlvrIndices);

msg_after_channelcoding = turboenc(msg_after_crc); %——>根据turbo码的编码情况3*（k+4），一个帧变为76800长度

%% 交织
%msg_after_interweave = lte_interweave(msg_after_channelcoding);
[xsize ysize] = size(msg_after_channelcoding);
depth=100;%交织深度
width=(xsize*ysize)/depth;
num=depth*width;
L =  length(msg_after_channelcoding)/num;
msg_after_interweave = [];
temp2 = [];
temp1 = reshape(msg_after_channelcoding,width,depth);
for j1 = 1:width
        temp2 = [temp2,temp1(j1,:)];
end
msg_after_interweave = temp2;
%% 符号调制
%msg_after_modulation = lte_modulation(msg_after_interweave);
QPSK=comm.PSKModulator(4,'BitInput',true,...
          'PhaseOffset',pi/4,'SymbolMapping','Gray');
msg_after_interweave = msg_after_interweave';
msg_after_modulation = QPSK(msg_after_interweave);
%constellation(QPSK)
% msg_after_modulation = msg_after_modulation.'; %——>根据QPSK的映射规则，一个帧变为38400
msg_after_modulation = msg_after_modulation.';
%% 插入导频
%msg_after_insertpilotsignal = lte_insertpilotsignal(msg_after_modulation);
N_fft = 128;   %  FFT点数
num_symbol = length(msg_after_modulation)/N_fft;
msg_compare1 = reshape(msg_after_modulation,N_fft,num_symbol);
P_f_inter = 4;      %导频间隔
rand('seed',2);
pilot_symbols = randi([0,1],M*N_fft,1);
pilot_symbols=QPSK(pilot_symbols);     %导频符号

num_pilot=floor(num_symbol/P_f_inter)+1;     %导频符号的数目
num_data=num_symbol+num_pilot;     %总的OFDM符号数目

pilot_Indx=zeros(1,num_pilot);
pilot_Indx=1:P_f_inter+1:num_data;     %导频的位置
Data_Indx=[];
for i2=1:num_data
    if mod(i2,P_f_inter+1)~=1
        Data_Indx=[Data_Indx,i2];
    end
end
msg_after_insertpilotsignal=zeros(N_fft,num_data);     %预设整个矩阵
msg_after_insertpilotsignal(:,Data_Indx)=reshape(msg_after_modulation,N_fft,num_symbol);     %插入数据
msg_after_insertpilotsignal(:,pilot_Indx)=repmat(pilot_symbols,1,num_pilot);     %插入导频

%% 成帧
%msg_framing = lte_framing(msg_after_insertpilotsignal);
% ifft
time_signal = ifft(msg_after_insertpilotsignal);
% 加循环前缀
N_cp = 12;
cyclic_signal = [time_signal(N_fft-N_cp+1:end,:);time_signal];%加入循环前缀
% 串并转换
msg_framing1 = reshape(cyclic_signal,1,[]);     %并串转换，待发送信号  插入导频和循环前缀后后一个帧由38400个符号变为47320个符号
msg_framing = normalize(msg_framing1);       %做了IFFT之后需要功率归一化——>保证后面加噪声是正确的

% 加上时间同步头和频率同步头
%生成pn码1和pn码2
load('pn1.mat');
load('pn2.mat');
synchronization = [zeros(1,20),pn1,zeros(1,20),zeros(1,20),pn2,pn2,zeros(1,20)];
%synchronization = [zeros(1,20),pn1,zeros(1,20)];
synchronization = QPSK(synchronization');
synchronization = synchronization.';
lengthsyn = length(synchronization);
msg_framing = [synchronization,msg_framing];%加入同步头后变为47744个符号 
msg_framing = msg_framing.';

%% 加频偏
fd = 1000; % 频偏
w = 2*pi*fd/(Rb*256*2*log2(M));
for i = 1:length(msg_framing)
    msg_framing(i) = msg_framing(i)*exp(1i*w*i);
end
%% 过信道
%msg_afterchannel = lte_afterchannel(msg_upfrequency);
%EVA信道
EVArayleighchan = comm.RayleighChannel( ...
    'SampleRate',          2e6,...
    'PathDelays',          [0 30 150 310 370 710 1090 1730 2510]./1e9, ...
    'AveragePathGains',    [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9], ...
    'NormalizePathGains',  true, ...
    'MaximumDopplerShift', 110, ...            %最大多普勒频移 =v*f/c
    'RandomStream',        'mt19937ar with seed', ...
    'Seed',                22, ...
    'PathGainsOutputPort', true);
rng(22,'twister');
% 高斯信道
EbN0 = 0;
SNR = EbN0 + 10*log10(M);
msg_afterchannel = EVArayleighchan(msg_framing);
msg_afterchannel = awgn(msg_afterchannel,SNR,'measured');
msg_downfrequency = msg_afterchannel;

 %% 时间同步  
% %msg_timesynchronization = lte_timesynchronization(msg_downfrequency);
msg_downfrequency = msg_downfrequency.';
localpn1 = QPSK(pn1');
localpn1 = localpn1.';
lengtht = length(msg_downfrequency);
beitax = zeros(1,lengtht);
beita = zeros(1,lengtht);
gama = zeros(1,lengtht);
Mx = zeros(1,lengtht); 
%计算门限
Tcc = 10;
for k1 = 1:lengtht-256/M
    for m = 1:256/M
        beitax(k1) = beitax(k1) + abs(msg_downfrequency(m+k1-1)*conj(msg_downfrequency(m+k1-1)));
    end
    beita(k1) = Tcc*beitax(k1);
end

%同步判断
for k2 = 1 : lengtht-256/M
    for m = 1 : 256/M
        gama(k2) = gama(k2) + conj(localpn1(m))*msg_downfrequency(m+k2-1);%%由于MATLAB数组的原因，此处k表示延时+1，即k=1表示无延时（见doc文档）
    end
   Mx(k2) = gama(k2)*conj(gama(k2));
   if Mx(k2) > beita(k2)
     timestartpoint = k2-7;
       break;
   end
end

%% 频率同步
PN_num_freq=2;  %pn2序列的重复个数
PN_num=128; %PN2序列调制后的位数
msg_downfrequency_conj = conj(msg_downfrequency);
p = zeros(1,(PN_num_freq-1)*PN_num);
for i = 1:(PN_num_freq-1)*PN_num
    p(i) = msg_downfrequency_conj(i+158)*msg_downfrequency(i +158+ PN_num);
end
P = sum(p);
w_es = abs(angle(P)/PN_num);
f_es = w_es*Rb*PN_num/(2*pi);
% 频偏校正
msg_correct_fre = zeros(1,length(msg_downfrequency));
for i = 1:length(msg_downfrequency)
   msg_correct_fre(i) = msg_downfrequency(i)*exp(-1i*w_es);
end
%% 数据分离
%msg_dataspliting = lte_dataspliting(msg_fresynchronization);
msg_dataspliting1 = msg_downfrequency(425:end);%%去除掉两个同步头，只剩下数据部分47320个符号
%% 信道估计
%msg_channelestimate = lte_channelestimate(msg_dataspliting);
%在此之前先去掉循环前缀，做好fft变换，提取导频和数据
%msg_channelestimate_interpolation = zeros(128,num_symbol);
%  去循环前缀
msg_dataspliting2 = reshape(msg_dataspliting1,(N_fft+N_cp),num_data);     %串并转换
msg_dataspliting3 = msg_dataspliting2(N_cp+1:end,:);     %去前缀

% FFT
msg_dataspliting4 = fft(msg_dataspliting3,128);
%msg_dataspliting4 = conj(msg_dataspliting4);
%msg_dataspliting4 = normalize(msg_dataspliting4);

% 导频和数据提取
msg_dataspliting_pliot = msg_dataspliting4(:,pilot_Indx);
msg_dataspliting_data = msg_dataspliting4(:,Data_Indx);

%信道估计
pilot_patt = repmat(pilot_symbols,1,num_pilot);
msg_channelestimate1 = msg_dataspliting_pliot./pilot_patt;

%线性插值
%msg_channelestimate_interpolation = [];
for ii=1:N_fft
        msg_channelestimate_interpolation(ii,:) = interp1(pilot_Indx,msg_channelestimate1(ii,1:(num_pilot)),Data_Indx,'spline');
end

%均衡（破零均衡 和 MMSE均衡）   破零均衡
est = conj(msg_channelestimate_interpolation)./(abs(msg_channelestimate_interpolation).^2);
msg_channelestimate=msg_dataspliting_data.*est;
msg_compare2 = msg_channelestimate;

%均衡前后对比
scatterplot(msg_dataspliting1);
msg_plot1 = reshape(msg_channelestimate,1,[]);
scatterplot(msg_plot1);

%msg_channelestimate = msg_dataspliting_data;
%% 符号检测
%msg_symboldetection = lte_symboldetection(msg_channelestimate); 
%QPSK解调
noisevar = 10^(-EbN0/10);
msg_symboldetectionbefore = reshape(msg_channelestimate,[],1);     %并串转换
DeQPSK = comm.PSKDemodulator(...
        'ModulationOrder',4,...
        'BitOutput', true, ...
        'PhaseOffset', pi/4, 'SymbolMapping','Gray',...
        'DecisionMethod','Approximate log-likelihood ratio',...
        'Variance',noisevar);
%constellation(DeQPSK)
msg_symboldetection =DeQPSK(msg_symboldetectionbefore);
%errorStats = hError(msg_symboldetection,msg_after_interweave);
%% 解交织
depth=100;
width = (xsize*ysize)/depth;
num = depth*width;
L =  length(msg_symboldetection)/num;
msg_deinterweave = [];
temp4 = [];
temp3 = reshape(msg_symboldetection,depth,width);
for j2 = 1:depth
     temp4 = [temp4,temp3(j2,:)];
end
msg_deinterweave  = temp4;
msg_deinterweave = msg_deinterweave';
%errorStats = hError(msg_deinterweave,msg_after_channelcoding)
%% 信道译码
%msg_channeldecoding = lte_channeldecoding(msg_deinterweave);
turbodec = comm.TurboDecoder('TrellisStructure', poly2trellis(4, [13 15], 13),'InterleaverIndices',intrlvrIndices);
msg_channeldecoding = turbodec(-msg_deinterweave);
%errorStats = hError(msg_channeldecoding,msg_after_crc);
%% 解CRC
%msg_decrc = lte_decrc(msg_channeldecoding);
crcdetector = comm.CRCDetector('Polynomial',[1 1 zeros(1, 16) 1 1 0 0 0 1 1]);     %采用的是CRC-24B
msg_decrc = crcdetector(msg_channeldecoding);

errorStats = hError(msg_source,msg_decrc);
fprintf('Bit error rate = %5.2e\nNumber of errors = %d\nTotal bits = %d\n', ...
    errorStats)




   
