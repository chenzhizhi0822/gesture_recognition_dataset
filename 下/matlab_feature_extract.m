function [featuresMatrix, V_Buffer, azimuth_Buffer, elevation_Buffer, meanBuffer] = matlab_feature_extract(Data_input, frame_idx, V_Buffer, azimuth_Buffer, elevation_Buffer, meanBuffer)
    %% 
    frame_idx = double(frame_idx + 1);
    %64*32*3
    Data = single(Data_input);%
%     n_frames = 1;%帧数
    
    fc = 60e9; % Center frequency
    % For CW, it is 0.
    %% 英飞凌手势采集数据参数
    NTS = 64; % Number of time samples per sweep采样点数
    n_chirps = 32;%chirp数
    % record_length=length(Data)/NTS*Tsweep; % length of recording in s
    Tsweep = 300e-6; % 扫描时间
    c = 3e8; % 光速
    lambda= c/fc;%波长
    N = 3; % 天线数量
    Vres = lambda/(2*Tsweep*n_chirps);

%     buffer1 = squeeze(meanBuffer(:,1,:));
%     buffer2 = squeeze(meanBuffer(:,2,:));
%     buffer3 = squeeze(meanBuffer(:,3,:));

    V = zeros(1);
    weighted_average_elevation = zeros(1);
    weighted_average_azimuth = zeros(1);
    detected_points = zeros(1);
    Range_points = zeros(1);
    doppler_azimuth_corr = zeros(1);
    doppler_elevation_corr = zeros(1);
    

    %% 预定义存储特征值的矩阵
    featuresMatrix = zeros(1, 7);
    
    persistent Data1;
    persistent Data2;
    persistent Data3;
%     Data1 = zeros(64,32);
%     Data2 = zeros(64,32);
%     Data3 = zeros(64,32);
    Data1 = Data(:,:,1);
    Data2 = Data(:,:,2);
    Data3 = Data(:,:,3);
%     
    Data1_mean = mean(Data1,2);
    Data2_mean = mean(Data2,2);
    Data3_mean = mean(Data3,2);
% 

    %短窗MTI 每个chirp减去均值
    data1_short = Data1-Data1_mean; %短窗MTI,1
    data2_short = Data2-Data2_mean; %短窗MTI,2
    data3_short = Data3-Data3_mean; %短窗MTI,3

    data1_short =data1_short/ max(max(abs(data1_short))); %%归一化,防止数值较大的特征主导模型的训练过程，提高模型的稳定性和性能。
    data2_short =data2_short/ max(max(abs(data2_short))); %%归一化
    data3_short =data3_short/ max(max(abs(data3_short))); %%归一化

    % Long window
    
%     buffer1(:,frame_idx) = Data1_mean;
%     buffer2(:,frame_idx) = Data2_mean;
%     buffer3(:,frame_idx) = Data3_mean; 

    meanBuffer(:,:,1) = ((frame_idx-1)*meanBuffer(:,:,1) + Data1)/frame_idx;
    meanBuffer(:,:,2) = ((frame_idx-1)*meanBuffer(:,:,2) + Data2)/frame_idx;
    meanBuffer(:,:,3) = ((frame_idx-1)*meanBuffer(:,:,3) + Data3)/frame_idx; %长窗，每帧更新
%  
% 
    %长窗MTI 减去15帧窗口的平均值
    if  frame_idx==1
        data1_long = Data1;
        data2_long = Data2;
        data3_long = Data3;
    else
        data1_long = Data1 - meanBuffer(:,:,1);
        data2_long = Data2 - meanBuffer(:,:,2);
        data3_long = Data3 - meanBuffer(:,:,3);
    end

    data1_long =data1_long/ max(max(abs(data1_long))); %%归一化
    data2_long =data2_long/ max(max(abs(data2_long)));
    data3_long =data3_long/ max(max(abs(data3_long)));

    Data1 = data1_short+data1_long;%加和长短窗值
    Data2 = data2_short+data2_long;
    Data3 = data3_short+data3_long;

%     Data1 = data1_short;
%     Data2 = data2_short;
%     Data3 = data3_short;

   
    %% 2DFFT
    Data1_fft = fftshift(fft2(single(Data1)), 2);
    Data2_fft = fftshift(fft2(single(Data2)), 2);
    Data3_fft = fftshift(fft2(single(Data3)), 2);
    
%     % 取单边谱
    Data1_fft = Data1_fft(1:64/2,:);
    Data2_fft = Data2_fft(1:64/2,:);
    Data3_fft = Data3_fft(1:64/2,:);

% %     clear data1_short data2_short data3_short;
% %     clear Data1 Data2 Data3 Data1_fft Data2_fft Data3_fft;
%     Data = abs(Data1_fft);
%     [flag] = start_judge(Data1_fft);
% 
%     if flag == true
    %     % 非相参积累
        Data = (abs(Data1_fft) + abs(Data2_fft) + abs(Data3_fft)) / 3;
    %     Data = abs(Data1_fft);
    %   
        
        %% 加权多普勒均值
        Zi = (sum(Data,1)).';%1x32 能量
        Di= ((-n_chirps/2):n_chirps/2-1)*Vres*100; %多普勒索引转换为真正的速度
        % Di=96:127;
        V = Di*Zi./sum(Zi);
        V_Buffer(frame_idx) = V;
     
        %% 选择热图中具有最高信号的N个像元
%         N_highest_pixels = 10;  % 选择最高信号的5个像元
  % 估计噪声水平并设置阈值
        noise_level = mean(Data(:)); %
        threshold = 6 * noise_level; % 将阈值设为噪声水平的3倍
        % 选择高于阈值的像元
        noncoherent_values = abs(Data);  % 获取非相干积累结果的绝对值
        
        % % 获取热图中信号强度最高的N个像元的索引
%         [~, sorted_indices] = sort(int16(noncoherent_values(:)), 'descend');  % 对信号强度排序
        
%         highest_indices = sorted_indices(1:N_highest_pixels);  % 选择排序后信号强度最高的N个像素的索引
        high_signal_indices =noncoherent_values > threshold;
        
        % 将一维索引转换为二维索引
        [range_bins, doppler_bins] = find(high_signal_indices);
        
        % 计算角度FFT并提取高程和方位角加权平均值
        elevation_angles = [];
        azimuth_angles = [];
        weights = [];
        
        for j = 1:length(range_bins)
            range_bin = range_bins(j);
            doppler_bin = doppler_bins(j);
        
            % 提取对应像元的数据
            pixel_data = [Data1_fft(range_bin, doppler_bin), ...
                          Data2_fft(range_bin, doppler_bin), ...
                          Data3_fft(range_bin, doppler_bin)];
        
            % 计算角度FFT
            angle_spectrum_3_1 = fftshift(fft([pixel_data(1), pixel_data(3)], N, 2), 2); % 水平对
            angle_spectrum_3_2 = fftshift(fft([pixel_data(2), pixel_data(3)], N, 2), 2); % 垂直对
        
            [~, max_index_3_1] = max(abs(angle_spectrum_3_1), [], 2); % 水平对最大值索引
            [~, max_index_3_2] = max(abs(angle_spectrum_3_2), [], 2); % 垂直对最大值索引
        
            azimuth = -90 + (max_index_3_1 - 1) * (180 / (N - 1)); % 水平角度
            elevation = -90 + (max_index_3_2 - 1) * (180 / (N - 1)); % 俯仰角度
        
            % 存储结果
            elevation_angles = [elevation_angles, elevation];
            azimuth_angles = [azimuth_angles,azimuth];
            weights = [weights; noncoherent_values(range_bin, doppler_bin)];
        end
        
         % 计算高程和方位角加权平均值
        
         weighted_average_elevation = (elevation_angles * weights) / sum(weights);
         weighted_average_azimuth = (azimuth_angles* weights) / sum(weights);

         elevation_Buffer(frame_idx) = weighted_average_elevation;
         azimuth_Buffer(frame_idx) = weighted_average_azimuth; 
        
    %      检测到的点数
         detected_points = sum(noncoherent_values(:) > threshold);
         
         % 计算范围平均值
         range_indices = (1:NTS/2).'; % 距离索引
         Range_points = sum(range_indices .* sum(Data, 2)) / sum(sum(Data));
         
         % 计算多普勒方位角相关性
         if frame_idx ==1
             doppler_azimuth_corr=0;
         else
             meanX = mean(V_Buffer);
             meanY = mean(azimuth_Buffer);
             numerator1 = sum((V_Buffer - meanX) .* (azimuth_Buffer - meanY));
             denominator1 = sqrt(sum((V_Buffer - meanX).^2) * sum((azimuth_Buffer - meanY).^2));
             doppler_azimuth_corr = numerator1 / denominator1;
%                  doppler_azimuth_corr = corr(V_Buffer.', azimuth_Buffer.');

             meanX = mean(V_Buffer);
             meanZ = mean(elevation_Buffer);
             numerator2 = sum((V_Buffer - meanX) .* (elevation_Buffer - meanZ));
             denominator2 = sqrt(sum((V_Buffer - meanX).^2) * sum((elevation_Buffer - meanZ).^2));
             doppler_elevation_corr = numerator2 / denominator2;

%                  doppler_elevation_corr = corr(V_Buffer.', elevation_Buffer.');
          
        end
            
        % 将特征存储在矩阵中，横轴为帧数纵轴为特征
        featuresMatrix = [V, weighted_average_elevation, weighted_average_azimuth, detected_points, Range_points, doppler_azimuth_corr, doppler_elevation_corr];

        
        
    end
% end


