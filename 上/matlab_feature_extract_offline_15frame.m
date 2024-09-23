features = zeros(300, 106);
for j = 1:300
    meanBuffer = single(zeros(64,32,3));
    V_Buffer = zeros(1,15);
    azimuth_Buffer = zeros(1,15);
    elevation_Buffer = zeros(1,15);
    frame_idx = (0);

    start_name = ['Data5','_'];
    read_name = [start_name,num2str(j),'.mat'];
    save_data = coder.load(read_name,'save_data');
    Data_input = uint16(save_data.save_data);
    
    adcDataRowRx1 = Data_input(:,[1601:2080],1);
    adcDataRowRx2 = Data_input(:,[1601:2080],2);
    adcDataRowRx3 = Data_input(:,[1601:2080],3);
    
    adcDataRowRx1 = reshape(adcDataRowRx1,[64 32 15]); %%
    adcDataRowRx2 = reshape(adcDataRowRx2,[64 32 15]);
    adcDataRowRx3 = reshape(adcDataRowRx3,[64 32 15]);
    
    newMatrix = uint16(zeros(64,32,3));
    
    featuresMatrix_result = zeros(1,105);
    % featuresMatrix_result2 = zeros(1,90);
    
    % featuresMatrix_result2 = frame_100_process(Data_input);
    
    for i = 1 : 15
        newMatrix(:,:,1) = adcDataRowRx1(:,:,i);
        newMatrix(:,:,2) = adcDataRowRx2(:,:,i);
        newMatrix(:,:,3) = adcDataRowRx3(:,:,i);
        [featuresMatrix, V_Buffer, azimuth_Buffer, elevation_Buffer, meanBuffer] = matlab_feature_extract(newMatrix, frame_idx, V_Buffer, azimuth_Buffer, elevation_Buffer, meanBuffer);
    
        featuresMatrix_result(i) = featuresMatrix(1);
        featuresMatrix_result(i+15) = featuresMatrix(2);
        featuresMatrix_result(i+30) = featuresMatrix(3);
        featuresMatrix_result(i+45) = featuresMatrix(4);
        featuresMatrix_result(i+60) = featuresMatrix(5);
        featuresMatrix_result(i+75) = featuresMatrix(6);
        featuresMatrix_result(i+90) = featuresMatrix(7);
    
        frame_idx = frame_idx + 1;
    end
    
   % 在向量末尾添加一个值为1的标签
    featuresMatrix_result(end + 1) = 5;
    features(j,:) = featuresMatrix_result;
end


save_name = ['train5.mat'];
save(save_name, 'features','-mat');
