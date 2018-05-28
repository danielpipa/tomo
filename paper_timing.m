clear; close all;

photons = 50;
readout_noise = 1;
noise_ph = (pi^2/(2*log(2)))*(1/photons);
noise_rd = (pi^2/3)*(readout_noise^2/photons^2)*(2^4/2^2);
sigma_2 = noise_ph+noise_rd;

for size_i = 9:10
    fprintf('###########################################\nBeginning size_i: %d\n\n',size_i);
    file_name = sprintf("paperInit_%.2d",size_i);
    load(file_name);

    %%
    for i = 1:length(sizes)
        tempX{i} = randn(1,sizes(i)^2)*L{i};
        tempY{i} = randn(1,sizes(i)^2)*L{i};
    end
    layersX = [tempX{:}]';     %X layers with random slopes
    layersY = [tempY{:}]';     %Y layers with random slopes




    %%**************************************************
    % simulate WFS
    %**************************************************
    [WFSX] = HHmex3(layersX,x_shift,y_shift,sizes);
    [WFSY] = HHmex3(layersY,x_shift,y_shift,sizes);


    %%**************************************************
    % RECONSTRUCTION OF THE LAYERS
    %****************************************************

    % Add some noise
    WFSX = WFSX + sqrt(sigma_2)*randn(size(WFSX));
    WFSY = WFSY + sqrt(sigma_2)*randn(size(WFSY));

    %%

    %Matrix
    % fprintf('\nReconstructing. Matrix\n');
    % HtH = Hmtx'*Hmtx + sigma_2*Cmtx_;
    % tic;
    % bb = Hmtx'*WFSX;
    % layersX_hat = pcg(HtH,bb,1e-6,1000);
    % bb = Hmtx'*WFSY;
    % layersY_hat = pcg(HtH,bb,1e-6,1000);
    % Time_Matrix = toc;
    % fprintf('Elapsed Time = %f\n',Time_Matrix);

    %Classic MVM
    fprintf('\nReconstructing. MVM\n');
    HtH = Hmtx'*Hmtx + sigma_2*Cmtx_;
    Rmat = inv(HtH)*Hmtx';
    layersX_hat = Rmat*WFSX; %Do one time to get function on memory
    tic;
    %layersX_hat = pcg(HtH,bb,1e-6,1000);
    layersX_hat = Rmat*WFSX;
    %layersY_hat = pcg(HtH,bb,1e-6,1000);
    layersY_hat = Rmat*WFSY;
    Time_MVM = toc;
    fprintf('Elapsed Time = %f\n',Time_MVM);

    %Matrix without inverting Cov
    % fprintf('\nReconstructing. Matrix without inverting Cov\n');
    % A = Hmtx*Cmtx*Hmtx';
    % A = A + sigma_2*eye(size(A));
    % tic;
    % layersX_hat2 = Cmtx*Hmtx'*pcg(A,WFSX,1e-6,20);
    % layersY_hat2 = Cmtx*Hmtx'*pcg(A,WFSY,1e-6,20);
    % Time_Matrix = toc;
    % fprintf('Elapsed Time = %f\n',Time_Matrix);

    %Filter without inverting Cov
    fprintf('\nReconstructing. Filter without inverting Cov\n');
    AA = @(x) HHmex3(Cparmex3(HHTmex3(x,x_shift,y_shift,sizes),sizes,frac,CovKernel),x_shift,y_shift,sizes) + sigma_2*x;
    temp = pcg(AA,WFSX,1e-4,200);%Do one time to get function on memory
    layersX_hat3 = Cparmex3(HHTmex3(temp,x_shift,y_shift,sizes),sizes,frac,CovKernel);
    tic;
    temp = pcg(AA,WFSX,1e-4,200);
    layersX_hat3 = Cparmex3(HHTmex3(temp,x_shift,y_shift,sizes),sizes,frac,CovKernel);
    temp = pcg(AA,WFSY,1e-4,200);
    layersY_hat3 = Cparmex3(HHTmex3(temp,x_shift,y_shift,sizes),sizes,frac,CovKernel);
    Time_Fil = toc;
    fprintf('Elapsed Time = %f\n',Time_Fil);


    %% **************************************************
    % SLOPE ERROR RESULTS
    %****************************************************

    TargetX = HH_projection3(layersX,x_shift,y_shift,sizes,alt_p,WFS_target);
    TargetY = HH_projection3(layersY,x_shift,y_shift,sizes,alt_p,WFS_target);

    %MVM
    TargetX_hat = HH_projection3(layersX_hat,x_shift,y_shift,sizes,alt_p,WFS_target);
    Error_TargetX_MVM = rms(TargetX-TargetX_hat)
    Error_layersX_MVM = rms(layersX-layersX_hat)
    TargetY_hat = HH_projection3(layersY_hat,x_shift,y_shift,sizes,alt_p,WFS_target);
    Error_TargetY_MVM = rms(TargetY-TargetY_hat);
    Error_layersY_MVM = rms(layersY-layersY_hat);

    %Matrix
%     TargetX_hat2 = HH_projection3(layersX_hat2,x_shift,y_shift,sizes,alt_p,WFS_target);
%     Error_TargetX_Mat = rms(TargetX-TargetX_hat2)
%     Error_layersX_Mat = rms(layersX-layersX_hat2)
%     TargetY_hat2 = HH_projection3(layersY_hat2,x_shift,y_shift,sizes,alt_p,WFS_target);
%     Error_TargetY_Mat = rms(TargetY-TargetY_hat2);
%     Error_layersY_Mat = rms(layersY-layersY_hat2);

    %Filter without inverting C
    TargetX_hat3 = HH_projection3(layersX_hat3,x_shift,y_shift,sizes,alt_p,WFS_target);
    Error_TargetX_FilterNI = rms(TargetX-TargetX_hat3)
    Error_layersX_FilterNI = rms(layersX-layersX_hat3)
    TargetY_hat3 = HH_projection3(layersY_hat3,x_shift,y_shift,sizes,alt_p,WFS_target);
    Error_TargetY_FilterNI = rms(TargetY-TargetY_hat3);
    Error_layersY_FilterNI = rms(layersY-layersY_hat3);


    %% **************************************************
    % WFE RESULTS
    %****************************************************
    %Based on the expected tilt removed WFE -> Eq. 9.11, pag 316 Hardy
%     expected_WFE = 0.141*(tel_diam/r0)^(5/3)*500/(2*pi);    %In Nanometers
%     d = tel_diam/WFS_size;
%     rec_target = WFreconWFS(TargetX,TargetY,sizes,d,0.005);
%     scaling_value = expected_WFE/rms(rec_target(:));
%     rec_target = rec_target*scaling_value;
% 
%     [temp1, temp2, temp3] = WFrecon(layersX,layersY,sizes,d,0.005);
%     rec_layers = [temp1; temp2; temp3]*scaling_value;
% 
%     %MVM
%     [temp1, temp2, temp3] = WFrecon(layersX_hat,layersY_hat,sizes,d,0.005);
%     rec_layers_hat = [temp1; temp2; temp3]*scaling_value;
%     WFE_layers_MVM = rms(rec_layers(:) - rec_layers_hat(:))
%     rec_target_hat = WFreconWFS(TargetX_hat,TargetY_hat,sizes,d,0.005)*scaling_value;
%     WFE_Target_MVM = rms(rec_target(:) - rec_target_hat(:))
% 
%     %Matrix without inverting Cov
% %     [temp1, temp2, temp3] = WFrecon(layersX_hat2,layersY_hat2,sizes,d,0.005);
% %     rec_layers_hat2 = [temp1; temp2; temp3]*scaling_value;
% %     WFE_layers_Mat = rms(rec_layers(:) - rec_layers_hat2(:))
% %     rec_target_hat2 = WFreconWFS(TargetX_hat2,TargetY_hat2,sizes,d,0.005)*scaling_value;
% %     WFE_Target_Mat = rms(rec_target(:) - rec_target_hat2(:))
% 
%     %Filter without inverting C
%     [temp1, temp2, temp3] = WFrecon(layersX_hat3,layersY_hat3,sizes,d,0.005);
%     rec_layers_hat3 = [temp1; temp2; temp3]*scaling_value;
%     WFE_layers_FilterNI = rms(rec_layers(:) - rec_layers_hat3(:))
%     rec_target_hat3 = WFreconWFS(TargetX_hat3,TargetY_hat3,sizes,d,0.005)*scaling_value;
%     WFE_Target_FilterNI = rms(rec_target(:) - rec_target_hat3(:))



    %%

    file_name = sprintf("paper_timing_%.2d",size_i);
    save(file_name);
end