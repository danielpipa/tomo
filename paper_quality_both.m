clear; close all;

%Parei no Size_i =6, trial_i = 39, salvo no temp.mat

photons = 50;
readout_noise = 1;
noise_ph = (pi^2/(2*log(2)))*(1/photons);
noise_rd = (pi^2/3)*(readout_noise^2/photons^2)*(2^4/2^2);
sigma_2 = noise_ph+noise_rd;
for size_i = 2:8
    fprintf('###########################################\nBeginning size_i: %d\n\n',size_i);
    file_name = sprintf("paperInit_%.2d",size_i);
	load(file_name);
    
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
        
        HtH = Hmtx'*Hmtx + sigma_2*Cmtx_;
        Rmat = inv(HtH)*Hmtx';


    for trial_i = 1:100
        fprintf('size_i: %d, trial_i: %d stop_i: %d\n',size_i,trial_i,0);
        
        % Add some noise
        WFSX = WFSX + sqrt(sigma_2)*randn(size(WFSX));
        WFSY = WFSY + sqrt(sigma_2)*randn(size(WFSY));
        
        %Classic MVM
        layersX_hat = Rmat*WFSX; %Do one time to get function on memory
        tic;
        %layersX_hat = pcg(HtH,bb,1e-6,1000);
        layersX_hat = Rmat*WFSX;
        %layersY_hat = pcg(HtH,bb,1e-6,1000);
        layersY_hat = Rmat*WFSY;
        Time_MVM(trial_i) = toc;
        
        %% **************************************************
        % SLOPE ERROR RESULTS
        %****************************************************

        TargetX = HH_projection3(layersX,x_shift,y_shift,sizes,alt_p,WFS_target);
        TargetY = HH_projection3(layersY,x_shift,y_shift,sizes,alt_p,WFS_target);
        
        %MVM
        TargetX_hat = HH_projection3(layersX_hat,x_shift,y_shift,sizes,alt_p,WFS_target);
        Error_TargetX_MVM(trial_i) = rms(TargetX-TargetX_hat);
        Error_layersX_MVM(trial_i) = rms(layersX-layersX_hat);
        TargetY_hat = HH_projection3(layersY_hat,x_shift,y_shift,sizes,alt_p,WFS_target);
        Error_TargetY_MVM(trial_i) = rms(TargetY-TargetY_hat);
        Error_layersY_MVM(trial_i) = rms(layersY-layersY_hat);


        %% **************************************************
        % WFE RESULTS
        %****************************************************
        %Based on the expected tilt removed WFE -> Eq. 9.11, pag 316 Hardy
        expected_WFE = 0.141*(tel_diam/r0)^(5/3)*500/(2*pi);    %In Nanometers
        d = tel_diam/WFS_size;
        rec_target = WFreconWFS(TargetX,TargetY,sizes,d,0.005);
        scaling_value = expected_WFE/rms(rec_target(:));
        rec_target = rec_target*scaling_value;

        [temp1, temp2, temp3] = WFrecon(layersX,layersY,sizes,d,0.005);
        rec_layers = [temp1; temp2; temp3]*scaling_value;

        
        %MVM
        [temp1, temp2, temp3] = WFrecon(layersX_hat,layersY_hat,sizes,d,0.005);
        rec_layers_hat = [temp1; temp2; temp3]*scaling_value;
        WFE_layers_MVM(trial_i) = rms(rec_layers(:) - rec_layers_hat(:));
        rec_target_hat = WFreconWFS(TargetX_hat,TargetY_hat,sizes,d,0.005)*scaling_value;
        WFE_Target_MVM(trial_i) = rms(rec_target(:) - rec_target_hat(:));
        
        for stop_i = 1:20
            fprintf('size_i: %d, trial_i: %d stop_i: %d\n',size_i,trial_i,stop_i);
            stop_itera = stop_i*5;

            %Filter without inverting Cov
            AA = @(x) HHmex3(Cparmex3(HHTmex3(x,x_shift,y_shift,sizes),sizes,frac,CovKernel),x_shift,y_shift,sizes) + sigma_2*x;
            [temp,~,~,~,~] = pcg(AA,WFSX,1e-15,stop_itera);%Do one time to get function on memory
            layersX_hat3 = Cparmex3(HHTmex3(temp,x_shift,y_shift,sizes),sizes,frac,CovKernel);
            tic;
            [temp,~,~,ITER_X(stop_i,trial_i),~] = pcg(AA,WFSX,1e-15,stop_itera);
            layersX_hat3 = Cparmex3(HHTmex3(temp,x_shift,y_shift,sizes),sizes,frac,CovKernel);
            [temp,~,~,ITER_Y(stop_i,trial_i),~] = pcg(AA,WFSY,1e-15,stop_itera);
            layersY_hat3 = Cparmex3(HHTmex3(temp,x_shift,y_shift,sizes),sizes,frac,CovKernel);
            Time_Fil(stop_i,trial_i) = toc;


            %% **************************************************
            % SLOPE ERROR RESULTS
            %****************************************************

            %Filter without inverting C
            TargetX_hat3 = HH_projection3(layersX_hat3,x_shift,y_shift,sizes,alt_p,WFS_target);
            Error_TargetX_FilterNI(stop_i,trial_i) = rms(TargetX-TargetX_hat3);
            Error_layersX_FilterNI(stop_i,trial_i) = rms(layersX-layersX_hat3);
            TargetY_hat3 = HH_projection3(layersY_hat3,x_shift,y_shift,sizes,alt_p,WFS_target);
            Error_TargetY_FilterNI(stop_i,trial_i) = rms(TargetY-TargetY_hat3);
            Error_layersY_FilterNI(stop_i,trial_i) = rms(layersY-layersY_hat3);


            %Filter without inverting C
            [temp1, temp2, temp3] = WFrecon(layersX_hat3,layersY_hat3,sizes,d,0.005);
            rec_layers_hat3 = [temp1; temp2; temp3]*scaling_value;
            WFE_layers_FilterNI(stop_i,trial_i) = rms(rec_layers(:) - rec_layers_hat3(:));
            rec_target_hat3 = WFreconWFS(TargetX_hat3,TargetY_hat3,sizes,d,0.005)*scaling_value;
            WFE_Target_FilterNI(stop_i,trial_i) = rms(rec_target(:) - rec_target_hat3(:));
        end

    end
    %%

    file_name = sprintf("paper_quality_noise%.3f_size%.2d.mat",sigma_2,size_i);
    save(file_name);
end