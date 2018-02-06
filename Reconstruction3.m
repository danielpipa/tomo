clear; close all;

load Reconstruction3_InitData

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
sigma_2 = var(WFSX)/1;
%sigma_2 = 0;
WFSX = WFSX + sqrt(sigma_2)*randn(size(WFSX));
WFSY = WFSY + sqrt(sigma_2)*randn(size(WFSY));

%%

%%C full Covariance
fprintf('\nReconstructing. C Full Covariance\n');
AA = @(x) HHTmex3(HHmex3(x,x_shift,y_shift,sizes),x_shift,y_shift,sizes) + sigma_2*Cmtx_*x;
tic;
bb = HHTmex3(WFSX,x_shift,y_shift,sizes);
layersX_hat = pcg(AA,bb,1e-6,1000);
bb = HHTmex3(WFSY,x_shift,y_shift,sizes);
layersY_hat = pcg(AA,bb,1e-6,1000);
Time_C = toc;
fprintf('Elapsed Time = %f\n',Time_C);

%C filter Covariance
fprintf('\nReconstructing. C Filter Covariance\n');
AA = @(x) HHTmex3(HHmex3(x,x_shift,y_shift,sizes),x_shift,y_shift,sizes) + sigma_2*Cmex3(x,sizes);
tic;
bb = HHTmex3(WFSX,x_shift,y_shift,sizes);
layersX_hat2 = pcg(AA,bb,1e-6,1000);
bb = HHTmex3(WFSY,x_shift,y_shift,sizes);
layersY_hat2 = pcg(AA,bb,1e-6,1000);
Time_C2 = toc;
fprintf('Elapsed Time = %f\n',Time_C2);

%Matrix
fprintf('\nReconstructing. Matrix\n');
HtH = Hmtx'*Hmtx + sigma_2*Cmtx_;
tic;
bb = Hmtx'*WFSX;
layersX_hat3 = pcg(HtH,bb,1e-6,1000);
bb = Hmtx'*WFSY;
layersY_hat3 = pcg(HtH,bb,1e-6,1000);
Time_Matrix = toc;
fprintf('Elapsed Time = %f\n',Time_Matrix);

%Matrix without inverting Cov
fprintf('\nReconstructing. Matrix without inverting Cov\n');
A = Hmtx*Cmtx*Hmtx';
A = A + sigma_2*eye(size(A));
tic;
b = WFSX;
layersX_hat4 = Cmtx*Hmtx'*pcg(A,b,1e-6,1000);
b = WFSY;
layersY_hat4 = Cmtx*Hmtx'*pcg(A,b,1e-6,1000);
Time_Matrix2 = toc;
fprintf('Elapsed Time = %f\n',Time_Matrix2);

%Filter without inverting Cov
fprintf('\nReconstructing. Filter without inverting Cov\n');
AA = @(x) HHmex3(Cparmex3(HHTmex3(x,x_shift,y_shift,sizes),sizes,CovKernel),x_shift,y_shift,sizes) + sigma_2*x;
tic;
bb = WFSX;
temp = pcg(AA,bb,1e-6,1000);
layersX_hat5 = Cparmex3(HHTmex3(temp,x_shift,y_shift,sizes),sizes,CovKernel);
bb = WFSY;
temp = pcg(AA,bb,1e-6,1000);
layersY_hat5 = Cparmex3(HHTmex3(temp,x_shift,y_shift,sizes),sizes,CovKernel);
Time_Fil2 = toc;
fprintf('Elapsed Time = %f\n',Time_Fil2);


%% **************************************************
% SLOPE ERROR RESULTS
%****************************************************

TargetX = HH_projection3(layersX,x_shift,y_shift,sizes,alt_p,WFS_target);
TargetY = HH_projection3(layersY,x_shift,y_shift,sizes,alt_p,WFS_target);

%C full Covariance
TargetX_hat = HH_projection3(layersX_hat,x_shift,y_shift,sizes,alt_p,WFS_target);
Error_TargetX_Full = rms(TargetX-TargetX_hat)
Error_layersX_Full = rms(layersX-layersX_hat)
TargetY_hat = HH_projection3(layersY_hat,x_shift,y_shift,sizes,alt_p,WFS_target);
Error_TargetY_Full = rms(TargetY-TargetY_hat);
Error_layersY_Full = rms(layersY-layersY_hat);

%C filter Covariance
TargetX_hat2 = HH_projection3(layersX_hat2,x_shift,y_shift,sizes,alt_p,WFS_target);
Error_TargetX_Filter = rms(TargetX-TargetX_hat2)
Error_layersX_Filter = rms(layersX-layersX_hat2)
TargetY_hat2 = HH_projection3(layersY_hat2,x_shift,y_shift,sizes,alt_p,WFS_target);
Error_TargetY_Filter = rms(TargetY-TargetY_hat2);
Error_layersY_Filter = rms(layersY-layersY_hat2);

%Matrix
TargetX_hat3 = HH_projection3(layersX_hat3,x_shift,y_shift,sizes,alt_p,WFS_target);
Error_TargetX_Mat = rms(TargetX-TargetX_hat3)
Error_layersX_Mat = rms(layersX-layersX_hat3)
TargetY_hat3 = HH_projection3(layersY_hat3,x_shift,y_shift,sizes,alt_p,WFS_target);
Error_TargetY_Mat = rms(TargetY-TargetY_hat3);
Error_layersY_Mat = rms(layersY-layersY_hat3);

%Matrix without inverting C
TargetX_hat4 = HH_projection3(layersX_hat4,x_shift,y_shift,sizes,alt_p,WFS_target);
Error_TargetX_MatNI = rms(TargetX-TargetX_hat4)
Error_layersX_MatNI = rms(layersX-layersX_hat4)
TargetY_hat4 = HH_projection3(layersY_hat4,x_shift,y_shift,sizes,alt_p,WFS_target);
Error_TargetY_MatNI = rms(TargetY-TargetY_hat4);
Error_layersY_MatNI = rms(layersY-layersY_hat4);

%Filter without inverting C
TargetX_hat5 = HH_projection3(layersX_hat5,x_shift,y_shift,sizes,alt_p,WFS_target);
Error_TargetX_FilterNI = rms(TargetX-TargetX_hat5)
Error_layersX_FilterNI = rms(layersX-layersX_hat5)
TargetY_hat5 = HH_projection3(layersY_hat5,x_shift,y_shift,sizes,alt_p,WFS_target);
Error_TargetY_FilterNI = rms(TargetY-TargetY_hat5);
Error_layersY_FilterNI = rms(layersY-layersY_hat5);


%% **************************************************
% WFE RESULTS
%****************************************************
%Based on the expected tilt removed WFE -> Eq. 9.11, pag 316 Hardy
expected_WFE = 0.141*(tel_diam/r0)^(5/3)*500/(2*pi);    %In Nanometers
d = tel_diam/WFS_size;
rec_target = WFrecon(TargetX,TargetY,sizes,d,0.005);
scaling_value = expected_WFE/rms(rec_target(:));
rec_target = rec_target*scaling_value;

rec_layers = WFrecon(layersX,layersY,sizes,d,0.005)*scaling_value;

%C full Covariance
rec_layers_hat = WFrecon(layersX_hat,layersY_hat,sizes,d,0.005)*scaling_value;
WFE_layers_Full = rms(rec_layers(:) - rec_layers_hat(:))
rec_target_hat = WFrecon(TargetX_hat,TargetY_hat,sizes,d,0.005)*scaling_value;
WFE_Target_Full = rms(rec_target(:) - rec_target_hat(:))

%C filter Covariance
rec_layers_hat2 = WFrecon(layersX_hat2,layersY_hat2,sizes,d,0.005)*scaling_value;
WFE_layers_Filter = rms(rec_layers(:) - rec_layers_hat2(:))
rec_target_hat2 = WFrecon(TargetX_hat2,TargetY_hat2,sizes,d,0.005)*scaling_value;
WFE_Target_Filter = rms(rec_target(:) - rec_target_hat2(:))

%Matrix
rec_layers_hat3 = WFrecon(layersX_hat3,layersY_hat3,sizes,d,0.005)*scaling_value;
WFE_layers_Mat = rms(rec_layers(:) - rec_layers_hat3(:))
rec_target_hat3 = WFrecon(TargetX_hat3,TargetY_hat3,sizes,d,0.005)*scaling_value;
WFE_Target_Mat = rms(rec_target(:) - rec_target_hat3(:))

%Matrix without inverting C
rec_layers_hat4 = WFrecon(layersX_hat4,layersY_hat4,sizes,d,0.005)*scaling_value;
WFE_layers_MatNI = rms(rec_layers(:) - rec_layers_hat4(:))
rec_target_hat4 = WFrecon(TargetX_hat4,TargetY_hat4,sizes,d,0.005)*scaling_value;
WFE_Target_MatNi = rms(rec_target(:) - rec_target_hat4(:))

%Filter without inverting C
rec_layers_hat5 = WFrecon(layersX_hat5,layersY_hat5,sizes,d,0.005)*scaling_value;
WFE_layers_FilterNI = rms(rec_layers(:) - rec_layers_hat5(:))
rec_target_hat5 = WFrecon(TargetX_hat5,TargetY_hat5,sizes,d,0.005)*scaling_value;
WFE_Target_FilterNi = rms(rec_target(:) - rec_target_hat5(:))


%%
%
%**************************************************
% VISUAL COMPARISON
%**************************************************
Target_vis = zeros(sizes(1));
Target_vis(subap_index) = TargetX;
figure;
imagesc(real(Target_vis));
title('Target original');
Target_hat_vis = zeros(sizes(1));
Target_hat_vis(subap_index) = TargetX_hat5;
figure;
imagesc(real(Target_hat_vis));
title('Target reconstructed');


%layer 1
figure;
imagesc(real(reshape(layersX(1:sizes(1)^2),sizes(1),[])))
title('layer 1 original');
figure;
imagesc(real(reshape(layersX_hat5(1:sizes(1)^2),sizes(1),[])))
title('layer 1 rec');


%layer 2
figure;
imagesc(real(reshape(layersX(sizes(1)^2+1:sizes(1)^2+sizes(2)^2),sizes(2),[])))
title('layer 2 original');
figure;
imagesc(real(reshape(layersX_hat5(sizes(1)^2+1:sizes(1)^2+sizes(2)^2),sizes(2),[])))
title('layer 2 rec');

%%

save Reconstruction3_Results