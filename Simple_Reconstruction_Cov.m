clear; close all;

load Simple_Reconstruction_Cov_InitData

%%
%layersXY = [randn(1,sizes(1)^2)*L1+1i*(randn(1,sizes(1)^2)*L1) randn(1,sizes(2)^2)*L2+1i*(randn(1,sizes(2)^2)*L2)]';     %Two layers with random slopes
layersXY = [randn(1,sizes(1)^2)*L1 randn(1,sizes(2)^2)*L2]';     %Two layers with random slopes

%layersXY = layersXY/sqrt(var(layersXY));  %Set variance to one ##There
%must be a better way##: Using now Cov_M1_norm

%clear Cov_M1 Cov_M2 L1 L2


%%**************************************************
% simulate WFS
%**************************************************
[WFS] = HHmex(layersXY,x_shift,y_shift,sizes);


%%**************************************************
% RECONSTRUCTION OF THE LAYERS
%****************************************************



% Add some noise
sigma_2 = var(WFS)/10;
WFS_noisy = WFS + sqrt(sigma_2)*randn(size(WFS));

%%
% %MATLAB
% AA = @(x) HHT(HH(x,x_shift,y_shift,sizes),x_shift,y_shift,sizes) + sigma_2*C(x,sizes);
% tic;
% bb = HHT(WFS_noisy,x_shift,y_shift,sizes);
% layersXY_hat1 = pcg(AA,bb,1e-15,1000);
% Time_Mat = toc

%%C full Covariance
AA = @(x) HHTmex(HHmex(x,x_shift,y_shift,sizes),x_shift,y_shift,sizes) + sigma_2*Cmtx_*x;
tic;
bb = HHTmex(WFS,x_shift,y_shift,sizes);
layersXY_hat2 = pcg(AA,bb,1e-6,1000);
Time_C2 = toc;

%C filter Covariance
AA = @(x) HHTmex(HHmex(x,x_shift,y_shift,sizes),x_shift,y_shift,sizes) + sigma_2*Cmex(x,sizes);
tic;
bb = HHTmex(WFS,x_shift,y_shift,sizes);
layersXY_hat = pcg(AA,bb,1e-6,1000);
Time_C = toc;

%Matrix
HtH = Hmtx'*Hmtx + sigma_2*Cmtx_;
tic;
bb = Hmtx'*WFS;
layersXY_hat3 = pcg(HtH,bb,1e-6,1000);
Time_Matrix = toc;

%Matrix without inverting Cov
A = Hmtx*Cmtx*Hmtx';
A = A + sigma_2*eye(size(A));
b = WFS;
tic;
layersXY_hat4 = Cmtx*Hmtx'*pcg(A,b,1e-6,1000);
Time_Matrix2 = toc;

%Filter without inverting Cov
AA = @(x) HHmex(Cparmex(HHTmex(x,x_shift,y_shift,sizes),sizes,CovKernel),x_shift,y_shift,sizes) + sigma_2*ones(size(WFS));
bb = WFS;
tic;
temp = pcg(AA,bb,1e-6,1000);
layersXY_hat5 = Cparmex(HHTmex(temp,x_shift,y_shift,sizes),sizes,CovKernel);
Time_Fil2 = toc;

%% **************************************************
% COMPARE
%**************************************************
%C full Covariance
Target_hat2 = HH_projection(layersXY_hat2,subap_index,sizes,alt_p,WFS_target);
Target = HH_projection(layersXY,subap_index,sizes,alt_p,WFS_target);
Error_Target_Full2 = mse(Target-Target_hat2)
Error_layers_Full2 = mse(layersXY-layersXY_hat2)

%C filter Covariance
Target_hat = HH_projection(layersXY_hat,subap_index,sizes,alt_p,WFS_target);
Target = HH_projection(layersXY,subap_index,sizes,alt_p,WFS_target);
Error_Target_Filter1 = mse(Target-Target_hat)
Error_layers_Filter1 = mse(layersXY-layersXY_hat)

%All full
Target_hat3 = HH_projection(layersXY_hat3,subap_index,sizes,alt_p,WFS_target);
Target = HH_projection(layersXY,subap_index,sizes,alt_p,WFS_target);
Error_Target_Full = mse(Target-Target_hat3)
Error_layers_Full = mse(layersXY-layersXY_hat3)

%All full without inverting C
Target_hat4 = HH_projection(layersXY_hat4,subap_index,sizes,alt_p,WFS_target);
Target = HH_projection(layersXY,subap_index,sizes,alt_p,WFS_target);
Error_Target_NIFull = mse(Target-Target_hat4)
Error_layers_NIFull = mse(layersXY-layersXY_hat4)

%Filter without inverting C
Target_hat5 = HH_projection(layersXY_hat5,subap_index,sizes,alt_p,WFS_target);
Target = HH_projection(layersXY,subap_index,sizes,alt_p,WFS_target);
Error_Target_NIFilter = mse(Target-Target_hat5)
Error_layers_NIFilter = mse(layersXY-layersXY_hat5)

%%
%
%**************************************************
% VISUAL COMPARISON
%**************************************************
% Target_vis = zeros(sizes(1));
% Target_vis(subap_index) = Target;
% figure;
% imagesc(real(Target_vis));
% title('Target original');
% Target_hat_vis = zeros(sizes(1));
% Target_hat_vis(subap_index) = Target_hat;
% figure;
% imagesc(real(Target_hat_vis));
% title('Target reconstructed');
% 
% 
% %layer 1
% figure;
% imagesc(real(reshape(layersXY(1:sizes(1)^2),sizes(1),[])))
% title('layer 1 original');
% figure;
% imagesc(real(reshape(layersXY_hat(1:sizes(1)^2),sizes(1),[])))
% title('layer 1 rec');
% 
% 
% %layer 2
% figure;
% imagesc(real(reshape(layersXY(sizes(1)^2+1:sizes(1)^2+sizes(2)^2),sizes(2),[])))
% title('layer 2 original');
% figure;
% imagesc(real(reshape(layersXY_hat(sizes(1)^2+1:sizes(1)^2+sizes(2)^2),sizes(2),[])))
% title('layer 2 rec');

%%

save Simple_Reconstruction_Cov_Results