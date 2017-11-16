clear; close all;

%**************************************************
% INITIALIZATION - WFS DATA
%**************************************************
r0 = 0.2;
L0 = 50;
f0 = 1/L0;
frac = [1 1];   %###This will need to be included on the HH_array, leaving at ones for now####
tel_diam = 8.0; %Telescope Diameter (m)
ang_dist = 30;  %Angular distance between outer guide stars (arcsec)
WFS_angles = ang_dist*[[sin(0),cos(0)]; [sin(72*pi/180),cos(72*pi/180)]; [sin(144*pi/180),cos(144*pi/180)]; [sin(216*pi/180),cos(216*pi/180)]; [sin(288*pi/180),cos(288*pi/180)]];   %for each WFS, in x,y
WFS_target = [0,0];
WFS_size = 41;   %Number of subapertures in diameter (Only odd numbers)
% subap_mask = [[0, 0, 1, 1, 1, 0, 0]
%               [0, 1, 1, 1, 1, 1, 0]
%               [1, 1, 1, 1, 1, 1, 1]
%               [1, 1, 1, 1, 1, 1, 1]
%               [1, 1, 1, 1, 1, 1, 1]
%               [0, 1, 1, 1, 1, 1, 0]
%               [0, 0, 1, 1, 1, 0, 0]];

subap_mask = ones(WFS_size);


subap_index = find(subap_mask); %index of valid subapertures
n_subaps = length(subap_index); %number of valid subapertures

%**************************************************
% INITIALIZATION - TURBULENCE DATA
%**************************************************
altitudes = [0.,4000.];     %layers altitudes (m)
layer_size = WFS_size;
alt_p = altitudes/(tel_diam/layer_size);  %layer altitudes relative to pixel size (pixels)
shifts = tan(ang_dist*pi/648000)*alt_p;        %maximum shift in each layer
extra_p = 2*ceil(shifts);        %extra pixels for layer movement
if mod(extra_p,2)
    extra_p = extra_p+1;    %must be even
end
sizes = [layer_size layer_size+extra_p(2)];


%**************************************************
% CREATE COVARIANCE MATRIX
%**************************************************
c = (24/5*gamma(6/5))^(5/6)*gamma(11/6)/(2^(5/6)*pi^(8/3))*(r0*f0)^(-5/3);
X_index = meshgrid(1:layer_size);
Y_index = X_index';
Indexes = [reshape(X_index,1,[]); reshape(Y_index,1,[])];
mindist = 1e-6; %Avoid zeros
Cov_M1 = zeros(layer_size^2);
for lin = 1:layer_size^2
    for col = 1:layer_size^2
        %Distance
        dist = sqrt(sum((Indexes(:,col) - Indexes(:,lin)).^2));
        dist = dist*tel_diam/layer_size + mindist;  %distance between voxels
        Cov_M1(lin,col) = c*frac(1)*(2*pi*dist*f0)^(5/6)*besselk(5/6,2*pi*f0*dist);
    end
end

X_index = meshgrid(1:sizes(2));
Y_index = X_index';
Indexes = [reshape(X_index,1,[]); reshape(Y_index,1,[])];
mindist = 1e-6; %Avoid zeros
Cov_M2 = zeros(sizes(2)^2);
for lin = 1:sizes(2)^2
    for col = 1:sizes(2)^2
        %Distance
        dist = sqrt(sum((Indexes(:,col) - Indexes(:,lin)).^2));
        dist = dist*tel_diam/layer_size + mindist;  %distance between voxels
        Cov_M2(lin,col) = c*frac(2)*(2*pi*dist*f0)^(5/6)*besselk(5/6,2*pi*f0*dist);
    end
end

%% **************************************************
% REGERATE RANDOM DISPLACEMENTS
%**************************************************


Cov_M1_norm = Cov_M1/Cov_M1(1,1);
Cov_M2_norm = Cov_M2/Cov_M2(1,1);

L1 = chol(Cov_M1_norm); % Cholesky decomposition
L2 = chol(Cov_M2_norm); % Cholesky decomposition
%layersXY = [randn(1,sizes(1)^2)*L1+1i*(randn(1,sizes(1)^2)*L1) randn(1,sizes(2)^2)*L2+1i*(randn(1,sizes(2)^2)*L2)]';     %Two layers with random slopes
layersXY = [randn(1,sizes(1)^2)*L1 randn(1,sizes(2)^2)*L2]';     %Two layers with random slopes

%layersXY = layersXY/sqrt(var(layersXY));  %Set variance to one ##There
%must be a better way##: Using now Cov_M1_norm

%clear Cov_M1 Cov_M2 L1 L2

%% **************************************************
% PRE CALCULATIONS FOR HH
%**************************************************
Tn_layers = length(sizes);     %Number of layers
Tn_wfs = length(WFS_angles);   %Number of WFS
n_subaps = sizes(1)^2;     %Number of subapertures
for j = 1:Tn_wfs
    for i = 2:Tn_layers  %skip ground layer
        x_shift(j,i-1) = tan(WFS_angles(j,1)*pi/648000)*alt_p(i);    %Shifts in higher layers
        y_shift(j,i-1) = tan(WFS_angles(j,2)*pi/648000)*alt_p(i);
    end
end


%% **************************************************
% simulate WFS
%**************************************************
[WFS] = HHmex(layersXY,x_shift,y_shift,sizes);


%% **************************************************
% RECONSTRUCTION OF THE LAYERS
%****************************************************

% Hmtx = zeros(numel(WFS),numel(layersXY));
% for i=1:numel(layersXY)
%     layers = zeros(size(layersXY));
%     layers(i) = 1;
%     Hmtx(:,i) = HHmex(layersXY,x_shift,y_shift,sizes);
%     i
% end

%

Cov_M1_norm_ = Cov_M1_norm^-1;
Cov_M2_norm_ = Cov_M2_norm^-1;
C_mtx = @(x,sizes) [Cov_M1_norm_*x(1:sizes(1)^2);Cov_M2_norm_*x(sizes(1)^2+1:end)];

% Add some noise

sigma_2 = var(WFS)/50;
WFS_noisy = WFS + sqrt(sigma_2)*randn(size(WFS));

%%
%AA = @(x) HHT(HH(x,x_shift,y_shift,sizes),x_shift,y_shift,sizes);
%AA = @(x) HHT(HH(x,x_shift,y_shift,sizes),x_shift,y_shift,sizes) + sigma_2*C(x,sizes);
AA = @(x) HHT(HH(x,x_shift,y_shift,sizes),x_shift,y_shift,sizes) + sigma_2*C_mtx(x,sizes);
tic;
bb = HHT(WFS_noisy,x_shift,y_shift,sizes);
layersXY_hat1 = pcg(AA,bb,1e-15,1000);
Time_Mat = toc

AA = @(x) HHTmex(HHmex(x,x_shift,y_shift,sizes),x_shift,y_shift,sizes) + sigma_2*Cmex(x,sizes);;
tic;
bb = HHTmex(WFS,x_shift,y_shift,sizes);
layersXY_hat = pcg(AA,bb,1e-6,1000); % without cov
Time_C = toc

% 
% HtH = Hmtx'*Hmtx;
% tic;
% bb = Hmtx'*WFS;
% layersXY_hat2 = pcg(HtH,bb,1e-6,1000); % without cov
% Time_Matrix = toc




%mse(layersXY_hat-layersXY_hat2)


%%**************************************************
% COMPARE
%**************************************************
Target_hat = HH_projection(layersXY_hat,subap_index,sizes,alt_p,WFS_target);
Target = HH_projection(layersXY,subap_index,sizes,alt_p,WFS_target);

Error_Target = mean(abs(Target-Target_hat).^2)

Error_layers = mean(abs(layersXY-layersXY_hat).^2)


%
%**************************************************
% VISUAL COMPARISON
%**************************************************
Target_vis = zeros(sizes(1));
Target_vis(subap_index) = Target;
figure;
imagesc(real(Target_vis));
title('Target original');
Target_hat_vis = zeros(sizes(1));
Target_hat_vis(subap_index) = Target_hat;
figure;
imagesc(real(Target_hat_vis));
title('Target reconstructed');


%layer 1
figure;
imagesc(real(reshape(layersXY(1:sizes(1)^2),sizes(1),[])))
title('layer 1 original');
figure;
imagesc(real(reshape(layersXY_hat(1:sizes(1)^2),sizes(1),[])))
title('layer 1 rec');


%layer 2
figure;
imagesc(real(reshape(layersXY(sizes(1)^2+1:sizes(1)^2+sizes(2)^2),sizes(2),[])))
title('layer 2 original');
figure;
imagesc(real(reshape(layersXY_hat(sizes(1)^2+1:sizes(1)^2+sizes(2)^2),sizes(2),[])))
title('layer 2 rec');

