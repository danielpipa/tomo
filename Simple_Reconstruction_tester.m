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
WFS_size = 5;   %Number of subapertures in diameter (Only odd numbers)
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
% CREATE KNOWN PHASE SCREEN
%**************************************************

layersXY = [zeros(1,sizes(1)^2) zeros(1,sizes(2)^2)]';     %Two layers with random slopes
%layersXY(7) = 1+1i;
layersXY(7) = 1;
%layersXY(sizes(1)^2+25) = 1+1i;
layersXY(sizes(1)^2+26) = 1;


%**************************************************
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


%**************************************************
% simulate WFS
%**************************************************
%[WFS] = HH(layersXY,x_shift,y_shift,sizes);
[WFS] = HHmex(layersXY,x_shift,y_shift,sizes);

[layer_teste] = HHT(WFS,x_shift,y_shift,sizes);
%[layer_teste] = HHTmex(WFS,x_shift,y_shift,sizes);



%% **************************************************
% RECONSTRUCTION OF THE LAYERS
%****************************************************

AA = @(x) HHTmex(HHmex(x,x_shift,y_shift,sizes),x_shift,y_shift,sizes);
tic;
bb = HHTmex(WFS,x_shift,y_shift,sizes);
layersXY_hat = pcg(AA,bb,1e-6,1000); % without cov
toc;


%**************************************************
% COMPARE
%**************************************************
Target_hat = HH_projection(layersXY_hat,subap_index,sizes,alt_p,WFS_target);
Target = HH_projection(layersXY,subap_index,sizes,alt_p,WFS_target);

Error = mean(abs(Target-Target_hat).^2)


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

