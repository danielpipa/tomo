clear; close all;

%**************************************************
% INITIALIZATION - WFS DATA
%**************************************************
r0 = 0.2;
L0 = 50;
f0 = 1/L0;
%frac = [0.74,0.02,0.02,0.10,0.12];  %-> Tokovinin2006
frac = [1,1,1,1,1];
tel_diam = 8.0; %Telescope Diameter (m)
ang_dist = 30;  %Angular distance between outer guide stars (arcsec)
WFS_angles = ang_dist*[[sin(0),cos(0)]; [sin(72*pi/180),cos(72*pi/180)]; [sin(144*pi/180),cos(144*pi/180)]; [sin(216*pi/180),cos(216*pi/180)]; [sin(288*pi/180),cos(288*pi/180)]];   %for each WFS, in x,y
WFS_target = [0,0];
WFS_size = 31;   %Number of subapertures in diameter (Only odd numbers)
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
altitudes = [0.,1000.,2000.,4000.,8000.];     %layers altitudes (m) -> Tokovinin2006
layer_size = WFS_size;
alt_p = altitudes/(tel_diam/layer_size);  %layer altitudes relative to pixel size (pixels)
shifts = tan(ang_dist*pi/648000)*alt_p;        %maximum shift in each layer
extra_p = 2*ceil(shifts);        %extra pixels for layer movement
extra_p = extra_p+mod(extra_p,2);    %must be even
sizes = extra_p+layer_size;


%**************************************************
% CREATE COVARIANCE MATRICES
%**************************************************
c = (24/5*gamma(6/5))^(5/6)*gamma(11/6)/(2^(5/6)*pi^(8/3))*(r0*f0)^(-5/3);

fprintf('Creating Covariance Matrices\nMatrix:                 ');
for i = 1:length(sizes)
    X_index = meshgrid(1:sizes(i));
    Y_index = X_index';
    Indexes = [reshape(X_index,1,[]); reshape(Y_index,1,[])];
    mindist = 1e-6; %Avoid zeros
    Cov_M{i} = zeros(sizes(i)^2);
    for lin = 1:sizes(i)^2
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%01d line %04d/%04d',i,lin,sizes(i)^2);
        for col = 1:sizes(i)^2
            %Distance
            dist = sqrt(sum((Indexes(:,col) - Indexes(:,lin)).^2));
            dist = dist*tel_diam/layer_size + mindist;  %distance between voxels
            Cov_M{i}(lin,col) = c*(2*pi*dist*f0)^(5/6)*besselk(5/6,2*pi*f0*dist);
        end
    end
end
fprintf('\n\n');



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
% PROCESS COVARIANCE MATRICES
%**************************************************

for i = 1:length(sizes)
    Cov_M_norm{i} = Cov_M{i}/Cov_M{i}(1,1)*frac(i);
    L{i} = chol(Cov_M_norm{i});    % Cholesky decomposition
end


% Full system matrix
Hmtx = zeros(WFS_size^2*Tn_wfs,sum(sizes.^2));
fprintf('Creating full system Matrix\nIteration:            ');
for i=1:sum(sizes.^2)
    fprintf('\b\b\b\b\b\b\b\b\b\b\b%05d/%05d',i,sum(sizes.^2));
    layers = zeros(sum(sizes.^2));
    layers(i) = 1;
    Hmtx(:,i) = HHmex5(layers,x_shift,y_shift,sizes);
end
fprintf('\n\n');

%Invert Covariances for use in reconstruction

for i = 1:length(sizes)
    Cov_M_norm_{i} = Cov_M_norm{i}^-1;
end

Cmtx_ = blkdiag(Cov_M_norm_{:});
Cmtx = blkdiag(Cov_M_norm{:});

%Create Kernel for Non Inverse Covariance
temp = zeros(layer_size^2,1);
temp((layer_size^2+1)/2) = 1;

[X, Y] = meshgrid(-sizes(end):sizes(end),-sizes(end):sizes(end));
Z = sqrt(X.^2 + Y.^2);
Z = Z.*tel_diam./layer_size + mindist;
for i = 1:sizes(end)*2
    for j = 1:sizes(end)*2
        CovKernel(j,i) = c*(2*pi*Z(j,i)*f0)^(5/6)*besselk(5/6,2*pi*f0*Z(j,i));
    end
end
CovKernel = CovKernel./(max(max(CovKernel)));

save Complete_Reconstruction_InitData
