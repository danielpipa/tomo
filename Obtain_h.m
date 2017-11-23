%Obtain h

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


%% **************************************************
% REGERATE RANDOM DISPLACEMENTS
%**************************************************


Cov_M1_norm = Cov_M1/Cov_M1(1,1);

Cov_M1_norm_ = Cov_M1_norm^-1;


%% **************************************************
% Get h From Cov
%**************************************************

temp = zeros(WFS_size);
temp((WFS_size+1)/2,(WFS_size+1)/2) = 1;

A = Cov_M1_norm_*temp(:);
A = reshape(A,WFS_size,[]);

center = (WFS_size+1)/2;
h = A(center-6:center+6,center-6:center+6);