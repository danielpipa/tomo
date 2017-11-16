function [ out ] = HH_projection( data, subap_index, sizes, alt_p, WFS_angle )
%data: If transpose-> Complex cell with WFS slopes. If not complex cell with layers slopes
%sizes: Array with layers sizes
%alt_p: Altitudes of layers in pixel lenght
%WFS_angles: Directions of each guide star in arcseconds. Each array line corresponds to a WFS
%t -> transpose or not


Tn_layers = length(sizes);
layersXY{1} = reshape(data(1:sizes(1)^2),sizes(1),[]);
for l_n = 2:Tn_layers
    layersXY{l_n} = reshape(data((l_n-1)*sizes(l_n-1)^2+1:(l_n-1)*sizes(l_n-1)^2+sizes(l_n)^2),sizes(l_n),[]);
end

out = P(S(layersXY,sizes,alt_p,WFS_angle));    %output is the WFSs slopes
out = reshape(out(subap_index),[],1);

end


function y = S(layersXY,sizes, alt_p,angles)
%Shift operator
%Receives layers and shifts them according with angles and altitudes
N_layers = length(alt_p);   %number of layers
sizeG = sizes(1);    %size of Ground Layer (and wfs)
for i = 2:N_layers  %skip ground layer
    x_shift = tan(angles(1)*pi/648000)*alt_p(i);
    y_shift = tan(angles(2)*pi/648000)*alt_p(i);
    temp = imtranslate(layersXY{i},[x_shift,y_shift]);   %shift
    crop = (length(temp)-sizeG)/2;
    layersXY{i} = temp(crop+1:end-crop,crop+1:end-crop);    %crop center
end

y = layersXY;
end

function y = St(layersXY,sizes,alt_p,angles)
%Shift operator
%Receives layers and shifts them according with angles and altitudes
N_layers = length(alt_p);   %number of layers
sizeG = sizes(1);    %size of Ground Layer (and wfs)
for i = 2:N_layers  %skip ground layer
    x_shift = tan(angles(1)*pi/648000)*alt_p(i);
    y_shift = tan(angles(2)*pi/648000)*alt_p(i);
    %zero padding to size of layer
    temp = padarray(layersXY{i},[(sizes(i)-sizeG)/2 (sizes(i)-sizeG)/2]);
    layersXY{i} = imtranslate(temp,[-x_shift,-y_shift]);   %shift
end

y = layersXY;
end

function y = P(layersXY)
%Projection Operator
%Sums all layers
N_layers = length(layersXY);   %number of layers
y = layersXY{1};
for i = 2:N_layers
    y = y+layersXY{i};
end

end

function y = Pt(layersXY,sizes)
%Projection Operator
%Sums all layers
N_layers = length(sizes);   %number of layers
for i = 1:N_layers
    y{i} = layersXY;
end

end