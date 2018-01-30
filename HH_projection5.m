function [ out ] = HH_projection5( layers, x_shift,y_shift,sizes, alt_p, WFS_target )


Tn_layers = length(sizes);     %Number of layers
n_subaps = sizes(1)^2;     %Number of subapertures
for j = 1:5
    for i = 2:Tn_layers  %skip ground layer
        x_shift(j,i-1) = tan(WFS_target(1)*pi/648000)*alt_p(i);    %Shifts in higher layers
        y_shift(j,i-1) = tan(WFS_target(2)*pi/648000)*alt_p(i);
    end
end

out = HHmex5(layers,x_shift,y_shift,sizes);

out = out(1:sizes(1)^2);