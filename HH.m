function [ out ] = HH( layers, x_shift, y_shift, sizes )
%layer0: Ground Layer
%layer1: Upper layer
%x_shifts: Layer shifts in pixels (X axis)
%y_shifts: Layer shifts in pixels (Y axis)
%sizes: Array with layers sizes

layer0 = reshape(layers(1:sizes(1)^2),[],sizes(1));
layer1 = reshape(layers(sizes(1)^2+1:end),[],sizes(2));

crop = floor((sizes(2)-sizes(1))/2);

% #WFS 1
i = 1;
temp = imtranslate(layer1,[x_shift(i),y_shift(i)]);   %shift
out1 = temp(crop+1:end-crop,crop+1:end-crop);    %crop center
out1 = out1 + layer0;   %Sum layers

% #WFS 2
i = 2;
temp = imtranslate(layer1,[x_shift(i),y_shift(i)]);   %shift
out2 = temp(crop+1:end-crop,crop+1:end-crop);    %crop center
out2 = out2 + layer0;   %Sum layers

% #WFS 3
i = 3;
temp = imtranslate(layer1,[x_shift(i),y_shift(i)]);   %shift
out3 = temp(crop+1:end-crop,crop+1:end-crop);    %crop center
out3 = out3 + layer0;   %Sum layers

% #WFS 4
i = 4;
temp = imtranslate(layer1,[x_shift(i),y_shift(i)]);   %shift
out4 = temp(crop+1:end-crop,crop+1:end-crop);    %crop center
out4 = out4 + layer0;   %Sum layers

% #WFS 5
i = 5;
temp = imtranslate(layer1,[x_shift(i),y_shift(i)]);   %shift
out5 = temp(crop+1:end-crop,crop+1:end-crop);    %crop center
out5 = out5 + layer0;   %Sum layers

out = [out1(:); out2(:); out3(:); out4(:); out5(:)];

end