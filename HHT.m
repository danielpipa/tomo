function [ out ] = HHT( WFS, x_shift, y_shift, sizes )
%WFS1-5: Wavefront sensor images, 1 to 5
%x_shifts: Layer shifts in pixels (X axis)
%y_shifts: Layer shifts in pixels (Y axis)
%sizes: Array with layers sizes


WFS1 = reshape(WFS(1:sizes(1)^2),[],sizes(1));
WFS2 = reshape(WFS(sizes(1)^2+1:2*sizes(1)^2),[],sizes(1));
WFS3 = reshape(WFS(2*sizes(1)^2+1:3*sizes(1)^2),[],sizes(1));
WFS4 = reshape(WFS(3*sizes(1)^2+1:4*sizes(1)^2),[],sizes(1));
WFS5 = reshape(WFS(4*sizes(1)^2+1:end),[],sizes(1));

crop = floor((sizes(2)-sizes(1))/2);

% #Layer 0
out0 = WFS1 + WFS2 + WFS3 + WFS4 + WFS5;

% #Layer 1
% #WFS 1
i = 1;
temp = padarray(WFS1,[crop crop]);
WFS11 = imtranslate(temp,[-x_shift(i),-y_shift(i)]);
% #WFS 2
i = 2;
temp = padarray(WFS2,[crop crop]);
WFS21 = imtranslate(temp,[-x_shift(i),-y_shift(i)]);
% #WFS 3
i = 3;
temp = padarray(WFS3,[crop crop]);
WFS31 = imtranslate(temp,[-x_shift(i),-y_shift(i)]);
% #WFS 4
i = 4;
temp = padarray(WFS4,[crop crop]);
WFS41 = imtranslate(temp,[-x_shift(i),-y_shift(i)]);
% #WFS 5
i = 5;
temp = padarray(WFS5,[crop crop]);
WFS51 = imtranslate(temp,[-x_shift(i),-y_shift(i)]);

out1 = WFS11 + WFS21 + WFS31 + WFS41 + WFS51;

out = [out0(:); out1(:)];
end