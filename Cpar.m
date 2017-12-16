function [ out ] = Cpar( x, sizes, h )
%The same as Cmex, but receiving an external kernel

layer0 = reshape(x(1:sizes(1)^2),[],sizes(1));
layer1 = reshape(x(sizes(1)^2+1:end),[],sizes(2));

layer0 = imfilter(layer0,h);
layer1 = imfilter(layer1,h);

% layer0 = imfilter(layer0,h,'symmetric');
% layer1 = imfilter(layer1,h,'symmetric');

% layer0 = imfilter(layer0,h,'replicate');
% layer1 = imfilter(layer1,h,'replicate');

% layer0 = imfilter(layer0,h,'circular');
% layer1 = imfilter(layer1,h,'circular');

% layer0 = imfilter(layer0,h,'conv');
% layer1 = imfilter(layer1,h,'conv');

out = [layer0(:); layer1(:)];

end

