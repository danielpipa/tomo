function [ out ] = C( x, sizes )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

layer0 = reshape(x(1:sizes(1)^2),[],sizes(1));
layer1 = reshape(x(sizes(1)^2+1:end),[],sizes(2));

% h = [0.0685,-0.3458,1];      %Empiracly obtained
% h = [h,fliplr(h(1:end-1))];
% h = conv2(h,h');

h = [0.0047   -0.0237    0.0685   -0.0237    0.0047;
   -0.0237    0.1196   -0.3458    0.1196   -0.0237;
    0.0685   -0.3458    1.0000   -0.3458    0.0685;
   -0.0237    0.1196   -0.3458    0.1196   -0.0237;
    0.0047   -0.0237    0.0685   -0.0237    0.0047];

layer0 = imfilter(layer0,h);
layer1 = imfilter(layer1,h);

out = [layer0(:); layer1(:)];

end

