function [  ] = View_images3( layers, WFS, sizes )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

figure
imagesc(real(reshape(layers(1:sizes(1)^2),sizes(1),[])))
title('layer 1')

figure
imagesc(real(reshape(layers(sizes(1)^2+1:sizes(1)^2+sizes(2)^2),sizes(2),[])))
title('layer 2')

figure
imagesc(real(reshape(layers(sizes(1)^2+sizes(2)^2+1:sizes(1)^2+sizes(2)^2+sizes(3)^2),sizes(3),[])))
title('layer 3')


figure
imagesc(real(reshape(WFS(1:sizes(1)^2),sizes(1),[])))
title('WFS 1')

figure
imagesc(real(reshape(WFS(sizes(1)^2+1:2*(sizes(1)^2)),sizes(1),[])))
title('WFS 2')

figure
imagesc(real(reshape(WFS(2*(sizes(1)^2)+1:3*(sizes(1)^2)),sizes(1),[])))
title('WFS 3')

figure
imagesc(real(reshape(WFS(3*(sizes(1)^2)+1:4*(sizes(1)^2)),sizes(1),[])))
title('WFS 4')

figure
imagesc(real(reshape(WFS(4*(sizes(1)^2)+1:5*(sizes(1)^2)),sizes(1),[])))
title('WFS 5')

end

