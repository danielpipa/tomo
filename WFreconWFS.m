function [layer1WF] = WFreconWFS(layersX,layersY,sizes,d,tol)
%layersX -> Vector with slopes in X direction
%layersX -> Vector with slopes in Y direction
%sizes -> Sizes of layers
%d -> size of subapertures, in m
%tol -> Tolerance for singular value truncation


layerX1 = layersX(1:sizes(1)^2);
layerY1 = layersY(1:sizes(1)^2);
layerXY1 = [layerX1';layerY1'];
layerXY1 = layerXY1(:);   %mix X and Y slopes

% Contruct Matrix B (n*n*2 rows by (n+1)*(n+1) columns)
n = sizes(1);
B = zeros(n*n*2,(n+1)*(n+1));
% Index matrix for the slopes interconections
index = reshape((1:(n+1)*(n+1)), n+1, n+1)';

% Go through every X and Y slope and scribe it's information into the
% matrix
row = 0;
for y = 1:n
    for x = 1:n
        % X Slope
        row = row + 1;
        B(row, index(y,x+1)) = 1;
        B(row, index(y,x)) = -1;
        B(row, index(y+1, x+1)) = 1;
        B(row, index(y+1, x)) = -1;
        % Y Slope
        row = row + 1;
        B(row, index(y,x+1)) = -1;
        B(row, index(y,x)) = -1;
        B(row, index(y+1, x+1)) = 1;
        B(row, index(y+1, x)) = 1;
    end
end
% % Add ones line to force 0 average
B(row+1, :) = ones(1, (n+1)*(n+1) );

B = B/(2*d);    %Scale matrix to subapertures sizes

[InvB,reduccion]=InvTLS(B,tol);           %Truncated inversion of Matrix B

[wn,wm]=size(InvB);
InvB=InvB(:,1:wm-1); % Subtract column corresponding to coefficients average

layer1WF = InvB*layerXY1;

%layer1WF = reshape(layer1WF,(sizes(1)+1),[]);



end

