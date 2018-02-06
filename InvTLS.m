% Funcion que invierte matriz rectangular, habiendo previamente reducido
% el orden de acuerdo a valores singulares muy pequenhos
function [Wi,red]=InvTLS(W,tol)
Sd = svd(W);    % calcula valores singulares, ordenados por valor
imax=sum(Sd>tol*Sd(1)); % numero de valores singulares sobre tolerancia
Smax=Sd(imax);  % tolerancia para valores singulares a considerar
Wi=pinv(W,Smax); % invierte matriz solo para val. singulares significativos
red=max(size(Sd))-imax; % numero de valores singulares truncados