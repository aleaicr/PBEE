function [EDP_cleared, IM_cleared] = clearNaNEDP(EDP,IM,IM_interp1)
% Contreras - Sanguinetti

% Crear dos nuevas matrices sin los valores de NaN que devuelve IIDAP

% Inputs
% EDP:      Data(i).EDP: Matriz de valores de EDP devueltos por IIDAP 
% IM:       Data(i).IM: Matriz de valores IM asociados a EDP devueltos por IIDAP
% IM_interp1: Matriz para interpolar datos

% Outputs
% EDP_cleared: Matriz de EDP limpiada de valores NaN e interpolada
% IM_cleared: Matriz de IM limpiada de valores NaN e interpolada

[m,n] = size(EDP);
EDP_rmm = rmmissing(EDP);
IM_rmm = rmmissing(IM);

EDP_cleared = interp1(IM_rmm,EDP_rmm,IM_interp1);
IM_cleared = IM_interp1(1:EDP_cleared);