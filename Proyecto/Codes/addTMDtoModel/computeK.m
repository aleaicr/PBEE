function [K] = computeK(k)
% P. Heresi THAMDOF

% Generar matriz K con las rigideces por piso (vector de rigideces)

% Input
% k: Vector de rigideces por piso

% Output
% K: matriz de rigideces de EDM por grado de libertad

if length(k) > 1
    k_aux = k(2:end);
    k_aux(end+1,1) = 0;
    K = diag(k+k_aux) - diag(k(2:end),1) - diag(k(2:end),-1);
else
    K = k;
end
end