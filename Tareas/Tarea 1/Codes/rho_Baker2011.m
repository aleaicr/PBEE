function [rho_val] = rho_Baker2011(Ti,T_ast)
% Determina rho(Tmin,Tmax) del Step 3 de Baker (2011)

% Determinamos cual es Tmin y Tmax
if Ti < T_ast
    Tmin = Ti;
    Tmax = T_ast;
else
    Tmin = T_ast;
    Tmax = Ti;
end

% Determinamos I
if Tmin < 0.189
    I = 1;
else
    I = 0;
end

% Aplicamos formulita
rho_val = 1-cos(pi/2 - (0.359 + 0.163*I*log(Tmin/0.189))*log(Tmax/Tmin));   % Retornamos rho_val
end