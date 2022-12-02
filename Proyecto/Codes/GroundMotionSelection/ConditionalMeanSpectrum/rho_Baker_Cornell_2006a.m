function [rho] = rho_Baker_2011(Ti,Tast)
% Determina rho(Tmin,Tmax) del Step 3 de Baker (2011)

% Determinamos cual es Tmin y Tmax
[Tmin,Tmax] = bounds([Ti,Tast]);

% Determinamos I
if Tmin < 0.189
    I = 1;
else
    I = 0;
end

% Aplicamos formulita
rho = 1-cos(pi/2 - (0.359 + 0.163*I*log(Tmin/0.189))*log(Tmax/Tmin));   % Retornamos rho_val

end