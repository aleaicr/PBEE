function [rho_val] = rho_Baker_Jayaram_2008(Ti,Tast)
% Determina rho(Tmin,Tmax) de Baker & Jayaram 2008
% Parte I.3.4 - Tarea 1 - Ingeniería Sísmica Avanzada
% Contreras - Sanguinetti

% Determinamos cual es Tmin y Tmax
[Tmin,Tmax] = bounds([Ti,Tast]);

% C1
C1 = 1 - cos(pi/2 - 0.366*log(Tmax/max([Tmin,0.109])));

% C2
if Tmax < 0.2
    C2 = 1 - 0.105*(1 - 1/(1+exp(100*Tmax-5)))*((Tmax-Tmin)/(Tmax-0.0099));
else
    C2 = 0;
end

% C3
if Tmax < 0.109
    C3 = C2;
else
    C3 = C1;
end

% C4
C4 = C1 + 0.5*(sqrt(C3) - C3)*(1 + cos(pi*Tmin/0.109));

% Calculamos rho
if Tmax < 0.109
    rho_val = C2;
elseif Tmin > 0.109
    rho_val = C1;
elseif Tmax < 0.2
    rho_val = min(C2,C4);
else
    rho_val = C4;
end

end