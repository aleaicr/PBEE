function [Sa_avg] = computeSaAvg(Sa,T_vect,T1,c1,c2,dT)
% Contreras - Sanguinetti
% Ingeniería Sísmica Avanzada - Proyecto de Investigación - USM 2022

% ComputeSaAvg con el método de Eads et al (2013-2015?, no me acuerdo)

% Inputs
% Sa:           Vector con todos los Sa (Aceleración espectral)
% T_vect:       Vector con todos los periodos asociados a Sa
% T:            Periodo fundamental
% c1            Coeficiente para periodo mínimo para Sa_avg(c1*T)
% c2:           Coeficiente para periodo mínimo para Sa_avg (c2*T)
% dT:           Paso del periodo para computar Sa_avg

% Outputs
% Sa_avg: Valor de aceleración espectral media computada con método de Eads
% et al 2013 (o 2015 no me acuerdo)
%%
% Sa_avg = pitatoria_i^N(Sa(c_i*T))^(1/N)
T_interp = (c1*T1:dT:c2*T1).';
Sa_interp = interp1(T_vect,Sa,T_interp,'linear','extrap');                                    % Sa interpolado
Sa_avg = exp(1/length(Sa_interp)*sum(log(Sa_interp)));                      % Sa_avg (Aceleracion espectral media)

end

