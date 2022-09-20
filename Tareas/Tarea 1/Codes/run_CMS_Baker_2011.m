%% Run CMS_baker_2011
% Script de prueba/ejemplo para ejecutar función CMS_Baker_2011.m
% (Parte I.3.5 Tarea 1 - Sísmica Avanzada)

%% Inicializar
clear variables
close all
clc

%% Inputs función
Ti_init = 0.05;                                                             % Mínimo TI que indica Baker (2011)
Ti_step = 0.01;                                                             % Paso para los periodos
Ti_final = 5;                                                               % Máximo Ti que indica Baker (2011)

Ti = (Ti_init:Ti_step:Ti_final)';                                           % sec, Vector de Periodos
Tast = 1;                                                                   % sec, Periodo Condicionante
Sa_Tast = 0.4319;                                                           % g,  Ordenada Espectral del periodo condicinante (UHS Parte I.2)
M = 7.48;                                                                   % Magnitud de momento (Mw) (Parte I.3.a)
R = 18.23;                                                                  % km (Parte I.3.a)
Vs30 = 537;                                                                 % Clasificación C NEHRP
mec_focal = 1;                                                              % Fault_Type de función BSSA_2014_nga.m (Strike-Slip)
region = 1;                                                                 % region de función BSSA_2014_nga.m (California)
z1 = 999;                                                                   % Basin depth (km) de función BSSA_2014_nga.m (Unspecified)

%% RUN
[median_CMS,sigma_CMS,mu_lnSa,sigma_lnSa,epsilon_Tast] = test_CMS_Baker_2011(Ti,Tast,Sa_Tast,M,R,Vs30,mec_focal,region,z1);

%% Graficamos
figure
plot(Ti, median_CMS)
xlabel('Periodo (T) [sec]')
ylabel('\mu_{lnSa(Ti)|lnSa(T*)}')
title('\mu CMS Baker 2011')
grid on

figure
plot(Ti, sigma_CMS)
xlabel('Periodo (T) [sec]')
ylabel('\sigma_{lnSa(Ti)|lnSa(T*)}')
title('\sigma CMS Baker 2011')
grid on

figure
plot(Ti, mu_lnSa)
hold on
plot(Ti, median_CMS)
plot(Ti, median_CMS+0.16*sigma_CMS)
plot(Ti, median_CMS+0.84*sigma_CMS)
hold off
xlabel('Periodo (T) [sec]')
ylabel('lnSa_{CMS}')
title('CMS Baker 2011')
legend('median','median_{CMS}','median_{CMS}+0.16sigma_{CMS}','median_{CMS}+0.84sigma_{CMS}')
grid on
