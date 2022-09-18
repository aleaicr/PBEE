%% Run CMS_baker_2011
% Script de prueba/ejemplo para ejecutar función CMS_Baker_2011.m

%% Inicializar
clear variables
close all
clc

%% Inputs función
Ti_init = 0.05;                                                             % Mínimo TI que indica Baker (2011)
Ti_step = 0.01;                                                             % Paso para los periodos
Ti_final = 5;                                                               % Máximo Ti que indica Baker (2011)

Ti = Ti_init:Ti_step:Ti_final; % sec
Tast = 1; % sec
SaTast = 0.4259; % g
M = 8; % Magnitud de momento (Mw)
R = 10; % km
Vs30 = 537; % Clasificación C NEHRP
mec_focal = 0; % Fault_Type de función BSSA_2014_nga.m (Unspecified)
region = 1; % region de función BSSA_2014_nga.m (California)
z1 = 999; % Basin depth (km) de función BSSA_2014_nga.m (Unspecified)

%% RUN

[median_CMS,sigma_CMS] = CMS_Baker_2011(Ti,Tast,SaTast,M,R,Vs30,mec_focal,region,z1);

%% Graficamos
figure
plot(Ti, median_CMS)
xlabel('Periodo (T) [sec]')
ylabel('\mu_{lnSa(Ti)|lnSa(T*)')
title('\mu CMS Baker 2011')
grid on

figure
plot(Ti, sigma_CMS)
xlabel('Periodo (T) [sec]')
ylabel('\sigma_{lnSa(Ti)|lnSa(T*)')
title('\sigma CMS Baker 2011')
