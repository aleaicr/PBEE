%% Conditional Mean Spectrum
% Contreras - Sanguinetti
% Ingeniería Sísmica Avanzada - Proyecto de Investigación - USM 2022

%% Inicializar
clear variables
close all
clc

%% Inputs
Ti_init = 0.05;
Ti_step = 0.01;
Ti_final = 5;

Ti = (Ti_init:Ti_step:Ti_final)';
Tast = 2.25;    % Sa_Tast hay que interpolarlo
M = 7.23;                                                                   % Magnitud de momento (Mw) (Parte I.3.a)
R = 19.67;                                                                  % km (Parte I.3.a)
mec_focal = 1;                                                              % Fault_Type de función BSSA_2014_nga.m (Strike-Slip)
region = 1;                                                                 % region de función BSSA_2014_nga.m (California)
z1 = 999;                                                                   % Basin depth (km) de función BSSA_2014_nga.m (Unspecified)
Vs30 = 537;                                                                 % Clasificación C NEHRP

%% UHS
% USGS --> Coords: 34.049663; -118.257513 and soil class C
UHS_periods = [0; 0.1; 0.2; 0.3; 0.5; 0.75; 1; 2; 3; 4; 5];
UHS_Spectrum = [0.4859; 0.9978; 1.1844; 1.0676; 0.7817; 0.5544; 0.4139; 0.1812; 0.1098; 0.0768; 0.0682];

% Interpolar Sa(T1) en T = 2.25
% x periodo, y Sa
[B,I] = mink(abs(UHS_periods-Tast),2);                                      % UHS_Periods(I(1)) = 2; UHS_Periods(I(2)) = 3;
SaTast = interp1([UHS_periods(I(1));UHS_periods(I(2))],[UHS_Spectrum(I(1));UHS_Spectrum(I(2))],Tast);

%% Conditional Mean Spectrum & Predicted Mean Spectrum
% returns the predicted too using BSSA_2014_nga.m
[median_CMS,sigma_CMS,mu_BSSA,sigma_BSSA,epsilon_Tast,rho] = CMS_Baker_2011(Ti,Tast,SaTast,M,R,Vs30,mec_focal,region,z1);

median_CMS = exp(median_CMS);
mu_BSSA = exp(mu_BSSA);

figure
plot(UHS_periods,UHS_Spectrum,'--','color','k')
hold on
plot(Ti,mu_BSSA,'-','color','k')
plot(Ti,median_CMS,'-','color','r')
plot(Ti,exp(log(median_CMS) + 1*sigma_CMS),'--','color','m')
plot(Ti,exp(log(median_CMS) - 1*sigma_CMS),'--','color','m')
hold off
xlabel('Periodos (T) [sec]')
ylabel('Aceleración Espectral [g]')
legend('UHS 10% en 50 Años', 'Predicted Mean BSSA\_2014\_nga.m (Mw = 7.48, R = 18.23km)','Conditional Mean Spectrum \mu_{CMS}','\mu_{CMS} +/- \sigma_{CMS}','')
grid on
title('Conditional Mean Spectrum - Baker 2011')
xlim([0 5])
ylim([0 1.4])

matObj = matfile('median_CMS.mat');
matObj.T = Ti;
matObj.median_CMS = median_CMS;