%% Run CMS_baker_2011
% Parte I.3.5 - Tarea 1 -  Ingeniería Sísmica Avanzada
% Contreras - Sanguinetti

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
SaTast = 0.4319;                                                            % g,  Ordenada Espectral del periodo condicinante (UHS Parte I.2)
M = 7.48;                                                                   % Magnitud de momento (Mw) (Parte I.3.a)
R = 18.23;                                                                  % km (Parte I.3.a)
mec_focal = 1;                                                              % Fault_Type de función BSSA_2014_nga.m (Strike-Slip)
region = 1;                                                                 % region de función BSSA_2014_nga.m (California)
z1 = 999;                                                                   % Basin depth (km) de función BSSA_2014_nga.m (Unspecified)
Vs30 = 537;                                                                 % Clasificación C NEHRP

%% Valores UHS (10% en 50 años)
% https://earthquake.usgs.gov/hazards/interactive/
% Para cierto nivel de amenaza que uno quiera (?)
UHS_periods = [0;0.1;0.2;0.3;0.5;0.75;1;2;3;4;5];
UHS_Spectrum = [0.4706;0.9631;1.1412;1.0322;0.779;0.5636;0.4319;0.2056;0.1340;0.0988;0.0774];
% UHS_Tast = 0.4319;                                                          % Sa(1 Sec)

%% Valores Design Spectrum
designSpectrum = load('designSpectrum.mat');
T_DS = designSpectrum.designSpectrum(:,1);
Sa_DS = designSpectrum.designSpectrum(:,2);

%% RUN
[median_CMS,sigma_CMS,mu_BSSA,sigma_BSSA,epsilon_Tast,rho] = CMS_Baker_2011(Ti,Tast,SaTast,M,R,Vs30,mec_focal,region,z1);
% Outputs: 
% median_CMS        -> Media del espectro condicionado [g]
% sigma_CMS         -> Desv.Est del espectro condicionado [g]
% mu_lnSa           -> Media del lnSa (BSSA_2014_nga.m) [g]
% sigma_lnSa        -> Desv. Est del lnSa (BSSA_2014_nga.m) [g]
% epsilon_Tast      -> Valor epsilon(T*)

%% Graficamos

% De lnSa a Ssa
median_CMS = exp(median_CMS);                                               % mu_lnSa       -> mu_Sa
% sigma_CMS = exp(sigma_CMS);                                                 % sigma_lnSa    -> sigma_Sa
mu_BSSA = exp(mu_BSSA);
% sigma_BSSA = exp(sigma_BSSA);

% figure
% plot(UHS_periods,UHS_Spectrum,'--','color','k')
% hold on
% plot(Ti,mu_BSSA,'-','color','k')
% hold off
% xlabel('Periodos (T) [sec]')
% ylabel('Aceleración Espectral [g]')
% legend('UHS 10% en 50 Años', 'Predicted Median BSSA\_2014\_nga.m (Mw = 7.48, R = 18.23km)')
% grid on
% title('UHS vs Predicted Median')

% figure
% plot(Ti,rho)
% xlabel('T*')
% ylabel('\rho(T_i,T*)')
% title('Correlación entre \epsilon(Ti,T*) / Backer & Jayaram (2008)')
% grid on
% legend('T* = 1 sec')

% figure
% plot(UHS_periods,UHS_Spectrum,'--','color','k')
% hold on
% plot(Ti,mu_BSSA,'-','color','k')
% plot(Ti,median_CMS,'-','color','r')
% plot(T_DS,Sa_DS,'color','#0072c3')
% hold off
% xlabel('Periodos (T) [sec]')
% ylabel('Aceleración Espectral [g]')
% legend('UHS 10% en 50 Años', 'Predicted Median BSSA\_2014\_nga.m (Mw = 7.48, R = 18.23km)','Conditional Mean Spectrum','Espectro de Diseño ASCE-7 (Chapter 11)')
% grid on
% title('CMS Baker 2011')
% xlim([0 5])
% ylim([0 1.4])

figure
plot(UHS_periods,UHS_Spectrum,'--','color','k')
hold on
plot(Ti,mu_BSSA,'-','color','k')
plot(Ti,median_CMS,'-','color','r')
plot(Ti,exp(log(median_CMS) + 1*sigma_CMS),'--','color','m')
plot(Ti,exp(log(median_CMS) - 1*sigma_CMS),'--','color','m')
plot(T_DS,Sa_DS,'color','#0072c3')
hold off
xlabel('Periodos (T) [sec]')
ylabel('Aceleración Espectral [g]')
legend('UHS 10% en 50 Años', 'Predicted Mean BSSA\_2014\_nga.m (Mw = 7.48, R = 18.23km)','Conditional Mean Spectrum \mu_{CMS}','\mu_{CMS} +/- \sigma_{CMS}','','Espectro de Diseño ASCE-7 (Chapter 11)')
grid on
title('Conditional Mean Spectrum - Baker 2011')
xlim([0 5])
ylim([0 1.4])

