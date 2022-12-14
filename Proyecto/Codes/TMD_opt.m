%% Optimización de TMD
% Contreras - Sanguinetti
% Ingeniería Sísmica Avanzada - Proyecto de Investigación - USM 2022

% Mostrar en figura las tasas anuales medias de colapso y compararla con la
% SIN TMD

%% Inicializar
clear variables
close all
clc

%% Inputs
CMAF = [1*10^-5; 3.9752*10^-07; 7.1555*10^-08; 5.3017*10^-08; 6.8389*10^-07; 3.3856*10^-06; 6.9115*10^-08; 4.5384*10^-07; 1.8939*10^-07; 6.9115*10^-08; 4.5384*10^-07; 9.758*10^-08; 8.8398*10^-08; 1.3315*10^-06; 1.3315*10^-07; 1.7605*10^-07; 5.3506*10^-07; 1.3315*10^-07; 1.7605*10^-07];                                                                  % Tasa anual media de colapso de todas las combinaciones
% [sin TMD; TMD1; TMD2; ... ; TMD17; TMD18];
T_TMD = [2.24; 2.24; 2.24; 2.18; 2.18; 2.18; 2.12; 2.12; 2.12; 2.35; 2.35; 2.35; 2.18; 2.18; 2.18; 2.06; 2.06; 2.06];
% [TMD1; TDM2; ... ; TMD18]
xi_TMD = [0.05; 0.10; 0.15; 0.05; 0.10; 0.15; 0.05; 0.10; 0.15; 0.05; 0.10; 0.15; 0.05; 0.10; 0.15; 0.05; 0.10; 0.15];

%% TMD óptimo
pos = find(CMAF == min(CMAF)) - 1;
T_TMD_opt = T_TMD(pos);
xi_TMD_opt = xi_TMD(pos);
CMAF_opt = CMAF(pos+1);


%% Generar figuras


figure
stem3(T_TMD,xi_TMD*100, CMAF(2:end,1),'LineWidth',1.7)
hold on
stem3(0,0,CMAF(1,1),'LineWidth',1.7)
stem3(T_TMD_opt,xi_TMD_opt*100,CMAF_opt,'color','r','LineWidth',1.7)
hold off
xlim([0 3])
ylim([0 30])
xlabel('T_{TMD} [sec]')
ylabel('\xi_{TMD} [%]')
zlabel('\lambda_c')
legend('Con TMD','Original','Con TMD óptimo')

figure
stem3(T_TMD,xi_TMD*100, CMAF(2:end,1),'LineWidth',1.7)
hold on
stem3(T_TMD_opt,xi_TMD_opt*100,CMAF_opt,'color','r','LineWidth',1.7)
hold off
xlim([0 3])
ylim([0 30])
xlabel('T_{TMD} [sec]')
ylabel('\xi_{TMD} [%]')
zlabel('\lambda_c')
legend('Con TMD','Con TMD óptimo')