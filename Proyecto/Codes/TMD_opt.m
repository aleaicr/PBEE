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
CMAF = [0.1; 0.05; 0.02; 0.09];                                                                  % Tasa anual media de colapso de todas las combinaciones
% [sin TMD; TMD1; TMD2; ... ; TMD17; TMD18];
T_TMD = [2; 1.5; 3];
% [TMD1; TDM2; ... ; TMD18]
xi_TMD = [0.1; 0.12; 0.15];

figure
stem3(T_TMD,xi_TMD*100, CMAF(2:end,1))
hold on
stem3(0,0,CMAF(1,1))
hold off
xlabel('T_{TMD} [sec]')
ylabel('\xi_{TMD} [%]')
zlabel('\lambda_c')


figure
surf(T_TMD,xi_TMD*100, CMAF(2:end,1))
hold on
surf(0,0,CMAF(1,1))
hold off
xlabel('T_{TMD} [sec]')
ylabel('\xi_{TMD} [%]')
zlabel('\lambda_c')
