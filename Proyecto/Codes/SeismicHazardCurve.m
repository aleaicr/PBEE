%% Curva de amenaza sísmica
% Alexis Contreras - Martina Sanguinetti
% Ingeniería sísmica avanzada - Proyecto Investigación - USM 2022

% Solo contiene los vectores de la amenaza sísmica

%% Inicializar
clear variables
close all
clc
%%
IM_Saavg = [0.01;0.013;0.015;0.019;0.023;0.028;0.033;0.041;0.049;0.06;0.073;0.089;0.108;0.131;0.159;0.193;0.234;0.285;0.346;0.421;0.511;0.621;0.755;0.917;1.115;1.354;1.646;2];
lambda_Saavg = [0.274529554;0.223874996;0.179741466;0.142766053;0.112165203;0.086762584;0.066803522;0.050979124;0.038715512;0.029373076;0.022110009;0.016878141;0.012882974;0.009983871;0.007564458;0.005724685;0.004229869;0.003054078;0.002095112;0.001332039;0.000780955;0.000424777;0.000215999;0.000114034;5.8321E-05;1.91771E-05;6.21444E-06;2.10781E-06];

figure
loglog(IM_Saavg,lambda_Saavg)
xlabel('IM: Sa_{avg} [g]')
ylabel('\lambda_{Sa_{avg}} [1/yr]')
title('Curva de amenaza sísmica')
grid on