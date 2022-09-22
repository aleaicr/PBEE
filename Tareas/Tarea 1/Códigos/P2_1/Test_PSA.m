%% Test PSa
% Ingeniería Sísmica Avanzada
% Contreras - Sanguinetti

% An Improved Method for Nonstationary Spectral Matching
% Linda Al Atik & Norman Abrahamson (2010)

%% Inicializar
clear varibales
close all
clc

%% Inputs
% Registro
load('RSN1227_CHICHI_CHY074-N.mat')
Reg = Acc; clear Acc Npts
% Reg   -> Registro de Aceleraciones
% dt    -> Paso temporal del registro de aceleraciones

% Periodos
Tn_init = 0;
Tn_step = 0.1;
Tn_final = 10;
Tn = Tn_init:Tn_step:Tn_final;

% beta (del Método de Newmark)
beta = 0.5;

% Razón de amortiguamiento
xi = 0.05;

% Condiciones iniciales
ui = 0;
udi = 0;

%% Run function
[~,~,~,~,PSa] = Newmark_Lineal(beta,xi,dt,ui,udi,Reg,Tn);

%% Gráficamos
figure
plot(Tn,PSa)
xlabel('Periodos T (sec)')
ylabel('Pseudo-espectro de aceleraciones PSa (g)')
grid on
title('Pseudo-Espectro de Aceleraciones')
legend('RSN1227 CHICHI CHY074-N')