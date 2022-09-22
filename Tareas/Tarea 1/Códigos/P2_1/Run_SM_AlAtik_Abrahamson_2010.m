%% Run of Spectral Matching Al Atik & Abrahamson 2010
% Ingeniería Sísmica Avanzada -  Tarea 1
% Contreras - Sanguinetti

% An Improved Method for Nonstationary Spectral Matching
% Linda Al Atik & Norman Abrahamson (2010)

% PSa se obtiene con el Método de Newmark (Ing.Sísmica - Contreras-Adasme)


%% Inicializar
clear variables
close all
clc
%% Load Registro
load('RSN1227_CHICHI_CHY074-N.mat')
Reg = Acc;
% dt: Paso temporal del registro
% Npts: Desconocido

%% Inputs
Tast = 1; % sec                                                             % Periodo de ajuste
xi = 0.05; % 5% de amortiguamiento                                          % Razón de amortiguamiento del espectro
val_Obj = 1;                                                                    % Valor objetivo


%% Run Function
[] = SM_AlAtik_Abrahamson_2010(Reg,dt,Tast,xi,val_Obj);

%% Gráficos

