%% Run of Spectral Matching Al Atik & Abrahamson 2010
% Ingeniería Sísmica Avanzada -  Tarea 1
% Contreras - Sanguinetti

% An Improved Method for Nonstationary Spectral Matching
% Linda Al Atik & Norman Abrahamson (2010)

%% Inicializar
clear variables
close all
clc

%% Inputs
Reg ;
dt = 1;
Tast = 1; % sec
xi = 0.05; % 5% de amortiguamiento
val_Obj;


%% Run Function
[] = SM_AlAtik_Abrahamson_2010(Reg,dt,Tast,xi,val_Obj);

%% Gráficos
