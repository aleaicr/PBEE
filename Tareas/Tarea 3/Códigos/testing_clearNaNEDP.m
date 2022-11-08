%% Inicializar
clear variables
close all
clc

%% Inputs
g = 9.81; % m/s2
cant_estr = 2;                                                              % Cantidad de estructutras
cant_variantes = 4;                                                         % Cantidad de variantes por cada estructura

ResultsFile = "est_1_A";
ResultsDir = "IIDAP_T3";
n_ests = length(ResultsFile);
IM_interp1 = (0.1:0.1:20).';

%% Load Data
[EDP,IM,IMc,Backbone] = getIdaCurves_v2(convertStringsToChars(ResultsDir), convertStringsToChars(ResultsFile));
EDP_rmm = rmmissing(EDP);
IM_rmm = rmmissing(IM);
EDP_interp1 = interp1(IM_rmm,EDP_rmm,IM_interp1);

figure
plot(IM_interp1,EDP_interp1)