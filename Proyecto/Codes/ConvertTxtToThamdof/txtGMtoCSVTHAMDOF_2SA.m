%% Leer Registro y Formato CSV THAMDOF
% Contreras - Sanguinetti
% Ingeniería Sísimica Avanzada, USM 2022

% Leer los registros en formato tareas (en txt) y transformarlo al formato
% THAMDOF
% https://github.com/pheresi/THAMDOF/blob/master/Sample%20Input%20Files/GMs.csv

% Comentarios
% - Cambiar la carpeta donde están los registros guardados

%% Inicializar
clear variables
close all
clc

%% Inputs
GMFolder = 'Registros';                                                     % Nombre de la carpeta donde se encuentran los registros (todos los registros deben estar en [g] -> las que da el profe están en [g])
cant_registros = 31;
GMDataName = 'GM Data.txt';                                                 % Nombre del archivo .txt con los dt de sampling de cada registro
IM_Sa_avg = [0.5; 1.5]; % g                                                 % Franjas para realizar el análisis
T = 2.18; % sec                                                             % Periodo fundamental de la estructura
beta_newmark = 1/4;                                                         % Beta del método de Newmark-beta (dejarlo como 1/4 para que sea incondicionalmente estable)
xi = 0.05;                                                                  % Amortiguamiento para el espectro
t_extra = 100*ones(cant_registros,1); % [filas]                             % Tiempo extra para el análisis, (n*ones() es el mismo tiempo extra para todos, especificar si se quiere t_extra distinto para un registro)
GMTHAMDOFName = "GMs_2SA.csv";
c1 = 0.2;
c2 = 3;

%% Generar archivos
% Registros con SF para que den todos las franjas (todos los registros para todas las franjas)
GMThamdofFormat(GMFolder,GMDataName,T,IM_Sa_avg,c1,c2,xi,beta_newmark,t_extra,GMTHAMDOFName);
