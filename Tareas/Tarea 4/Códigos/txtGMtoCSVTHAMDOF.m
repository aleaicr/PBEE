%% Leer Registro y Formato CSV THAMDOF
% Contreras - Sanguinetti
% Ingeniería Sísimica Avanzada, USM 2022

% Leer los registros en formato tareas (en txt) y transformarlo al formato
% THAMDOF
% https://github.com/pheresi/THAMDOF/blob/master/Sample%20Input%20Files/GMs.csv

%% Inicializar
clear variables
close all
clc

%% Inputs
GMFolder = '../Registros';                                                  % Nombre de la carpeta donde se encuentran los registros (todos los registros deben estar en [g] -> las que da el profe están en [g])
unit_GM = 'g';  % solo 'g' o 'm/s2', nada mas puede ir acá, si no va a tirar error, los registros de las tareas vienen en 'g'
cant_registros = 20;                                                        % Cantidad de Registros
GMDataName = 'GM Data.txt';                                                 % Nombre del archivo .txt con los dt de sampling de cada registro
IM_SaT1 = [0.1; 0.4; 0.6; 0.8]; % g                                         % Franjas para realizar el análisis
T = 2; % sec                                                                % Periodo fundamental de la estructura
beta_newmark = 1/4;                                                         % Beta del método de Newmark-beta (dejarlo como 1/4 para que sea incondicionalmente estable)
xi = 0.05;                                                                  % Amortiguamiento para el espectro
t_extra = 3000*ones(cant_registros,1); % [filas]                            % Tiempo extra para el análisis, (n*ones() es el mismo tiempo extra para todos, especificar si se quiere t_extra distinto para un registro)
GMTHAMDOFName = "THAMDOF_data/GMs.csv";
% GMTHAMDOFName_SF1 = "THAMDOF_data/GMS_SF1.csv";

%% Generar archivos
% Registros con SF para que den todos las franjas (todos los registros para todas las franjas)
GMThamdofFormat(GMFolder,GMDataName,T,IM_SaT1,xi,beta_newmark,t_extra,GMTHAMDOFName,unit_GM);
% La aceleración queda en 'g' en el .csv, independiente la unidad con la
% que estén en el registro, para eso era definir la unidad, para poder
% cambiarla.

% Fijarse muy bien, si las unidades del registro están en g o en m/s2,
% cambiar función si están en m/s2, dejar como está si están en g

% Registros con SF = 1 (solo los que están en la carpeta, sin escalar para todas las franjas)
% GMThamdofFormat_SF1(GMFolder,GMDataName,T,xi,beta_newmark,t_extra,GMTHAMDOFName_SF1);