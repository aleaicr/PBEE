%% Tarea 2 - Ingeniería Sísmica Avanzada
% Alexis Contreras - Martina Sanguinetti
% Ingeniería Sísmica Avanzada - USM

% Parte 2 Estructura 1

%% Inicializar
clear variables
close all
clc

%% Estructura
% W = 550; % tonf                                                             % Peso de la estructura
% T = 1.5; % sec                                                              % Periodo fundamental
% zeta = 0.05; % %                                                            % Fracción de amortiguamiento
% Cy = 0.1; %                                                                 % Coeficiente sísmico inelástico
% alfa = 0.02; %                                                              % Coeficiente de endurecimiento post-fluencia
g = 9.81;   % m/s2

%% Parte 2 - Estructura 1
%% 1.
% Currvas IDA utilizando el desplazamiento máximo como EDP y Sa(T1,xi=5%)
% como IM.
% No debería haber colapso ya que no estamos considerando efectos P-Delta y
% modelo no considera degradación.
% Indicar punto de fluencia en las curvas IDA
% incr = 0.1; % g                                                             % Incremento de IM
% IMmax = 3; % g                                                              % IM máxima

% Nombre de carpeta donde están los archivos
ResultsDir = "Resultados1\1\";

% Escribir el nombre de todos los archivos
ResultsNames = dir(ResultsDir);                                                 % Trucazo para leer nombres de una carpeta

% Generar curvas IDA
figure
hold on
for i = 3:length(ResultsNames)
    [EDP,IM] = getIdaCurves2(ResultsDir, string(ResultsNames(i).name));
    % EDP(nFranjas,nRegistros)
    % IM(nFranjas,nRegistros)
    plot(EDP,IM/g,'color','k')
end
xlabel('EDP: Desplazamiento Máximo [m]')
ylabel('IM: Sa(T_1,\xi) [g]')                                               % En IIDAP, se pusieron unidades de metros y segundos => aceleraciones en m/s2
title('Multi-Record IDA Curves', ' Estructura T_1 = 1.5s; \xi = 5%')
grid on

%% 2.

