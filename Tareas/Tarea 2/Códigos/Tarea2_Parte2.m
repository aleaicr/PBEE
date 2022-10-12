%% Tarea 2 - Ingeniería Sísmica Avanzada
% Alexis Contreras - Martina Sanguinetti
% Ingeniería Sísmica Avanzada - USM

% Parte 2 Estructura 1

%% Inicializar
clear variables
close all
clc

%% Estructura
W = 550; % tonf                                                             % Peso de la estructura
T = 1.5; % sec                                                              % Periodo fundamental
zeta = 0.05; % %                                                            % Fracción de amortiguamiento
Cy = 0.1; %                                                                 % Coeficiente sísmico inelástico
alfa = 0.02; %                                                              % Coeficiente de endurecimiento post-fluencia


%% Parte 2 - Estructura 1

%% 1.
% Currvas IDA utilizando el desplazamiento máximo como EDP y Sa(T1,xi=5%)
% como IM

% No debería haber colapso ya que no estamos considerando efectos P-Delta y
% modelo no considera degradación

% Indicar punto de fluencia en las curvas IDA

incr = 0.1; % g                                                             % Incremento de IM
IMmax = 3; % g                                                              % IM máxima

ResultsFile = [];                                                           % Ingresar todos los nombres de los archivos resultado
figure
hold on
for i = 1:length(ResultsFile)
    [EDP,IM] = getIdaCurves(ResultsFile(i));
    % EDP(nFranjas,nRegistros)
    % IM(nFranjas,nRegistros)
    plot(EDP,IM)
    leyenda = [leyenda]
end
xlabel('EDP: Desplazamiento Máximo')
ylabel('IM: Sa(T1,\xi)')
title('Multi-Record IDA Curves')