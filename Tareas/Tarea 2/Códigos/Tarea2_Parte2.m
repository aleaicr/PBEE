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
ResultsDir = 'Resultados1';

% Escribir el nombre de todos los archivos
ResultsName = 'estructura1';                                                 % Trucazo para leer nombres de una carpeta

% Generar curvas IDA
[EDP,IM] = getIdaCurves(ResultsDir, ResultsName);
% EDP(nFranjas,nRegistros)
% IM(nFranjas,nRegistros)
figure
plot(EDP,IM/g,'.','color','#909090')
xlabel('EDP: Desplazamiento Máximo [m]')
ylabel('IM: Sa(T_1,\xi) [g]')                                               % En IIDAP, se pusieron unidades de metros y segundos => aceleraciones en m/s2
title('Multi-Record IDA Curves', ' Estructura T_1 = 1.5s; \xi = 5%')
grid on

figure
plot(EDP,IM/g,'k')
xlabel('EDP: Desplazamiento Máximo [m]')
ylabel('IM: Sa(T_1,\xi) [g]')                                               % En IIDAP, se pusieron unidades de metros y segundos => aceleraciones en m/s2
title('Multi-Record IDA Curves', ' Estructura T_1 = 1.5s; \xi = 5%')
grid on

%% 2.
% Graficar curva IDA mediana y la desviación estándar logarítmica de los
% desplazamientos como función de Sa(T1)

% mu_ln = median(log(x))
% median = geomean(x)
% std_ln = std(log(x))

EDP_muln = mean(log(EDP'));                                                 % Media
EDP_median = geomean(EDP');                                                 % Mediana
EDP_stdln = std(log(EDP'));                                                 % Desviación estandar

% Gráficos

% Mediana logarítmica
figure
plot(IM(:,1)/g,EDP_median','.')
ylabel('Mediana logarítmica de EDP, geomean(EDP)')
xlabel('IM: Sa(T_1,\xi) [g]') 
grid on
title('Mediana logarítmica EDP')

figure
plot(IM(:,1)/g,EDP_median','color','k')
ylabel('Mediana logarítmica de EDP, geomean(EDP)')
xlabel('IM: Sa(T_1,\xi) [g]') 
grid on
title('Mediana logarítmica EDP')

% Desviación Estándar Logarítmica
figure
plot(IM(:,1)/g,EDP_stdln','.')
ylabel('Desviación estandar logarítmica de EDP, std(log(EDP))')
xlabel('IM: Sa(T_1,\xi) [g]') 
grid on
title('Desviación estándar EDP')

figure
plot(IM(:,1)/g,EDP_stdln','color','k')
ylabel('Desviación estandar logarítmica de EDP, std(log(EDP))')
xlabel('IM: Sa(T_1,\xi) [g]') 
grid on
title('Desviación estándar EDP')

% Media
figure
plot(IM(:,1)/g,EDP_muln','.')
ylabel('Media logarítmica de EDP, mean(log(EDP))')
xlabel('IM: Sa(T_1,\xi) [g]') 
grid on
title('Media logarítmica EDP')

figure
plot(IM(:,1)/g,EDP_muln','color','k')
ylabel('Media logarítmica de EDP, mean(log(EDP))')
xlabel('IM: Sa(T_1,\xi) [g]') 
grid on
title('Media logarítmica EDP')

%% 3
% Copiar gráficos media vs Sa, pero agregar curva de desplazamientos
% máximos que se obtendrían asumiendo principio de igualdad de
% desplazamientos

% Se cumpliría este principio solo en la primera estructura ya tiene periodo
% mayor a 1[s], pero no para la segunda estructura

%% 4
% Suponiendo distribución lognormal para EDP dado IM, calcular probabilidad
% de que los desplazamientos máximos excedan 25cm para una ordenada
% espectral de 1g.

% Desplazamiento máximo objetivo
desplMaxObj = 25/100; % 25  cm = 25/100 m

% Posición IM = 1g
posIM1g = find(IM(:,1)/g == 1);

% Columna de EDP que necesito para este análisis
EDP_col = EDP(:,posIM1g);
EDP_logncdf = logncdf(EDP_col);
EDP_valuesInterp = interp1(EDP_col,EDP_logncdf,desplMaxObj);
% Probabilidad
figure
plot(EDP_col,EDP_logncdf)
hold on
plot(desplMaxObj,EDP_valuesInterp,'^r')
hold off
xlabel('EDP|IM = 1g')
ylabel('logncdf(EDP)|IM = 1g')
text(desplMaxObj+0.1,EDP_valuesInterp,['P(EDP > 25cm | IM = 1g) =' string(EDP_valuesInterp)])
grid on

%% 5
% Repetir proceso pero ahora con Sa_avg como IM, para definir Sa_avg
% utilizar 20 periodos entre 0.2T y 3.0T





