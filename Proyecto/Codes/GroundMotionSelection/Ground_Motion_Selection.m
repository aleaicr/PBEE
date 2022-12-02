%% Ground_Motion_Selection
% Contreras - Sanguinetti
% Ingeniería Sísmica Avanzada - Proyecto de Investigación - USM 2022

% Eads, Miranda & Lignos (2015)
% Selección de Registros por Sa_avg

% Parámetros
% T:            Periodo fundamental, c1*T,c2*T,....cN*T
% c1, cN:       Coeficientes para definir el rango para calcular Sa_avg
% Tn:           Vector de periodos para el análisis (para generar el espectro de aceleraciones)
% n_ interpol:  Cantidad de interpolaciones (para mejorar presición)

%% Inicializar
clear variables
close all
clc

%% Inputs
Tast = 2.18; % sec                                                          % Periodo del primer modo de la estructura
c1 = 0.2;                                                                   % Factor para determinar el T_min para calcular Sa_avg
cN = 3.0;                                                                   % Factor para determinar el T_max para calcular Sa_avg
n_reg = 20;                                                                 % Cantidad de registros a seleccionar

% Periodos para espectros
Ti_init = 0.05;
Ti_step = 0.01;
Ti_final = 6.6;

%% Load Database

% Espectro de media condicionada
load('ConditionalMeanSpectrum/median_CMS.mat')                              % Cargar espectro de media condicionada
% median_CMS y T

% Data del NGA
load('NGa_Data/NGA_W2_meta_data.mat','Sa_1','Sa_2','Periods','magnitude','Rjb','soil_Vs30','station_name')               % Cargar base da datos de registros
% Sa_1 y Sa_2 son los espectros de mismos registros pero en dos direcciones,
% se pueden tomar como registros distintos y dejarlos en una sola matriz de
% forma consecutiva
Periods = Periods.';

%% Cálculo de Sa_avg del CMS (Conditional Mean Spectrum)
% Rango para obtener Sa_avg
Tn_inicial = c1*Tast;
Tn_final = cN*Tast;
Tn_interp = sort([(Tn_inicial:Ti_step:Tn_final).'; Tast]);

% Interpolar CMS
median_CMS_interp = interp1(T,median_CMS,Tn_interp,'linear','extrap');
Sa_avg_CMS = exp((1/length(median_CMS_interp))*sum(median_CMS_interp));

figure
plot(Tn_interp,median_CMS_interp,'.')
hold on
plot(Tast,median_CMS_interp(Tn_interp == Tast,1),'o','color','r','LineWidth',2)
hold off
grid on
xlabel('Periodo [s]')
ylabel('Sa(T_1) [g]')
title('Rango de CMS')
legend('Puntos para obtener Sa_{avg}','Sa(T_1)')

fprintf('Sa_avg_CMS = %.4f [g]\n\n',Sa_avg_CMS)

%% Ground Motion Selection
% Los datos están ordenados de la siguiente manera:
% Sa_1 y Sa_2 (rows: Ordenada Espectral, columns: Periodos)
% Periods = (rows: Periodos, columns: 1)

% Inicializar vector Sa_avg para cada registro
Sa_avg = zeros(2*length(Sa_1.'),1);

% Espectros ordenados
Sa = [Sa_1; Sa_2];
cant_reg_select = 0;
while cant_reg_select ~= 20
    % Calcular Sa_avg para cada registro
    for i = 1:2*length(Sa_1.')
        % Interpolar valores de T para el espectro del registro i
        Sa_interp = interp1(Periods,Sa(i,:).',Tn_interp);
        Sa_avg(i,1) = exp((1/length(Sa_interp))*sum(Sa_interp));
    end
    
    % Seleccionar aquel que está más cerca
    [Error,Index] = mink(abs(Sa_avg - Sa_avg_CMS),n_reg);                   % length(Index) --> cantidad de registros
    if length(unique(Index)) < 20                                           % Si 
       n_reg = n_reg + 1;
    end
    cant_reg_select = length(unique(Index));
end

% Detecar cual registro es
componente = ones(n_reg,1);
for i = 1:length(Index)
    if Index(i,1) > length(Sa_1.')
        componente(i,1) = 2;
        Index(i,1) = Index(i,1) - length(Sa_1.');
    end
end

vs30_reg = soil_Vs30(Index);
rjb_reg = Rjb(Index);
magnitude_reg = magnitude(Index);
station_name_reg = station_name(Index);

tabla = table();
tabla.Reg = (1:1:n_reg)';
tabla.StationName = station_name_reg;
tabla.Sequence = Index;
tabla.Sa_avg = Sa_avg(Index);
tabla.Component = componente;
tabla.Magnitude = magnitude_reg;
tabla.rjb = rjb_reg;
tabla.vs30 = vs30_reg;
disp(tabla)
clear tabla


