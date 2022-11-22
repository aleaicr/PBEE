%% Collapse EDP
% Contreras - Sanguinetti
% Ingeniería Sísimca Avanzada - USM 2022

% Comentarios:
% - Determinar los parámetros de respuesta estrutural que me generan el
% colapso para cada piso de mi estructura

%% Inicializar
clear variables
close all
clc

%% Inputs
BuildingData = "THAMDOF_data/BuildingData.csv";

% Los datos en formato THAMDOF están estructurados de la siguienta forma
% Story | h | Wi | P | Ki | Fy | as | ac | dcdy | xi | do | Vo
%    1  | 2 | 3  | 4 | 5  | 6  | 7  | 8  | 9    | 10 | 11 | 12 
%% Cargar datos
B = readmatrix(BuildingData);                                               % B de Building Data

%% Determinar defomarción de colapso
dy = zeros(B(end,1),1);
dc = zeros(B(end,1),1);
du = zeros(B(end,1),1);
for i = 1:B(end,1)
    dy(i) = B(i,6)/B(i,5);                                                  % dy = Fy/Ke;
    dc(i) = B(i,9)*dy(i);                                                   % dc = dc/dy * dy; c de Capping
    du(i) = dc(i)-B(i,8)*B(i,5)*(B(i,6)+B(i,7)*B(i,5)*(dc(i)-dy(i)));       % en unidad de [long] de Ke, ej: Ke en [tonf/cm] => dc en [cm]
end

% du es el desplazamiento de colapso
