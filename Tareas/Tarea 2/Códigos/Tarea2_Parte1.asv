%% Tarea 2 Parte 1 - Ingeniería Sísmica Avanzada
% Alexis Contreras - Martina Sanguinetti
% Ingeniería Sísmica Avanzada - USM

%% Inicializar
clear variables
close all
clc

%% Registros a utilizar
% Incluir nuevo registro
AT2fileDir = '';                                                            % Carpeta en la que se encuentra el archivo AT2
AT2fileName = 'RSN1101_KOBE_AMA000.AT2';                                    % Nombre del archivo AT2
newTxtFileDir = '..\Registros';                                             % Carpeta en la que se quiere escribir el nuevo archivo
newTxtFileName = 'RSN1101_KOBE.txt';                                        % Nombre del nuevo archivo

[Reg,dt] = readAT2(AT2fileDir,AT2fileName,newTxtFileDir,newTxtFileName);

% Escribir nuevo registro en GM Data.txt
GMDataDir = '..\Registros';
GMDataName = 'GM Data.txt';

writeGMData(newTxtFileName,dt,GMDataDir,GMDataName)                         % Escribir 'newTxtFileName   dt' en archivo GM Data.txt
% No pude modificar la cantidad de veces que me lo escribe, cuidado con
% correr varias veces esta línea ya que  volvería a escribir varias veces
% en archivo GM Data.txt

