%% Leer Archivos AT2 y escribit txt
% Alexis Contreras - Martina Sanguinetti
% Ingeniería Sísmica Avanzada - Proyecto Investigación - USM 2022

%% Inicializar
clear variables
close all
clc

%% Registros a utilizar
AT2fileDir = 'GMSelected_AT2';
GMDataDir = 'GMSelected_txt';
GMDataName = 'GM Data.txt';
newTxtFileDir = 'GMSelected_txt';                                             % Carpeta en la que se quiere escribir el nuevo archivo
files = dir(AT2fileDir);
files(1:1:2) = [];
for i = 1:length(files)
    % Incluir nuevo registro
    AT2fileName = files(i).name;                                    % Nombre del archivo AT2
    newTxtFileName = convertStringsToChars(string(files(i).name(1:end-4)) + ".txt");                                        % Nombre del nuevo archivo
    [Reg,dt] = readAT2writeTxt(AT2fileDir,AT2fileName,newTxtFileDir,newTxtFileName);
    
    % Escribir nuevo registro en GM Data.txt    
%     writeGMData(newTxtFileName,dt,GMDataDir,GMDataName)                         % Escribir 'newTxtFileName   dt' en archivo GM Data.txt
    % No pude modificar la cantidad de veces que me lo escribe, cuidado con
    % correr varias veces esta línea ya que  volvería a escribir varias veces
    % en archivo GM Data.txt
end