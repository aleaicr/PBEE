function [] = writeGMData(newRegName,newReg_dt,GMDataDir,GMDataName)
% Escribir nuevo registro en arhchivo GM Data.txt con formato:
% newRegName.txt dt

% Alexis Contreras - Martina Sanguinetti
% Ingeniería Sísmica Avanzada - USM

% Inputs:
% newRegName        -> Nombre del nuevo registro (con .txt al final)
% newReg_dt         -> Paso temporal (dt) del nuevo registro
% GMADataDir        -> Dirección del arhicov GM Data.txt
% GMDataName        -> Por defecto poner lo que entrega PEER, GM Data.txt pero especificar si es otro

% Outputs:
%

%%
% Abrir para leer
fileID = fopen([GMDataDir '\' GMDataName],'r');                             % Abrir archivo para lectura

% Contar la cantidad de líneas
fLine = fgetl(fileID);
j = 0;                                                                      % Supongamos que no está
while ischar(fLine)
    fLine = fgetl(fileID);
    if isequal(fLine,strcat([char(newRegName) '    ' char(string(round(newReg_dt,3)))])) % Contar las que no sean iguales
        j = 1;                                                              % Si está entonces no debemos escribirlo otra vez
    end
end
fclose(fileID);                                                             % Se cierra para lectura

% Abrir para escribir
fileID = fopen([GMDataDir '\' GMDataName],'a+');                            % Abrir archivo para escritura (escribir al final)

% Si ya existe (j = 1), entonces no se añade, si no existe (else) entonces se
% añade
if j == 1
    % no se añade
elseif j == 0
    fprintf(fileID,'%s\t%.2f\r',newRegName,newReg_dt);                    % Escribir al final
end

% Cerrar archivo
fclose(fileID);                                                             % Cerrar archivo
end
