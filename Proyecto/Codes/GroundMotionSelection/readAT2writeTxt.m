function [Reg,dt] = readAT2writeTxt(AT2fileDir,AT2fileName,newTxtFileDir,newTxtFileName)
% Leer formato AT2, dejarlo como vector columna y crear archivo .txt que
% guarde el registro, si se añade un nuevo 'path', no es necesario incluir
% direcciones, solo el nombre y fijarse que no se repita

% Alexis Contreras - Martina Sanguinetti
% Ingeniería Sísmica Avanzada - USM

% Inputs
% AT2fileDir        -> Dirección del arhicov 
% AT2fileName       -> Nombre del archivo AT2 de PEER (char)
% newTxtFileDir     -> Dirección donde se quiere guardar nuevo archivo txt
% newTxtFileName    -> Nombre para guardar el archivo en txt

% Outputs
% Reg               -> Vector con aceleraciones del registro
% dt                -> Valor del paso temporal del muestreo

%% Leer dt
if isempty(AT2fileDir)                                                      % Está en la misma carpeta donde se corre la función
    fileID = fopen(AT2fileName,'r');                                        % Abrir modo lectura
else
    fileID = fopen([AT2fileDir '\' AT2fileName],'r');                       % Abrir modo lectura pero el archivo .AT2 está en otra carpeta
end
unoAUno = textscan(fileID, '%s');                                           % Escan a cada elemento como string
dt = str2double(unoAUno{1}(22));                                            % El elemento 22 tiene el valor de dt y se convierte en double
fclose(fileID);                                                             % Cerrar archivo

%% Leer registro, reordenar como vector y escribir txt
if isempty(AT2fileDir)                                                      % Está en la misma carpeta donde se corre la función
    A = dlmread(AT2fileName,'',4,0);                                        % Leer documento .AT2 si está en la misma carpeta
else
    A = dlmread([AT2fileDir '\' AT2fileName],'',4,0);                       % Leer documento .AT2 pero está en otra carpeta
end

Reg = reshape(A',size(A,1)*size(A,2),1);                                    % Reordenar matriz como vector  

if isempty(newTxtFileDir)   
    fID = fopen(newTxtFileName,'w');
    fprintf(fID,'%.10f\r\n',Reg);                                           % Escribir registro en nuevo .txt
    fclose(fID);
else
    fID = fopen([newTxtFileDir '\' newTxtFileName],'w');
    fprintf(fID,'%.10f\r\n',Reg);                                           % Escribir registro en nuevo .txt
    fclose(fID);
end
fprintf('Archivo AT2 tiene dt = %f \n',dt)
fprintf('Se ha creado el archivo %s \n',newTxtFileName)
end
