%% Leer Registro y Formato CSV THAMDOF
% Contreras - Sanguinetti
% Ingeniería Sísimica Avanzada, USM 2022

% Leer los registros en formato tareas (en txt) y transformarlo al formato
% THAMDOF
% https://github.com/pheresi/THAMDOF/blob/master/Sample%20Input%20Files/GMs.csv

% Comentarios
% - Cambiar la carpeta donde están los registros guardados

%% Inicializar
clear variables
close all
clc

%% Inputs
GMFolder = 'Registros';                                                  % Nombre de la carpeta donde se encuentran los registros (todos los registros deben estar en [g] -> las que da el profe están en [g])
cant_registros = 20;
GMDataName = 'GM Data.txt';                                                 % Nombre del archivo .txt con los dt de sampling de cada registro
IM_Sa_avg = [0.1; 0.5; 0.6; 0.7; 0.8; 0.9; 1; 1.1; 1.2; 1.4; 1.5; 1.9]; % g                                         % Franjas para realizar el análisis
T = 2.18; % sec                                                                % Periodo fundamental de la estructura
beta_newmark = 1/4;                                                         % Beta del método de Newmark-beta (dejarlo como 1/4 para que sea incondicionalmente estable)
xi = 0.05;                                                                  % Amortiguamiento para el espectro
t_extra = 3000*ones(cant_registros,1); % [filas]                            % Tiempo extra para el análisis, (n*ones() es el mismo tiempo extra para todos, especificar si se quiere t_extra distinto para un registro)
GMTHAMDOFName = "GMs.csv";
c1 = 0.2;
c2 = 3;

%% Generar archivos
% Registros con SF para que den todos las franjas (todos los registros para todas las franjas)
GMThamdofFormat(GMFolder,GMDataName,T,IM_Sa_avg,c1,c2,xi,beta_newmark,t_extra,GMTHAMDOFName);

% %% Leer archivos de la carpeta
% files = dir(GMFolder);
% % los nombres de los archivos parten desde el 4 
% % => files(4).name = nombre del primer archivo
% % files(5).name = nombre del segundo archivo ...
% a = struct();
% max_length_GM = 0;                                                          % Para guardar el registro con el máximo largo
% for i = 4:length(files)
%     a(i-3).name = files(i).name;
%     a(i-3).GM = load([GMFolder '\' files(i).name]);                         % Guardar registro 
%     a(i-3).L = length(a(i-3).GM);                                           % Guardar Largo del registro
%     if a(i-3).L > max_length_GM
%         max_length_GM = a(i-3).L;                                           % Guardar máximo largo de registro
%     end
% end
% 
% %% Leer Archivo GM Data.txt
% fileID = fopen([GMFolder '\' GMDataName],'r');                              % Abrir archivo para lectura
% fLine = fgetl(fileID);                                                      % Char de la primera linea del GM Data.txt
% while ischar(fLine)
%     string_data_line = split(string(fLine));                                % Separamos en dos (por el tab)
%     for i = 4:length(files)
%         if string_data_line(1) == string(files(i).name)                     % Solo si es el mismo registro
%             a(i-3).dt = double(string_data_line(2));                        % Guardamos dt
%             % Aprovechar de obtener la aceleración espectral
%             Periods = (c1*T:a(i-3).dt:c2*T).';
%             a(i-3).SaT1 = Newmark_Lineal_Sa(beta_newmark,xi,a(i-3).dt,0,0,a(i-3).GM,Periods); % Guardar aceleración espectral del primer modo, GM está en m/s
%             a(i-3).Sa_avg = computeSaAvg(a(i-3).SaT1,Periods,T,c1,c2,a(i-3).dt);
%         end
%     end
%     fLine = fgetl(fileID);                                                  % Siguiente línea del GM Data.txt
% end
% fclose(fileID);
% 
% %% Determinar factores de escala
% IM_length = length(IM_Sa_avg);                                                % Cantidad de franjas
% SF = zeros(IM_length,length(files)-3);                                      % (sf, registro)
% for i = 4:length(files)
%     for im = 1:IM_length
%         SF(im,i-3) = IM_Sa_avg(im,1)/a(i-3).Sa_avg;                             % SF = Sa_obj/Sa_original  (Sa espectro de aceleraciones lineal)
%     end    
% end
% 
% %% Formato THAMDOF
% % GM; length; dt; factor de escala; registro
% matrix_csv_ = zeros(max_length_GM+3+max(t_extra),(length(files)-3)*IM_length);
% for i = 4:length(files)                                                     % Cantidad de registros
%     for j = 1:IM_length                                                     % Cantidad de franjas que uno quiere
%         for r = 1:4                                                         % Cuatro primeras columnas tienen cosas distintas
%             if r == 1
% %                 matrix_csv(r,IM_length*((i-3)-1)+j) = i-3;                  % Nombre del registro 1
%             elseif r == 2
%                 matrix_csv_(r,IM_length*((i-3)-1)+j) = a(i-3).L + t_extra(i-3);            % Largo del registro
%             elseif r == 3
%                 matrix_csv_(r,IM_length*((i-3)-1)+j) = a(i-3).dt;           % Paso temporal del sampling del registro (muestreo?)
%             elseif r == 4
%                 matrix_csv_(r,IM_length*((i-3)-1)+j) = SF(j,i-3);           % Factor de Escala para las IM que se quieren
%             end
%         end
%         for k = 5:length(a(i-3).GM)+4
%             matrix_csv_(k,IM_length*((i-3)-1)+j) = a(i-3).GM(k-4,1);        % Registro
%         end
%     end
% end
% 
% %% Escribir matriz en CSV
% matrix_csv_(1,:) = [];
% titles = convertStringsToChars("GM" + string(1:1:(length(files)-3)) + " - IM" + string(1:1:IM_length).').';   % El título de la columna. Ej: "GM1 - IM3", donde GM1 es el primer registro e IM3 es la tercera franja
% A = cell(1,IM_length*(length(files)-3));
% % Quedan ordenados de la siguiente forma, cada registro, en su franja por columnas:
% % GM1 - IM1 | GM1 - IM2 | GM1 - IM3| GM1 - IM4 | GM2 -IM1 | GM2 - IM2 | ... | GM20 - IM3 | GM20 - IM4|
% for i = 4:length(files)
%     for j = 1:IM_length
%         A{IM_length*(i-3-1)+j} = titles(i-3,j);
%     end
% end
% table_csv = array2table(matrix_csv_);
% table_csv.Properties.VariableNames(1:IM_length*(length(files)-3)) = string(A);
% writetable(table_csv,GMTHAMDOFName)
% fprintf('Se ha creado el archivo %s\n\n',GMTHAMDOFName)