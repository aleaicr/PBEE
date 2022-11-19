function []= SingleGMThamdofFormat(GMFolder,GMDataName,T,xi,beta_newmark,t_extra,GMTHAMDOFName_collapse)
% Contreras - Sanguinetti
% Ingeniería Sísmica Avanzada - USM 2022

% Leer registros en txt y traspasarlos a formato THAMDOF (AT2? y en .csv)

% Inputs
% GMFolder: Carpeta donde se encuentran los registros (y GM Data.txt)
%           Los registros deben venir en [g], si vienen en m/s2, en linea
%           54, dividir todo por g = 9.81
% GMDataName: Nombre del archivo "GM Data.txt" donde se guardan los pasos
%             temporales de cada registro (debe estar ordenado)
% T: Periodo fundamental de la estructura [sec]
% beta_newmark: factor del método de Newmark (ej: beta_newmark = 1/4)
% t_extra: Tiempo extra de simulación (cantidad de filas extra en matriz matrix_csv_
% GMTHAMDOFName_collapse: Nombre del archivo THAMDOF para análisis del riesgo de colapso (requisito que sea con .csv si o si)

% Outputs
% Ninguno, pero Genera un archivo csv con Todos los registros de la carpeta
% en formato THAMDOF

% Comentarios: 
% - Hasta el momento, solo se consideran IM = Sa(T1) aceleración
% espectral del primer modo (newmarkbeta)
% - Los -3 y +3 que se ven son por la posición con la que se ven los
% archivos al utilizar la función dir(), los registros parten del 3


%% Leer archivos de la carpeta
files = dir(GMFolder);
% los nombres de los archivos parten desde el 4 
% => files(4).name = nombre del primer archivo
% files(5).name = nombre del segundo archivo
a = struct();
max_length_GM = 0;
for i = 4:length(files)
    a(i-3).name = files(i).name;
    a(i-3).GM = load([GMFolder '\' files(i).name]);                         % Guardar registro 
    a(i-3).L = length(a(i-3).GM);                                           % Guardar Largo del registro
    if a(i-3).L > max_length_GM
        max_length_GM = a(i-3).L;                                           % Guardar máximo largo de registro
    end
end

%% Leer Archivo GM Data.txt
fileID = fopen([GMFolder '\' GMDataName],'r');                              % Abrir archivo para lectura
fLine = fgetl(fileID);                                                      % Char de la primera linea del GM Data.txt
while ischar(fLine)
    string_data_line = split(string(fLine));                                % Separamos en dos (por el tab)
    for i = 4:length(files)
        if string_data_line(1) == string(files(i).name)                     % Solo si es el mismo registro
            a(i-3).dt = double(string_data_line(2));                        % Guardamos dt
            % Aprovechar de obtener la aceleración espectral
            a(i-3).SaT1 = Newmark_Lineal_Sa(beta_newmark,xi,a(i-3).dt,0,0,a(i-3).GM,T); % Guardar aceleración espectral del primer modo, GM está en m/s
        end
    end
    fLine = fgetl(fileID);                                                  % Siguiente línea del GM Data.txt
end
fclose(fileID);

%% Formato THAMDOF
% GM; length; dt; factor de escala; registro
IM_length = 1;
matrix_csv_ = zeros(max_length_GM+3+max(t_extra),(length(files)-3)*IM_length);
for i = 4:length(files)                                                     % Cantidad de registros
    for j = 1:IM_length                                                     % Cantidad de franjas que uno quiere
        for r = 1:4                                                         % Cuatro primeras columnas tienen cosas distintas
            if r == 1
%                 matrix_csv(r,IM_length*((i-3)-1)+j) = i-3;                  % Nombre del registro 1
            elseif r == 2
                matrix_csv_(r,IM_length*((i-3)-1)+j) = a(i-3).L + t_extra(i-3);            % Largo del registro
            elseif r == 3
                matrix_csv_(r,IM_length*((i-3)-1)+j) = a(i-3).dt;           % Paso temporal del sampling del registro (muestreo?)
            elseif r == 4
                matrix_csv_(r,IM_length*((i-3)-1)+j) = 1;           % Factor de Escala para las IM que se quieren
            end
        end
        for k = 5:length(a(i-3).GM)+4
            matrix_csv_(k,IM_length*((i-3)-1)+j) = a(i-3).GM(k-4,1);        % Registro
        end
    end
end

%% Escribir matriz en CSV
matrix_csv_(1,:) = [];
titles = convertStringsToChars(["GM" + string([1:1:(length(files)-3)])]).';
A = cell(1,IM_length*(length(files)-3));
for i = 4:length(files)
    for j = 1:IM_length
        A{IM_length*(i-3-1)+j} = titles(i-3,j);
    end
end
table_csv = array2table(matrix_csv_);
table_csv.Properties.VariableNames(1:IM_length*(length(files)-3)) = string(A);
writetable(table_csv,GMTHAMDOFName_collapse)
fprintf('Se ha creado el archivo %s\n\n',GMTHAMDOFName_collapse)
end

%% Método de Newmark-beta utilizado para determinar el espectro, descomentar ctrl+shift+R si no se tiene programado

% function [Sa] = Newmark_Lineal_Sa(beta,xi,dt,ui,udi,uddg,Tn)
% % Respuesta estructural por el método de Newmkar-beta
% 
% % Contreras - Adasme 
% % Ingeniería Sísmica - USM 2022
% 
% % Inputs
% % beta  -> factor del método de Newmark (ej beta = 0 o 1/4)
% % xi    -> Razón de amortiguamiento de cada modo a analizar  (ej: 0.05, no 5)
% % dt    -> [sec] Rango de partición temporal del registro
% % ui    -> [m] Desplazamiento inicial (Condición inicial)
% % udi   -> [m/s] Velocidad Inicial (Condición inicial)
% % uddg  -> [m/sw] Registro de aceleración del suelo
% % Tn    -> [sec] Vector de Periodos de cada modo a analizar
% 
% % Outputs
% % Sd    -> [m] Espectro de desplazamientos
% % Sv    -> [m/s] Espectro de velocidades
% % Sa    -> [m/s2] Espectro de aceleraciones
% % PSv   -> [m/s] Pseudo-Espectro de velocidades
% % PSa   -> [m/s2] Pseudo-Espectro de aceleraciones
% 
% % Comentario
% % Si beta = 0 -> verificar estabilidad del método omega*delta_t < 2
% % si beta = 1/4 -> método es incondicionalmente establa
% 
% gamma = 0.5;                                                                % Método según cómo se toma la aceleración
% uddg_length = length(uddg);
% Tn_length = length(Tn);
% 
% %% Inicializar vectores
% % Para no reescribir sobre memoria cada vez que se modifica el tamaño
% 
% % Inicialización de vector de desplazamientos, velocidad y aceleración
% u = zeros(uddg_length,1);
% ud = zeros(uddg_length,1);
% udd = zeros(uddg_length,1);
% 
% % Inicialización de espectro de aeleración
% Sa = zeros(Tn_length,1);
% 
% % Condiciones Iniciales
% u(1,1) = ui;
% ud(1,1) = udi;
% 
% for j = 1:Tn_length                                                         % j es el periodo a evaluar
%     if Tn(1,j) == 0
%         for i = 1:size(uddg,1)                                              % i es el número de dato del registro
%             udd(i,1) = uddg(i,1);
%         end
%     else
%         wn = 2*pi/Tn(1,j);                                                  % Frecuencia natural del sistema
%         udd(1,1) = uddg(1,1) - 2*xi*wn*udi - wn^2*ui;
%         a1 = 1/(beta*dt^2) + 2*xi*wn*gamma/(beta*dt);
%         a2 = 1/(beta*dt) + 2*xi*wn*(gamma/beta-1);
%         a3 = (1/(2*beta)-1) + 2*xi*wn*dt*(gamma/(2*beta)-1);
%         k_ton = a1 + wn^2;
%         for i = 2:size(uddg,1)
%             p_ton = -uddg(i,1) + a1*u(i-1,1) + a2*ud(i-1,1) + a3*udd(i-1,1);
%             u(i,1) = p_ton/k_ton;
%             ud(i,1) = gamma/(beta*dt)*(u(i,1)-u(i-1,1)) + (1-gamma/beta)*ud(i-1,1) + dt*(1-gamma/(2*beta))*udd(i-1,1);
%             udd(i,1) = (u(i,1)-u(i-1,1))/(beta*dt^2) - ud(i-1,1)/(beta*dt) - (1/(2*beta)-1)*udd(i-1,1);
%         end
%         udd_tot = udd + uddg;
%         udd_tot(1,1) = uddg(1,1);
%         Sa(j,1) = max(abs(udd_tot(:,1)));
%     end
%     u(:,1) = zeros(size(u,1),1);
%     ud(:,1) = zeros(size(ud,1),1);
%     udd(:,1) = zeros(size(udd,1),1);
%     u(1,1) = ui;
%     ud(1,1) = udi;
% end
% end