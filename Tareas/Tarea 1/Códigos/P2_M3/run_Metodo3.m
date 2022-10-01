%% Método 3: Ajuste Utilizando el IM alternativo Sa_{avg}
% Eads, Miranda & Lignos (2015)
% Selección de Registros por Sa_avg

% Ingeniería Sísmica Avanzada
% Contreras Sanguinetti

% Parámetros
% T1            -> Periodo T para c1*T,c2*T....cN*T
% c1, cN        -> Coeficientes para el rango para calcular Sa_avg (inicial y final)
% Tn            -> Vector de Periodos para el análisis
% n_interpol    -> Cantidad de interpolaciones (para mejorar presición) para Sa1, Sa2 y Periods


%% Inicializar
clear variables
close all
clc

%% Load Databases
load('../PI_3_4y5/median_CMS.mat')
load('../P2_M2/CS_Selection-master/CS_Selection-master/Databases/NGA_W2_meta_data.mat')
Periods = Periods';                                                         % Para que quede en vector

%% Calculo Savg
T1 = 1;                                                                     % sec
c1 = 0.2;
cN = 3.0;
Tn_inicial = c1*T1;
Tn_final = cN*T1;

Ti_init = 0.05;                                                             % Mínimo Ti
Ti_step = 0.01;                                                             % Paso para los periodos
Ti_final = 5;                                                               % Máximo Ti

Tn = (Ti_init:Ti_step:Ti_final)';                                           % Rango de periodos

Sa_avg_vect = zeros(length(Tn));
median_CMS_Tn = [Tn median_CMS];

for i = 1:length(Tn)
    if Tn(i) >= Tn_inicial  && Tn(i) <= Tn_final
        Sa_avg_vect(i) = median_CMS(i);
    end
end

% Los valores que entran dentro del rango de periodos
Sa_avg_vect = nonzeros(Sa_avg_vect);        % Quitamos todos los ceros

% Por lo tanto Sa_avg
Sa_avg = prod(Sa_avg_vect)^(1/length(Sa_avg_vect));

% Sa_1 -> Componente 1 de 21539 registros
% Sa_2 -> Componente 2 de los mismos 21539 registros
% Periods -> Vector de periodos

Sa = [Sa_1;Sa_2];
[m,n] = size(Sa);
%% 5) Calcular Sa_avg(T=1sec) de cada espectro
n_interpol = 1;                                                             % Creo que función interp no funciona como creo que funciona

Periods_new = interp(Periods,n_interpol);
Sa_avg_registros = zeros(m,1);
Sa_positions = zeros(length(Periods_new),1);

for j = 1:length(Periods_new)
    if Periods_new(j) >= Tn_inicial && Periods_new(j) <= Tn_final
        Sa_positions(j) = j;
    end
end
Sa_positions = nonzeros(Sa_positions);
Sa_avg_5 = zeros(m,1);
for j = 1:m
    Sa_vals = interp(Sa(j,:),n_interpol)';
    Sa_avg_5(j) = prod(Sa_vals(min(Sa_positions):max(Sa_positions)))^(1/length(Sa_vals(min(Sa_positions):max(Sa_positions))));
end

%% Seleccionar 20 registros
n_reg = 20;                                                                 % Cantidad de registros a seleccionar
c = 0;
n_reg_original = n_reg;
length_registros = length(Sa_1);
componente = ones(n_reg,1);                                                 % Todos son del 1, a menos que corresponda al dos

while c ~= 1
    [Error,Reg_index] = mink(abs(Sa_avg_5-Sa_avg),n_reg);
    for i = 1:n_reg
        if Reg_index(i) > length_registros                                  % Si se pasa del largo de registro, entonces es la segunda componente
            Reg_index(i) = Reg_index(i) - length_registros; 
            componente(i) = 2;                                              % Todos son del 1, a menos que corresponda al dos
        end
    end
%     Reg_index = unique(Reg_index);
    % Si se repiten algunos, unique los borra y abría que buscar otros registros
    if length(Reg_index) < n_reg_original
        n_reg = n_reg + 1; 
    elseif length(Reg_index) == n_reg_original
        c = 1;
    end
end
fprintf('Los registros seleccionados son: %.0f\npasando por %.0f \n\n',n_reg_original,n_reg)
Reg_index = sort(Reg_index);

%% Mostrar resultados
vs30_reg = soil_Vs30(Reg_index);
rjb_reg = Rjb(Reg_index);
magnitude_reg = magnitude(Reg_index);

tabla = table();
tabla.Registro = (1:1:n_reg)';
tabla.Secuencia = Reg_index;
tabla.Componente = componente;
tabla.Magnitud = magnitude_reg;
tabla.Distancia = rjb_reg;
tabla.vs30 = vs30_reg;
disp(tabla)
clear tabla

disp('Creo que no funciona bien, ya que no estoy seguro si interp() genera "n" puntos')
disp('entre dos valores del vector, creo que hay que usar interp1() (como dice la tarea xd)')
disp('interp(Periods,Sa(j,:),0:0.01:Periods(end)')


