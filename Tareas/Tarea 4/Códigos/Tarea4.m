%% Tarea 4: Análisis Dinámico Incremental de Sistemas de MGDL
% Contreras - Sanguinetti
% Ingeniería Sísmica Avanzada - USM 2022

% Comentarios
% - Al obtener los resultados desde THAMDOF, guardarlos en una carpeta con
% los mismos nombres que muestra THAMDOF (i.e: Export -> Displacement =>
% guardar como Displacement.xlsx)
% - En los inputs de texto usar chars, no strings  (i.e: usar apostrofe
% ('), no comillas ("))
% - Revisar todos los vectores y matrices definidas como inputs en las preguntas
% correspondientes
% - También revisar el límite para considerar qué es colapso (10^20 en
% general debería funcionar, revisar todos, depende la unidad de medida)

%% Inicializar
clear variables
close all
clc

%% Inputs
ResultsDir = 'THAMDOF_results';                                             % Carpeta donde están los resultados
cant_pisos = 9;                                                             % Cantidad de pisos
cant_registros = 20;                                                        % Cantidad de registros
IMs = [0.1; 0.4; 0.6; 0.8]; % g
cant_franjas = length(IMs);                                                           % Cantidad de franjas
g = 9.81;

% Nombre de archivos
nA_IDR = 'IDR.xlsx';                                                        % Interstory Drift Ratio
nA_PFA = 'Total Acceleration.xlsx';                                         % Peak Floor Aceleration
nA_rIDR = 'Residual IDR.xlsx';                                              % Residual Interstory Drift Ratio
nA_Disp = 'Displacement.xlsx';
% Alphabet
alphabet = string('A':'Z');
alphab = strings(cant_registros*cant_franjas,1);
counter = 0;
a_l = 26;           % largo del abecedario
for i = 1:cant_registros*cant_franjas    % llega hasta AZ no mas
    if i <= 1*a_l
        alphab(i) = alphabet{1}(i);
    end
    if and(i > 1*a_l, i <= 2*a_l)
        alphab(i) = string(alphabet{1}(1)) + string(alphabet{1}(i-1*a_l));
    end
    if and(i > 2*a_l, i <= 3*a_l)
        alphab(i) = string(alphabet{1}(2)) + string(alphabet{1}(i-2*a_l));
    end
    if and(i > 3*a_l, i <= 4*a_l)
        alphab(i) = string(alphabet{1}(3)) + string(alphabet{1}(i-3*a_l));
    end
end

%% Load Data
files = dir(ResultsDir);                                                    % Ver carpeta
files(1:1:2) = [];

for i = 1:length(files)
    if isequal(files(i).name, nA_IDR)
        IDRid = i;
        continue
    elseif isequal(files(i).name, nA_PFA)
        PFAid = i;
        continue
    elseif isequal(files(i).name,nA_rIDR)
        rIDRid = i;
        continue
    elseif isequal(files(i).name,nA_Disp)
        Dispid = i;
        continue
    end
end

%% P1
% Generar un gráfico de la mediana del EDP como función de IM condicionado
% a que no hay colapso para:

% a. La máxima razón de derivas de entrepiso delta del primer piso
% b. La máxima razón de derivas de entrepiso delta del noveno piso
% c. La aceleración máxima de piso (PFA) en el techo
% d. El máximo valor de delta en cualquier piso
% e. El máximo valor residual de delta en cualquier piso


%% P1 a.  (y b)
% Importar EDP = IDR [%] desde THAMDOF para todos los pisos
% range = convertStringsToChars(alphab(1) + "2:" + alphab(end) + string(cant_pisos+1)); 
IDR_data = importdata(convertStringsToChars(string(ResultsDir) + "\" + string(files(IDRid).name)));

% Obtener EDP para todos los pisos
F_EDP_pos = IDR_data.data.IDR(1:cant_pisos,:);                              % Floor EDP_positivos
F_EDP_neg = IDR_data.data.IDR((cant_pisos+4):(2*cant_pisos+3),:);           % Floor EDP negativos
F_EDP = max(abs(F_EDP_pos),abs(F_EDP_neg));                                 % Floor EDP (el máximo entre los positivos y los negativos)

% Mediana dado no colapso para todos los pisos
EDP_median = zeros(length(IMs),cant_pisos);
EDP_stdln = zeros(length(IMs),cant_pisos);

for i = 1:cant_pisos
    figure
    hold on
    leyenda_string = strings(length(IMs),1);
    for j = 1:length(IMs)
        vect = j+1 + 0:cant_franjas:cant_registros*cant_franjas;                       % j para ir por cada IM, el +1 es para empezar desde la segunda y el 
        F_EDP_vect = F_EDP(i,vect);
        F_EDP_vect(F_EDP_vect > 10^20) = [];                                % Límite para definir colapso para TODOS LOS pisos de 10^20, puede ser menos y hasta se puede poner como input
        EDP_median_vect = F_EDP(i,vect);
        EDP_median_vect(EDP_median_vect > 10^2) = [];
        EDP_median(j,i) = geomean(EDP_median_vect);
        EDP_stdln(j,i) = std(log(EDP_median_vect));
        plot(IMs(j,1), F_EDP_vect.'*100,'.','color','#909090')
%         leyenda_string(j) = "Datos";
    end
    plot([0; IMs], [0 ; EDP_median(:,i)*100],'-o','color','r','LineWidth',1.5)
    plot([0; IMs], [0 ; exp(log(EDP_median(:,i)) + EDP_stdln(:,i))*100],'--','color','r')
    plot([0; IMs], [0 ; exp(log(EDP_median(:,i)) - EDP_stdln(:,i))*100],'--','color','r')
%     legend([leyenda_string;"Mediana"])
    hold off
    grid on
    title('Median IDR',convertStringsToChars("Piso " + string(i)));
    xlabel('IM = Sa(T1) [g]')
    ylabel('EDP = IDR [%]')
end

% Todos juntos
% Extraer Data para otras preguntas
EDP_stdln_IDR = EDP_stdln;                                                  % por piso

%% P1 b.
% Listo en a

%% P1 c.
% Importar EDP = PFA [g] desde THAMDOF para todos los pisos
% range = convertStringsToChars(alphab(1) + "2:" + alphab(end) + string(cant_pisos+1)); 
PFA_data = importdata(convertStringsToChars(string(ResultsDir) + "\" + string(files(PFAid).name)));

% Obtener EDP para todos los pisos
F_EDP_pos = PFA_data.data.a_t(1:cant_pisos,:);                              % Floor EDP_positivos
F_EDP_neg = PFA_data.data.a_t((cant_pisos+4):(2*cant_pisos+3),:);           % Floor EDP negativos
F_EDP = max(abs(F_EDP_pos),abs(F_EDP_neg));                                 % Floor EDP (el máximo entre los positivos y los negativos), matriz de EDP(piso,RegistroEscalado), piso = 9, RegistroEscalado = 80

% Mediana dado no colapso para todos los pisos
EDP_median = zeros(length(IMs),cant_pisos);
EDP_stdln = zeros(length(IMs),cant_pisos);

for i = 1:cant_pisos
    figure
    hold on
%     leyenda_string = strings(length(IMs),1);
    for j = 1:length(IMs)
        vect = j+1 + 0:cant_franjas:cant_registros*cant_franjas;
        F_EDP_vect = F_EDP(i,vect);
        F_EDP_vect(F_EDP_vect > 10^20) = [];                                % Límite para definir colapso para TODOS LOS pisos de 10^20, puede ser menos y hasta se puede poner como input
        EDP_median_vect = F_EDP(i,vect);
        EDP_median_vect(EDP_median_vect > 10^2) = [];
        EDP_median(j,i) = geomean(EDP_median_vect);
        EDP_stdln(j,i) = std(log(EDP_median_vect));
        plot(IMs(j,1), F_EDP_vect.','.','color','#909090')
%         leyenda_string(j) = "Datos";
    end
    plot([0; IMs], [0 ; EDP_median(:,i)],'-o','color','r','LineWidth',1.5)
    plot([0; IMs], [0 ; exp(log(EDP_median(:,i)) + EDP_stdln(:,i))],'--','color','r')
    plot([0; IMs], [0 ; exp(log(EDP_median(:,i)) - EDP_stdln(:,i))],'--','color','r')
%     legend([leyenda_string;"Mediana"])
    hold off
    grid on
    title('Median PFA',convertStringsToChars("Piso " + string(i)));
    xlabel('IM = Sa(T1) [g]')
    ylabel('EDP = PFA [g]')
end

% Extraer Data para otras preguntas
EDP_stdln_PFA = EDP_stdln;                                                  % por piso

%% P1 d.
% Importar EDP = IDR [%] desde THAMDOF para todos los pisos
% range = convertStringsToChars(alphab(1) + "2:" + alphab(end) + string(cant_pisos+1)); 
IDR_data = importdata(convertStringsToChars(string(ResultsDir) + "\" + string(files(IDRid).name)));

% Obtener EDP para todos los pisos
F_EDP_pos = IDR_data.data.IDR(1:cant_pisos,(2:(cant_registros*cant_franjas+1)));                              % Floor EDP_positivos
F_EDP_neg = IDR_data.data.IDR((cant_pisos+4):(2*cant_pisos+3),(2:(cant_registros*cant_franjas+1)));           % Floor EDP negativos
F_EDP = max(abs(F_EDP_pos),abs(F_EDP_neg));                                 % Floor EDP (el máximo entre los positivos y los negativos), matriz de EDP(piso,RegistroEscalado), piso = 9, RegistroEscalado = 80

% Mediana dado no colapso para el máximo IDR de todos los pisos
EDP_median = zeros(length(IMs),1);
EDP_stdln = zeros(length(IMs),1);

figure
hold on
for j = 1:length(IMs)
    vect = j + 0:cant_franjas:cant_registros*cant_franjas;                         
    F_EDP_vect = F_EDP(:,vect);
    F_EDP_vect_new = zeros(cant_registros,1);
    for r = 1:cant_registros                                                % Para cada registro
        clear F_EDP_nv % new_vect
        F_EDP_nv = F_EDP_vect(:,r);                                         % Guardo el vector específico del registro
        F_EDP_nv(F_EDP_nv > 10^2) = [];                                     % Le saco el colapso
        F_EDP_max_floor = max(F_EDP_nv);                                    % Innecesario xd, pero ya lo hice
        if isempty(F_EDP_max_floor)
            continue
        else
            F_EDP_vect_new(r,1) = F_EDP_max_floor;
        end
    end
    F_EDP_vect_new(F_EDP_vect_new == 0) = [];
%     F_EDP_vect(F_EDP_vect > 10^30) = [];                                % Límite para definir colapso para TODOS LOS pisos de 10^20, puede ser menos y hasta se puede poner como input
%     F_EDP_vect = max(F_EDP_vect);
    EDP_median(j,1) = geomean(F_EDP_vect_new);
    EDP_stdln(j,1) = std(log(F_EDP_vect_new));

    plot(IMs(j,1), F_EDP_vect_new.'*100,'.','color','#909090')
end
plot([0; IMs], [0 ; EDP_median(:,1)*100],'-o','color','r','LineWidth',1.5)
plot([0; IMs], [0 ; exp(log(EDP_median(:,1)) + EDP_stdln(:,1))*100],'--','color','r')
plot([0; IMs], [0 ; exp(log(EDP_median(:,1)) - EDP_stdln(:,1))*100],'--','color','r')
hold off
grid on
title('Median maxFloor IDR');
xlabel('IM = Sa(T1) [g]')
ylabel('EDP = maxFloor IDR [%]')

% Extraer Data para otras preguntas
EDP_stdln_maxFloorIDR = EDP_stdln;                                          % máximo del piso

%% P1 e.
% Importar EDP = ResIDR [%] desde THAMDOF para todos los pisos
% range = convertStringsToChars(alphab(1) + "2:" + alphab(end) + string(cant_pisos+1));
rIDR_data = importdata(convertStringsToChars(string(ResultsDir) + "\" + string(files(rIDRid).name)));

% Obtener EDP para todos los pisos
F_EDP = abs(rIDR_data.data.ResIDR(1:cant_pisos,(2:(cant_registros*cant_franjas+1))));                              % Floor EDP_positivos
% F_EDP_neg = IDR_data.data.ResIDR((cant_pisos+4):(2*cant_pisos+3),(2:(cant_registros*cant_franjas+1)));           % Floor EDP negativos
% F_EDP = max(abs(F_EDP_pos),abs(F_EDP_neg));                                 % Floor EDP (el máximo entre los positivos y los negativos), matriz de EDP(piso,RegistroEscalado), piso = 9, RegistroEscalado = 80

% Mediana dado no colapso para el máximo IDR de todos los pisos
EDP_median = zeros(length(IMs),1);
EDP_stdln = zeros(length(IMs),1);

figure
hold on
for j = 1:length(IMs)
    vect = j + 0:cant_franjas:cant_registros*cant_franjas;                         
    F_EDP_vect = F_EDP(:,vect);
    F_EDP_vect_new = zeros(cant_registros,1);
    for r = 1:cant_registros                                                % Para cada registro
        clear F_EDP_nv % new_vect
        F_EDP_nv = F_EDP_vect(:,r);                                         % Guardo el vector específico del registro
        F_EDP_nv(F_EDP_nv > 10^10) = [];                                     % Le saco el colapso
        F_EDP_max_floor = max(F_EDP_nv);                                    % Innecesario xd, pero ya lo hice
        if isempty(F_EDP_max_floor)
            continue
        else
            F_EDP_vect_new(r,1) = F_EDP_max_floor;
        end
    end
    F_EDP_vect_new(F_EDP_vect_new == 0) = [];
%     F_EDP_vect(F_EDP_vect > 10^30) = [];                                % Límite para definir colapso para TODOS LOS pisos de 10^20, puede ser menos y hasta se puede poner como input
%     F_EDP_vect = max(F_EDP_vect);
    EDP_median(j,1) = geomean(F_EDP_vect_new);
    EDP_stdln(j,1) = std(log(F_EDP_vect_new));
    plot(IMs(j,1), F_EDP_vect_new.'*100,'.','color','#909090')
end
plot([0; IMs], [0 ; EDP_median(:,1)*100],'-o','color','r','LineWidth',1.5)
plot([0; IMs], [0 ; exp(log(EDP_median(:,1)) + EDP_stdln(:,1))*100],'--','color','r')
plot([0; IMs], [0 ; exp(log(EDP_median(:,1)) - EDP_stdln(:,1))*100],'--','color','r')
hold off
grid on
title('Median maxFloor ResIDR');
xlabel('IM = Sa(T1) [g]')
ylabel('EDP = maxFloor resIDR [%]')

% Extraer Data para otras preguntas
EDP_stdln_maxFloorresIDR = EDP_stdln;                                       % máximo del piso

%% P2
% Desviación de EDP dado no colapso para los mismos casos que en P1
% a
EDP_stdln = EDP_stdln_IDR;
figure
plot([0; IMs], [0 ; EDP_stdln(:,1)],'-o','color','r')
xlabel('IM = Sa(T1) [g]')
ylabel('EDP stdln')
title('Desviación estándar logarítmica','Máximo IDR Piso 1')
grid on

%b
EDP_stdln = EDP_stdln_IDR;
figure
plot([0; IMs], [0 ; EDP_stdln(:,9)],'-o','color','r')
xlabel('IM = Sa(T1) [g]')
ylabel('EDP stdln')
title('Desviación estándar logarítmica','Máximo IDR Piso 9')
grid on

% c
EDP_stdln = EDP_stdln_PFA;
figure
plot([0; IMs], [0 ; EDP_stdln(:,9)],'-o','color','r')
xlabel('IM = Sa(T1) [g]')
ylabel('EDP stdln')
title('Desviación estándar logarítmica','Máximo PFA Techo')
grid on

% d
EDP_stdln = EDP_stdln_maxFloorIDR;
figure
plot([0; IMs], [0 ; EDP_stdln(:,1)],'-o','color','r')
xlabel('IM = Sa(T1) [g]')
ylabel('EDP stdln')
title('Desviación estándar logarítmica','Máximo IDR en cualquier piso')
grid on

% e
EDP_stdln = EDP_stdln_maxFloorresIDR;
figure
plot([0; IMs], [0 ; EDP_stdln(:,1)],'-o','color','r')
xlabel('IM = Sa(T1) [g]')
ylabel('EDP stdln')
title('Desviación estándar logarítmica','Máximo ResIDR en cualquier piso')
grid on

%% P3   (Martina)
% Generar gráficos que muestren la distribución en altura de la media
% geométrica de los siguientes EDP, dado no colapso

% Eje x -> Media geométrica para cada piso para cada IM
% EJe y -> Número del piso

pisos = 1:1:cant_pisos;

% a. La máxima razón de derivas de piso (delta)




% b. El máximo valor residual de delta

% c. PFA

%% P4
% Determinar Curva de Fragilidad de Colapso (CFC) con método de mínimos
% cuadrados y método de máxima verosimilitud y comparar

% El colapso quedará definido cuando el desplazamiento, de cualquier piso,
% exceda el desplazamiento de resistencia 0 (du - desplazamiento último)

Disp_data = importdata(convertStringsToChars(string(ResultsDir) + "\" + string(files(Dispid).name)));
F_EDP_pos = Disp_data.data.u(1:cant_pisos,(2:(cant_registros*cant_franjas+1)));                              % Floor EDP_positivos
F_EDP_neg = Disp_data.data.u((cant_pisos+4):(2*cant_pisos+3),(2:(cant_registros*cant_franjas+1)));           % Floor EDP negativos
F_EDP = max(abs(F_EDP_pos),abs(F_EDP_neg)); 

n_collapse = zeros(cant_franjas,1);
fraccion = zeros(cant_franjas,1);
for i = 1:cant_franjas
    vect = i + 0:cant_franjas:cant_franjas*cant_registros;
    F_EDP_vect = F_EDP(:,vect);
    for r = 1:cant_registros
        for j = 1:cant_pisos
            if F_EDP_vect(j,r) > 10^20
                n_collapse(i,1) = n_collapse(i,1) + 1;
                break
            end
        end
    end
    fraccion(i,1) = n_collapse(i,1)./cant_registros;
end

% Rango de Mu y Sigma para realizar la "optimización" del fittind de la
% curva CFC por método de minimos cuadrados y 

% Inputs
mus = -5:0.01:5;
sigmas = 0.001:0.001:1;

% Método de los mínimos cuadrados
E = zeros(cant_franjas,1);
mc = +inf;
for m = 1:length(mus)
    for s = 1:length(sigmas)
        for i = 1:cant_franjas
            E(i) = abs(fraccion(i) - normcdf((log(IMs(i)) - mus(m))./sigmas(s)))^2;
        end
        if sum(E) < mc
            mc = sum(E);
            sigma_mc = sigmas(s);
            mu_mc = mus(m);
        end
    end
end

% Método de máxima verosimilitud
mv = -inf;
L = zeros(cant_franjas,1);
for m = 1:length(mus)
    for s = 1:length(sigmas)
        for i = 1:cant_franjas
            combinatoria = nchoosek(cant_registros,n_collapse(i));          % n: cantidad registros, k: cantidad de colapsos
            pc = normcdf((log(IMs(i)) - mus(m))/sigmas(s))^n_collapse(i);
            pnc = (1-normcdf((log(IMs(i)) - mus(m))/sigmas(s)))^(cant_registros-n_collapse(i));
            L(i,1) = combinatoria*pc*pnc;
        end
        if exp(sum(L)) > mv
            mv = exp(sum(L));
            mu_mv = mus(m);
            sigma_mv = sigmas(s);
        end
    end
end
tabla = table();
tabla.par = ["mu";"sigma"];
tabla.mc = [mu_mc; sigma_mc];
tabla.mv = [mu_mv; sigma_mv];
disp(tabla)
clear tabla

% Generar figura para corrobarar método
% Inputs
IM_range = 0.01:0.01:3; % g
% Figura
figure
plot(IMs,fraccion,'o','Color','r','linewidth',1.5)
hold on
plot(IM_range,normcdf((log(IM_range)-mu_mc)/sigma_mc),'color','b','Linewidth',1.5)
plot(IM_range,normcdf((log(IM_range)-mu_mv)/sigma_mv),'color','g','Linewidth',1.5)
hold off
xlabel('IM: Sa(T1) [g]')
ylabel('P(C|IM=im)')
title('Curva de Fragilidad de Colapso')
legend('Data','Mínimos Cuadrados','Máxima Verosimilitud')
grid on

%% P5
% Obtener la curva de amenaza sísmica desde el sitio web de USGS, suelo
% clase D, coord(37.785,-122.44), realizar ajuster polinomial de tercer
% grado, graficar e incluir R2

% Data USGS
IM_USGS = [ 0.0025;0.0045;0.0075;0.0113;0.0169;0.0253;0.038;0.057;0.0854;0.128;0.192;0.288;0.432;0.649;0.973;1.46;2.19;3.28;4.92;7.38];
lambda_USGS = [0.56700011;0.3744033;0.24882876;0.17413627;0.12000236;0.081077782;0.053541436;0.03454773;0.02167907;0.013171584;0.0076841789;0.0042117965;0.0020870034;8.8886321E-4;3.1228863E-4;8.4600432E-5;1.5386503E-5;1.3194825E-6;1.5856989E-8;6.4405085E-13];

% Ajuste polinomial

%% P6
% Calcular valores

% a


% b


% c


% d


% e

