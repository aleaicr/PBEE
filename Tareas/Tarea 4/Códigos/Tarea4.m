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
cant_pisos = 9;                                                             % Cantidad de pisos del edificio
cant_registros = 20;                                                        % Cantidad de registros porcada IM (se considera que todas las IMs tienen la misma cantidad de registros)
IMs = [0.1; 0.4; 0.6; 0.8]; % g                                             % Valor IMs de las franjas del "IDA"
cant_franjas = length(IMs);                                                 % Cantidad de franjas
g = 9.81;

% Nombre de archivos a utilizar en la tarea
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
EDP_muln_PFA = EDP_median;                                                  % Media
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
EDP_muln_maxFloorIDR = EDP_median;                                          % Máximo del piso
EDP_stdln_maxFloorIDR = EDP_stdln;                                          % Máximo del piso

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

close all
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
        if exp(sum(log(L))) > mv
            mv = exp(sum(log(L)));  % máxima verosimilitud
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
more_IM = (0.057:0.01:1.5).';                                               % Mismo rango que los datos para el ajuste polinomial
IM_range = sort([more_IM; IMs]);                          % Vector para el cual se va a realizar la interpolación, notar que se agregaron los IMs de Interés

% Figuras
figure
plot(IMs,fraccion,'o','Color','r','linewidth',1.5)
xlabel('IM: Sa(T1) [g]')
ylabel('P(C|IM=im)')
xlim([0 IM_range(end)])
ylim([0 1])
title('Fracciones zj/nj')
legend('Data')
grid on

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
% Datos que otorga USGS para la amenaza sísmica

original_IM_USGS = [0.0025; 0.0045; 0.0075; 0.0113; 0.0169; 0.0253; 0.038; 0.057; 0.0854;0.128; 0.192; 0.288; 0.432; 0.649; 0.973; 1.46; 2.19; 3.28; 4.92; 7.38];
original_lambda_USGS = [0.56700011; 0.3744033; 0.24882876; 0.17413627; 0.12000236; 0.081077782; 0.053541436; 0.03454773; 0.02167907; 0.013171584; 0.0076841789; 0.0042117965; 0.0020870034; 8.8886321E-4; 3.1228863E-4; 8.4600432E-5; 1.5386503E-5; 1.3194825E-6; 1.5856989E-8; 6.4405085E-13];

% Modificando valores originales, para obtener un mejor ajuste
IM_USGS = [0.057; 0.0854;0.128; 0.192; 0.288; 0.432; 0.649; 0.973; 1.46; 2.19; 3.28];
lambda_USGS = [0.03454773; 0.02167907; 0.013171584; 0.0076841789; 0.0042117965; 0.0020870034; 8.8886321E-4; 3.1228863E-4; 8.4600432E-5; 1.5386503E-5; 1.3194825E-6];

% Ajuste polinomial ------
% Determinar coeficientes
[P,S] = polyfit(log(IM_USGS),log(lambda_USGS),3);

% Ajustar curva
% El IM_range a utilizar para calcular los puntos serán los mismos que en
% la parte anterior
more_IM = (IM_USGS(1,1):0.01:3).';                                               % Mismo rango que los datos para el ajuste polinomial
IM_range = sort([more_IM; IMs]);                          % Vector para el cual se va a realizar la interpolación, notar que se agregaron los IMs de Interés
lambda_poly = exp(P(4)*ones(length(IM_range),1) + P(3)*log(IM_range) + P(2)*log(IM_range).^2 + P(1)*log(IM_range).^3); 
R_square = 1 - S.normr^2 / norm(log(lambda_USGS)-mean(log(lambda_USGS)))^2; % Valor de R^2

figure
loglog(original_IM_USGS,original_lambda_USGS,'o','color','#076F51')
hold on
loglog(IM_USGS,lambda_USGS,'o','color','r')
loglog(IM_range,lambda_poly,'.-','color','b')
hold off
xlabel('IM: Sa(T1) [g]')
ylabel('\lambda_{IM} [1/yr]')
grid on
legend('Datos USGS','Datos USGS Utilizados para el ajuste','Ajuste polinomial tercer orden')
text(10^0,10^-2,"R2 = " + string(R_square))
title('Curva de amenaza sísimca')

% Derivada
parte1 = (P(3) + 2*P(2)*log(IM_range) + 3*P(1)*log(IM_range).^2)./IM_range;
parte2 = exp(P(4)*ones(length(IM_range),1) + P(3)*log(IM_range) + P(2)*log(IM_range).^2 + P(1)*log(IM_range).^3); % lambda_
dlim_poly = abs(parte1.*parte2);

figure
loglog(IM_range,dlim_poly,'.-','color','b')
xlabel('IM: Sa(T1) [g]')
ylabel('|d\lambda_{IM}(IM=im)/d(IM)|')
grid on
legend('Derivada ajuste polinomial tercer orden')
title('Derivada de curva de amenaza sísimca')

% Derivada en puntos de interés (IMs)

dlim_IMs = zeros(length(IMs),1);
for i = 1:length(IMs)
    [im_row, ~] = find(IM_range == IMs(i));
    dlim_IMs(i,1) = abs(dlim_poly(im_row,1));
end


%% P6 a
fprintf('Parte 6 a) \n \n \n')
% Calcular valores

% Obtener EDP para todos los pisos

% a - calcular frecuencia anual media (lambda_EDP) EDP: delta de todos los
% pisos

% Teniendo la media y la desviación de delta_max para no colapso se puede
% calcular la probabilidad de obtener EDP > edp como:
% P(EDP > edp | IM = im, nonCollapse) = 1 - normcdf((log(edp)-mu)/sigma)

EDP_median = EDP_muln_maxFloorIDR;
EDP_stdln = EDP_stdln_maxFloorIDR;

% Interpolamos estos valores para un rango de IM (IM_range, aprovechando)
EDP_median_interp = interp1(IMs,EDP_median,IM_range,'linear','extrap');
EDP_stdln_interp = interp1(IMs,EDP_stdln,IM_range,'linear','extrap');

% Graficamos
figure
plot(IMs,EDP_median,'o','color','r')
hold on
plot(IM_range,EDP_median_interp,'color','b')
hold off
xlabel('IM: Sa(T_1) [g]')
ylabel('\mu_{ln} (%)')
legend('Media para franjas','Media interpolada')
grid on

figure
plot(IMs,EDP_stdln,'o','color','r')
hold on
plot(IM_range,EDP_stdln_interp,'color','b')
hold off
xlabel('IM: Sa(T_1) [g]')
ylabel('\sigma_{ln} (%)')
legend('Desv. Estándar para franjas','Desv. Estandar interpolada')
grid on

% Ahora calculamos la probabilidad
edp1 = 0.01;                                                                % Drift máximo para los pisos
edp2 = 0.04;                                                                % Drift máximo para los pisos
Pexc_EDP_nc_1 = ones(length(IM_range),1) - normcdf(ones(length(IM_range),1)*log(edp1), log(EDP_median_interp), EDP_stdln_interp);
Pexc_EDP_nc_2 = ones(length(IM_range),1) - normcdf(ones(length(IM_range),1)*log(edp2), log(EDP_median_interp), EDP_stdln_interp);

% Graficamos
figure
plot(IM_range,Pexc_EDP_nc_1,'color','b')
xlabel('IM: Sa(T_1) [g]')
ylabel('P(EDP > 0.01 | IM = im, NC')
title('Probabilidad de Excedencia')
legend('Prob. Excedencia \delta = 0.01')
grid on

figure
plot(IM_range,Pexc_EDP_nc_2,'color','b')
xlabel('IM: Sa(T_1) [g]')
ylabel('P(EDP > 0.04 | IM = im, NC')
title('Probabilidad de Excedencia')
legend('Prob. Excedencia \delta = 0.04')
grid on

% Ahora multiplicamos Probabilida*|dlambda(IM=im)/delta(im))|
des_lambda_1 = Pexc_EDP_nc_1.*dlim_poly;                                    % des_lambda de desagregación de lambda
des_lambda_2 = Pexc_EDP_nc_2.*dlim_poly;

figure
plot(IM_range,des_lambda_1,'color','b')
xlabel('IM: Sa(T_1) [g]')
ylabel('P(EDP > 0.01 | IM = im, NC)*|\lambda_{IM}/dIM|')
title('Desagregación \lambda_{EDP > 0.01}')
grid on

figure
plot(IM_range,des_lambda_2,'color','b')
xlabel('IM: Sa(T_1) [g]')
ylabel('P(EDP > 0.04 | IM = im, NC)*|\lambda_{IM}/dIM|')
title('Desagregación \lambda_{EDP > 0.04}')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SI CALCULAMOS LAS "FRACCIONES"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ESAS SERÍAN NUESTRAS PROBABILIDADES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LUEGO LAS MULTIPLICAMOS POR DLAMBDA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RESPECTIVO en su IMs
d1 = 0.01; % interstory drift ratio (delta)
d2 = 0.02;
IDR_data = importdata(convertStringsToChars(string(ResultsDir) + "\" + string(files(IDRid).name)));
EDP_pos = IDR_data.data.IDR(1:cant_pisos,:);                                % Floor EDP_positivos
EDP_neg = IDR_data.data.IDR((cant_pisos+4):(2*cant_pisos+3),:);             % Floor EDP negativos
EDP = max(abs(EDP_pos),abs(EDP_neg));                                       % Floor EDP (el máximo entre los positivos y los negativos)

frac_d1 = zeros(length(IMs),1);
frac_d2 = zeros(length(IMs),1);
for i = 1:length(IMs)
    % Fijar en una franja
    vect = i+1 + 0:cant_franjas:cant_registros*cant_franjas;                % Vector que
    EDP_franja = EDP(:,vect);
    
    % Quitar franjas que den colapsos (en cualquier piso)
    [rows,cols] = find(EDP_franja > 100);
    EDP_franja(:,cols) = [];

    % Largo real del im (cant registros nueva)
    EDP_im_length_new = length(EDP_franja);                                     % Cantidad de registros que no dan colapso para IM = im_i

    % Inicializar Contadores de "cuantas veces da mayor a d1(o d2,d3,d4...)"
    count_d1 = 0;
    count_d2 = 0;

    % Fijar en un registro
    for reg = 1:EDP_im_length_new
        EDP_franja_reg = EDP_franja(reg);                                   % Vector que fija un registro de una franja

        % EDP en verdad es el máximo delta entre todos los pisos
        EDP_maxIDR = max(EDP_franja_reg);                                   % EDP de ese registro (EDP = max(IDR_todosLosPisos))
        
        % ¿Es edp mayor a d1?
        if EDP_maxIDR >= d1
            count_d1 = count_d1 + 1;
        end

        % ¿Es edp mayor a d2?
        if EDP_maxIDR >= d2
            count_d2 = count_d2 + 1;
        end

        % Se pueden realizar preguntas cuantas veces se quiera, o cuantos
        % lambda_EDPs se quieran realizar (en este caso son solo dos)
    end

    % fracción (Probabilidad EDP > edp (edp = d1 o d2 o d3...) | IM, NC)
    frac_d1(i,1) = count_d1/EDP_im_length_new;
    frac_d2(i,1) = count_d2/EDP_im_length_new;
end
des_lambda_frac_1 = frac_d1.*dlim_IMs;
des_lambda_frac_2 = frac_d2.*dlim_IMs;

figure
plot(IM_range,des_lambda_1,'color','b')
hold on
plot(IMs, des_lambda_frac_1,'-o','color','r')
hold off
xlabel('IM: Sa(T_1) [g]')
ylabel('P(EDP > 0.01 | IM = im, NC)*|\lambda_{IM}/dIM|')
title('Desagregación \lambda_{EDP > 0.01}')
legend('Interpolado','Franjas')
grid on

figure
plot(IM_range,des_lambda_2,'color','b')
hold on
plot(IMs, des_lambda_frac_2,'-o','color','r')
hold off
xlabel('IM: Sa(T_1) [g]')
ylabel('P(EDP > 0.04 | IM = im, NC)*|\lambda_{IM}/dIM|')
title('Desagregación \lambda_{EDP > 0.04}')
legend('Interpolado','Franjas')
grid on

lambda_EDP_1 = trapz(IM_range,des_lambda_1);
lambda_EDP_2 = trapz(IM_range,des_lambda_2);

%% P6 b)
fprintf('Parte 6 b) \n \n \n')
% Calcular lambda_EDP para EDP = maxFloor_delta = 0.01 y 0.04 pero ahora
% considerando la posibilidad de colapso
% Ver reporte para saber porqué se hace esto.

% P(C|IM=im) ya se tiene desde máxima verosimilitud
PColIM = normcdf((log(IM_range)-mu_mv)/sigma_mv);

% Cambio Pexc
Pexc_EDP_nc_1 = Pexc_EDP_nc_1.*(ones(length(IM_range),1)-PColIM) + PColIM;
Pexc_EDP_nc_2 = Pexc_EDP_nc_2.*(ones(length(IM_range),1)-PColIM) + PColIM;

% des_lambda
des_lambda_1_C = Pexc_EDP_nc_1.*dlim_poly;                                    % des_lambda de desagregación de lambda
des_lambda_2_C = Pexc_EDP_nc_2.*dlim_poly;

% Figuras
figure
plot(IM_range,des_lambda_1,'color','r')
hold on
plot(IM_range,des_lambda_1_C,'color','b')
hold off
xlabel('IM: Sa(T_1) [g]')
ylabel('P(EDP > 0.01 | IM = im, NC)*|\lambda_{IM}/dIM|')
title('Desagregación \lambda_{EDP}(EDP>0.01)')
legend('\lambda_{EDP}(EDP>0.01|NC)','\lambda_{EDP}(EDP>0.01)')
grid on

figure
plot(IM_range,des_lambda_2,'color','r')
hold on
plot(IM_range,des_lambda_2_C,'color','b')
hold off
xlabel('IM: Sa(T_1) [g]')
ylabel('P(EDP > 0.04 | IM = im, NC)*|\lambda_{IM}/dIM|')
title('Desagregación \lambda_{EDP}(EDP>0.04)')
legend('\lambda_{EDP}(EDP>0.04|NC)','\lambda_{EDP}(EDP>0.04)')
grid on

lambda_EDP_1_C = trapz(IM_range,des_lambda_1_C);
lambda_EDP_2_C = trapz(IM_range,des_lambda_2_C);

% Comparación dado NC y dado NC y C
tabla = table();
tabla.delta_max = [0.01; 0.04];
tabla.lambda_NC = [lambda_EDP_1; lambda_EDP_2];
tabla.lambda = [lambda_EDP_1_C; lambda_EDP_2_C];
tabla.diferencia = [lambda_EDP_1-lambda_EDP_1_C; lambda_EDP_2-lambda_EDP_2_C];
tabla.diferencia_porcentual = [(lambda_EDP_1- lambda_EDP_1_C)/lambda_EDP_1*100 ; (lambda_EDP_2- lambda_EDP_2_C)/lambda_EDP_2*100];
disp(tabla)
clear tabla

%% P6 c}
fprintf('Parte 6 d) \n \n \n')
% Calcular lambda_EDP, EDP PFA_roof para 0.5g y 1g, i.e.:
% lambda_EDP(PFA_roof> 0.5g) y lambda_EDP(PFA_roof>1g), suponer que
% PFA_roof> 1 g si la estructura colapsa

% Tenemos que hacer lo mismo que la a y b al mismo tiempo.

% Teniendo la media y la desviación de PFA_roof para no colapso se puede
% calcular la probabilidad de obtener EDP > edp | NC como:
% P(EDP > edp | IM = im, nonCollapse) = 1 - normcdf((log(edp)-mu)/sigma)

EDP_median = EDP_muln_PFA(:,cant_pisos);                                    % Obtenemos solo la del techo, el resto no nos importa
EDP_stdln = EDP_stdln_PFA(:,cant_pisos);

% Interpolamos estos valores para un rango de IM (IM_range, aprovechando)
EDP_median_interp = interp1(IMs,EDP_median,IM_range,'linear','extrap');
EDP_stdln_interp = interp1(IMs,EDP_stdln,IM_range,'linear','extrap');

% Graficamos
figure
plot(IMs,EDP_median,'o','color','r')
hold on
plot(IM_range,EDP_median_interp,'color','b')
hold off
xlabel('IM: Sa(T_1) [g]')
ylabel('\mu_{ln} (%)')
legend('Media para franjas','Media interpolada')
grid on

figure
plot(IMs,EDP_stdln,'o','color','r')
hold on
plot(IM_range,EDP_stdln_interp,'color','b')
hold off
xlabel('IM: Sa(T_1) [g]')
ylabel('\sigma_{ln} (%)')
legend('Desv. Estándar para franjas','Desv. Estandar interpolada')
grid on

% Ahora calculamos la probabilidad
edp1 = 0.5;                                                                % Drift máximo para los pisos
edp2 = 1;                                                                % Drift máximo para los pisos
Pexc_EDP_nc_1 = ones(length(IM_range),1) - normcdf(ones(length(IM_range),1)*log(edp1), log(EDP_median_interp), EDP_stdln_interp);
Pexc_EDP_nc_2 = ones(length(IM_range),1) - normcdf(ones(length(IM_range),1)*log(edp2), log(EDP_median_interp), EDP_stdln_interp);

% Graficamos
figure
plot(IM_range,Pexc_EDP_nc_1,'color','b')
xlabel('IM: Sa(T_1) [g]')
ylabel('P(PFA_r > 0.5 [g] | IM = im, NC')
title('Probabilidad de Excedencia')
legend('Prob. Excedencia PFA_{roof} = 0.5[g]')
grid on

% Graficamos
figure
plot(IM_range,Pexc_EDP_nc_2,'color','b')
xlabel('IM: Sa(T_1) [g]')
ylabel('P(PFA_r > 1 [g] | IM = im, NC')
title('Probabilidad de Excedencia')
legend('Prob. Excedencia PFA_{roof} = 1 [g]')
grid on

% Ahora multiplicamos Probabilida*|dlambda(IM=im)/delta(im))|
des_lambda_1 = Pexc_EDP_nc_1.*dlim_poly;                                    % des_lambda de desagregación de lambda
des_lambda_2 = Pexc_EDP_nc_2.*dlim_poly;

figure
plot(IM_range,des_lambda_1,'color','b')
xlabel('IM: Sa(T_1) [g]')
ylabel('P(PFA_{roof} > 0.5 [g] | IM = im, NC)*|\lambda_{IM}/dIM|')
title('Desagregación \lambda_{EDP}(PFA_{roof} > 0.5 [g])')
grid on

figure
plot(IM_range,des_lambda_2,'color','b')
xlabel('IM: Sa(T_1) [g]')
ylabel('P(PFA_{roof} > 1 [g] | IM = im, NC)*|\lambda_{IM}/dIM|')
title('Desagregación \lambda_{EDP}(PFA_{roof} > 1 [g])')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SI CALCULAMOS LAS "FRACCIONES"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ESAS SERÍAN NUESTRAS PROBABILIDADES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LUEGO LAS MULTIPLICAMOS POR DLAMBDA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RESPECTIVO EN SU IMs

d1 = 0.5; % Peak floor acceleration in roof (PFA_roof)
d2 = 1;

PFA_data = importdata(convertStringsToChars(string(ResultsDir) + "\" + string(files(PFAid).name)));

% Obtener EDP para todos los pisos
EDP_pos = PFA_data.data.a_t(1:cant_pisos,:);                                % Floor EDP_positivos
EDP_neg = PFA_data.data.a_t((cant_pisos+4):(2*cant_pisos+3),:);             % Floor EDP negativos
EDP = max(abs(EDP_pos),abs(EDP_neg));                                       % Floor EDP (el máximo entre los positivos y los negativos), matriz de EDP(piso,RegistroEscalado), piso = 9, RegistroEscalado = 80

% Seleccionar el del último piso
EDP = EDP(cant_pisos,:);

frac_d1 = zeros(length(IMs),1);
frac_d2 = zeros(length(IMs),1);
for i = 1:length(IMs)
 
    % Fijar en una franja
    vect = i+1 + 0:cant_franjas:cant_registros*cant_franjas;                % Vector que
    EDP_franja = EDP(:,vect);
    
    % Quitar franjas que den colapsos (en cualquier piso)
    [rows,cols] = find(EDP_franja > 100);
    EDP_franja(:,cols) = [];

    % Largo real del im (cant registros nueva)
    EDP_im_length_new = length(EDP_franja);                                     % Cantidad de registros que no dan colapso para IM = im_i

    % Inicializar Contadores de "cuantas veces da mayor a d1(o d2,d3,d4...)"
    count_d1 = 0;
    count_d2 = 0;

    % Fijar en un registro
    for reg = 1:EDP_im_length_new
        EDP_franja_reg = EDP_franja(reg);                                   % Vector que fija un registro de una franja

        % EDP en verdad es el máximo delta entre todos los pisos
        EDP_maxIDR = max(EDP_franja_reg);                                   % EDP de ese registro (EDP = max(IDR_todosLosPisos))
        
        % ¿Es edp mayor a d1?
        if EDP_maxIDR >= d1
            count_d1 = count_d1 + 1;
        end

        % ¿Es edp mayor a d2?
        if EDP_maxIDR >= d2
            count_d2 = count_d2 + 1;
        end

        % Se pueden realizar preguntas cuantas veces se quiera, o cuantos
        % lambda_EDPs se quieran realizar (en este caso son solo dos)
    end

    % fracción (Probabilidad EDP > edp (edp = d1 o d2 o d3...) | IM, NC)
    frac_d1(i,1) = count_d1/EDP_im_length_new;
    frac_d2(i,1) = count_d2/EDP_im_length_new;
end
des_lambda_frac_1 = frac_d1.*dlim_IMs;
des_lambda_frac_2 = frac_d2.*dlim_IMs;

figure
plot(IM_range,des_lambda_1,'color','b')
hold on
plot(IMs, des_lambda_frac_1,'-o','color','r')
hold off
xlabel('IM: Sa(T_1) [g]')
ylabel('P(EDP > 0.5 [g] | IM = im, NC)*|\lambda_{IM}/dIM|')
title('Desagregación \lambda_{EDP}(PFA_r > 0.5 [g])')
legend('Interpolado','Franjas')
grid on

figure
plot(IM_range,des_lambda_2,'color','b')
hold on
plot(IMs, des_lambda_frac_2,'-o','color','r')
hold off
xlabel('IM: Sa(T_1) [g]')
ylabel('P(EDP > 1 [g] | IM = im, NC)*|\lambda_{IM}/dIM|')
title('Desagregación \lambda_{EDP}(PFA_r > 1 [g])')
legend('Interpolado','Franjas')
grid on

lambda_EDP_1 = trapz(IM_range,des_lambda_1);
lambda_EDP_2 = trapz(IM_range,des_lambda_2);

% P(C|IM=im) ya se tiene desde máxima verosimilitud
PColIM = normcdf((log(IM_range)-mu_mv)/sigma_mv);                           % Desde curva de fragilidad

% Cambio Pexc
Pexc_EDP_nc_1 = Pexc_EDP_nc_1.*(ones(length(IM_range),1)-PColIM) + PColIM;
Pexc_EDP_nc_2 = Pexc_EDP_nc_2.*(ones(length(IM_range),1)-PColIM) + PColIM;

% des_lambda
des_lambda_1_C = Pexc_EDP_nc_1.*dlim_poly;                                    % des_lambda de desagregación de lambda
des_lambda_2_C = Pexc_EDP_nc_2.*dlim_poly;

% Figuras
figure
plot(IM_range,des_lambda_1,'color','r')
hold on
plot(IM_range,des_lambda_1_C,'color','b')
hold off
xlabel('IM: Sa(T_1) [g]')
ylabel('P(PFA_{roof} > 0.5[g] | IM = im, NC)*|\lambda_{IM}/dIM|')
title('Desagregación \lambda_{EDP}(PFA_r > 0.5 [g])')
legend('\lambda_{EDP}(PFA_r > 0.5|NC)','\lambda_{EDP}(PFA_r > 0.5)')
grid on

% Figuras
figure
plot(IM_range,des_lambda_2,'color','r')
hold on
plot(IM_range,des_lambda_2_C,'color','b')
hold off
xlabel('IM: Sa(T_1) [g]')
ylabel('P(PFA_{roof} > 1 [g] | IM = im, NC)*|\lambda_{IM}/dIM|')
title('Desagregación \lambda_{EDP}(PFA_r > 1 [g])')
legend('\lambda_{EDP}(PFA_r > 1|NC)','\lambda_{EDP}(PFA_r > 1)')
grid on

lambda_EDP_1_C = trapz(IM_range,des_lambda_1_C);
lambda_EDP_2_C = trapz(IM_range,des_lambda_2_C);

% Comparación dado NC y dado NC y C
tabla = table();
tabla.delta_max = [0.01; 0.04];
tabla.lambda_NC = [lambda_EDP_1; lambda_EDP_2];
tabla.lambda = [lambda_EDP_1_C; lambda_EDP_2_C];
tabla.diferencia = [lambda_EDP_1-lambda_EDP_1_C; lambda_EDP_2-lambda_EDP_2_C];
tabla.diferencia_porcentual = [(lambda_EDP_1- lambda_EDP_1_C)/lambda_EDP_1*100 ; (lambda_EDP_2- lambda_EDP_2_C)/lambda_EDP_2*100];
disp(tabla)
clear tabla

%% P6 d
% Encontrar P(C|IM_10%50años) y P(C|IM_2%50Años)
fprintf('Parte 6 d) \n \n \n')
% Para ello, primer se requiere encontrar la tasa anual media de IM_%
% sabiendo que la probabilidad distribuye como poisson y posteriormente
% evaluarlo para la curva de fragilidad

% lambda_IM = -log(1-prob)/#años
lambda_im_10_50 = -log(1-0.1)/50;
lambda_im_2_50 = -log(1-0.02)/50;

% Con ello encontramos qué valor de lambda cumple con el IM según la curva
% de amenaza utilizando el polinomio de Miranda
% log(lambda) = a0 + a1*log(IM) + a2*log(IM)^2 + a3*log(IM)^3
syms im
IM_1 = solve(   log(lambda_im_10_50) == P(4) + P(3)*log(im) + P(2)*log(im)^2 + P(1)*log(im)^3,im   );
IM_2 = solve(   log(lambda_im_2_50) == P(4) + P(3)*log(im) + P(2)*log(im)^2 + P(1)*log(im)^3 ,im   );

% Evaluar los IMs en la curva de fragilidad
% con Máxima verosimilitud
PC_IM_1_mv = normcdf((log(IM_1)-mu_mv)/sigma_mv);
PC_IM_2_mv = normcdf((log(IM_2)-mu_mv)/sigma_mv);
% con Mínimos cuadrados
PC_IM_1_mc = normcdf((log(IM_1)-mu_mc)/sigma_mc);
PC_IM_2_mc = normcdf((log(IM_2)-mu_mc)/sigma_mc);

% figuras
figure
plot(IM_range,normcdf((log(IM_range)-mu_mv)/sigma_mv),'color','g')
hold on
plot(IM_1,PC_IM_1_mv,'O','color','r')
plot(IM_2,PC_IM_2_mv,'^','color','r')
plot(IM_range,normcdf((log(IM_range)-mu_mc)/sigma_mc),'color','b')
plot(IM_1,PC_IM_1_mc,'O','color','#05295c')
plot(IM_2,PC_IM_2_mc,'^','color','#05295c')
xlabel('IM: Sa(T_1) [g]')
ylabel('P(C|IM)')
legend('CFC Máx.Verosimulitud','IM 10% 50 años MV','IM 2% 50 años MV','CFC Mínimos Cuadrados','IM 10% 50 años MC','IM 2% 50 años MC')
title('Curva de Fragilidad')
grid on

tabla = table();
tabla.Parametro = ["lambda_IM";"IM=im";"P(C|IM=im)_maxVer";"P(C|IM=im)_minCuad"];
tabla.pPR_10porc_50anios = [double(lambda_im_10_50);double(IM_1);double(PC_IM_1_mv);double(PC_IM_1_mc)];
tabla.pPR_2porc_50anios = [double(lambda_im_2_50);double(IM_2);double(PC_IM_2_mv);double(PC_IM_2_mc)];
disp(tabla)

%% P6 e)
% Calcular desagregación lambda_c y lambda(c|t=50)
fprintf('Parte 6 e) \n \n \n')

% lambda_c(im) = P(C|IM=im)|dlambda/dIM|dim

% PCollapse|IM
PCim = normcdf((log(IM_range)-mu_mv)/sigma_mv);

% Derivada
% desde parte 5: dlim_poly

% lambda_c(im)
lambda_c_im = PCim.*dlim_poly;

figure
plot(IM_range,lambda_c_im)
colororder('r')
grid on
xlabel('IM: Sa(T_1) [g]')
ylabel('P(C|IM=im)|d\lambda_{IM}/dIM|')
title('Desagregación de MAF colapso')

% lambda_C
lambda_c = trapz(IM_range,lambda_c_im);

fprintf('lambda_c = %.10f\n',lambda_c)

% probabilidad de que estructura colapse de acá en 50 años
t = 50; % años
Pcol_t = 1 - exp(-lambda_c*50);

fprintf('Probabilidad de que la estructura colapse de acá en %.0f años es %.5f',t,Pcol_t)











