%% Evaluación del Riesgo de Colapso
% Alexis Contreras - Martina Sanguinetti
% Ingeniería Sísmica Avanzada - Proyecto de Investigación - USM 2022

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
ResultsDir = 'Results/Results_TMD18';                                             % Carpeta donde están los resultados
cant_pisos = 10;                                                            % Cantidad de pisos del edificio
cant_registros = 31;                                                        % Cantidad de registros por cada IM (se considera que todas las IMs tienen la misma cantidad de registros)
% IMs = [0.2; 0.5; 1; 1.2; 2];
IMs = [2.5; 2.9]; % g                                                       % Valor IMs de las franjas del "IDA"
cant_franjas = length(IMs);                                                 % Cantidad de franjas
% g = 9.81;
IM_limit = 7;
% Nombre de archivos a utilizar para determinar colapso
nA_Disp = 'Displacement.xlsx';                                              % Nombre archivo obtenido desde THAMDOF de desplazamiento
% Inputs en Curva de Fragilidad
% Inputs en Estimación lambda_c

%% Load Data
files = dir(ResultsDir);                                                    % Ver carpeta
files(1:1:2) = [];

for i = 1:length(files)
    if isequal(files(i).name,nA_Disp)
        Dispid = i;
        continue
    end
end

%% Collapse Fragiility Curve
% Determinar Curva de Fragilidad de Colapso (CFC) con método de mínimos
% cuadrados y método de máxima verosimilitud y comparar

% El colapso quedará definido cuando el desplazamiento, de cualquier piso,
% exceda el desplazamiento de resistencia 0 (du - desplazamiento último),
% hasta el momento, es solo hasta 10^5 cualquier piso, pero debería ser el
% desplazamiento de resistencia 0

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
            if F_EDP_vect(j,r) > 10^5
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
sigmas = 0.01:0.01:1;

% % Método de los mínimos cuadrados
% E = zeros(cant_franjas,1);
% mc = +inf;
% for m = 1:length(mus)
%     for s = 1:length(sigmas)
%         for i = 1:cant_franjas
%             E(i) = abs(fraccion(i) - normcdf((log(IMs(i)) - mus(m))./sigmas(s)))^2;
%         end
%         if sum(E) < mc
%             mc = sum(E);
%             sigma_mc = sigmas(s);
%             mu_mc = mus(m);
%         end
%     end
% end

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
% tabla.mc = [mu_mc; sigma_mc];
tabla.mv = [mu_mv; sigma_mv];
disp(tabla)
clear tabla

% Generar figura para corrobarar método

more_IM = (0.1:0.01:IM_limit).';                                               % Mismo rango que los datos para el ajuste polinomial
% IM_range = sort([more_IM; IMs]);                                            % Vector para el cual se va a realizar la interpolación, notar que se agregaron los IMs de Interés
IM_range = more_IM;
% Fracciones
figure
plot(IMs,fraccion,'o','Color','r','linewidth',1.5)
xlabel('IM: Sa(T1) [g]')
ylabel('P(C|IM=im)')
xlim([0 IM_range(end)])
ylim([0 1])
title('Fracciones zj/nj')
legend('Data')
grid on

% Curva de Fragilidad
Pc_im = normcdf((log(IM_range)-mu_mv)/sigma_mv);

figure
plot(IMs,fraccion,'o','Color','r','linewidth',1.5)
hold on
% plot(IM_range,normcdf((log(IM_range)-mu_mc)/sigma_mc),'color','b','Linewidth',1.5)
plot(IM_range,normcdf((log(IM_range)-mu_mv)/sigma_mv),'color','b','Linewidth',1.5)
hold off
xlabel('IM: Sa(T1) [g]')
ylabel('P(C|IM=im)')
title('Curva de Fragilidad de Colapso')
legend('Data','Máxima Verosimilitud')
grid on

% original_more_IM = more_IM;
original_IM_range = IM_range;
%% Seismic Hazard Curve
% Obtener la curva de amenaza sísmica desde el sitio web SHC (nshcweb2)
% suelo % clase D, coord(34.05,-118.25), realizar el ajuste polinomial de
% cuarto orden

% Data USGS
% Datos que otorga USGS para la amenaza sísmica

original_IM_Saavg = [0.01;0.013;0.015;0.019;0.023;0.028;0.033;0.041;0.049;0.06;0.073;0.089;0.108;0.131;0.159;0.193;0.234;0.285;0.346;0.421;0.511;0.621;0.755;0.917;1.115;1.354;1.646;2];
original_lambda_Saavg = [0.274529554;0.223874996;0.179741466;0.142766053;0.112165203;0.086762584;0.066803522;0.050979124;0.038715512;0.029373076;0.022110009;0.016878141;0.012882974;0.009983871;0.007564458;0.005724685;0.004229869;0.003054078;0.002095112;0.001332039;0.000780955;0.000424777;0.000215999;0.000114034;5.8321E-05;1.91771E-05;6.21444E-06;2.10781E-06];

% Modificando valores originales, para obtener un mejor ajuste
IM_Saavg = [0.028;0.033;0.041;0.049;0.06;0.073;0.089;0.108;0.131;0.159;0.193;0.234;0.285;0.346;0.421;0.511;0.621;0.755;0.917;1.115;1.354;1.646;2];
lambda_Saavg = [0.086762584;0.066803522;0.050979124;0.038715512;0.029373076;0.022110009;0.016878141;0.012882974;0.009983871;0.007564458;0.005724685;0.004229869;0.003054078;0.002095112;0.001332039;0.000780955;0.000424777;0.000215999;0.000114034;5.8321E-05;1.91771E-05;6.21444E-06;2.10781E-06];

% Ajuste polinomial ------
% Determinar coeficientes
[P,S] = polyfit(log(IM_Saavg),log(lambda_Saavg),4);

% Ajustar curva de polinomio cuarto orden
% El IM_range a utilizar para calcular los puntos serán los mismos que en
% la parte anterior
% more_IM = (0.1:0.01:3).';                                                   % Mismo rango que los datos para el ajuste polinomial
% IM_range = sort([more_IM; IMs]);                                            % Vector para el cual se va a realizar la interpolación, notar que se agregaron los IMs de Interés

% lambda_poly = exp(P(5)*ones(length(IM_range),1) + P(4)*log(IM_range) + P(3)*log(IM_range).^2 + P(2)*log(IM_range).^3 + P(1)*log(IM_range).^4); 
% R_square = 1 - S.normr^2 / norm(log(lambda_Saavg)-mean(log(lambda_Saavg)))^2; % Valor de R^2
% 
% figure
% loglog(original_IM_Saavg,original_lambda_Saavg,'o','color','#076F51','LineWidth',1.5)
% hold on
% loglog(IM_Saavg,lambda_Saavg,'o','color','r','LineWidth',1.5)
% loglog(IM_range,lambda_poly,'.-','color','b','LineWidth',1.5)
% hold off
% xlabel('IM: Sa_{avg} [g]')
% ylabel('\lambda_{IM} [1/yr]')
% grid on
% legend('Datos USGS','Datos USGS Utilizados para el ajuste','Ajuste polinomial tercer orden')
% text(10^0,10^-2,"R2 = " + string(R_square))
% title('Curva de amenaza sísimca')
% 
% % Derivada de la amenaza sísmica
% parte1 = (P(4) + 2*P(3)*log(IM_range) + 3*P(2)*log(IM_range).^2 + 4*P(1)*log(IM_range).^3)./IM_range;
% parte2 = exp(P(5)*ones(length(IM_range),1) + P(4)*log(IM_range) + P(3)*log(IM_range).^2 + P(2)*log(IM_range).^3 + P(1)*log(IM_range).^4); % lambda_
% dlim_poly = abs(parte1.*parte2);
% 
% figure
% loglog(IM_range,dlim_poly,'.-','color','b','LineWidth',1.5)
% xlabel('IM: Sa_{avg} [g]')
% ylabel('|d\lambda_{IM}(IM=im)/d(IM)|')
% grid on
% legend('Derivada ajuste polinomial tercer orden')
% title('Derivada de curva de amenaza sísimca')


% Calculando para IM range original (el de toda la curva de fragilidad,
% para que tengan el mismo largo de elementos)
IM_range = original_IM_range;
lambda_poly = exp(P(5)*ones(length(IM_range),1) + P(4)*log(IM_range) + P(3)*log(IM_range).^2 + P(2)*log(IM_range).^3 + P(1)*log(IM_range).^4); 
parte1 = (P(4) + 2*P(3)*log(IM_range) + 3*P(2)*log(IM_range).^2 + 4*P(1)*log(IM_range).^3)./IM_range;
parte2 = exp(P(5)*ones(length(IM_range),1) + P(4)*log(IM_range) + P(3)*log(IM_range).^2 + P(2)*log(IM_range).^3 + P(1)*log(IM_range).^4); % lambda_
dlim_poly = abs(parte1.*parte2);


%% Desagregación lambda_c
des_lambda_c = Pc_im.*dlim_poly;

% figure
% plot(IM_range,des_lambda_c)
% xlabel('IM: Sa_{avg}')
% ylabel('P(C|IM=im)|d\lambda_{IM}(IM=im)/d(IM)|')
% title('Desagregación de \lambda_c')
% grid on

%% lambda_c
lambda_c = trapz(IM_range,des_lambda_c);
disp(lambda_c)

%% Figuras para mostrar en reporte

figure
subplot(3,1,1) 
plot(IM_range,Pc_im,'color','b','LineWidth',1.5)
hold on
plot(IMs,fraccion,'o','color','r','LineWidth',1.5)
legend('Curva de Fragilidad de colapso')
ylabel('P(C|IM=im)')
xlim([0 IM_limit])
grid on
subplot(3,1,2)
plot(IM_range,dlim_poly,'color','r','LineWidth',1.5)
legend('Derivada de amenaza sísmica')
ylabel('|d\lambda_{IM}/d(IM)|')
xlim([0 IM_limit])
grid on
subplot(3,1,3)
plot(IM_range,des_lambda_c,'LineWidth',1.5)
xlabel('IM: Sa_{avg}')
ylabel('P(C|IM)|d\lambda_{IM}/d(IM)|')
grid on
xlim([0 IM_limit])
legend(convertStringsToChars("lambda_c = " + string(lambda_c)))

% PC_old = Pc_im;