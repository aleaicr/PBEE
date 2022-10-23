%% Tarea 2 - Ingeniería Sísmica Avanzada
% Alexis Contreras - Martina Sanguinetti
% Ingeniería Sísmica Avanzada - USM

% Parte 2 Estructura 1

%% Inicializar
clear variables
close all
clc

%% Estructura 1
W = 550; % tonf                                                             % Peso de la estructura
T = 1.5; % sec                                                              % Periodo fundamental
zeta = 0.05; % %                                                            % Fracción de amortiguamiento
Cy = 0.1; %                                                                 % Coeficiente sísmico inelástico
% alfa = 0.02; %                                                              % Coeficiente de endurecimiento post-fluencia
g = 9.81;   % m/s2

%% Parte 2 - Estructura 1
%% 1.
% Currvas IDA utilizando el desplazamiento máximo como EDP y Sa(T1,xi=5%)
% como IM.
% No debería haber colapso ya que no estamos considerando efectos P-Delta y
% modelo no considera degradación.
% Indicar punto de fluencia en las curvas IDA
dIM = 0.1; % g                                                             % Incremento de IM
IMmax = 3; % g                                                              % IM máxima

% Nombre de carpeta donde están los archivos
ResultsDir = 'Resultados1';

% Escribir el nombre de todos los archivos
ResultsName = 'estructura1';                                                 % Trucazo para leer nombres de una carpeta

% Generar curvas IDA
[EDP,IM] = getIdaCurves(ResultsDir, ResultsName);
% EDP(nFranjas,nRegistros)
% IM(nFranjas,nRegistros)

ay = Cy;            % ay = Cy*g, pero en unidad g hay que volver a dividirlo
dy = Cy*W/(4*pi^2*W/g/T^2);

% Figuras 
figure
plot(EDP,IM/g,'.-','color','#606060')
hold on
plot(dy,ay,'^','Color','r')
xline(dy,'Color','r')
hold off
text(2,0.5,['Fluencia (dy,ay)' string([dy ay])])
ylim([0 IMmax])
xlabel('EDP: Desplazamiento Máximo [m]')
ylabel('IM: Sa(T_1,\xi) [g]')                                               % En IIDAP, se pusieron unidades de metros y segundos => aceleraciones en m/s2
title('Multi-Record IDA Curves', ' Estructura T_1 = 1.5s; \xi = 5%')
grid on


% figure
% plot(IM/g,EDP,'.-','Color','#606060')
% ylabel('EDP: Desplazamiento Máximo [m]')
% xlabel('IM: Sa(T_1,\xi) [g]')                                               % En IIDAP, se pusieron unidades de metros y segundos => aceleraciones en m/s2
% title('Multi-Record IDA Curves', ' Estructura T_1 = 1.5s; \xi = 5%')
% grid on

%% 2.
% Graficar curva IDA mediana y la desviación estándar logarítmica de los
% desplazamientos como función de Sa(T1)

% mu_ln = median(log(x))
% median = geomean(x)
% std_ln = std(log(x))

% EDP_muln = mean(log(EDP')); % me da negativo usando lo del ppt
EDP_median = geomean(EDP');                                                 % Mediana
EDP_stdln = std(log(EDP'));                                                 % Desviación estandar
EDP_muln = exp(log(EDP_median) + 0.5*EDP_stdln.^2);                         % Estimación de la Media

% Gráficos

% Mediana logarítmica + Curvas IDA
figure
plot(EDP_median',IM(:,1)/g,'color','r','LineWidth',2)
hold on
plot(EDP,IM/g,'.-','color','#606060')
xline(dy,'Color','b')
ylim([0 IMmax])
hold off
xlabel('EDP y geomean(EDP)')
ylabel('IM: Sa(T_1,\xi) [g]') 
grid on
title('Mediana logarítmica EDP')


% Mediana logarítmica
figure
plot(IM(:,1)/g,EDP_median','.-','color','k')
ylabel('Mediana logarítmica de EDP, geomean(EDP)')
xlabel('IM: Sa(T_1,\xi) [g]') 
grid on
title('Mediana logarítmica EDP')

% Desviación Estándar Logarítmica
figure
plot(IM(:,1)/g,EDP_stdln','.-','color','k')
ylabel('Desviación estandar logarítmica de EDP, std(log(EDP))')
xlabel('IM: Sa(T_1,\xi) [g]') 
grid on
title('Desviación estándar EDP')

% Media
figure
plot(IM(:,1)/g,EDP_muln','.-','color','k')
ylabel('Media logarítmica de EDP, mean(log(EDP))')
xlabel('IM: Sa(T_1,\xi) [g]') 
grid on
title('Media logarítmica EDP')

% Mediana +- sigma
figure
plot(IM(:,1)/g,EDP_median','.-','color','r')
hold on
plot(IM(:,1)/g,exp(log(EDP_median')+EDP_stdln'),'.-','color','k')
plot(IM(:,1)/g,exp(log(EDP_median')-EDP_stdln'),'.-','color','k')
hold off
ylabel('Mediana logarítmica de EDP, geomean(EDP)')
xlabel('IM: Sa(T_1,\xi) [g]') 
grid on
title('Mediana logarítmica EDP')

figure
plot(EDP_median',IM(:,1)/g,'.-','color','r')
hold on
plot(exp(log(EDP_median')+EDP_stdln'),IM(:,1)/g,'.-','color','k')
plot(exp(log(EDP_median')-EDP_stdln'),IM(:,1)/g,'.-','color','k')
hold off
xlabel('Mediana logarítmica de EDP, geomean(EDP)')
ylabel('IM: Sa(T_1,\xi) [g]') 
grid on
title('Mediana logarítmica EDP')

%% 3
% Copiar gráficos media vs Sa, pero agregar curva de desplazamientos
% máximos que se obtendrían asumiendo principio de igualdad de
% desplazamientos

% Se cumpliría este principio solo en la primera estructura ya tiene periodo
% mayor a 1[s], pero no para la segunda estructura

% Principio de grandes desplazamientos
omega = 2*pi/T;
Sd = IM(:,1)./omega^2;  % IM m/s2 omega en rad2/s2 -> Sd = metros

% Media
figure
plot(IM(:,1)/g,EDP_muln','.-','color','k')
hold on
plot(IM(:,1)/g,Sd)
hold off
ylabel('Media logarítmica de EDP, mean(log(EDP))')
xlabel('IM: Sa(T_1,\xi) [g]') 
grid on
title('Media logarítmica EDP')
legend('EDP_muln','Sd')

% Media
figure
plot(EDP_muln',IM(:,1)/g,'.-','color','k')
hold on
plot(Sd,IM(:,1)/g)                                                        % Dejamos Sd en centimetros
hold off
xlabel('Media logarítmica de EDP, mean(log(EDP))')
ylabel('IM: Sa(T_1,\xi) [g]') 
grid on
title('Media logarítmica EDP')
legend('EDP_muln','Sd')


%% 4
% Suponiendo distribución lognormal para EDP dado IM, calcular probabilidad
% de que los desplazamientos máximos excedan 25cm para una ordenada
% espectral de 1g.

% % Desplazamiento máximo objetivo
% desplMaxObj = 25/100; % 25  cm = 25/100 m
% 
% % Posición IM = 1g
% posIM1g = find(IM(:,1)/g == 1);
% 
% % Columna de EDP que necesito para este análisis
% EDP_col = EDP(:,posIM1g);
% EDP_logncdf = logncdf(EDP_col);
% EDP_valuesInterp = interp1(EDP_col,EDP_logncdf,desplMaxObj);                % Interpolando, el valor de la probabilidad es este
% 
% % Distribución
% figure
% plot(EDP_col,EDP_logncdf)
% hold on
% plot(desplMaxObj,EDP_valuesInterp,'^r')
% hold off
% xlabel('EDP|IM = 1g')
% ylabel('logncdf(EDP)|IM = 1g')
% text(desplMaxObj+0.1,EDP_valuesInterp,['P(EDP > 25cm | IM = 1g) =' string(EDP_valuesInterp)])
% grid on
% 
% fprintf('P4. P(EDP > 25cm | IM = 1g) = %.4f\n',EDP_valuesInterp);

% P(EDP>edp | IM = im) = 1 - cdf((ln(edp)-ln(EDP_median|IM=im))/EDP_sigmaln|IM=im)
% Ya tenemos EDP_median|IM = 1g
% Ya tenemos EDP_stdln|IM = 1g

% Desplazamiento máximo objetivo
desplMaxObj = 25/100; % 25  cm = 25/100 m

% Posición IM = 1g
posIM1g = find(IM(:,1)/g == 1);

% Probabilidad
probabilidad = 1-normcdf(log(desplMaxObj),log(EDP_median(posIM1g)),EDP_stdln(posIM1g));
fprintf('%f\n',probabilidad)


%% 5
% Repetir proceso pero ahora con Sa_avg como IM, para definir Sa_avg
% utilizar 20 periodos entre 0.2T y 3.0T

% Se corre un nuevo análisis en II-DAP pero con IM alternativa Sa_avg en
% 'analysis type'

% Guardar datos anteriores
EDP_Sa = EDP;                                                               % EDP cuando IM es Sa(T1)
IM_Sa = IM;                                                                 % IM  = Sa(T1)

% Nombre para importar datos
ResultsDir_Sa_avg = 'Resultados1';
ResultsName_Sa_avg = 'estructura1_sa_avg';

% getIdaCurves
[EDP,IM] = getIdaCurves(ResultsDir_Sa_avg, ResultsName_Sa_avg);

% Repetir paso 2 para IM = Sa_avg
% Graficar curva IDA mediana y la desviación estándar logarítmica de los
% desplazamientos como función de Sa(T1)

EDP_median = geomean(EDP');                                                 % Mediana
EDP_stdln = std(log(EDP'));                                                 % Desviación estandar
EDP_muln = exp(log(EDP_median) + 0.5*EDP_stdln.^2);                         % Estimación de la Media

% Gráficos
% EDP vs IM

figure
plot(EDP,IM/g,'.-','color','#909090')
ylabel('EDP: Desplazamiento Máximo [m]')
xlabel('IM: Sa(T_1,\xi) [g]')                                               % En IIDAP, se pusieron unidades de metros y segundos => aceleraciones en m/s2
title('Multi-Record IDA Curves', ' Estructura T_1 = 1.5s; \xi = 5%')
grid on

% Gráficos
% Mediana logarítmica
figure
plot(IM(:,1)/g,EDP_median','.-','color','k')
ylabel('Mediana logarítmica de EDP, geomean(EDP)')
xlabel('IM: Sa(T_1,\xi) [g]') 
grid on
title('Mediana logarítmica EDP')

% Desviación Estándar Logarítmica
figure
plot(IM(:,1)/g,EDP_stdln','.-','color','k')
ylabel('Desviación estandar logarítmica de EDP, std(log(EDP))')
xlabel('IM: Sa(T_1,\xi) [g]') 
grid on
title('Desviación estándar EDP')

% Media
figure
plot(IM(:,1)/g,EDP_muln','.-','color','k')
ylabel('Media logarítmica de EDP, mean(log(EDP))')
xlabel('IM: Sa(T_1,\xi) [g]') 
grid on
title('Media logarítmica EDP')

% Mediana +- sigma
figure
plot(IM(:,1)/g,EDP_median','.-','color','r')
hold on
plot(IM(:,1)/g,exp(log(EDP_median')+EDP_stdln'),'.-','color','k')
plot(IM(:,1)/g,exp(log(EDP_median')-EDP_stdln'),'.-','color','k')
hold off
xlabel('Mediana logarítmica de EDP, geomean(EDP)')
ylabel('IM: Sa(T_1,\xi) [g]') 
grid on
title('Mediana logarítmica EDP')

figure
plot(EDP_median',IM(:,1)/g,'.-','color','r')
hold on
plot(exp(log(EDP_median')+EDP_stdln'),IM(:,1)/g,'.-','color','k')
plot(exp(log(EDP_median')-EDP_stdln'),IM(:,1)/g,'.-','color','k')
hold off
ylabel('Mediana logarítmica de EDP, geomean(EDP)')
xlabel('IM: Sa(T_1,\xi) [g]') 
grid on
title('Mediana logarítmica EDP')


% Volvemos a dejar como estaba antes
EDP = EDP_Sa;                                                               % EDP cuando IM es Sa(T1)
IM = IM_Sa;                                                                 % IM  = Sa(T1)

%% 6
% Extraer curvas de amenaza sísmica desde USGS para periodos en estudio
% Realizar interpolación lineal de los datos  de 0.1g
% Notar que no hay para T = 1.5sec -> Interpolar de forma hiperbólica en
% periodo (Sa prop. 1/T)

% USGS Data
USGS_data_Dir = 'Resultados1';
USGS_data_Name = 'USGS_data_est1';

USGS_data = readmatrix([USGS_data_Dir '\' USGS_data_Name]);

% Hazazrd Curve (T = 1sec)
IM_1 = USGS_data(:,1);      % [g]
lambda_1 = USGS_data(:,2);

% Hazard Curve (T = 2sec)
IM_2 = USGS_data(:,3);      % [g]
lambda_2 = USGS_data(:,4);

% Figura
figure
loglog(IM_1,lambda_1)
hold on
loglog(IM_2,lambda_2)
xlabel('Ground Motion [g] IM = Sa(T_1)')
ylabel('Annual Frequency of Excedence \lambda_{IM}')
title('Hazard Curves')
grid on
legend('T=1[sec]','T=2[sec]')

% % Interpolamos ambas para IM
% IM_1_interp1 = 0.1:0.1:IMmax;
% IM_2_interp1 = 0.1:0.1:IMmax;
% 
% lambda_1_interp1 = interp1(IM_1,lambda_1,IM_1_interp1);
% lambda_2_interp1 = interp1(IM_2,lambda_2,IM_2_interp1);
% 
% figure
% loglog(IM_1_interp1,lambda_1_interp1,'.-','Color','r')
% hold on
% loglog(IM_2_interp1,lambda_2_interp1,'.-','Color','b')
% hold off
% xlabel('IM = Sa(T1) [g]')
% ylabel('Frecuencia Anual de excedencia \lambda_{IM}')
% title('Curvas de amenaza interpoladas')
% legend('T=1[sec]','T=2[sec]')
% grid on

% Trucazo, haciendo trampa

% Interpolamos para tener mismas lambda_IM
lambda_IM_interp1 = logspace(-5,-1,20)';    % Vector de 200 elementos entre 0.00001 y 0.1 en escala logaritmica

IM_1sec_interp1 = interp1(lambda_1,IM_1,lambda_IM_interp1);
IM_2sec_interp1 = interp1(lambda_2,IM_2,lambda_IM_interp1);

% figure
% loglog(IM_1sec_interp1,lambda_IM_interp1,'.-','Color','r')
% hold on
% loglog(IM_2sec_interp1,lambda_IM_interp1,'.-','Color','b')
% hold off
% xlabel('IM = Sa(T1) [g]')
% ylabel('Frecuencia Anual de excedencia \lambda_{IM}')
% title('Curvas de amenaza interpoladas')
% legend('T=1[sec]','T=2[sec]')
% grid on

% Interpolación hiperbólica para obtener lambda vs IM para T = 1.5sec
% IM1 = alfa*(1/T1) + beta
% IM2 = alfa*(1/T2) + beta
% Solucionar para alfa y beta

T1 = 1; % sec
T2 = 2; % sec
IM_1_5_sec = zeros(length(IM_1sec_interp1),1);

syms alfa beta
for i = 1:length(lambda_IM_interp1)                                              % IM_1_interp1 y IM_2_interp1 son los mismos
    sol = solve(IM_1sec_interp1(i)-alfa*(1/T1)-beta,IM_2sec_interp1(i)-alfa*(1/T2)-beta);
    alfa_val = double(sol.alfa);
    beta_val = double(sol.beta);
    IM_1_5_sec(i,1) = alfa_val*(1/T)+beta_val;
end

figure
loglog(IM_1sec_interp1,lambda_IM_interp1,'.-','Color','#007d79')
hold on
loglog(IM_2sec_interp1,lambda_IM_interp1,'.-','Color','#6929c4')
loglog(IM_1_5_sec,lambda_IM_interp1,'.-','Color','r')
loglog(IM_1,lambda_1,'Color','#007d79')
loglog(IM_2,lambda_2,'Color','#6929c4')
hold off
xlabel('IM = Sa(T1) [g]')
ylabel('Frecuencia Anual de excedencia')
title('Curvas de amenaza interpoladas')
legend('T = 1[sec]','T = 2[sec]','T = 1.5[sec]','T = 1[sec] USGS', 'T = 2[sec] USGS')
grid on


% Terminando trucazo vuelvo a interpolar en IM para mi curva de amenaza sismica de T = 1.5s
IM_interp1 = 0.1:0.1:IMmax;
lambda_interp1 = interp1(IM_1_5_sec,lambda_IM_interp1,IM_interp1);

figure
loglog(IM_interp1,lambda_interp1,'.-','Color','r')
xlabel('IM = Sa(T1) [g]')
ylabel('Frecuencia Anual de excedencia')
title('Curva de amenaza interpoladas')
legend('T = 1.5s')
grid on


%% 6.2 Análisis probabilístico de la demanda sísmica
% Relizar un análisis probabilístico de la demanda sísmica para caluclar
% lambda_EDP, con IM = Sa(T1) y despl_max como EDP

% Ya tenemos la curva de amenaza sísmica -> podemos obtener d/dIM
% (lambda(IM))

% Hasta la fecha no me dan los vectores de IM_T_15sec y lambdaIM_15sec (1.5 sec)
% De la forma como las pide el profe, pero trabajemos con los que tenemos 
% como si estuvieran correctos

% Supongo que numéricamente, al no contar con dIM, se calcula como una
% sumatoria en vez de una integral

% Renombramos parámetros
% IM = IM_1_5_sec;

% Desde amenaza
IM_amenaza = IM_interp1;                                                    % 0.1:0.1:3   [g]
lambda_amenaza = lambda_interp1;                                            % Interpolación
dIM = 0.1; % g

% Desde EDP
IM_IDA = IM_Sa;                                                             % Es el mismo que la interpolación de IM_amenaza (0.1:0.1:3 [g]) 
EDP_IDA = EDP_Sa;

% Obtención de abs(d/dIM (lambda_IM(x)))
dlambdadIM = abs(diff(lambda_amenaza)/dIM)';

% Rango de edp para lambda_EDP
edps = (0.001:0.001:2)';

lambda_EDP = zeros(length(edps),1);
for e = 1:length(edps)
    edp = edps(e);
    multiplicacion = zeros(length(IM),1);
    for i = 1:length(IM)
        im = IM(i)/g;
         % Posición IM = 1g en curva de 
        posIM = find(IM_IDA(:,1)/g == im);
        
        % EDP_median
        EDP_median_val = EDP_median(posIM); 
        
        % EDP_stdln
        EDP_stdln_val = EDP_stdln(posIM);
        
        % Parte EDP|IM
        probabilidad = 1-normcdf(log(edp),log(EDP_median_val),EDP_stdln_val);
    
        % Parte IM
        dlambdadIM_val = dlambdadIM(im);

        % dIM
        % dIM = 0.1; % g

        % Multiplicacion
        multiplicacion(im) = probabilidad*dlambdadIM_val*dIM;
    end
    lambda_EDP(edp) = sum(multiplicacion);
end











