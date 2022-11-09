%% Tarea 2 - Ingeniería Sísmica Avanzada
% Alexis Contreras - Martina Sanguinetti
% Ingeniería Sísmica Avanzada - USM

% Parte 2 Estructura 2

%% Inicializar
clear variables
close all
clc

%% Estructura 2
W = 125; % tonf                                                             % Peso de la estructura
T = 0.3; % sec                                                              % Periodo fundamental
zeta = 0.05; % %                                                            % Fracción de amortiguamiento
Cy = 0.2; %                                                                 % Coeficiente sísmico inelástico
% alfa = 0.02; %                                                              % Coeficiente de endurecimiento post-fluencia
g = 9.81;   % m/s2

%% Parte 2 - Estructura 1
%% 1.
% Currvas IDA utilizando el desplazamiento máximo como EDP y Sa(T1,xi=5%)
% como IM.
% No debería haber colapso ya que no estamos considerando efectos P-Delta y
% modelo no considera degradación.
% Indicar punto de fluencia en las curvas IDA
dIM = 0.1; % g                                                              % Incremento de IM
IM_max = 3; % g                                                              % IM máxima
IM_step = 0.1;
IM_min = 0.1;

% Nombre de carpeta donde están los archivos
ResultsDir = 'Resultados2';

% Escribir el nombre de todos los archivos
ResultsName = 'estructura2';                                                 % Trucazo para leer nombres de una carpeta

% Generar curvas IDA
[EDP,IM] = getIdaCurves(ResultsDir, ResultsName);
% EDP(nFranjas,nRegistros)
% IM(nFranjas,nRegistros)

ay = Cy;            % ay = Cy*g, pero en unidad g hay que volver a dividirlo
dy = Cy*W/(4*pi^2*W/g/T^2);
EDP_median = geomean(EDP');   
% Figuras 
figure
plot(EDP_median',IM(:,1)/g,'color','r','LineWidth',2)
hold on
plot(EDP,IM/g,'.-','color','#606060')
hold on
plot(dy,ay,'^','Color','r')
xline(dy,'Color','b')
hold off
text(2,0.5,['Fluencia (dy,ay)' string([dy ay])])
ylim([0 IM_max])
xlabel('EDP: Desplazamiento Máximo [m]')
ylabel('IM: Sa(T_1,\xi) [g]')                                               % En IIDAP, se pusieron unidades de metros y segundos => aceleraciones en m/s2
title('Multi-Record IDA Curves', ' Estructura T_1 = 0.3s; \xi = 5%')
grid on

%% 2.
% Graficar curva IDA mediana y la desviación estándar logarítmica de los
% desplazamientos como función de Sa(T1)

% mu_ln = median(log(x))
% median = geomean(x)
% std_ln = std(log(x))

% EDP_muln = mean(log(EDP')); % me da negativo usando lo del ppt
% EDP_median = geomean(EDP');                                                 % Mediana
EDP_stdln = std(log(EDP'));                                                 % Desviación estandar
EDP_muln = exp(log(EDP_median) + 0.5*EDP_stdln.^2);                         % Estimación de la Media

% Gráficos

% % Mediana logarítmica + Curvas IDA
% figure
% plot(EDP_median',IM(:,1)/g,'color','r','LineWidth',2)
% hold on
% plot(EDP,IM/g,'.-','color','#606060')
% xline(dy,'Color','b')
% ylim([0 IM_max])
% hold off
% xlabel('EDP y geomean(EDP)')
% ylabel('IM: Sa(T_1,\xi) [g]') 
% grid on
% title('Mediana logarítmica EDP')

% % Mediana logarítmica
% figure
% plot(IM(:,1)/g,EDP_median','.-','color','k')
% ylabel('Mediana logarítmica de EDP, geomean(EDP)')
% xlabel('IM: Sa(T_1,\xi) [g]') 
% grid on
% title('Mediana logarítmica EDP')

% Desviación Estándar Logarítmica
figure
plot(EDP_stdln',IM(:,1)/g,'.-','color','k')
xlabel('Desviación estandar logarítmica de EDP, std(log(EDP))')
ylabel('IM: Sa(T_1,\xi) [g]') 
grid on
title('Estructura 2')
legend('Desviación Estaándar logarítmica de EDP)')

% % Media
% figure
% plot(IM(:,1)/g,EDP_muln','.-','color','k')
% ylabel('Media logarítmica de EDP, mean(log(EDP))')
% xlabel('IM: Sa(T_1,\xi) [g]') 
% grid on
% title('Media logarítmica EDP')

% % Mediana +- sigma
% figure
% plot(IM(:,1)/g,EDP_median','.-','color','r')
% hold on
% plot(IM(:,1)/g,exp(log(EDP_median')+EDP_stdln'),'.-','color','k')
% plot(IM(:,1)/g,exp(log(EDP_median')-EDP_stdln'),'.-','color','k')
% hold off
% ylabel('Mediana logarítmica de EDP, geomean(EDP)')
% xlabel('IM: Sa(T_1,\xi) [g]') 
% grid on
% title('Mediana logarítmica EDP')

figure
plot(EDP_median',IM(:,1)/g,'.-','color','r')
hold on
plot(exp(log(EDP_median')+EDP_stdln'),IM(:,1)/g,'.-','color','k')
plot(exp(log(EDP_median')-EDP_stdln'),IM(:,1)/g,'.-','color','k')
plot(EDP_muln',IM(:,1)/g,'.-','color','b')
hold off
xlabel('Mediana y media logarítmica de EDP, geomean(EDP)')
ylabel('IM: Sa(T_1,\xi) [g]') 
legend('Mediana','Mediana+-\sigma','','Media')
grid on
title('Estructura 2')

%% 3
% Copiar gráficos media vs Sa, pero agregar curva de desplazamientos
% máximos que se obtendrían asumiendo principio de igualdad de
% desplazamientos

% Se cumpliría este principio solo en la primera estructura ya tiene periodo
% mayor a 1[s], pero no para la segunda estructura

% Principio de grandes desplazamientos
omega = 2*pi/T;
Sd = IM(:,1)./omega^2;  % IM m/s2 omega en rad2/s2 -> Sd = metros

% % Media
% figure
% plot(IM(:,1)/g,EDP_muln','.-','color','k')
% hold on
% plot(IM(:,1)/g,Sd)
% hold off
% ylabel('Media logarítmica de EDP, mean(log(EDP))')
% xlabel('IM: Sa(T_1,\xi) [g]') 
% grid on
% title('Media logarítmica EDP')
% legend('EDP_muln','Sd')

% Media
figure
plot(EDP_muln',IM(:,1)/g,'color','r')
hold on
plot(Sd,IM(:,1)/g,'.-','color','k')                                                        % Dejamos Sd en centimetros
hold off
xlabel('Desplazamiento [m]')
ylabel('IM: Sa(T_1,\xi) [g]') 
grid on
title('Estructura 2')
legend('Media','Sd(T_1)')


%% 4
% Suponiendo distribución lognormal para EDP dado IM, calcular probabilidad
% de que los desplazamientos máximos excedan 25cm para una ordenada
% espectral de 1g.

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
ResultsDir_Sa_avg = 'Resultados2';
ResultsName_Sa_avg = 'estructura2_sa_avg';

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
xlabel('EDP: Desplazamiento Máximo [m]')
ylabel('IM: Sa_avg [g]')                                                    % En IIDAP, se pusieron unidades de metros y segundos => aceleraciones en m/s2
title('Multi-Record IDA Curves', ' Estructura T_1 = 0.3s; \xi = 5%')
grid on

% Gráficos
% % Mediana logarítmica
% figure
% plot(IM(:,1)/g,EDP_median','.-','color','k')
% ylabel('Mediana logarítmica de EDP, geomean(EDP)')
% xlabel('IM: Sa_avg [g]') 
% grid on
% title('Mediana logarítmica EDP')

% Desviación Estándar Logarítmica
figure
plot(EDP_stdln',IM(:,1)/g,'.-','color','k')
xlabel('Desviación estandar logarítmica de EDP, std(log(EDP))')
ylabel('IM: Sa_avg [g]') 
grid on
title('Estructura 2')
legend('Desviación estándar logarítmica')

% % Media
% figure
% plot(IM(:,1)/g,EDP_muln','.-','color','k')
% ylabel('Media logarítmica de EDP, mean(log(EDP))')
% xlabel('IM: Sa_avg [g]') 
% grid on
% title('Media logarítmica EDP')

% % Mediana +- sigma
% figure
% plot(IM(:,1)/g,EDP_median','.-','color','r')
% hold on
% plot(IM(:,1)/g,exp(log(EDP_median')+EDP_stdln'),'.-','color','k')
% plot(IM(:,1)/g,exp(log(EDP_median')-EDP_stdln'),'.-','color','k')
% hold off
% ylabel('Mediana logarítmica de EDP, geomean(EDP)')
% xlabel('IM: Sa_avg [g]') 
% grid on
% title('Mediana logarítmica EDP')

figure
plot(EDP_median',IM(:,1)/g,'.-','color','r')
hold on
plot(exp(log(EDP_median')+EDP_stdln'),IM(:,1)/g,'.-','color','k')
plot(exp(log(EDP_median')-EDP_stdln'),IM(:,1)/g,'.-','color','k')
plot(EDP_muln',IM(:,1)/g,'.-','color','b')
hold off
xlabel('Mediana y media logarítmica de EDP, geomean(EDP)')
ylabel('IM: Sa_avg [g]') 
grid on
title('Estructura 2 pregunta 5')
legend('Mediana','Mediana+-\sigma','','Media')

% Volvemos a dejar como estaba antes
EDP = EDP_Sa;                                                               % EDP cuando IM es Sa(T1)
IM = IM_Sa;                                                                 % IM  = Sa(T1)

%% 6
% Extraer curvas de amenaza sísmica desde USGS para periodos en estudio
% Realizar interpolación lineal de los datos  de 0.1g
% Notar que no hay para T = 1.5sec -> Interpolar de forma hiperbólica en
% periodo (Sa prop. 1/T)

% USGS Data
USGS_data_Dir = 'Resultados2';
USGS_data_Name = 'USGS_data_est2';

USGS_data = readmatrix([USGS_data_Dir '\' USGS_data_Name]);
% Col1: Ground Motion, Col2: Lambda                                                                    % Periodo de las dos segundas columnas

% Hazazrd Curve (T = 1sec)
IM_03 = USGS_data(:,1);      % [g]
lambda_03 = USGS_data(:,2);

% Figura
figure
loglog(IM_03,lambda_03,'Color','r')
xlabel('Ground Motion [g] IM = Sa(T_1)')
ylabel('Annual Frequency of Excedence \lambda_{IM}')
title('Hazard Curves USGS')
grid on
legend('T = 0.3[s]')

% Interpolamos ambas para IM
IM_obj = IM_Sa(:,1)/g;
lambda_03_int = zeros(length(IM_obj),1);

for i = 1:length(IM_obj)
    lambda_03_int(i) = interp1(IM_03,lambda_03,IM_obj(i));
end

figure
loglog(IM_obj,lambda_03_int,'.-','color','r')
hold off
xlabel('Ground Motion [g] IM = Sa(T_1)')
ylabel('Annual Frequency of Excende \lambda_{IM}')
title('Hazard Curves USGS')
legend('T = 0.3[s]')
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
% IM_amenaza = IM_interp1(1:IM_max/IM_step);                                                    % 0.1:0.1:3   [g]
IM_amenaza = IM_obj;
IM_amenaza(1,:) = [];
lambda_amenaza = lambda_03_int;                                            % Interpolación
lambda_amenaza(1,:) = [];

% Desde Curvas IDA, las arreglamos para que tengan mismo paso que amenaza
% sismica

IM_IDA = IM_amenaza;                                                        % IM_IDA = IM_Sa(:,1)/g; IM_IDA(31:33)=[]; IM_IDA(0) = [];
EDP_IDA = EDP_Sa;
EDP_IDA(1,:) = [];

% Parámetros distribución lognormal
EDP_median = geomean(EDP_IDA')';                                            % Mediana
EDP_stdln = std(log(EDP_IDA'))';                                            % Desviación estandar
EDP_muln = exp(log(EDP_median) + 0.5*EDP_stdln.^2);                        % Estimación de la Media

% Obtención de abs(d/dIM (lambda_IM(x)))

% dlambdadIM = zeros(length(lambda_amenaza),1);
% for i = 1:length(lambda_amenaza)-1
%     dlambdadIM(i) = abs((lambda_amenaza(i+1)-lambda_amenaza(i))/(IM_amenaza(i+1)-IM_amenaza(i)));
% end

dlambdadIM = abs(diff(lambda_amenaza)./diff(IM_amenaza));

% figure
% loglog(IM_amenaza(2:end,1),dlambdadIM)
% % loglog(IM_amenaza,dlambdadIM)
% xlabel('IM = Sa(T_1) [g]')
% ylabel('$|\frac{d\lambda(IM=im)}{dIM}|$','interpreter','latex')
% title('Derivada de la curva de amenaza sísmica')
% grid on

% Rango de edp para lambda_EDP
edps = (0.001:0.0001:0.3)';

% Obtención de lambda_EDP
lambda_EDP = zeros(length(edps),1);
multiplicacion = zeros(length(IM_IDA),length(edps));
probabilidad = zeros(length(IM_IDA),length(edps));

for e = 1:length(edps)                                                      % Recorro los edps para obtener lambda_EDP(edp)
    edp = edps(e);
    for i = 1:length(IM_amenaza)-1                                          % Recorro los IM para sumarlos
        % Valor de im
        im = IM_amenaza(i);
        
        % Posición IM = 1g en curvas IDA 
        posIM = find(IM_IDA == im);                                         % En verdad es igual a "i"
        
        % EDP_median
        EDP_median_val = EDP_median(posIM);                                 % EDP_median obtenido para todas las IM=im en punto 2
        
        % EDP_stdln
        EDP_stdln_val = EDP_stdln(posIM);                                   % EDP_stdln_val otebnido para todas las IM=im en punto 2
        
        % EDP_muln
        EDP_muln_val = EDP_muln(posIM);
        % Parte EDP|IM
        % Asumir que EDP distribuye como lognormal en función de IM
        probabilidad(i,e) = 1 - normcdf(log(edp),log(EDP_median_val),EDP_stdln_val);

        % Parte IM
        dlambdadIM_val = dlambdadIM(i);

        % dIM
        % dIM = 0.1; % g (en inputs arriba)

        % Multiplicacion
        multiplicacion(i,e) = probabilidad(i,e)*dlambdadIM_val*dIM;
    end
    lambda_EDP(e) = sum(multiplicacion(:,e));
end

figure
semilogy(edps,lambda_EDP)
xlabel('EDP: Desplazamiento máximo [m]')
ylabel('\lambda_{edp}: Frecuencia anual de excedencia')
title('Curva de demanda sísmica (P6)')
grid on

% figure
% loglog(edps,lambda_EDP)
% xlabel('EDP: Desplazamiento máximo [m]')
% ylabel('\lambda_{edp}: Frecuencia anual de excedencia')
% title('Curva de demanda sísmica')
% grid on

edps_lineal = edps;
lambda_EDP_lineal = lambda_EDP;

%% 7
% Solo para estructura de T = 0.3s (segunda estructura)

%% 7 a
% Realizar interpolación lineal de la curva de amenaza sísmica lambda_IM en
% el espacio log-log. Graficar curva de amenaza sismica (lambda_IM vs IM)
% para ambas interpolaciones (lineal de la parte 6 y loglog de esta parte)
% comparar ambas curvas. repetir integración numérica del paso 6 con esta
% curva y grafique la curva de amenaza de EDP resultante

% Interpolacion loglog
IM_03_intlog = IM_Sa(:,1)/g;
IM_03_intlog(1) = [];
lambda_03_intlog = zeros(length(IM_03_intlog),1);
for i = 1:length(IM_03_intlog)
    lambda_03_intlog(i) = exp(interp1(log(IM_03),log(lambda_03),log(IM_03_intlog(i))));
end

figure
loglog(IM_obj,lambda_03_int,'.-')
hold on
loglog(IM_03_intlog,lambda_03_intlog,'.-')
hold off
xlabel('Ground Motion: IM = Sa(T_1) [g]')
ylabel('Annual Frequency of Excedence \lambda_{IM}')
title('Curva de Amenaza sísmica (P7a)')
grid on
legend('Interpolación Lineal','Interpolación Log-Log')

% Amenaza
IM_amenaza = log(IM_03_intlog);
lambda_amenaza = lambda_03_intlog;

% IDA
IM_IDA = IM_amenaza;
EDP_IDA = EDP_Sa;
EDP_IDA(1,:) = [];

% Parámetros distribución lognormal
EDP_median = geomean(EDP_IDA')';                                            % Mediana
EDP_stdln = std(log(EDP_IDA'))';                                            % Desviación estandar
EDP_muln = exp(log(EDP_median) + 0.5*EDP_stdln.^2);                         % Estimación de la Media

% abs(dlambda_IM=im/dIM)
dlambdadIM = abs(diff(lambda_amenaza)./diff(IM_amenaza));

% figure
% loglog(IM_amenaza(2:end,1),dlambdadIM)
% % loglog(IM_amenaza,dlambdadIM)
% xlabel('IM = Sa(T_1) [g]')
% ylabel('$|\frac{d\lambda(IM=im)}{dIM}|$','interpreter','latex')
% title('Derivada Curva Amenaza Sismica (P7a)')
% grid on

% Rango de edp para lambda_EDP
edps = (0.001:0.0001:0.3)';

% Obtención de lambda_EDP
lambda_EDP = zeros(length(edps),1);
multiplicacion = zeros(length(IM_IDA),length(edps));
probabilidad = zeros(length(IM_IDA),length(edps));

for e = 1:length(edps)                                                      % Recorro los edps para obtener lambda_EDP(edp)
    edp = edps(e);
    for i = 1:length(IM_amenaza)-1                                          % Recorro los IM para sumarlos
        % Valor de im
        im = IM_amenaza(i);
        
        % Posición IM = 1g en curvas IDA 
        posIM = find(IM_IDA == im);                                         % En verdad es igual a "i"
        
        % EDP_median
        EDP_median_val = EDP_median(posIM);                                 % EDP_median obtenido para todas las IM=im en punto 2
        
        % EDP_stdln
        EDP_stdln_val = EDP_stdln(posIM);                                   % EDP_stdln_val otebnido para todas las IM=im en punto 2
        
        % EDP_muln
        EDP_muln_val = EDP_muln(posIM);

        % Parte EDP|IM
        % Asumir que EDP distribuye como lognormal en función de IM
        probabilidad(i,e) = 1 - normcdf(log(edp),log(EDP_median_val),EDP_stdln_val);

        % Parte IM
        dlambdadIM_val = dlambdadIM(i);

        % dIM
        % dIM = 0.1; % g (en inputs arriba)

        % Multiplicacion
        multiplicacion(i,e) = probabilidad(i,e)*dlambdadIM_val*dIM;
    end
    lambda_EDP(e) = sum(multiplicacion(:,e));
end

figure
semilogy(edps_lineal,lambda_EDP_lineal)
hold on
semilogy(edps,lambda_EDP)
hold off
xlabel('EDP: Desplazamiento máximo [m]')
ylabel('\lambda_{edp}: Frecuencia anual de excedencia')
title('Curva de demanda sísmica (P7a)')
legend('Interpolación Lineal','Interpolación LogLog')
grid on

% figure
% loglog(edps,lambda_EDP)
% xlabel('EDP: Desplazamiento máximo [m]')
% ylabel('\lambda_{edp}: Frecuencia anual de excedencia')
% title('Curva de demanda sísmica')
% grid on

% Save Data
edps_loglog = edps;
lambda_EDP_loglog = lambda_EDP;

%% 7 b
% Ajustar polinomio de cuarto orden obtenida en 7.a
% Improved Analytical Representation of Seismic Hazard Curves (Miranda)

% Ajuste polinomio de interpolación loglog anterior (4 orden)
[P,S] = polyfit(log(IM_03_intlog),log(lambda_03_intlog),4);

% Guardamos
IM_Sa = IM;

% Simplificamos la varibale
IM = IM_03_intlog;

% lambda(IM=im) y |dmafimdim| de Miranda
lambda_03_7b = zeros(length(IM),1);
dmafimdim = zeros(length(IM),1);
for im = 1:length(IM)
    lambda_03_7b(im) = exp(P(5)+P(4)*log(IM(im))+P(3)*log(IM(im))^2+P(2)*log(IM(im))^3+P(1)*log(IM(im))^4);
    parte1 = (P(4)+2*P(3)*log(IM(im))+3*P(2)*log(IM(im))^2+4*P(1)*log(IM(im))^3)/IM(im);
    exponente = exp(lambda_03_7b(im));
    dmafimdim(im) = abs(parte1*exponente);
end

figure
loglog(IM_obj,lambda_03_int,'.-')
hold on
loglog(IM_03_intlog,lambda_03_intlog,'.-')
loglog(IM,lambda_03_7b)
hold off
xlabel('Ground Motion: IM = Sa(T_1) [g]')
ylabel('Annual Frequency of Excedence \lambda_{IM}')
title('Curva de Amenaza sísmica (P7b)')
legend('Interpolación Lineal','Interpolación Log-Log','Polinomio cuarto orden (Miranda)')
grid on

figure
loglog(IM,lambda_03_7b,'.-')
hold on
loglog(IM_03_intlog,lambda_03_intlog,'o')
hold off
xlabel('Ground Motion: IM = Sa(T_1) [g]')
ylabel('Annual Frequency of Excedence \lambda_{IM}')
title('Curva de Amenaza sísmica (P7b)')
grid on
legend('Polinomio Cuarto Orden (Miranda)','Interpolación Log-Log')
% Amenaza
IM_amenaza = IM;
lambda_amenaza = lambda_03_7b;

% IDA
IM_IDA = IM_amenaza;
EDP_IDA = EDP_Sa;
EDP_IDA(1,:) = [];

% Parámetros distribución lognormal
EDP_median = geomean(EDP_IDA')';                                            % Mediana
EDP_stdln = std(log(EDP_IDA'))';                                            % Desviación estandar
EDP_muln = exp(log(EDP_median) + 0.5*EDP_stdln.^2);                        % Estimación de la Media

% Rango de edp para lambda_EDP
edps = (0.001:0.0001:0.3)';

% Obtención de lambda_EDP
lambda_EDP = zeros(length(edps),1);
multiplicacion = zeros(length(IM_IDA),length(edps));
probabilidad = zeros(length(IM_IDA),length(edps));

for e = 1:length(edps)                                                      % Recorro los edps para obtener lambda_EDP(edp)
    edp = edps(e);
    for i = 1:length(IM_amenaza)-1                                          % Recorro los IM para sumarlos
        % Valor de im
        im = IM_amenaza(i);
        
        % Posición IM = 1g en curvas IDA 
        posIM = find(IM_IDA == im);                                         % En verdad es igual a "i"
        
        % EDP_median
        EDP_median_val = EDP_median(posIM);                                 % EDP_median obtenido para todas las IM=im en punto 2
        
        % EDP_stdln
        EDP_stdln_val = EDP_stdln(posIM);                                   % EDP_stdln_val otebnido para todas las IM=im en punto 2
        
        % EDP_muln
        EDP_muln_val = EDP_muln(posIM);

        % Parte EDP|IM
        % Asumir que EDP distribuye como lognormal en función de IM
        probabilidad(i,e) = 1 - normcdf(log(edp),log(EDP_median_val),EDP_stdln_val);

        % Parte IM
        dlambdadIM_val = dlambdadIM(i);

        % dIM
        % dIM = 0.1; % g (en inputs arriba)

        % Multiplicacion
        multiplicacion(i,e) = probabilidad(i,e)*dlambdadIM_val*dIM;
    end
    lambda_EDP(e) = sum(multiplicacion(:,e));
end

figure
semilogy(edps_lineal,lambda_EDP_lineal)
hold on
semilogy(edps_loglog,lambda_EDP_loglog)
semilogy(edps,lambda_EDP)
hold off
xlabel('EDP: Desplazamiento máximo [m]')
ylabel('\lambda_{edp}: Frecuencia anual de excedencia')
title('Curva de demanda sísmica (P7b)')
legend('Interpolación Lineal','Interpolación Log-Log','Polinomio cuarto orden (Miranda)')
grid on

% figure
% loglog(edps,lambda_EDP)
% xlabel('EDP: Desplazamiento máximo [m]')
% ylabel('\lambda_{edp}: Frecuencia anual de excedencia')
% title('Curva de demanda sísmica')
% grid on



