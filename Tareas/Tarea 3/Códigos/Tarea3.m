%% Tarea 3 - Ingeniería Sísmica Avanzada
% Contreras - Sanguinetti

% Archivos desde IIDAP

%% Inicializar
clear variables
close all
clc

%% Inputs
g = 9.81; % m/s2
cant_estructuras = 2;                                                              % Cantidad de estructutras
cant_variantes = 4;                                                         % Cantidad de variantes de capacidad de cada estructura (toda estructura debe tener la misma cantidad de variantes)
cant_registros = 20;                                                        % Cantidad de registros para curvas IDA (ingresado a IIDAP)
IM_maximo = 20; % g                                                         % IM máximo considerado para curvas IDA
n_interp = 100;                                                             % Cantidad de puntos en linspace para EDP interp1 para P2
ResultsFiles = ["est_1_A";"est_1_B";"est_1_C";"est_1_D";"est_2_A";"est_2_B";"est_2_C";"est_2_D"];   % Nombre de archivos a analizar
ResultsDir = "IIDAP_T3";                                                    % Dirección de carpeta donde están los archivos
n_modelos = length(ResultsFiles);                                              % Número total de 

Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';                                    % Letra para cada estructura (para la leyenda de gráficos)
ResultsFilesString= ["estructura 1A";"estructura 1B";"estructura 1C";"estructura 1D";"estructura 2A";"estructura 2B";"estructura 2C";"estructura 2D"];
Colors = ["#0072BD";"#D95319";"#EDB120";"#7E2F8E";"#77AC30";"#4DBEEE";"#A2142F"]; % hay 7 colores, si se quieren más, fijarlos

IM_range = (0.1:0.1:6).';
sequence = zeros(cant_estructuras,cant_variantes);
for i = 1:cant_estructuras
    sequence(i,:) = (i-1)*cant_variantes + (1:1:cant_variantes);
end

%% Load Data
Data = struct();
for i = 1:n_modelos
    [Data(i).EDP,Data(i).IM,Data(i).IMc,Data(i).Backbone,Data(i).nGM] = getIdaCurves_v2_mod(convertStringsToChars(ResultsDir), convertStringsToChars(ResultsFiles(i)));
end

% Tenemos todos los datos para cada una de las 2 estructuras 4 variantes enumeradas del 1 al 8 
% (1,2,3,4 son est1 A,B,C,D y 5,6,7,8 son est2 A,B,C,D respectivamente)
% De haber una tercera estructura entonces 9,10,11,12 de est3 A,B,C,D
% Ej: Data(1).EDP = matriz EDP que retorna getIdaCurves_v2 para estructura 1 A

%% P1
% Graficar curvas de capacidad monotónica de las 8 figuras (2 estructuras, 4 variantes de IMK bilineal)

% Para graficar se necesita dy dp dpc du; Fy Fmax Fres
A = cell(cant_variantes,cant_estructuras);
for i = 1:cant_estructuras
    figure
    hold on
    for j = 1:cant_variantes
        dy = Data((i-1)*cant_variantes+j).Backbone.Uy_pos0;
        dp = Data((i-1)*cant_variantes+j).Backbone.Up_pos0;
%         dpc = Data((i-1)*cant_variantes+j).Backbone.Upc_pos0;
        du = Data((i-1)*cant_variantes+j).Backbone.Uu_pos0;
        dres = Data((i-1)*cant_variantes+j).Backbone.Ures_pos0;
        fy = Data((i-1)*cant_variantes+j).Backbone.Fy_pos0;
        fmax = Data((i-1)*cant_variantes+j).Backbone.FmaxFy_pos0*fy;
        fres = Data((i-1)*cant_variantes+j).Backbone.FresFy_pos0*fy; 
        X = [0; dy; dp; dres; du; du];
        Y = [0; fy; fmax; fres; fres; 0];
        plot(X,Y)
        A{j,i} = Alphabet(j);
    end
    hold off
    xlabel('\delta [m]')
    ylabel('F [kN]')
    title('Backbone IMK Bilinear Model')
    legend(A{:,i})
    grid on
end
clear dy dp dpc du dres fy fmax fres X Y

%% P2 INTERPOLANDO
% Graficar curvas IDA y la mediana (análisis entre 0.1g y 20g y límite de
% desplazamiento muy grande) comparar curvas de mediana de desplazamiento
n_franjas = 200;                                                            % Cantidad de franjas para interpolación de IM (0.1g hasta el máximo del "escalamiento" del registro)
for i = 1:n_modelos
    % Interpolación
    for n = 1:Data(i).nGM
        EDP_rmm = rmmissing(Data(i).EDP(:,n));
        IM_rmm = rmmissing(Data(i).IM(:,n));
        Data(i).IM_interp1(:,n) = linspace(0.1,max(Data(i).IM(:,n)),n_franjas); % g
        Data(i).EDP_interp1(:,n) = interp1(IM_rmm,EDP_rmm,Data(i).IM_interp1(:,n),'linear','extrap');
%         Data(i).EDP_interp1(:,n) = interp1(Data(i).IM(:,n),Data(i).EDP(:,n),Data(i).IM_interp1(:,n),'linear','extrap');
%         Data(i).EDP_interp1_pchip(:,n) = interp1(Data(i).IM(:,n),Data(i).EDP(:,n),Data(i).IM_interp1(:,n),'pchip');
    end

    % Distribución
%     EDP_stdln = std(log(Data(i).EDP_interp1.'));                           % Sin dispersión porq no me da razonable, IDKW
    Data(i).EDP_median(:,1) = geomean(Data(i).EDP_interp1.','omitnan').';                  % Creo que no es necesario omitnan ahora ya que interpolamos con rmmissing
    
    % Figura
    figure
    plot(Data(i).EDP,Data(i).IM,'.-','color','#606060')
    hold on
%     plot(Data(i).EDP_interp1,Data(i).IM_interp1,'color','#606060');
    plot(Data(i).EDP_median(:,1),Data(i).IM_interp1(:,1),'color','r','LineWidth',2)
%     plot(exp(log(EDP_median.')+EDP_stdln.'),Data(i).IM_interp1(:,1),'--','color','r','LineWidth',2)
%     plot(exp(log(EDP_median.')-EDP_stdln.'),Data(i).IM_interp1(:,1),'--','color','r','LineWidth',2)
    hold off
    xlabel('EDP: desplazamiento latereal [m]')
    ylabel('IM: Sa(T_1) [g]')
    title('Curvas IDA', ResultsFilesString(i))
    grid on
    legend(['Curvas IDA'; convertStringsToChars(repmat(string( ),Data(i).nGM-1,1));'Mediana'])
end

% % Todas las medianas juntas
% figure
% plot(Data.EDP_median,Data.IM_interp1(:,1))
% xlabel('EDP: Desplazamiento [m]')
% ylabel('IM: Sa(T_1) [g]')
% title('Medianas')

% Separadas por estructura
for i = 1:cant_estructuras
    figure
    hold on
    for j = 1:cant_variantes
        plot(Data((i-1)*cant_variantes+j).EDP_median,Data((i-1)*cant_variantes+j).IM_interp1(:,1))
    end
    hold off
    xlabel('EDP: Desplazamiento [m]')
    ylabel('IM: Sa(T_1) [g]')
    title('Medianas')
    legend(A{:,i})
    grid on
end


 %% P2 SIN INTERPOLAR
% for i = 1:n_ests
%    Data(i).EDP_median2(:,1) = geomean(Data(i).EDP.','omitnan').';                  % Creo que no es necesario omitnan ahora ya que interpolamos con rmmissing
%     figure
%     plot(Data(i).EDP,Data(i).IM,'.-','color','#606060')
%     hold on
%     plot(Data(i).EDP_median2(:,1),Data(i).IM(:,14),'color','r','LineWidth',2)
%     hold off
%     xlabel('EDP: desplazamiento latereal [m]')
%     ylabel('IM: Sa(T_1) [g]')
%     title('Curvas IDA', ResultsFilesString(i))
%     grid on
%     legend(['Curvas IDA'; convertStringsToChars(repmat(string( ),cant_registros-1,1));'Mediana'])
% end
% 
% for i = 1:cant_estr
%     figure
%     hold on
%     for j = 1:cant_variantes
%         plot(Data((i-1)*cant_variantes+j).EDP_median2,Data((i-1)*cant_variantes+j).IM(:,14))
%     end
%     hold off
%     xlabel('EDP: Desplazamiento [m]')
%     ylabel('IM: Sa(T_1) [g]')
%     title('Medianas')
%     legend(A{:,i})
%     grid on
% end


%% P3
% Graficar las curvas de fragilidad de colapso para cada estructura


for i = 1:cant_estructuras     % 2
    % Datos de curva de fragilidad de colapso
    figure
    hold on
    for j = 1:cant_variantes        % 4
        Data((i-1)*cant_variantes+j).muln = mean(log(Data((i-1)*cant_variantes+j).IMc).');
        Data((i-1)*cant_variantes+j).stdln = std(log(Data((i-1)*cant_variantes+j).IMc).');
        Data((i-1)*cant_variantes+j).PC = normcdf((log(Data((i-1)*cant_variantes+j).IMc).'-Data((i-1)*cant_variantes+j).muln)./Data((i-1)*cant_variantes+j).stdln);
        % Contar
        Data((i-1)*cant_variantes+j).nj = Data(i).nGM;
        Data((i-1)*cant_variantes+j).zj = 0;
        for r = 1:length(Data((i-1)*cant_variantes+j))
            for p = 1:length(Data(i).IMc)
                if Data(i).IMc(r,1) >= max(Data(i).IM(:,p))
                    Data((i-1)*cant_variantes+j).zj = Data((i-1)*cant_variantes+j).zj + 1;
                end
            end
        end
        Data((i-1)*cant_variantes+j).fraccion = Data((i-1)*cant_variantes+j).zj/Data((i-1)*cant_variantes+j).nj;
        plot(Data((i-1)*cant_variantes+j).IMc,Data((i-1)*cant_variantes+j).PC,'o','color',Colors(j))
    end
    hold off
    xlabel('IM: Sa(T_1) [g]')
    ylabel('P(C|IM=im)')
    legend(A{:,i})
    title('Datos de Probabilidad de Colapso')
    grid on
    
    % Curva de fragilidad de colapso (lognormal)
    figure
    hold on
    for j = 1:cant_variantes
        frag((i-1)*cant_variantes+j).curvFrag = logncdf(IM_range, Data((i-1)*cant_variantes+j).muln,Data((i-1)*cant_variantes+j).stdln);
        plot(Data((i-1)*cant_variantes+j).IMc,Data((i-1)*cant_variantes+j).PC,'o','color',Colors(j))
        plot(IM_range,frag((i-1)*cant_variantes+j).curvFrag,'color',Colors(j))
    end
    hold off
    xlabel('IM: Sa(T_1) [g]')
    ylabel('P(C|IM=im)')
    legend(A{:,i})
    title('Curva de Fragilidad de Colapso')
    grid on
end


%% P4
% Nombre de archivos
USGSdataName_est1 = 'USGS_data_est1.xlsx';
USGSdataName_est2 = 'USGS_data_est2.xlsx';

% Importar matrices
USGS_est1 = readmatrix(USGSdataName_est1);
USGS_est2 = readmatrix(USGSdataName_est2);
% Como USGS_est1(end,2) = 0, lo eliminamos
USGS_est1(end,:) = []; 

hazard = struct();
hazard(1).IM = USGS_est1(:,1);
hazard(2).IM = USGS_est2(:,1);
hazard(1).lambda = USGS_est1(:,2);
hazard(2).lambda = USGS_est2(:,2);

% Gráficos
figure
loglog(USGS_est1(:,1),USGS_est1(:,2),'.-')
hold on
loglog(USGS_est2(:,1),USGS_est2(:,2),'.-')
hold off
xlabel('IM: Sa(T_1) [g]')
ylabel('Annual Frequency of Exceedence \lambda_c')
title('Hazard Curves USGS')
legend('USGS T = 3.00 [sec]','USGS T = 0.20 [sec]')
grid on

% Interpolación para generar más datos
IM_length = 100;                                                            % Cantidad de valores a interpolar
USGSinterplog_est1(:,1) = linspace(0.0025,3.28,IM_length).';                   % Mejor modificarlo a meno, si no, puede que el ajuste polinomial de Miranda salga mal
USGSinterplog_est2(:,1) = linspace(0.0025,7.3,IM_length).';                   % Los límites deben ser de la parte curva

USGSinterplog_est1(:,2) = exp(interp1(log(USGS_est1(:,1)),log(USGS_est1(:,2)),log(USGSinterplog_est1(:,1))));
USGSinterplog_est2(:,2) = exp(interp1(log(USGS_est2(:,1)),log(USGS_est2(:,2)),log(USGSinterplog_est2(:,1))));

% % Borrar NaN
USGSinterplog_est1 = rmmissing(USGSinterplog_est1);
USGSinterplog_est2 = rmmissing(USGSinterplog_est2);

% 2% en 50 años
lambda_2_50_annos = 1/2475;
IM_2_50_annos_est1 = interp1(USGSinterplog_est1(:,2),USGSinterplog_est1(:,1),lambda_2_50_annos);
IM_2_50_annos_est2 = interp1(USGSinterplog_est2(:,2),USGSinterplog_est2(:,1),lambda_2_50_annos);

% Gráficos
figure
loglog(USGSinterplog_est1(:,1),USGSinterplog_est1(:,2),'.-')
hold on
loglog(USGSinterplog_est2(:,1),USGSinterplog_est2(:,2),'.-')
loglog(IM_2_50_annos_est1,lambda_2_50_annos,'o','color','r')
loglog(IM_2_50_annos_est2,lambda_2_50_annos,'o','color','r')
text(3*10^-1,10^-4,string(IM_2_50_annos_est1))
text(3,10^-4,string(IM_2_50_annos_est2))
hold off
xlim([10^-3 10^1])
xlabel('IM: Sa(T_1) [g]')
ylabel('Annual Frequency of Exceedence \lambda_im')
title('Hazard Curves USGS', 'Interpolado')
legend('USGS_Interp T = 3.00 [sec]','USGS_Interp T = 0.20 [sec]','2% en 50 años (\lambda = 1/2475)')
grid on

USGS = struct();
USGS(1).interpvals = USGSinterplog_est1;                                       % primera columna son IM, segunda columna son lambda (MAF)
USGS(2).interpvals = USGSinterplog_est2;                                       % Primera columna son IM, segunda columna MAF

%% P5
% Coeficientes del ajuste polinomial
[P1,~] = polyfit(log(USGSinterplog_est1(:,1)),log(USGSinterplog_est1(:,2)),4);   % Busqueda de constantes del ajuste polinomial para estructura 1
[P2,~] = polyfit(log(USGSinterplog_est2(:,1)),log(USGSinterplog_est2(:,2)),4);   % Busqueda de constantes del ajuste polinomial para estructura 2

% Generación de los valores debido al ajuste polinomial de cuarto orden
lambda1_pf = zeros(IM_length,1);                                            % Inicializar lambda del polinomio est 1 (lambda = MAF: Mean Annual Frequency)
lambda2_pf = zeros(IM_length,1);                                            % Inicializar almbda del polinomio est 2
dldim1 = zeros(IM_length,1);                                                % Inicializar derivada
dldim2 = zeros(IM_length,1);                                                % Inicializar derivada

for i = 1:IM_length
    % Ajuste polinomial primera estructura
    im = USGS(1).interpvals(i,1);
    lambda1_pf(i) = exp(P1(5) + P1(4)*log(im) + P1(3)*log(im)^2 + P1(2)*log(im)^3 + P1(1)*log(im)^4); % Eq.3 Miranda Polinomio
    parte1 = (P1(4)+2*P1(3)*log(im)+3*P1(2)*log(im)^2+4*P1(1)*log(im)^3)/im;
    dldim1(i) = abs(parte1*exp(lambda1_pf(i)));
    
    % Ajuste polinomial segunda estructura
    im = USGS(2).interpvals(i,1);
    lambda2_pf(i) = exp(P2(5) + P2(4)*log(im) + P2(3)*log(im)^2 + P2(2)*log(im)^3 + P2(1)*log(im)^4); % Eq.3 Miranda Polinomio
    parte1 = (P2(4)+2*P2(3)*log(im)+3*P2(2)*log(im)^2+4*P2(1)*log(im)^3)/im;
    dldim2(i) = abs(parte1*exp(lambda2_pf(i)));
end

% Esto que escribí lo tira al revés xd
% IM_vals = logspace(-3,0,IM_length).';                                       % Valores de IM para interpolar
% lambda1_pf = exp(P1(5)*ones(IM_length,1) + P1(4).*log(IM_vals) + P1(3)*log(IM_vals).^2 + P1(4)*log(IM_vals).^3 + P1(5)*log(IM_vals).^4);
% lambda2_pf = exp(P2(5)*ones(IM_length,1) + P2(4).*log(IM_vals) + P2(3)*log(IM_vals).^2 + P2(4)*log(IM_vals).^3 + P2(5)*log(IM_vals).^4);

figure
loglog(USGS_est1(:,1),USGS_est1(:,2),'.-')
hold on
loglog(USGS_est2(:,1),USGS_est2(:,2),'.-')
loglog(USGS(1).interpvals(:,1),lambda1_pf)
loglog(USGS(2).interpvals(:,1),lambda2_pf)
hold off
xlabel('IM: Sa(T_1) [g]')
ylabel('Annual Frequency of Exceedence \lambda_im (MAF)')
title('Hazard Curves USGS')
legend('USGS T = 3.00 [sec]','USGS T = 0.20 [sec]','Polinomial Fit T = 3.00 [sec]','Polyfit T = 0.20 [sec]')
grid on

figure
loglog(USGS(1).interpvals(:,1),dldim1);
hold on
loglog(USGS(2).interpvals(:,1),dldim2);
hold off
xlabel('IM: Sa(T_1) [g]')
ylabel('Derivada lambda_im')
title('Derivada De curva de Amenaza sísmica')
legend('T = 3.00 [sec]', 'T = 0.20 [Sec]')
grid on


%% P6
% Caluclar y graficar la desagregación (por IM=im) de la tasa anual media
% de colapso de las estructuras 1A y 2A, considerando los siguientes casos
% para la deriva de la curva de amenaza sísmica

%% P6 a Interpolación lineal

% Inputs iniciales
lambda_col = zeros(n_modelos,1);                                            % Lambdas de Colapso
desCol = struct();                                                          % Desagregación del colapso
mod_analizar = [1 5];                                                       % Número del modelo que se quiere analizar (primeras variantes de cada estructura)

for r = 1:cant_estructuras                                                  % Cantidad de curvas de amenaza (1 para cada estructura, 1 estructura es un periodo => dos estructuras)
    j = mod_analizar(r);                                                    % Número del modelo que quiero analizar
    [i,~] = find(sequence == j);                                            % Buscar qué estructura es dependiendo del número del modelo
    hazard(i).IMlineal = IM_range;
    hazard(i).lambdalineal = interp1(hazard(i).IM,hazard(i).lambda,hazard(i).IMlineal);
    hazard(i).absdiff = abs(diff(hazard(i).lambdalineal)./diff(hazard(i).IMlineal));
end

for r = 1:numel(mod_analizar)                                               % Para los modelos que quiero analizar
    j = mod_analizar(r);                                                    % Número del modelo que quiero analizar
    [i,~] = find(sequence == j);                                            % Buscar qué estructura es dependiendo del número del modelo
    figure
    % Curva de Fragilidad P(C|IM=im)
    sp1 = subplot(3,1,1);
    plot(IM_range,frag(j).curvFrag,'color','r')
    grid on

    % abs de Derivada de amenaza sísmica |dlambda_IM(IM=im)/dIM=im|
    sp2 = subplot(3,1,2);
    plot(hazard(i).IMlineal, [hazard(i).absdiff; hazard(i).absdiff(end,1)]) % Repito el último número de la derivada para que sean 70 elementos
    grid on

    % Desagregación del Riesgo de Colapso P(C|IM=im)*|dlambda_IM(IM=im)/dIM=im|
    sp3 = subplot(3,1,3);
    desCol(j).Desagregacion_lineal = frag(j).curvFrag.*[hazard(i).absdiff; hazard(i).absdiff(end,1)];
    plot(hazard(i).IMlineal,desCol(j).Desagregacion_lineal)
    grid on

    % Tasa anual media de colapso
    lambda_col(j,1) = trapz(hazard(i).IMlineal,frag(j).curvFrag.*[hazard(i).absdiff; hazard(i).absdiff(end,1)]);
    xlabel('IM = Sa(T_1)')
    ylabel(sp1,'P(C|IM=im)')
    ylabel(sp2,'|d\lambda_IM/dIM|')
    ylabel(sp3,'P(C|IM=im)*|d\lambda_IM/dIM|')
    sgtitle(['Desagregación Tasa Anual Media de Colapso', 'Interpolación Lineal'])
end


%% P6. b Interpolación loglog

lambda_col_log = zeros(n_modelos,1);                                            % Lambdas de Colapso
for r = 1:cant_estructuras
    j = mod_analizar(r);                                                    % Número del modelo que quiero analizar
    [i,~] = find(sequence == j);                                            % Buscar qué estructura es dependiendo del número del modelo
    hazard(i).IMlog = IM_range;
    hazard(i).lambdalog =  exp(interp1(log(hazard(i).IM),log(hazard(i).lambda),log(hazard(i).IMlog)));
    hazard(i).absdifflog = abs(diff(hazard(i).lambdalog)./diff(hazard(i).IMlog));
end

for r = 1:numel(mod_analizar)                                               % Para los modelos que quiero analizar
    j = mod_analizar(r);                                                    % Número del modelo que quiero analizar
    [i,~] = find(sequence == j);                                            % Buscar qué estructura es dependiendo del número del modelo
    figure
    % Curva de Fragilidad P(C|IM=im)
    sp1 = subplot(3,1,1);
    plot(IM_range,frag(j).curvFrag,'color','r')
    grid on

    % abs de Derivada de amenaza sísmica |dlambda_IM(IM=im)/dIM=im|
    sp2 = subplot(3,1,2);
    plot(hazard(i).IMlog, [hazard(i).absdifflog; hazard(i).absdifflog(end,1)]) % Repito el último número de la derivada para que sean 70 elementos
    grid on

    % Desagregación del Riesgo de Colapso P(C|IM=im)*|dlambda_IM(IM=im)/dIM=im|
    sp3 = subplot(3,1,3);
    desCol(j).Desagregacion_log = frag(j).curvFrag.*[hazard(i).absdifflog; hazard(i).absdifflog(end,1)];
    plot(hazard(i).IMlog,desCol(j).Desagregacion_log)
    grid on

    % Tasa anual media de colapso
    lambda_col_log(j,1) = trapz(hazard(i).IMlog,frag(j).curvFrag.*[hazard(i).absdifflog; hazard(i).absdifflog(end,1)]);
    xlabel('IM = Sa(T_1)')
    ylabel(sp1,'P(C|IM=im)')
    ylabel(sp2,'|d\lambda_IM/dIM|')
    ylabel(sp3,'P(C|IM=im)*|d\lambda_IM/dIM|')
    sgtitle(['Desagregación Tasa Anual Media de Colapso','Interpolación LogLog'])
end

%% P6.c Ajuste Polinomial
lambda_col_poli = zeros(n_modelos,1);

% [P1,~] = polyfit(log(hazard(1).IMlog),log(hazard(1).lambdalog),4);   % Busqueda de constantes del ajuste polinomial para estructura 1
% [P2,~] = polyfit(log(hazard(2).IMlog),log(hazard(2).lambdalog),4);   % Busqueda de constantes del ajuste polinomial para estructura 2

for r = 1:cant_estructuras
    j = mod_analizar(r);                                                    % Número del modelo que quiero analizar
    [i,~] = find(sequence == j);                                            % Buscar qué estructura es dependiendo del número del modelo
    hazard(i).IMpoli = IM_range;
    for iv = 1:length(IM_range)
        im = IM_range(iv,1);
        if i == 1
            P = P1;
        elseif i == 2
            P = P2;
        end
        hazard(i).lambdapoli(iv,1) = exp(P(5) + P(4)*log(im) + P(3)*log(im)^2 + P(2)*log(im)^3 + P(1)*log(im)^4); % Eq.3 Miranda Polinomio
        parte1 = (P(4) + 2*P(3)*log(im) + 3*P(2)*log(im)^2 + 4*P(1)*log(im)^3)/im;
        hazard(i).absdiffpoli(iv,1) = abs(parte1*exp(log(hazard(i).lambdapoli(iv,1))));
    end
end

for r = 1:numel(mod_analizar)                                               % Para los modelos que quiero analizar
    j = mod_analizar(r);                                                    % Número del modelo que quiero analizar
    [i,~] = find(sequence == j);                                            % Buscar qué estructura es dependiendo del número del modelo
    figure
    % Curva de Fragilidad P(C|IM=im)
    sp1 = subplot(3,1,1);
    plot(IM_range,frag(j).curvFrag,'color','r')
    grid on

    % abs de Derivada de amenaza sísmica |dlambda_IM(IM=im)/dIM=im|
    sp2 = subplot(3,1,2);
    plot(hazard(i).IMpoli, hazard(i).absdiffpoli)
    grid on

    % Desagregación del Riesgo de Colapso P(C|IM=im)*|dlambda_IM(IM=im)/dIM=im|
    sp3 = subplot(3,1,3);
    desCol(j).Desagregacion_poli = frag(j).curvFrag.*hazard(i).absdiffpoli;
    plot(hazard(i).IMpoli,desCol(j).Desagregacion_poli)
    grid on

    % Tasa anual media de colapso
    lambda_col_poli(j,1) = trapz(hazard(i).IMpoli,frag(j).curvFrag.*hazard(i).absdiffpoli);
    xlabel('IM = Sa(T_1)')
    ylabel(sp1,'P(C|IM=im)')
    ylabel(sp2,'|d\lambda_IM/dIM|')
    ylabel(sp3,'P(C|IM=im)*|d\lambda_IM/dIM|')
    sgtitle(['Desagregación Tasa Anual Media de Colapso','Ajuste polinomio cuarto orden'])
end


% Plot ambos juntos
% No parametrizado, a la rápida no mas

figure
plot(hazard(1).IMlineal,desCol(1).Desagregacion_lineal)
hold on
plot(hazard(i).IMlog, desCol(1).Desagregacion_log)
plot(hazard(1).IMlog, desCol(1).Desagregacion_poli)
hold off
ylabel(sp3,'P(C|IM=im)*|d\lambda_IM/dIM|')
sgtitle(['Desagregación Tasa Anual Media de Colapso','Ajuste polinomio cuarto orden'])
legend('Interpolación Lineal', 'Interpolación LogLog', 'Ajuste Polinomial Miranda')
grid on
