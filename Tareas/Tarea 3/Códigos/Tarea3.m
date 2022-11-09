%% Tarea 3 - Ingeniería Sísmica Avanzada
% Contreras - Sanguinetti

%% Inicializar
clear variables
close all
clc

%% Inputs
g = 9.81; % m/s2
cant_estr = 2;                                                              % Cantidad de estructutras
cant_variantes = 4;                                                         % Cantidad de variantes por cada estructura (toda estructura debe tener la misma cantidad de variantes)
cant_registros = 20;
IM_maximo = 20; % g
ResultsFiles = ["est_1_A";"est_1_B";"est_1_C";"est_1_D";"est_2_A";"est_2_B";"est_2_C";"est_2_D"];
ResultsDir = "IIDAP_T3";
n_ests = length(ResultsFiles);
IM_interpolacion = (0.1:0.1:IM_maximo).'; % g

ResultsFilesString= ["estructura 1A";"estructura 1B";"estructura 1C";"estructura 1D";"estructura 2A";"estructura 2B";"estructura 2C";"estructura 2D"];
n_interp = 100;                                                             % Cantidad de puntos en linspace para EDP interp1 para P2

%% Load Data
Data = struct();
for i = 1:n_ests
    [Data(i).EDP,Data(i).IM,Data(i).IMc,Data(i).Backbone] = getIdaCurves_v2_mod(convertStringsToChars(ResultsDir), convertStringsToChars(ResultsFiles(i)));
end

% Tenemos todos los datos para cada una de las 2 estructuras 4 variantes enumeradas del 1 al 8 
% (1,2,3,4 son est1 A,B,C,D y 5,6,7,8 son est2 A,B,C,D respectivamente)
% De haber una tercera estructura entonces 9,10,11,12 de est3 A,B,C,D
% Ej: Data(1).EDP = matriz EDP que retorna getIdaCurves_v2 para estructura 1 A

%% P1
% Graficar 
% Para graficar se necesita dy dp dpc du; Fy Fmax Fres

for i = 1:cant_estr
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
    end
    hold off
    xlabel('\delta [m]')
    ylabel('F [kN]')
    title('Backbone IMK Bilinear Model')
    legend('A','B','C','D')
    grid on
end
clear dy dp dpc du dres fy fmax fres X Y

%% P2
% Graficar curvas IDA
n_franjas = 200;                                                            % Cantidad de franjas para interpolación de IM (0.1g hasta el máximo del "escalamiento" del registro)
for i = 1:n_ests
    % Interpolación
    for n = 1:cant_registros
        EDP_rmm = rmmissing(Data(i).EDP(:,n));
        IM_rmm = rmmissing(Data(i).IM(:,n));
        Data(i).IM_interp1(:,n) = linspace(0.1,max(Data(i).IM(:,n)),n_franjas); % g
        Data(i).EDP_interp1(:,n) = interp1(IM_rmm,EDP_rmm,Data(i).IM_interp1(:,n),'linear','extrap');
%         Data(i).EDP_interp1(:,n) = interp1(Data(i).IM(:,n),Data(i).EDP(:,n),Data(i).IM_interp1(:,n),'linear','extrap');
%         Data(i).EDP_interp1_pchip(:,n) = interp1(Data(i).IM(:,n),Data(i).EDP(:,n),Data(i).IM_interp1(:,n),'pchip');
    end

    % Distribución
%     EDP_stdln = std(log(Data(i).EDP_interp1.'));                           % Sin dispersión porq no me da razonable, IDKW
    EDP_median = geomean(Data(i).EDP_interp1.','omitnan');                  % Creo que no es necesario omitnan ahora ya que interpolamos con rmmissing
    
    % Figura
    figure
    plot(Data(i).EDP,Data(i).IM,'.-','color','#606060')
    hold on
%     plot(Data(i).EDP_interp1,Data(i).IM_interp1,'color','#606060');
    plot(EDP_median.',Data(i).IM_interp1(:,1),'color','r','LineWidth',2)
%     plot(exp(log(EDP_median.')+EDP_stdln.'),Data(i).IM_interp1(:,1),'--','color','r','LineWidth',2)
%     plot(exp(log(EDP_median.')-EDP_stdln.'),Data(i).IM_interp1(:,1),'--','color','r','LineWidth',2)
    hold off
    xlabel('EDP: desplazamiento latereal [m]')
    ylabel('IM: Sa(T_1) [g]')
    title('Curvas IDA', ResultsFilesString(i))
    grid on
    legend(['Curvas IDA'; convertStringsToChars(repmat(string( ),cant_registros-1,1));'Mediana'])
end

%% P3
% Graficar las curvas de fragilidad de colapso para cada estructura
% (modelo)

% % Hay que identificar en qué EDP=edpc ocurre el colapso
% for i = 1:n_ests
%     edp_c(i) = identifyCollapseEDP(Data(i).EDP_interp1);
%     for n = 1:n_registros
%         if Data(i).EDP > edp_c 
% end


%% P4

% Nombre de archivos
USGSdataName_est1 = 'USGS_data_est1.xlsx';
USGSdataName_est2 = 'USGS_data_est2.xlsx';

% Importar matrices
USGS_est1 = readmatrix(USGSdataName_est1);
USGS_est2 = readmatrix(USGSdataName_est2);

% Gráficos
figure
loglog(USGS_est1(:,1),USGS_est1(:,2),'.-')
hold on
loglog(USGS_est2(:,1),USGS_est2(:,2),'.-')
hold off
xlabel('IM: Sa(T_1) [g]')
ylabel('Annual Frequency of Exceedence \lambda_c')
title('Hazard Curves USGS')
legend('T = 3.00 [sec]','T = 0.2 [sec]')
grid on

% Interpolación para generar más datos
IM_length = 100;                                                            % Cantidad de valores a interpolar
USGSinterp_est1(:,1) = logspace(-3,0,IM_length).';
USGSinterp_est2(:,1) = logspace(-3,0,IM_length).';

USGSinterp_est1(:,2) = exp(interp1(log(USGS_est1(:,1)),log(USGS_est1(:,2)),log(USGSinterp_est1(:,1))));
USGSinterp_est2(:,2) = exp(interp1(log(USGS_est2(:,1)),log(USGS_est2(:,2)),log(USGSinterp_est2(:,1))));

% % Borrar NaN
% USGSinterp_est1 = rmmissing(USGSinterp_est1);
% USGSinterp_est2 = rmmissing(USGSinterp_est2);

% Gráficos
figure
loglog(USGSinterp_est1(:,1),USGSinterp_est1(:,2),'.-')
hold on
loglog(USGSinterp_est2(:,1),USGSinterp_est2(:,2),'.-')
hold off
xlim([10^-3 10^1])
xlabel('IM: Sa(T_1) [g]')
ylabel('Annual Frequency of Exceedence \lambda_c')
title('Hazard Curves USGS', 'Interpolado')
legend('T = 3.00 [sec]','T = 0.2 [sec]')
grid on

%% P5

% IM_length = 100;                                                            % Cantidad de valores a interpolar

% Ajuste polinomio
[P1,S1] = polyfit(log(USGSinterp_est1(:,1)),log(USGSinterp_est1(:,2)),4);   % Busqueda de constantes del ajuste polinomial para estructura 1
[P2,S2] = polyfit(log(USGSinterp_est2(:,1)),log(USGSinterp_est2(:,2)),4);   % Busqueda de constantes del ajuste polinomial para estructura 2

% Generación de los valores debido al ajuste polinomial de cuarto orden
IM_vals = logspace(-3,0,IM_length).';                                       % Valores de IM para interpolar
lambda1_pf = zeros(IM_length,1);                                            % Inicializar lambda del polinomio est 1
lambda2_pf = zeros(IM_length,1);                                            % Inicializar almbda del polinomio est 2
for i = length(IM_vals)
    im = IM_vals(i);
    lambda1_pf(i) = exp(P1(5) + P1(4)*log(im) + P1(3)*log(im)^2 + P1(2)*log(im)^3 + P1(1)*log(im)^4); % Eq.3 Miranda Polinomio
    lambda2_pf(i) = exp(P2(5) + P2(4)*log(im) + P2(3)*log(im)^2 + P2(2)*log(im)^3 + P2(1)*log(im)^4); % Eq.3 Miranda Polinomio
end

figure
loglog(IM_vals,lambda1_pf)
hold on
loglog(IM_vals,lambda2_pf)
hold off
xlabel('IM: Sa(T_1) [g]')
ylabel('Annual Frequency of Exceedence \lambda_c')
title('Hazard Curves USGS')
legend('T = 3.00 [sec]','T = 0.2 [sec]')
grid on
