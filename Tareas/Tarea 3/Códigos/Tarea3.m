%% Tarea 3 - Ingeniería Sísmica Avanzada
% Contreras - Sanguinetti

%% Inicializar
clear variables
close all
clc

%% Inputs
g = 9.81; % m/s2
cant_estr = 2;                                                              % Cantidad de estructutras
cant_variantes = 4;                                                         % Cantidad de variantes por cada estructura
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
    [Data(i).EDP,Data(i).IM,Data(i).IMc,Data(i).Backbone] = getIdaCurves_v2(convertStringsToChars(ResultsDir), convertStringsToChars(ResultsFiles(i)));
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
n_franjas = 200;
for i = 1:n_ests
    for n = 1:cant_registros
        EDP_rmm = rmmissing(Data(i).EDP(:,n));
        IM_rmm = rmmissing(Data(i).IM(:,n));
        Data(i).IM_interp1(:,n) = linspace(0.1,max(Data(i).IM(:,n)),n_franjas); % g
        Data(i).EDP_interp1(:,n) = interp1(IM_rmm,EDP_rmm,Data(i).IM_interp1(:,n),'linear','extrap');
%         Data(i).EDP_interp1(:,n) = interp1(Data(i).IM(:,n),Data(i).EDP(:,n),Data(i).IM_interp1(:,n),'linear','extrap');
%         Data(i).EDP_interp1_pchip(:,n) = interp1(Data(i).IM(:,n),Data(i).EDP(:,n),Data(i).IM_interp1(:,n),'pchip');
    end
%     EDP_stdln = std(log(Data(i).EDP_interp1.'));                           % Sin dispersión porq no me da razonable, IDKW
    EDP_median = geomean(Data(i).EDP_interp1.','omitnan');
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
    legend('Curvas IDA','Mediana','','Mediana +- \sigma')
end

%% P3
% Graficar las curvas de fragilidad de colapso para cada estructura
% (modelo)

% Hay que identificar en qué EDP=edpc ocurre el colapso
for i = 1:n_ests
    edp_c(i) = identifyCollapseEDP(Data(i).EDP_interp1);
    for n = 1:n_registros
        if Data(i).EDP > edp_c 
end