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

ResultsFile = ["est_1_A";"est_1_B";"est_1_C";"est_1_D";"est_2_A";"est_2_B";"est_2_C";"est_2_D"];
ResultsDir = "IIDAP_T3";
n_ests = length(ResultsFile);

%% Load Data
Data = struct();
for i = 1:n_ests
    [Data(i).EDP,Data(i).IM,Data(i).IMc,Data(i).Backbone] = getIdaCurves_v2(convertStringsToChars(ResultsDir), convertStringsToChars(ResultsFile(i)));
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
        dpc = Data((i-1)*cant_variantes+j).Backbone.Upc_pos0;
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

%% P2
% Graficar curvas IDA
for i = 1:n_ests
    figure
    plot(Data(i).EDP,Data(i).IM./g,'.-','color','#606060')
    xlabel('EDP: desplazamiento latereal [m]')
    ylabel('IM: Sa(T_1) [g]')
    title('Curvas IDA')
    grid on
end

%% P3
% Graficar las curvas de fragilidad de colapso para cada estructura
% (modelo)

% Por separado
for i = 1:n_estr
    figure
    plot()
end

% Todas juntas
for i = 1:cant_estr
    figure
    for j = 1:cant_variantes

    end
end