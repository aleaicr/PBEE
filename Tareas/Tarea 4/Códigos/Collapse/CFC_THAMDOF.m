%% Curva de Fragilidad de Colapso
%%   NO FUNCIONA  %%
%%   NO FUNCIONA  %%
%%   NO FUNCIONA  %%
%%   NO FUNCIONA  %%

% Contreras - Sanguinetti
% Ingeniería Sísmica Avanzada - USM 2022

% Determinar la curva de fragilidad de colapso por:
%   - Método de Mínimos cuadrados
%   - Método de máxima verosimilitud

% Comentario:
% Desde THAMDOF, para cada IM, ver cuantos provocaron colapso.

%% Inicializar
clear variables
close all
clc

%% Inputs
IM_stripes = [0.1; 0.4; 0.6; 0.8]; % g                                      % Medidas de intensidad de las franjas de análisis (Sa_T1)
nj = [0; 0; 5; 11];                                                         % Cantidad de registros que provocaron el colapso por cada IM (observaciones de colapso)
zj = [20; 20; 20; 20];                                                      % Cantidad de registro por cada IM (simulaciones)
muln_min = - 2;
muln_step = 0.01;
muln_max = 2;
sigma_min = -1;
sigma_step = 0.01;
sigma_max = 1;
IM_plot_range = (0: 0.1: 3).';                                              % Rango de IM para las figuras al final
%% Calculos previos
frac_PCol = nj./zj;
IM_length = length(IM_stripes);
muln_range = (muln_min:muln_step:muln_max).';
sigma_range = (sigma_min:sigma_step:sigma_max).';

%% Método mínimos cuadrados
E = zeros(IM_length,1);
mc = +inf;
for m = 1:length(muln_range)
    muln = muln_range(m);
    for s = 1:length(sigma_range)
        sigma = sigma_range(s);
        for i = 1:IM_length
            E(i) = abs(frac_PCol(i) - normcdf((log(IM_stripes(i)) - muln)./sigma))^2;
            if sum(E) < mc
                mc = sum(E);
                sigma_mc = sigma;
                muln_mc = muln;
            end
        end
    end
end

%% Método de máxima verosimilitud
ver = zeros(IM_length,1);
ver_ = zeros(IM_length,1);
mv = -inf;
for m = 1:length(muln_range)
    muln = muln_range(m);
    for s = 1:length(sigma_range)
        sigma = sigma_range(s);
        for i = 1:IM_length
            parte1 = nchoosek(zj(i),nj(i));
            parte2 = normcdf((log(IM_stripes(i))-muln)./sigma)^zj(i);
            parte3 = (1-normcdf((log(IM_stripes(i))-muln)./sigma))^(nj(i)-zj(i));
            ver_(i) = parte1*parte2*parte3;
            ver = exp(sum(ver_));
            if ver > mv                                           % Pitatoria() = exp(log(Pitatoria()) = exp(sumatoria)
                mv = ver;
                sigma_mv = sigma;
                muln_mv = muln;
            end
        end
    end
end


%% Mostrar resultados
tabla = table();
tabla.metodo = ["Mínimos cuadrados"; "Máxima Verosimilitud"];
tabla.valor = [mc; mv];
tabla.muln = [muln_mc; muln_mv];
tabla.sigma = [sigma_mc; sigma_mv];
disp(tabla)
clear tabla

%% Plot
figure
plot(IM_stripes,frac_PCol,'o')
hold on
plot(IM_plot_range,normcdf(IM_plot_range,muln_mc,sigma_mc),'color','r')
plot(IM_plot_range,normcdf(IM_plot_range,muln_mv,sigma_mv),'color','k')
hold off
grid on
legend('Data','Mínimos cuadrados','Máxima verosimilitud')
title('Curva de fragilidad de colapso')
ylabel('P(C|IM = im)')
xlabel('IM = Sa(T_1) [g]')

%% MALO
% %% Inputs
% collapse_SF = 'Collapse SF.xlsx';                                           % Nombre de archivo que contiene los SF para llegar al colapso
% GMFolder = '../Registros';                                                  % Nombre de la carpeta donde se encuentran los registros (todos los registros deben estar en [g] -> las que da el profe están en [g])
% GMDataName = 'GM Data.txt';                                                 % Nombre del archivo .txt con los dt de sampling de cada registro
% cant_registros = 20;
% T = 2; % sec                                                                % Periodo fundamental de la estructura
% beta_newmark = 1/4;                                                         % Beta del método de Newmark-beta (dejarlo como 1/4 para que sea incondicionalmente estable)
% xi = 0.05;                                                                  % Amortiguamiento para el espectro
% Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
% 
% % Aceleración espectral para Registros con SF = 1 (solo los que están en la carpeta, sin escalar para todas las franjas)
% Sa = Sa_or(GMFolder,GMDataName,T,xi,beta_newmark);
% 
% % SF Colapso
% SF = readmatrix(collapse_SF,'Sheet','colSF','Range',convertStringsToChars("A3:" + Alphabet(20) + "3"));
% 
% % IM colapso
% Sa_collapse = Sa.*SF.'; % g
% 
% % Median
% collapse_median = geomean(Sa_collapse);
% 
% % std
% collapse_stdln = std(log(Sa_collapse));
% 
% % plot
% figure
% plot(sort(Sa_collapse),logncdf(log(sort(Sa_collapse)),collapse_median,collapse_stdln),'.-')


