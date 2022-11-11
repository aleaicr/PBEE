%% Tarea 3 - Ingeniería Sísmica Avanzada
% Contreras - Sanguinetti

% Archivos desde IIDAP

%% Inicializar
clear variables
close all
clc

%% Inputs
g = 9.81; % m/s2
cant_estr = 2;                                                              % Cantidad de estructutras
cant_variantes = 4;                                                         % Cantidad de variantes de capacidad de cada estructura (toda estructura debe tener la misma cantidad de variantes)
cant_registros = 20;                                                        % Cantidad de registros para curvas IDA (ingresado a IIDAP)
IM_maximo = 20; % g                                                         % IM máximo considerado para curvas IDA
ResultsFiles = ["est_1_A";"est_1_B";"est_1_C";"est_1_D";"est_2_A";"est_2_B";"est_2_C";"est_2_D"];   % Nombre de archivos a analizar
ResultsDir = "IIDAP_T3";                                                    % Dirección de carpeta donde están los archivos
n_ests = length(ResultsFiles);
% IM_interpolacion = (0.1:0.1:IM_maximo).'; % g
Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';                                    % Letra para cada estructura (para la leyenda de gráficos)
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
% Graficar curvas de capacidad monotónica de las 8 figuras (2 estructuras, 4 variantes de IMK bilineal)

% Para graficar se necesita dy dp dpc du; Fy Fmax Fres
A = cell(cant_variantes,cant_estr);
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
    legend(['Curvas IDA'; convertStringsToChars(repmat(string( ),cant_registros-1,1));'Mediana'])
end

% % Todas las medianas juntas
% figure
% plot(Data.EDP_median,Data.IM_interp1(:,1))
% xlabel('EDP: Desplazamiento [m]')
% ylabel('IM: Sa(T_1) [g]')
% title('Medianas')

% Separadas por estructura
for i = 1:cant_estr
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
for i = 1:n_ests
   Data(i).EDP_median2(:,1) = geomean(Data(i).EDP.','omitnan').';                  % Creo que no es necesario omitnan ahora ya que interpolamos con rmmissing
    figure
    plot(Data(i).EDP,Data(i).IM,'.-','color','#606060')
    hold on
    plot(Data(i).EDP_median2(:,1),Data(i).IM(:,14),'color','r','LineWidth',2)
    hold off
    xlabel('EDP: desplazamiento latereal [m]')
    ylabel('IM: Sa(T_1) [g]')
    title('Curvas IDA', ResultsFilesString(i))
    grid on
    legend(['Curvas IDA'; convertStringsToChars(repmat(string( ),cant_registros-1,1));'Mediana'])
end

for i = 1:cant_estr
    figure
    hold on
    for j = 1:cant_variantes
        plot(Data((i-1)*cant_variantes+j).EDP_median2,Data((i-1)*cant_variantes+j).IM(:,14))
    end
    hold off
    xlabel('EDP: Desplazamiento [m]')
    ylabel('IM: Sa(T_1) [g]')
    title('Medianas')
    legend(A{:,i})
    grid on
end


%% P3
% Graficar las curvas de fragilidad de colapso para cada estructura

frag = struct();                                                            % Struct de la curva de fragilidad
for i = 1:n_ests                                                            % Para cada estructura i
    % Contar cantidad de simulaciones que producen colapso para cada IM
    all_IM = rmmissing(Data(i).IM(:));                                      % Todos los IM que están en esa estructura, pero con ceros
    all_IM(all_IM == 0) = [];                                               % Eliminar todos los ceros del vector
    frag(i).all_IM = sort(unique(all_IM));                                  % Todos los IM que están en esa estructura, No se me ocurrió cómo hacerlo sin trabajar antes con un vector fuera del struct
    frag(i).zj_col = zeros(length(frag(i).all_IM),1);                       % Inicializar vector de zj observaciones de colapso para la IM j
    frag(i).nj_sim = zeros(length(frag(i).all_IM),1);                       % Inicializar vector nj simulaciones para la IM j
    for n = 1:cant_registros                                                % Ahora vemos para cada registro n de la estructura i
        frag(i).edp_c(n,1) = max(rmmissing(Data(i).EDP(:,n)));              % Buscamos el edp que corresponde al colapso para este registro
        for r = 1:length(rmmissing(Data(i).EDP(:,n)))                       % Para cada EDP = edp
            pos = find(frag(i).all_IM == Data(i).IM(r,n));                  % Definimos la posición del IM(edp)
            frag(i).nj_sim(pos,1) = frag(i).nj_sim(pos,1) + 1;              % Sumamos un 1 ya que observamos una simulación en el IM j
            if Data(i).EDP(r,n) >= frag(i).edp_c(n,1)                       % pero Si mi EDP = edp > edp_c entonces
                frag(i).zj_col(pos,1) = frag(i).zj_col(pos,1) + 1;          % Sumamos un 1 ya que obsevmos un colapso en el IM j
            end
        end                                                                 % Siguiente EDP=edp 
    end                                                                     % Siguiente registro n
    % Fracción
    frag(i).fraccion(:,1) = frag(i).zj_col./frag(i).nj_sim;
    
    % mu_ln
    frag(i).muln = log(mean(frag(i).fraccion(:,1)));
    
    % stdln
    frag(i).stdln = log(std(frag(i).fraccion(:,1)));
    for im = 1:length(frag(i).all_IM(:,1))
        frag(i).Pc(im,1) = normcdf((log(frag(i).all_IM(im,1)) - frag(i).muln)./frag(i).stdln);
        frag(i).pbinom(im,1) = nchoosek(frag(i).nj_sim(im,1),frag(i).zj_col(im,1))*(frag(i).Pc(im,1).^frag(i).zj_col(im,1))*(1-frag(i).Pc(im,1)).^(frag(i).nj_sim(im,1) - frag(i).zj_col(im,1)); % P Binomial
    end
    frag(i).ver = exp(sum(frag(i).pbinom(:,1)));                            % Pitatoria() = exp(log(pitatoria)) = exp(sumatoria)
    
    



    figure
    plot(frag(i).all_IM,logncdf(frag(i).all_IM,frag(i).muln,frag(i).stdln))
    xlabel('IM: Sa(T_1) [g]')
    ylabel('Probabilidad de Colapso')
    title('Curva de fragilidad de colapso')
    legend(ResultsFilesString(i))
    grid on
end                                                                         % Ya que vimos todos los registros, vamos a la siguiente estructura

% Descripción del Algoritmo
% Vemos una estructura, luego vemos registro por registro, para el primer
% registro se guarda la máxima EDP (que corresponde al colapso
% coincidentemente porque el IIDAP para en el colapso) y luego se recorre
% cada EDP = edp, vemos el primero, comparamos que no es el colapso y
% guardamos que observamos un registro (nj) en la intensidad IMj, luego 
% vemos el segundo, el tercero.... hasta llegar al último, ese si 
% corresponde al colapso por lo que lo contamos y guardamos que observamos
% un registro, también que observamos un colapso y 




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
legend('USGS T = 3.00 [sec]','USGS T = 0.20 [sec]')
grid on

% Interpolación para generar más datos
IM_length = 100;                                                            % Cantidad de valores a interpolar
USGSinterp_est1(:,1) = linspace(0.0025,3.28,IM_length).';                   % Mejor modificarlo a meno, si no, puede que el ajuste polinomial de Miranda salga mal
USGSinterp_est2(:,1) = linspace(0.0025,7.3,IM_length).';                   % Los límites deben ser de la parte curva

USGSinterp_est1(:,2) = exp(interp1(log(USGS_est1(:,1)),log(USGS_est1(:,2)),log(USGSinterp_est1(:,1))));
USGSinterp_est2(:,2) = exp(interp1(log(USGS_est2(:,1)),log(USGS_est2(:,2)),log(USGSinterp_est2(:,1))));

% % Borrar NaN
USGSinterp_est1 = rmmissing(USGSinterp_est1);
USGSinterp_est2 = rmmissing(USGSinterp_est2);

% Gráficos
figure
loglog(USGSinterp_est1(:,1),USGSinterp_est1(:,2),'.-')
hold on
loglog(USGSinterp_est2(:,1),USGSinterp_est2(:,2),'.-')
hold off
xlim([10^-3 10^1])
xlabel('IM: Sa(T_1) [g]')
ylabel('Annual Frequency of Exceedence \lambda_im')
title('Hazard Curves USGS', 'Interpolado')
legend('USGS_Interp T = 3.00 [sec]','USGS_Interp T = 0.20 [sec]')
grid on

USGS = struct();
USGS(1).interpvals = USGSinterp_est1;                                       % primera columna son IM, segunda columna son lambda (MAF)
USGS(2).interpvals = USGSinterp_est2;                                       % Primera columna son IM, segunda columna MAF

%% P5
% Coeficientes del ajuste polinomial
[P1,S1] = polyfit(log(USGSinterp_est1(:,1)),log(USGSinterp_est1(:,2)),4);   % Busqueda de constantes del ajuste polinomial para estructura 1
[P2,S2] = polyfit(log(USGSinterp_est2(:,1)),log(USGSinterp_est2(:,2)),4);   % Busqueda de constantes del ajuste polinomial para estructura 2

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
