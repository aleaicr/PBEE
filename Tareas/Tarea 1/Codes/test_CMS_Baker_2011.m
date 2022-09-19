function [median_CMS,sigma_CMS] = test_CMS_Baker_2011(Ti,Tast,Sa_Tast,M,R,Vs30,mec_focal,region,z1)
% EJECUTAR FUNCIÓN CON ARCHIVO "run_CMS_Baker_2011.m"
% CMS de Baker 2011 (Conditional Mean Spectru: Tool for Ground Motion Selection)

% GMPE utilizada Boore et al (2014) (Archivo: BSSA_2014_nga.m)
% Correlación entre epsilon en diferentes periodos Baker & Jayaram (2008)

% Inputs
% Ti        -> Vector o lista de periodos para el análisis
% Tast      -> Periodo condicionante
% SaTast    -> Ordenada epsectral en el periodo condicionante
% M         -> Magnitud de momento para estimación en GMPE
% R         -> Distancia para estimación en GMPE (km)
% Vs30      -> Velocidad de ondas de corte del suelo (Clasificación) en km/h
% mec_focal -> Mecanismo focal de la falla (Fault_Type de BSSA_2014_nga.m)
% region    -> Región (region de BSSA_2014_nga.m)

% Outputs
% median_CMS -> Vector de medianas condicionadas
% sigma_CMS -> Vector de desviaciones estándars condicionadas
% (ambas condicionadas a Tast)

%% Step 1:
% Determine Target Sa at a Given Period Sa(T*) and the Associates M, R and epsilon

% Tast es input
% Sa_Tast es input
% M es input
% R es input
% epsilon_Tast se obtiene en Step 2

%% Step 2:
% Compute the Mean and Standard Deviation of the Response Spectrum, Given M
% and R

% Usamos BSSA_2014_nga (Descripción de Inputs en script de función)
Rjb = R;
Fault_Type = mec_focal;

Ti_length = length(Ti);                                                     % Tamaño del vector de periodos
mu_lnSa = zeros(Ti_length,1);                                               % Inicializar vectores para no reasignar en memoria
sigma_lnSa = zeros(Ti_length,1);


% Media y Sigma para todos los periodos de análisis
for i = 1:Ti_length
    [median_BSSA, sigma_BSSA, ~] = BSSA_2014_nga(M, Ti(i), Rjb, Fault_Type, region, z1, Vs30);    % No nos importa period1
    mu_lnSa(i) = log(median_BSSA);
    sigma_lnSa(i) = log(sigma_BSSA);
    if Ti(i) == Tast
        disp('lol xd')
        mu_Tast = log(median_BSSA);
        sigma_Tast = log(sigma_BSSA);
    end
end
lnSa_Tast  = log(Sa_Tast);
epsilon_Tast = (lnSa_Tast - mu_Tast)/sigma_Tast;

%% Desde ahora se realizan Step 3 y Step 4 para todos los Ti que se quieren analizar
% Inicializar vectores (para no reasignar espacio en memoria)
rho = zeros(Ti_length,1);
mu_eTi_eTast = zeros(Ti_length,1);
mu_lnSaTi_lnSaTast = zeros(Ti_length,1);
sigma_lnSaTi_lnSaTast = zeros(Ti_length,1);

for Ti_pos = 1:Ti_length
    Ti_val = Ti(Ti_pos);
    %% Step 3:
    % Compute epsilon at Other Periods
    % Calculamos rho para cada periodo
    rho(Ti_pos) = rho_Baker_2011(Ti_val,Tast);                               % rho(Ti) = rho(Ti,Tast)
    mu_eTi_eTast(Ti_pos) = rho(Ti_pos)*epsilon_Tast;                           % de qué me sirve?

    %% Step 4:
    % Compute CMS
    % Resolviendo Eq.1 para lnSa(T) reproduce la siguiente eq para cada Ti
    % dado lnSa(T*)
    mu_lnSaTi_lnSaTast(Ti_pos) = mu_lnSa(Ti_pos) + rho(Ti_pos)*epsilon_Tast*sigma_lnSa(Ti_pos);
    sigma_lnSaTi_lnSaTast(Ti_pos) = sigma_lnSa(Ti_pos)*sqrt(1-(rho(Ti_pos))^2);
end
median_CMS = mu_lnSaTi_lnSaTast;
sigma_CMS = sigma_lnSaTi_lnSaTast;
end
