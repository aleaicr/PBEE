%% Run of Spectral Matching Al Atik & Abrahamson 2010
% Ingeniería Sísmica Avanzada -  Tarea 1
% Contreras - Sanguinetti

% An Improved Method for Nonstationary Spectral Matching
% Linda Al Atik & Norman Abrahamson (2010)

% PSa se obtiene con el Método de Newmark (Contreras - Adasme (Ing.Sísimca 2022))

%% Inicializar
clear variables
close all
clc

%% Inputs
g = 980; % g = 980cm/s2

% Registro
load('RSN1227_CHICHI_CHY074-N.mat')
Reg = Acc; clear Acc Npts
% Reg   -> Registro de Aceleraciones [g]
% dt    -> Paso temporal del registro de aceleraciones [sec]
t_vect = (0:dt:(length(Reg)-1)*dt)';                                        % Vector con tiempos asociados a cada elemento del vector registro
Reg = Reg*g;    % Pasamos el registro de [g] a cm/s2

% CMS
load('median_CMS.mat');

% Periodos
Tn_init = 0;                                                                % La media condicionada no puede empezar desde antes de 0.05 según Baker(2011)
Tn_step = 0.1;
Tn_final = 3;
Tn = Tn_init:Tn_step:Tn_final;
Tast = 1; % sec                                                             % Periodo de ajuste (Condicionante?)
Tast_pos = find(Tn == Tast);

% Razón de amortiguamiento
xi = 0.05;

% Condiciones iniciales
ui = 0;
udi = 0;

% PSa Objetivo
PSaObj_Tast = median_CMS(Tast_pos); % g                                     % Valor objetivo del espectro de pseudo-aceleración condicionado en T*
PSa_CMS_Tast = PSaObj_Tast; % g                                             

% Error admisible
error_adm = 10^-10;

% Beta para método de Newmark
beta_Newmark  = 1/4;                                                        % Aceleración promedio constante
% gamma_Newmark = 0.5;

% Parámetro de relajación
gamma_relax = 1;                                                            % Enunciado

%% Run Function
PSaObj_Tast = PSaObj_Tast*g;

% Escalar
[~,~,~,~,PSa_reg,~,~,~] = PSa_peak(beta_Newmark,xi,dt,ui,udi,Reg,Tn');      % cm/s2 PSa del Registro
fact_escala = PSaObj_Tast/PSa_reg(Tast_pos);                                % Factor de escala
fprintf('El factor de escala es %.3f \n \n',fact_escala)                    % Mostrar factor de escala
Reg = 1*Reg;                                                      % Escalamos el registro

% An Improved Method for Nonstationary Spectral Matching
% Linda Al Atik & Norman Abrahamson (2010)

% Inputs
% Reg               -> Registro de aceleraciones a modificar
% dt                -> Espaciamiento de tiempo de muetro
% Tast              -> Periodo(s) al cual el espectro será modificado (Periodo(s) de ajuste T*) (de la Estructura 1GDL)
% xi                -> Razón de amortiguamiento del espectro a modificar (de la estructura 1GDL)
% val_Obj           -> Valor objetivo del espectro de pseudo-aceleración en T*
% error_objetivo    -> Error admisible para la diferencia entre PSa_Tast obtenido con onda modificada y el val_Obj_Tast
% beta_Newmark      -> Beta del método de Newmark
% ui, udi           -> Condiciones iniciales de desplazamiento y velocidad respectivamente
% gamma_relax       -> Parámetro de relajación

% Ouputs
% fj                -> Registro de aceleración de la onda modificada
% da                -> La onda modificada (Eq.2 del artículo)
% a1                -> Registro de aceleraciones modificado (Eq. 11 del artículo, gamma = 1)

% El método de Newmark ocupa gamma_Newmark = 0.5

% Parámetros útiles
N = length(Tast);                                                           % N = 1,según enunciado, Número total de puntos espectrales (pares: frecuencia,amortiguamiento), podría ser length(Tast) para considerar varios y hacer la suma, luego se tendría que modificar el código
Reg_length = length(Reg);                                                   % Largo del registro
t_length = Reg_length;                                                      % Cantidad de tiempos, evidentemente t_length = Reg_length

% Inicializar valores y vectores
hi = zeros(t_length,1);
fj = zeros(t_length,1);

% Pseudo-Espectro de aceleraciones y tiempo al peak para el periodo de la estructura

[~,~,~,~,PSa_Tast,~,tj,signPeak] = PSa_peak(beta_Newmark,xi,dt,ui,udi,Reg,Tast);   % Obtenemos solo el PSa en el periodo de interés

for j = 1:N   % (j: spectral point)                                         % Recorremos periodo-amortiguamiento (Solo 1 en este caso)
    % Como N = 1, no se van a ocupar vectores, pero cada uno debería
    % depender del vector periodo-amortiguamiento (i), el vector de tiempo (t)
    % y del vector de periodos a ajustar (j)
    wi = 2*pi/Tast;                                                         % Frecuencia circular de estructura de 1GDL                                                        
    betai = xi;                                                             % Razón de amortiguamiento de estructura 1GDL
    wi_prima = wi*sqrt(1-betai^2);                                          % (Eq.6) Frecuencia circular amortiguada de estructura de 1GDL
    freq = wi_prima/(2*pi);                                                 % (Frecuencia)
    gamma_f = 1.178*freq^-0.93;                                             % (Eq.17) Coeficiente que depende de la frecuencia (para Eq.16)
    Dtj = atan(sqrt(1-betai^2)/betai)/wi_prima;                             % (Eq.14)
    for t = 1:t_length  % (t
        hi(t) = -wi/(sqrt(1-betai^2))*exp(-wi*betai*(t-1)*dt)*sin(wi_prima*(t-1)*dt);    % (Eq. 18)
        fj(t) = cos(wi_prima*((t-1)*dt-tj+Dtj))*exp(-(((t-1)*dt - tj + Dtj)/gamma_f)^2);   % (Eq. 16)
    end
end

% Formar vector hi(ti-tau) para integrar en (Eq.7)
n_ti = round(tj/dt,0) + 1;
hi_titau = [flip(hi(1:n_ti)); zeros(length(hi)-n_ti,1)];

% cij
cij = trapz(dt,fj.*hi_titau);                                               % (Eq.7)

% dR
% Eq.9 es útil solo si dRi es el error espectral de ajuste (spectral misfit)
dRi = (PSaObj_Tast-PSa_Tast)*signPeak;                                    % (Eq.1)

% b
b = (cij^-1)*dRi;                                                           % (Eq.10)

% da: Onda
da = b.*fj;                                                                 % (Eq.2)

% a1: Registro de aceleraciones modificado
a1 = Reg + gamma_relax*da;                                                  % (Eq.11)

% Verificación error, si no cumple, iteramos otra vez
[~,~,~,~,PSamod,~,~,~] = PSa_peak(beta_Newmark,xi,dt,ui,udi,a1,Tast);       % Obtenemos PSa del registro modificado  

if abs(PSamod - PSaObj_Tast) > error_adm
    Reg = a1;   % reg = a(t) de (Eq.11)
    disp(1)
    [fj,a1,da,cij,b,dRi] = Copy_of_SM_AlAtik_Abrahamson_2010(Reg,dt,Tast,xi,PSaObj_Tast,error_adm,beta_Newmark,ui,udi,gamma_relax);
end
[~,~,~,~,PSa_a1,~,~,~] = PSa_peak(beta_Newmark,xi,dt,ui,udi,a1,Tn');        % Obtenemos PSa_a1 en cm/s2

a1 = a1./g;
PSa_a1 = PSa_a1./g;                                                         % Pasamos de cm/s2 a g
Reg = Reg./g;                                                               % Pasamos de cm/s2 a g
PSa_reg = PSa_reg./g;                                                       % Pasamos de cm/s2 a g

%% Gráficamos
% Registro
figure('Name','Registro Sísmico')
plot(t_vect,Reg)
xlabel('Tiempo t (sec)')
ylabel('Aceleración [g]')
grid on
title('Registro de aceleraciones')
legend('RSN1227 CHICHI CHY074-N')

% PSa
figure('Name', 'PSa')
plot(Tn,PSa_reg)
xlabel('Periodos T (sec)')
ylabel('Pseudo-espectro de aceleraciones PSa (g)')
grid on
title('Pseudo-Espectro de Aceleraciones')
legend('RSN1227 CHICHI CHY074-N')

% fj
figure('Name','fj')
plot(t_vect,fj)
xlabel('Tiempo t [sec]')
ylabel('Función de ajuste f_j')
title('Función de ajuste f_j')
legend('Tapered Cosine Wavelet')
grid on

% da
figure('Name','da')
plot(t_vect,da)
xlabel('Tiempo t [sec]')
ylabel('Onda Modificadora \delta a')
title('Onda Modificadora \delta a')
legend('Registro con Spectral-Matching')
grid on

% a1
figure('Name','a1')
plot(t_vect,a1)
xlabel('Tiempo t [sec]')
ylabel('Registro de Aceleraciones Modificado a_1 [g]')
title('Registro de Aceleraciones Modificado a_1')
legend('Registro con Spectral-Matching')
grid on

% PSa_a1
figure('Name', 'PSa a_1')
plot(Tn,PSa_a1)
xlabel('Periodos T (sec)')
ylabel('Pseudo-espectro de aceleraciones PSa a_1 [g]')
grid on
title('Pseudo-Espectro de Aceleraciones a_1 Al Atik & Abrahamson (2010)')
legend('Registro con Spectral-Matching')

% Comparación Registro
figure('Name', 'Comparación Registro')
plot(t_vect,Reg)
hold on
plot(t_vect,a1)
hold off
grid on
xlabel('tiempo (t) [sec]')
ylabel('Aceleración [g]')
legend('RSN1227 CHICHI CHY074-N','Registro con Spectral-Matching')

% Comparación PSa
figure('Name', 'Comparación PSa')
plot(Tn,PSa_reg)
hold on
plot(Tn,PSa_a1)
hold off
xlabel('Periodo (T) [sec]')
ylabel('PSa [g]')
legend('RSN1227 CHICHI CHY074-N','Registro con Spectral-Matching')
grid on


