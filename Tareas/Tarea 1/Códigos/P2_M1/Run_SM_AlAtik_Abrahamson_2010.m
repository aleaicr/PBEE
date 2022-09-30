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
Tn_init = 0.05;                                                             % Mínimo TI que indica Baker (2011)
Tn_step = 0.01;                                                             % Paso para los periodos
Tn_final = 5;
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

% Error admisible
error_adm = 10^-5;

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
Reg = 1.3*Reg;                                                      % Escalamos el registro

[fj,da,a1] = SM_AlAtik_Abrahamson_2010(Reg,dt,Tast,xi,PSaObj_Tast,error_adm,beta_Newmark,ui,udi,gamma_relax);  % Realizamos SM
[~,~,~,~,PSa_a1,~,~,~] = PSa_peak(beta_Newmark,xi,dt,ui,udi,a1,Tn');        % Obtenemos PSa_a1 en cm/s2

a1 = a1./g;
PSa_a1 = PSa_a1./g;                                                         % Pasamos de cm/s2 a g
Reg = Reg./g;                                                               % Pasamos de cm/s2 a g
PSa_reg = PSa_reg./g;                                                       % Pasamos de cm/s2 a g

%% Gráficamos
% % Registro
% figure('Name','Registro Sísmico')
% plot(t_vect,Reg)
% xlabel('Tiempo t (sec)')
% ylabel('Aceleración [g]')
% grid on
% title('Registro de aceleraciones')
% legend('RSN1227 CHICHI CHY074-N')

% % PSa
% figure('Name', 'PSa')
% plot(Tn,PSa_reg)
% xlabel('Periodos T (sec)')
% ylabel('Pseudo-espectro de aceleraciones PSa (g)')
% grid on
% title('Pseudo-Espectro de Aceleraciones')
% legend('RSN1227 CHICHI CHY074-N')

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

% % a1
% figure('Name','a1')
% plot(t_vect,a1)
% xlabel('Tiempo t [sec]')
% ylabel('Registro de Aceleraciones Modificado a_1 [g]')
% title('Registro de Aceleraciones Modificado a_1')
% legend('Registro con Spectral-Matching')
% grid on

% % PSa_a1
% figure('Name', 'PSa a_1')
% plot(Tn,PSa_a1)
% xlabel('Periodos T (sec)')
% ylabel('Pseudo-espectro de aceleraciones PSa a_1 [g]')
% grid on
% title('Pseudo-Espectro de Aceleraciones a_1 Al Atik & Abrahamson (2010)')
% legend('Registro con Spectral-Matching')

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
title('Comparación de Registros')

% Comparación PSa
figure('Name', 'Comparación PSa')
plot(Tn,PSa_reg)
hold on
plot(Tn,PSa_a1)
hold off
xlabel('Periodo (T) [sec]')
ylabel('PSa [g]')
legend('RSN1227 CHICHI CHY074-N','Registro con Spectral-Matching')
title('Comparación de Espectros de Respuesta')
grid on


