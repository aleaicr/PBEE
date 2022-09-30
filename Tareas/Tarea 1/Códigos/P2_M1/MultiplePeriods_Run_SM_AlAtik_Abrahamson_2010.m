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

% Periodos (mismo que del median_CMS.mat)
Tn_init = 0.05;                                                             % Mínimo TI que indica Baker (2011)
Tn_step = 0.01;                                                             % Paso para los periodos
Tn_final = 5;
Tn = Tn_init:Tn_step:Tn_final;

%% MULTIPLES PERIODOS, definir vector
Tast = [0.8;1]; % sec                                                             % Periodo de ajuste (Condicionante?)
Tast_length = length(Tast);
Tast_posi = zeros(Tast_length,1);

% Identificar posición de todos los Tast para cada 
for i = 1:length(Tn)
    for j = 1:length(Tast)
        if Tn(i) == Tast(j)
            Tast_posi(j) = i;
        end
    end
end

% Razón de amortiguamiento
xi = 0.05;

% Condiciones iniciales
ui = 0;
udi = 0;

% PSa Objetivo para cada periodo de ajuste igual al del CMS
PSaObj_Tast = zeros(Tast_length,1);
clear j
for j_val = 1:Tast_length
    Tastposicion = Tast_posi(j_val);
    medianCMS_val = median_CMS(Tastposicion);
    PSaObj_Tast(j_val) = medianCMS_val;
end

% Error admisible
error_adm = 10^-5;

% Beta para método de Newmark
beta_Newmark  = 1/4;                                                        % Aceleración promedio constante
% gamma_Newmark = 0.5;

% Parámetro de relajación
gamma_relax = 1;                                                            % Enunciado

%% Run Function
PSaObj_Tast = PSaObj_Tast*g;                                                % en cm/s2

% Escalar
RegOriginal = Reg;                                                          % Registro Original 
% [~,~,~,~,PSa_regOriginal,~,~,~] = PSa_peak(beta_Newmark,xi,dt,ui,udi,Reg,Tn');      % cm/s2 PSa del Registro
% fact_escala = PSaObj_Tast/PSa_regOriginal(Tast_pos);                                % Factor de escala
% fprintf('El factor de escala es %.3f \n \n',fact_escala)                    % Mostrar factor de escala
Reg = 1.3*Reg;                                                              % Escalamos el registro
RegEscalado = Reg;                                                          % Registro Escalado
[~,~,~,~,PSa_regEsc,~,~,~] = PSa_peak(beta_Newmark,xi,dt,ui,udi,Reg,Tn');      % cm/s2 PSa del Registro

[fj,da,a1] = MultiplePeriods_SM_AlAtik_Abrahamson_2010(Reg,dt,Tast,xi,PSaObj_Tast,error_adm,beta_Newmark,ui,udi,gamma_relax);  % Realizamos SM
[~,~,~,~,PSa_a1,~,~,~] = PSa_peak(beta_Newmark,xi,dt,ui,udi,a1,Tn');        % Obtenemos PSa_a1 en cm/s2

a1 = a1./g;
PSa_a1 = PSa_a1./g;                                                         % Pasamos de cm/s2 a g
RegEscalado = RegEscalado./g;                                                               % Pasamos de cm/s2 a g
RegOriginal = RegOriginal./g;
PSa_regOriginal = PSa_regOriginal./g;                                                       % Pasamos de cm/s2 a g
PSa_regEsc = PSa_regEsc./g;

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
plot(t_vect,RegEscalado,'color','r')
hold on
plot(t_vect,a1,'color','#007d79')
hold off
grid on
xlabel('tiempo (t) [sec]')
ylabel('Aceleración [g]')
legend('RSN1227 CHICHI CHY074-N Escalado','Registro con Spectral-Matching a_{1,final}')
title('Comparación de Registros')

% Comparación PSa
figure('Name', 'Comparación PSa')
plot(Tn,median_CMS,'color','r')
hold on
plot(Tn,PSa_regOriginal,'color','c')
plot(Tn,PSa_a1,'color','#007d79')
hold off
xlabel('Periodo (T) [sec]')
ylabel('Pseudo-Aceleración Espectral [g]')
legend('Conditional Mean Spectrum','RSN1227 CHICHI CHY074-N','Registro con Spectral-Matching a_{1,final}')
title('Comparación de Pseudo-Espectros de Respuesta')
grid on
xlim([0 3])

figure
plot(Tn,median_CMS,'color','r')
hold on
plot(Tn,PSa_regOriginal,'color','c')
plot(Tn,PSa_regEsc,'color','b')
hold off
grid on
xlabel('Periodos')
ylabel('Pseudo-Aceleración Espectral [g]')
legend('Conditional Mean Spectrum','PSa Registro sin esacalar','PSa Registro escalado (x1.3)')
xlim([0 3])