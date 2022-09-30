%% Test PSa
% Ingeniería Sísmica Avanzada
% Contreras - Sanguinetti

% An Improved Method for Nonstationary Spectral Matching
% Linda Al Atik & Norman Abrahamson (2010)

% Para calcular el pseudo-espectro de aceleraciones se utiliza la función
% PSa_peak.m

%% Inicializar
clear variabales
close all
clc

%% Inputs
g = 980; % cm/s2

% Registro
load('RSN1227_CHICHI_CHY074-N.mat')
Reg = Acc*g;                                                                % Se ingresa registro en cm/s2
clear Acc Npts
% Reg   -> Registro de Aceleraciones
% dt    -> Paso temporal del registro de aceleraciones
t_vect = (0:dt:(length(Reg)-1)*dt)';
fprintf('Duración del registro sísmico es de t = %.2f sec\n\n',(t_vect(end)))

% Periodos
Tn_init = 0;
Tn_step = 0.1;
Tn_final = 10;
Tn = (Tn_init:Tn_step:Tn_final)';                                           % Vector columna

% beta (del Método de Newmark)
beta_Newmark = 1/4;

% Razón de amortiguamiento
xi = 0.05;

% Condiciones iniciales
ui = 0;
udi = 0;

%% Run function
[~,~,~,~,PSa,uPeak,tPeak,signPeak] = PSa_peak(beta_Newmark,xi,dt,ui,udi,Reg,Tn);

PSa = PSa./g;                                                                % Pasamos de cm/s2 a g
% uPeak [cm], tPeak [s]

%% Gráficamos
% Registro
figure('Name','Registro Sísmico')
plot((dt:dt:length(Reg)*dt)',Reg)
xlabel('Tiempo t (sec)')
ylabel('Aceleración [g]')
grid on
title('Registro de aceleraciones')
legend('RSN1227 CHICHI CHY074-N')

% PSa
figure('Name', 'PSa')
plot(Tn,PSa)
xlabel('Periodos T (sec)')
ylabel('Pseudo-espectro de aceleraciones PSa (g)')
grid on
title('Pseudo-Espectro de Aceleraciones')
legend('RSN1227 CHICHI CHY074-N')

% u_peak
figure('Name','u_peak')
plot(Tn',uPeak)
xlabel('Periodo T (sec)')
ylabel('Desplazamiento peak [cm]')
title('Desplazamiento peak')
grid on
legend('RSN1227 CHICHI CHY074-N')

% t_peak
figure('Name','t_peak')
plot(Tn,tPeak)
xlabel('Periodo T (sec)')
ylabel('Tiempo al peak (t@p) [sec]')
title('Tiempo al peak')
grid on
legend('RSN1227 CHICHI CHY074-N')

% sign_peak
figure('Name','sign_peak')
plot(Tn,signPeak)
xlabel('Periodo T (sec)')
ylabel('signo del peak del desplazamiento')
title('Signo del desplazamiento peak')
grid on
legend('RSN1227 CHICHI CHY074-N')
ylim([-2 2])

%% Tablas en ventana de comandos
tabla = table();
tabla.Tn = Tn;
tabla.PSa = PSa;
tabla.uPeak = uPeak;
tabla.tPeak = tPeak;
tabla.signPeak = signPeak;
disp(tabla)
clear tabla

