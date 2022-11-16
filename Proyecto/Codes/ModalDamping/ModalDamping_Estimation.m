%% Estimación de los amortiguamientos modales
% Optimización de TMD para mejora de desempeño sísmico a través de análisis
% del riesgo de colapso

%% Inicializar
clear variables
close all
clc

%% Inputs Parámetros Edificio                                               % Benchmark Ohtori et al 2004 (9 pisos)
k = [];                                                                     % Rigidez de cada piso (abajo hacia arriba)
m = [];                                                                     % Masa de cada piso (abajo hacia arriba)
xi_ = [2/100; 2/100];                                                                    % Amortiguamiento modal de dos modos conocidos
xi_modes = [1; 5];                                                              % Número del modo de los amortiguamientos modales (de arriba)

%% Inputs Parámetros TMD
C_tmd = [];                                                                 % Amortiguamientos del TMD


%% Estimar xi_n con TMDs

% ComputeK & M
K = computek(k);                                                            % Matriz de rigidez de la estructura
M = diag(m);                                                                % Matriz de masas de la estructura

% eigenvector
[Phi, lambda] = eig(K,M);                                                   % Phi: Forma modal, lambda = wn^2
wn = diag(sqrt(lambda));                                                          % Frecuencia circular de cada modo

% Amortiguamiento modal
    % Método de Rayleigh C = alfa*M + beta*K
    % xi_i = 1/2(alfa/omega_i + beta*omega_i)
syms alfa beta
R_eq1 = xi_(1) == 1/2*(alfa/wn(xi_modes(1)) + beta*wn(xi_modes(1)));
R_eq2 = xi_(2) == 1/2*(alfa/wn(xi_modes(2)) + beta*wn(xi_modes(2)));
sol = solve(R_eq1,R_eq2,[alfa beta]);

C = double(sol.alfa)*M + double(sol.beta)*K;                                % Matriz de amortiguamiento diagonalizable

% Modal Damping
Cn = Phi.'*C*Phi;
