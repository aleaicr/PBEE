%% Estimación de los amortiguamientos modales
% Optimización de TMD para mejora de desempeño sísmico a través de análisis
% del riesgo de colapso

% Comentarios:
% - Ingresar rigidez desde el pushover
% - Masa debe ser LA MASA SÍSMICA, no la masa de la estructura
%               Masa Sísmica es desde D + 0.25L (Ohtori et al 2004 ya
%               entrega la masa sísmica)
%% Inicializar
clear variables
close all
clc

%% Inputs Parámetros Edificio                                               % Benchmark Ohtori et al 2004 (9 pisos)
g = 980;                                                                    % cm/s2, aceleración de gravedad

% Propiedades edificio
k = [80;200;228;220;234;205;200;120;115];                                   % tonf/cm, Rigidez de cada piso (abajo hacia arriba)
w = [535;494.5;494.5;494.5;494.5;494.5;494.5;494.5;494.5];                  % tonf, Peso de cada piso (abajo hacia arriba)
xi_ = [2/100; 2/100];                                                       % Amortiguamiento modal de dos modos conocidos
xi_modes = [1; 5];                                                          % Número del modo de los amortiguamientos modales (de arriba)
% [T_heresi,G_heresi] = EigAnalysis(w,k,g);                                 Solo para verificar
% fprintf('Periodo Heresi = %.2f[sec]\n\n',T_heresi(1))

% Propiedades TMD
xi_TMD = [0.1; 0.15; 0.2];                                                % Fracciones de amortiguamiento de cada TMD
k_TMD = [1.2; 1.5; 1.7];                                                    % Rigideces del TMD (para ajustarse al periodo del primer modo, aproximadamente) 
m_TMD_per = 5/100;                                                          % Porcentaje de peso que (masaTMD/masaEddificio)

%% Estimar xi_n SIN TMD
% Cantidad de grados de libertad
nDOF = length(k);

% ComputeK & M
m = w/g;
K = computeK(k);                                                            % Matriz de rigidez de la estructura
M = diag(m);                                                                % Matriz de masas de la estructura

% eigenvector
[Phi, lambda] = eig(K,M);                                                   % Phi: Forma modal, lambda = wn^2
wn = diag(sqrt(lambda));                                                    % Frecuencia circular de cada modo
Kn = Phi'*K*Phi;
Mn = Phi'*M*Phi;

% Amortiguamiento proporcional
% Método de Rayleigh: C = alfa*M + beta*K
% xi_i = 1/2(alfa/omega_i + beta*omega_i)
syms alfa beta
R_eq1 = xi_(1) == 1/2*(alfa/wn(xi_modes(1)) + beta*wn(xi_modes(1)));
R_eq2 = xi_(2) == 1/2*(alfa/wn(xi_modes(2)) + beta*wn(xi_modes(2)));
sol = solve(R_eq1,R_eq2,[alfa beta]);

% Compute xi
xi = zeros(nDOF,1);
for i = 1:nDOF
    xi(i) = 1/2*(double(sol.alfa)/wn(i) + double(sol.beta)*wn(i));          % Fracción de amortiguamiento para cada modo (thamdof->piso)
end

% Compute Proportional Damping Matrix C
C = double(sol.alfa)*Mn + double(sol.beta)*Kn;                                % Matriz de amortiguamiento diagonalizable

% Modal Damping
Cn = Phi.'*C*Phi;

figure
plot(wn,xi*100,'-o','color','r')
xlabel('\omega_n [rad/s]')
ylabel('\zeta_n (%)')
grid on
title('Modal Damping')

structure1(1) = struct();
structure1(1).Kn = Kn;
structure1(1).Mn = Mn;
structure1(1).Cn = Cn;
structure1(1).wn = wn;
structure1(1).xi = xi;
structure1(1).T = 2*pi./wn;
structure1(1).T1 = 2*pi/wn(1);

%% Estimar xi_n CON TMD
% Se agrega un nuevo grado de libertad

% Para la optimización, se realiza la 
m_TMD = m_TMD_per*sum(m);                                                   % Masa del TMD

% Nueva matriz de masa -> incluyendo TMD como nuevo piso
m_new = [m; m_TMD];
M_new = diag(m_new);                                                        % Nueva matriz M considerando TMD

% valores de C_TMD para cada combinación de TMD
omega_TMD = zeros(length(k_TMD),1);                                         % Frecuencias del TMD
C_TMD_vals = zeros(length(k_TMD),length(xi_TMD));                           % Valores de C del TMD (para expandir matriz)
TMDs = struct();                                                            % struct() del TMD para guardar las propiedades de los distintos TMDs
structures = struct();                                                      % struct() de la estructura para guardar las propiedades de las etructuras con distintos TMDs

for i = 1:length(k_TMD)                                                     % Para cada una de las distintas rigideces
    omega_TMD(i) = sqrt(k_TMD(i)/m_TMD);                                    % omega_TMD = w1_TMD = sqrt(k_tmd,m_tmd) (frecuencia del TMD)
    K_new = computeK([k; k_TMD(i)]);                                        % Nueva matriz K considerando Ktmd(i)
    % Valor y vectores propios                                              
    [Phi_new,lambda_new] = eig(K_new,M_new);
    wn_new = sqrt(diag(lambda_new));
    Kn_new = Phi_new'*K_new*Phi_new;
    Mn_new = Phi_new'*M_new*Phi_new;
    for j = 1:length(xi_TMD)
        p = length(xi_TMD)*(i-1)+j;                                         % posición
        % Valor de C_TMD
        C_TMD_vals(i,j) = 2*m_TMD*xi_TMD(j)*omega_TMD(i);                   % C_TMD = 2 m_TMD xi_TMD w1_TMD
        
        % Mezclar C y CTMD
        structures(p).C_new = [C zeros(length(C),1); zeros(1,length(C)+1)];
        structures(p).C_new(end,end) = C_TMD_vals(i,j);
        structures(p).C_new(end-1,end-1) = structures(p).C_new(end-1,end-1)+ C_TMD_vals(i,j);
        structures(p).C_new(end-1,end) = -C_TMD_vals(i,j);
        structures(p).C_new(end,end-1) = -C_TMD_vals(i,j);
        
        % Guardar propiedades TMD
        TMDs(p).k_TMD = k_TMD(i);                                           % Guardar rigidez del TMD
        TMDs(p).M_tmd = m_TMD;                                              % Guardar la masa del TMD (igual para todos)
        TMDs(p).T_TMD = 2*pi/sqrt(k_TMD(i)/m_TMD);                          % Guardar periodo del TMD
        TMDs(p).wn_TMD = sqrt(k_TMD(i)/m_TMD);                              % Guardar frecuencia del TMD (wn = 2pi/T)
        TMDs(p).xi_TMD = xi_TMD(j);                                         % Guardar el amortiguamiento del TMD
        TMDs(p).C_TMD_val = C_TMD_vals(i,j);                                % Guardar el valor de C_TMD del TMD
        
        % Cn_new
        structures(p).Cn_new = Phi_new.'*structures(p).C_new*Phi_new; % 
        
        % xi_n_new
        structures(p).xi_new = diag(structures(p).Cn_new)./(2.*diag(Mn_new).*wn_new);

        % Guardar otros datos
        structures(p).Kn = Kn_new;
        structures(p).Mn = Mn_new;
        structures(p).K = K_new;
        structures(p).M = M_new;
    end
end
%% Crear datos para analisis
xi_new = zeros(length(structures(1).xi_new),length(structures));  % fila: amortiguamiento del modo, columna: j
for i = 1:length(structures)
    xi_new(:,i) = structures(i).xi_new;
end
k_TMD_vect = zeros(length(TMDs),1);
xi_TMD_vect = zeros(length(TMDs),1);
T_TMD_vect = zeros(length(TMDs),1);
for i = 1:length(TMDs)
    k_TMD_vect(i,1) = TMDs(i).k_TMD;
    xi_TMD_vect(i,1) = TMDs(i).xi_TMD;
    T_TMD_vect(i,1) = TMDs(i).T_TMD;
end
% 
% %% Guardar datos
% matObj = matfile('ModalDamping_Data.mat');
% matObj.Properties.Writable = true;
% matObj.xi_modal = xi_new;
% matObj.k_TMD_vect = k_TMD_vect;
% matObj.xi_TMD_vect = xi_TMD_vect;
% clear matObj

fprintf('Revisar struct() para ver las propiedades de la estructura original\n')
disp(structure1)
fprintf('El amortiguamiento de la estructura SIN TMD es:\n')
disp(xi)
fprintf('Revisar struct() para ver las propiedades de las estructuras de %.0f pisos (%.0f pisos + 1 TMD arriba)\n',nDOF+1,nDOF)
disp(structures)
fprintf('Revisar structures.xi_modal para ver las razones de amortiguamiento modal simplificadas por método explicado por P.Heresi\n\n\n')

fprintf('El periodo fundamental de la estructura SIN TMD: T1 = %.2f [sec]\n\n',structure1(1).T1)

fprintf('Los periodos de cada TMD son los siguientes: \n')
tabla = table();
tabla.combinacion = (1:1:length(k_TMD)*length(xi_TMD)).';
tabla.k_TMD_tonf_cm = k_TMD_vect;
tabla.xi_TMD = xi_TMD_vect;
tabla.Periodo_TMD_sec = T_TMD_vect; 
disp(tabla)
clear tabla
fprintf('Los amortiguamientos modales obtenidos para cada combinacion T_TMD,xi_TMD\n')

tabla = table();
% tabla.k_TMD = k_TMD_vect;
% tabla.xi_TMD = xi_TMD_vect;
tabla.xi_modal1 = structures(1).xi_new;
tabla.xi_modal2 = structures(2).xi_new;
tabla.xi_modal3 = structures(3).xi_new;
tabla.xi_modal4 = structures(4).xi_new;
tabla.xi_modal5 = structures(5).xi_new;
tabla.xi_modal6 = structures(6).xi_new;
tabla.xi_modal7 = structures(7).xi_new;
tabla.xi_modal8 = structures(8).xi_new;
tabla.xi_modal9 = structures(9).xi_new;
disp(tabla)
clear tabla
fprintf('Donde el número de xi_modal representa el número de la combinación\n')
fprintf('Peso del TMD w_TMD = %.2f [tonf]\n\n', m_TMD*g)
fprintf('Dejar P(columna P-Delta) como 0')
fprintf('')