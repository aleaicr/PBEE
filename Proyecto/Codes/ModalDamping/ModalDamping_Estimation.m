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
k = [80;200;228;220;234;205;200;120;115];                                   % tonf/cm, Rigidez de cada piso (abajo hacia arriba)
m = [535;494.5;494.5;494.5;494.5;494.5;494.5;494.5;494.5]*1016;             % tonf,  Masa de cada piso (abajo hacia arriba)
xi_ = [2/100; 2/100];                                                       % Amortiguamiento modal de dos modos conocidos
xi_modes = [1; 5];                                                          % Número del modo de los amortiguamientos modales (de arriba)
nDOF = length(k);

% Con TMD
xi_TMD = [0.05; 0.1; 0.2];                                                  % Fracciones de amortiguamiento de cada TMD
k_TMD = [1;2;3];                                                            % Rigideces del TMD (para ajustarse al periodo del primer modo, aproximadamente) 
m_TMD_per = 5/100;                                                          % Todos los TMDs tendrán una masa del 5% de la amsa del edificio

%% Estimar xi_n SIN TMD

% ComputeK & M
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


%% Estimar xi_n CON TMD
% Se agrega un nuevo grado de libertad

% Para la optimización, se realiza la 
m_TMD = m_TMD_per*sum(m);                                                   % Masa del TMD

% Nueva matriz de masa -> incluyendo TMD como nuevo piso
m_new = [m; m_TMD];
M_new = diag(m_new);                                                        % Nueva matriz M considerando TMD

% valores de C_TMD para cada combinación de TMD
omega_TMD = zeros(length(k_TMD),1);
C_TMD_vals = zeros(length(k_TMD),length(xi_TMD));
TMD = struct();
structure = struct();

for i = 1:length(k_TMD)
    omega_TMD(i) = sqrt(k_TMD(i)/m_TMD);                                    % omga = sqrt(K/M)
    K_new = computeK([k; k_TMD(i)]);                                          % Nueva matriz K considerando Ktmd(i)
    [Phi_new,lambda_new] = eig(K_new,M_new);
    wn_new = sqrt(diag(lambda_new));
    Kn_new = Phi_new'*K_new*Phi_new;
    Mn_new = Phi_new'*M_new*Phi_new;
    for j = 1:length(xi_TMD)
        C_TMD_vals(i,j) = 2*m_TMD*xi(j)*omega_TMD(i);                       % C = 2 M xi Omega
        
        % Mezclar C y CTMD
        TMD(length(xi_TMD)*(i-1)+j).C_TMD = [C zeros(length(C),1); zeros(1,length(C)+1)];
        TMD(length(xi_TMD)*(i-1)+j).C_TMD(end,end) = C_TMD_vals(i,j);
        TMD(length(xi_TMD)*(i-1)+j).C_TMD(end-1,end-1) = TMD(length(xi_TMD)*(i-1)+j).C_TMD(end-1,end-1)+ C_TMD_vals(i,j);
        TMD(length(xi_TMD)*(i-1)+j).C_TMD(end-1,end) = - C_TMD_vals(i,j);
        TMD(length(xi_TMD)*(i-1)+j).C_TMD(end,end-1) = - C_TMD_vals(i,j);
        
        % Guardar a qué T (k) y xi está asociada este C_TMD
        TMD(length(xi_TMD)*(i-1)+j).k_TMD = k_TMD(i);
        TMD(length(xi_TMD)*(i-1)+j).T_TMD = sqrt(k_TMD(i)/m_TMD);
        TMD(length(xi_TMD)*(i-1)+j).xi_TMD = xi_TMD(j);
        TMD(length(xi_TMD)*(i-1)+j).C_TMD_val = C_TMD_vals(i,j);
        
        % Cn_new
        structure(length(xi_TMD)*(i-1)+j).Cn_new = Phi_new.'* TMD(length(xi_TMD)*(i-1)+j).C_TMD*Phi_new;
        
        % xi_n_new
        structure(length(xi_TMD)*(i-1)+j).xi_modal = diag(structure(length(xi_TMD)*(i-1)+j).Cn_new)./(2*Mn_new(length(xi_TMD)*(i-1)+j,length(xi_TMD)*(i-1)+j)*wn_new(length(xi_TMD)*(i-1)+j));
        structure(length(xi_TMD)*(i-1)+j).Kn = Kn_new;
        structure(length(xi_TMD)*(i-1)+j).Mn = Mn_new;
    end
end
%% Crear datos para analisis
xi_modal = zeros(length(structure(1).xi_modal),length(structure));  % fila: amortiguamiento del modo, columna: j
for i = 1:length(structure)
    xi_modal(:,i) = structure(i).xi_modal;
end
k_TMD_vect = zeros(1,length(TMD));
xi_TMD_vect = zeros(1,length(TMD));
for i = 1:length(TMD)
    k_TMD_vect(1,i) = TMD(i).k_TMD;
    xi_TMD_vect(1,i) = TMD(i).xi_TMD;
end

%% Guardar datos
matObj = matfile('ModalDamping_Data.mat');
matObj.Properties.Writable = true;
matObj.xi_modal = xi_modal;
matObj.k_TMD_vect = k_TMD_vect;
matObj.xi_TMD_vect = xi_TMD_vect;
clear matObj

fprintf('Revisar struct: "structure" para obtener las fracciones de amortiguamiento modales para ingresar a THAMDOF\n')
disp(structure)
fprintf('El de la estructura SIN TMD es:\n')
disp(xi)




