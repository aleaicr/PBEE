function [Sa] = Newmark_Lineal_Sa(beta,xi,dt,ui,udi,uddg,Tn)
% Respuesta estructural por el método de Newmkar-beta

% Contreras - Adasme 
% Ingeniería Sísmica - USM 2022

% Inputs
% beta  -> factor del método de Newmark (ej beta = 0 o 1/4)
% xi    -> Razón de amortiguamiento de cada modo a analizar  (ej: 0.05, no 5)
% dt    -> [sec] Rango de partición temporal del registro
% ui    -> [m] Desplazamiento inicial (Condición inicial)
% udi   -> [m/s] Velocidad Inicial (Condición inicial)
% uddg  -> [m/sw] Registro de aceleración del suelo
% Tn    -> [sec] Vector de Periodos de cada modo a analizar

% Outputs
% Sd    -> [m] Espectro de desplazamientos
% Sv    -> [m/s] Espectro de velocidades
% Sa    -> [m/s2] Espectro de aceleraciones
% PSv   -> [m/s] Pseudo-Espectro de velocidades
% PSa   -> [m/s2] Pseudo-Espectro de aceleraciones

% Comentario
% Si beta = 0 -> verificar estabilidad del método omega*delta_t < 2
% si beta = 1/4 -> método es incondicionalmente establa

gamma = 0.5;                                                                % Método según cómo se toma la aceleración
uddg_length = length(uddg);
Tn_length = length(Tn);
Tn = Tn.';
%% Inicializar vectores
% Para no reescribir sobre memoria cada vez que se modifica el tamaño

% Inicialización de vector de desplazamientos, velocidad y aceleración
u = zeros(uddg_length,1);
ud = zeros(uddg_length,1);
udd = zeros(uddg_length,1);

% Inicialización de espectro de aeleración
Sa = zeros(Tn_length,1);

% Condiciones Iniciales
u(1,1) = ui;
ud(1,1) = udi;

for j = 1:Tn_length                                                         % j es el periodo a evaluar
    if Tn(1,j) == 0
        for i = 1:size(uddg,1)                                              % i es el número de dato del registro
            udd(i,1) = uddg(i,1);
        end
    else
        wn = 2*pi/Tn(1,j);                                                  % Frecuencia natural del sistema
        udd(1,1) = uddg(1,1) - 2*xi*wn*udi - wn^2*ui;
        a1 = 1/(beta*dt^2) + 2*xi*wn*gamma/(beta*dt);
        a2 = 1/(beta*dt) + 2*xi*wn*(gamma/beta-1);
        a3 = (1/(2*beta)-1) + 2*xi*wn*dt*(gamma/(2*beta)-1);
        k_ton = a1 + wn^2;
        for i = 2:size(uddg,1)
            p_ton = -uddg(i,1) + a1*u(i-1,1) + a2*ud(i-1,1) + a3*udd(i-1,1);
            u(i,1) = p_ton/k_ton;
            ud(i,1) = gamma/(beta*dt)*(u(i,1)-u(i-1,1)) + (1-gamma/beta)*ud(i-1,1) + dt*(1-gamma/(2*beta))*udd(i-1,1);
            udd(i,1) = (u(i,1)-u(i-1,1))/(beta*dt^2) - ud(i-1,1)/(beta*dt) - (1/(2*beta)-1)*udd(i-1,1);
        end
        udd_tot = udd + uddg;
        udd_tot(1,1) = uddg(1,1);
        Sa(j,1) = max(abs(udd_tot(:,1)));
    end
    u(:,1) = zeros(size(u,1),1);
    ud(:,1) = zeros(size(ud,1),1);
    udd(:,1) = zeros(size(udd,1),1);
    u(1,1) = ui;
    ud(1,1) = udi;
end
end