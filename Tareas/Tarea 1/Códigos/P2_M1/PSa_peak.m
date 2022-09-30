function [Sd,Sv,Sa,PSv,PSa,uPeak,tPeak,signPeak] = PSa_peak(beta,xi,dt,ui,udi,Reg,Tn)
% Con Método de Newmark, calcular espectros y pseudo-espectros de un
% registro de aceleraciones para un vector de periodos

% Contreras - Adasme (Ing.Sísimca 2022) Modificado
% Modificación entre líneas 83 y 86; y línea 30

% Inputs
% beta      -> factor del método de Newmark
% xi        -> Razón de amortiguamiento
% dt        -> Paso temporal del registro
% ui        -> Desplazamiento inicial
% udi       -> Velocidad Inicial
% Reg       -> Registro de aceleración del suelo [cm/s2]
% Tn        -> Vector columna de Periodos a analizar [s]

% Outputs
% Sd        -> Espectro de desplazamientos
% Sv        -> Espectro de velocidades
% Sa        -> Espectro de aceleraciones
% PSv       -> Pseudo-Espectro de velocidades
% PSa       -> Pseudo-Espectro de aceleraciones
% uPeak     -> Desplazamiento máximo
% tPeak     -> Tiempo al peak
% signPeak  -> Signo del peak

gamma = 0.5;                                                                % Método según cómo se toma la aceleración
Reg_length = length(Reg);
Tn_length = length(Tn);
Tn = Tn';  % Vector fila                                                    % Damos vuelta el vector (modificación)
t_vect = (0:dt:(length(Reg)-1)*dt)';                                        % Vector de tiempos del Registro

%% Inicializar vectores
% Para no reescribir sobre memoria cada vez que se modifica el tamaño

% Inicialización de vector de desplazamientos, velocidad y aceleración
u = zeros(Reg_length,1);
ud = zeros(Reg_length,1);
udd = zeros(Reg_length,1);

% Inicialización de espectro de desplazamiento, velocidad y aceleración
Sd = zeros(Tn_length,1);
Sv = zeros(Tn_length,1);
Sa = zeros(Tn_length,1);

% Inicialización de pseudo-espectros de velocidad y aceleración
PSv = zeros(Tn_length,1);
PSa = zeros(Tn_length,1);

% Condiciones Iniciales
u(1,1) = ui;
ud(1,1) = udi;

% Inicializar outputs
uPeak = zeros(Tn_length,1);
tPeak = zeros(Tn_length,1);
signPeak = zeros(Tn_length,1);
Pos_tPeak = zeros(Tn_length,1);

for j = 1:Tn_length                                                         % j es el periodo a evaluar
    u(:,1) = zeros(size(u,1),1);
    ud(:,1) = zeros(size(ud,1),1);
    udd(:,1) = zeros(size(udd,1),1);
    u(1,1) = ui;
    ud(1,1) = udi;
    if Tn(1,j) == 0
        for i = 1:Reg_length               % i es el número de dato del registro
            udd(i,1) = Reg(i,1);
        end
    else
        wn = 2*pi/Tn(1,j);                                                  % Frecuencia natural del sistema
        udd(1,1) = Reg(1,1) - 2*xi*wn*udi - wn^2*ui;
        a1 = 1/(beta*dt^2) + 2*xi*wn*gamma/(beta*dt);
        a2 = 1/(beta*dt) + 2*xi*wn*(gamma/beta-1);
        a3 = (1/(2*beta)-1) + 2*xi*wn*dt*(gamma/(2*beta)-1);
        k_ton = a1 + wn^2;
        for i = 2:Reg_length
            p_ton = -Reg(i,1) + a1*u(i-1,1) + a2*ud(i-1,1) + a3*udd(i-1,1);
            u(i,1) = p_ton/k_ton;
            ud(i,1) = gamma/(beta*dt)*(u(i,1)-u(i-1,1)) + (1-gamma/beta)*ud(i-1,1) + dt*(1-gamma/(2*beta))*udd(i-1,1);
            udd(i,1) = (u(i,1)-u(i-1,1))/(beta*dt^2) - ud(i-1,1)/(beta*dt) - (1/(2*beta)-1)*udd(i-1,1);
        end
        udd_tot = udd + Reg;
        udd_tot(1,1) = Reg(1,1);
        Sd(j,1) = max(abs(u(:,1)));                                         % cm   (esto también es uPeak)
        Sv(j,1) = max(abs(ud(:,1)));                                        % cm/s
        Sa(j,1) = max(abs(udd_tot(:,1)));                                   % cm/s2
        PSv(j,1) = wn*Sd(j,1);                                              % cm/s
        PSa(j,1) = wn^2*Sd(j,1);                                            % PSa (cm/s2)
    end
    % Inicio Modificación
    [uPeak(j,1),Pos_tPeak(j,1)] = max(abs(u));                               % Peak de desplazamiento (cm) y posición del tiempo al peak
    signPeak(j,1) = sign(u(Pos_tPeak(j,1)));                                 % Signo del peak
    tPeak(j,1) = t_vect(Pos_tPeak(j,1));                                     % Tiempo al peak
    % Fin modificación
end
end