function [fj,da,a1] = MultiplePeriods_SM_AlAtik_Abrahamson_2010(Reg,dt,Tast,xi,PSaObj_Tast,error_adm,beta_Newmark,ui,udi,gamma_relax)
% An Improved Method for Nonstationary Spectral Matching
% Linda Al Atik & Norman Abrahamson (2010)

% FUNCIÓN AJUSTADA PARA UNA SOLA RAZÓN DE AMORTIGUAMIENTO, NO PONER UN
% VECTOR, solo un valor escalar


% Inputs
% Reg               -> Registro de aceleraciones a modificar
% dt                -> Espaciamiento de tiempo de muetro
% Tast              -> Periodos al cual el espectro será modificado (Periodo(s) de ajuste T*) (de la Estructura 1GDL)
% xi                -> Razón de amortiguamiento del espectro a modificar (de la estructura 1GDL),
% PSaObj_Tast       -> Valor objetivo del espectro de pseudo-aceleración en T*
% error_adm    -> Error admisible para la diferencia entre PSa_Tast obtenido con onda modificada y el val_Obj_Tast
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
hi = zeros(t_length,N);
fj = zeros(t_length,N);
wi = zeros(N,1);
wi_prima = zeros(N,1);
freq = zeros(N,1);
gamma_f = zeros(N,1);
Dtj = zeros(N,1);
n_ti = zeros(N,1);
hi_titau = zeros(t_length,N);
cij = zeros(N,1);
dRi = zeros(N,1);
tj = zeros(N,1);
signPeak = zeros(N,1);
PSa_Tast = zeros(N,1);
b = zeros(N,1);
da = zeros(t_length,1);

for j = 1:N   % (j: spectral point)                                         % Recorremos periodo-amortiguamiento (Solo 1 en este caso)
    [~,~,~,~,PSa_aster,~,tPeak,signoPeak] = PSa_peak(beta_Newmark,xi,dt,ui,udi,Reg,Tast(j));   % Obtenemos solo el PSa en el periodo de interés
    
    tj(j) = tPeak;
    signPeak(j) = signoPeak;
    PSa_Tast(j) = PSa_aster;

    % Como N = 1, no se van a ocupar vectores, pero cada uno debería
    % depender del vector periodo-amortiguamiento (i), el vector de tiempo (t)
    % y del vector de periodos a ajustar (j)
    wi(j) = 2*pi/Tast(j);                                                         % Frecuencia circular de estructura de 1GDL                                                        
    betai = xi;                                                             % Razón de amortiguamiento de estructura 1GDL
    wi_prima(j) = wi(j)*sqrt(1-betai^2);                                          % (Eq.6) Frecuencia circular amortiguada de estructura de 1GDL
    freq(j) = wi_prima(j)/(2*pi);                                                 % (Frecuencia)
    gamma_f(j) = 1.178*freq(j)^-0.93;                                             % (Eq.17) Coeficiente que depende de la frecuencia (para Eq.16)
    Dtj(j) = atan(sqrt(1-betai^2)/betai)/wi_prima(j);                             % (Eq.14)
    
    % Funciones hi y fj
    for t = 1:t_length  % (t
        hi(t,j) = -wi(j)/(sqrt(1 - betai^2))*exp(-wi(j)*betai*(t-1)*dt)*sin(wi_prima(j)*(t-1)*dt);    % (Eq. 18)
        fj(t,j) = cos(wi_prima(j)*((t-1)*dt - tj(j) + Dtj(j)))*exp(-(((t-1)*dt - tj(j) + Dtj(j)/gamma_f(j))^2));   % (Eq. 16)
    end

    % Formar vector hi(ti-tau) para integrar en (Eq.7)
    n_ti(j) = round(tj(j)/dt,0) + 1;
    hi_titau(:,j) = [flip(hi(1:n_ti,j)); zeros(length(hi(:,j))-n_ti(j),1)];
end

for j = 1:N
    % cij
    cij(j) = trapz(dt,fj(j).*hi_titau(:,j));                                % (Eq.7) % cjj

    % dR
    % Eq.9 es útil solo si dRi es el error espectral de ajuste (spectral misfit)
    dRi(j) = (PSaObj_Tast(j)-PSa_Tast(j))*signPeak(j);                      % (Eq.1)

    % b
    b(j) = (cij(j)^-1)*dRi(j);                                              % (Eq.10)

    % da: Onda
    da(:,1) = da(:,1) + b(j).*fj(:,j);                                      % (Eq.2)
end

% a1: Registro de aceleraciones modificado
a1 = Reg + gamma_relax*da;                                                  % (Eq.11)

% Verificación error, si no cumple, iteramos otra vez
[~,~,~,~,PSamod,~,~,~] = PSa_peak(beta_Newmark,xi,dt,ui,udi,a1,Tast);       % Obtenemos PSa del registro modificado  

if abs(PSamod - PSaObj_Tast) > error_adm
    Reg = a1;   % reg = a(t) de (Eq.11)
    disp(1)
    [fj,a1,da] = MultiplePeriods_SM_AlAtik_Abrahamson_2010(Reg,dt,Tast,xi,PSaObj_Tast,error_adm,beta_Newmark,ui,udi,gamma_relax);
end

end