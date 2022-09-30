function [fj,da,a1] = SM_AlAtik_Abrahamson_2010(Reg,dt,Tast,xi,PSaObj_Tast,error_adm,beta_Newmark,ui,udi,gamma_relax)
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
        hi(t) = -wi/(sqrt(1 - betai^2))*exp(-wi*betai*(t-1)*dt)*sin(wi_prima*(t-1)*dt);    % (Eq. 18)
        fj(t) = cos(wi_prima*((t-1)*dt - tj + Dtj))*exp(-(((t-1)*dt - tj + Dtj)/gamma_f)^2);   % (Eq. 16)
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
b = ((cij(end))^-1)*dRi;                                                           % (Eq.10)

% da: Onda
da = b.*fj;                                                                 % (Eq.2)

% a1: Registro de aceleraciones modificado
a1 = Reg + gamma_relax*da;                                                  % (Eq.11)

% Verificación error, si no cumple, iteramos otra vez
[~,~,~,~,PSamod,~,~,~] = PSa_peak(beta_Newmark,xi,dt,ui,udi,a1,Tast);       % Obtenemos PSa del registro modificado  

if abs(PSamod - PSaObj_Tast) > error_adm
    Reg = a1;   % reg = a(t) de (Eq.11)
    disp(1)
    [fj,a1,da] = SM_AlAtik_Abrahamson_2010(Reg,dt,Tast,xi,PSaObj_Tast,error_adm,beta_Newmark,ui,udi,gamma_relax);
end

end