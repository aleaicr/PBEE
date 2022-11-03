function [EDP,IM,IMc] = getIdaCurves_v2(ResultsFile)

% Version    : 2.0
% Creado por : Cristian Cruz (CC) y Pablo Heresi (PH)
% Fecha      : 11/12/2020 (m/d/y)
% Rev log    : [PH] Editado para poder leer curvas EDP vs IM de registros
%              con diferentes largos (importante cuando hay registros que
%              producen colapso).
%              [PH] Editado para que IM sea entregado en unidades de [g].
%              [PH] Nuevo output: IMc
%-------------------------------------------------------------------
% DESCRIPCION
%-------------------------------------------------------------------
% [EDP,IM] = getIdaCurves(nombreArchivo)
%
% Obtiene las curvas de EDP vs IM para todos los registros a partir de un
% analisis IDA realizado en II-DAP v1.3 o superior.
%
%-------------------------------------------------------------------
% FUNCIONES ADICIONALES LLAMADAS
%-------------------------------------------------------------------
% < ninguna >
%
%-------------------------------------------------------------------
% INPUTS
%-------------------------------------------------------------------
% nombreArchivo : String con el nombre de archivo (y extensión) de los
%                 resultados entregados por IIDAP. Por ejemplo:
%                 'ResultadosIIDAP.mat'
% 
%-------------------------------------------------------------------
% OUTPUTS
%-------------------------------------------------------------------
% EDP(nFranjas,nRegistros): Matriz con los resplazamientos máximos del
%                           sistema de 1GDL. Las filas corresponden a los
%                           resultados de cada franja mientras que las
%                           columnas corresponden a los resultados de cada
%                           registro.
% 
% IM(nFranjas,nRegistros):  Matriz con el IM utilizado en cada analisis
%                           puede ser Sa(T1) o bien Sa_avg). Las filas
%                           corresponden a los resultados de cada franja
%                           mientras que las columnas corresponden a los
%                           resultados de cada registro. Unidades: [g].
%
% IMc(1,nRegistros):        Vector con el IM de colapso de cada registro.
%                           Unidades: [g].
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cargar resultados
load(ResultsFile)

% variable nGM indica el numero de registros
% variable IDA contiene resultados del analisis

%% Calcular numero de franjas máximo de todos los registros
nStripes = 0;
for i = 1:nGM
    nStripes = max(nStripes, length(IDA.(['Sa' num2str(i)]))-1);
end

%% Cargar EDP y IM en matrices:

g    = 9.81;                   % Ac. gravedad

% Inicializar variables:
EDP = nan(nStripes,nGM);
IM  = nan(nStripes,nGM);
IMc = zeros(1,nGM);

% Cargar datos:
for i = 1:nGM
    nStripes_i = length(IDA.(['Sa' num2str(i)]))-1;
    EDP(1:nStripes_i,i) = IDA.(['U' num2str(i)])(1:nStripes_i);
    IM(1:nStripes_i,i) =  IDA.(['Sa' num2str(i)])(1:nStripes_i)/g;
    IMc(i) = IDA.(['Sa' num2str(i)])(end)/g;
end


