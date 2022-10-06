function [EDP,IM] = getIdaCurves(ResultsFile)

% Version    : 1.0
% Creado por : Cristian Cruz (CC) 
% Fecha      : 11/12/2020 (m/d/y)
% Rev log    : -
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
% nombreArchivo :           String con el nombre de archivo (y extensión) de los
%                           resultados entregados por IIDAP. Por ejemplo:
%                           'ResultadosIIDAP.mat'
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
%                           resultados de cada registro.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cargar resultados
load(ResultsFile)

% variable nGM indica el numero de registros
% variable IDA contiene resultados del analisis

%% Calcular numero de franjas:
nStripes = length(IDA.Sa1)-1; 

%% Cargar EDP y IM en matrices:
% Inicializar variables:
EDP = zeros(nStripes,nGM);
IM  = zeros(nStripes,nGM);

% Cargar datos:
for i = 1:nGM
    EDP(:,i) = IDA.(['U' num2str(i)])(1:nStripes);
    IM(:,i) =  IDA.(['Sa' num2str(i)])(1:nStripes);
end


