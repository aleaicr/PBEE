function [EDP,IM,IMc,Backbone_Data] = getIdaCurves_v2_mod(ResultsDir, ResultsFile)
% Modificada por: Alexis Contreras en Tarea 3 para obtener struct backbone_data

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
% Backbone_Data:                 Vector columna con los datos para graficar la
%                           Backbone (aprovechando que se abren los datos
%                           dentro de esta función).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cargar resultados
A = load([ResultsDir '\' ResultsFile]);

% variable nGM indica el numero de registros
% variable IDA contiene resultados del analisis

%% Calcular numero de franjas máximo de todos los registros
nStripes = 0;
for i = 1:A.nGM
    nStripes = max(nStripes, length(A.IDA.(['Sa' num2str(i)]))-1);
end

%% Cargar EDP y IM en matrices:

g = 9.81;                   % Ac. gravedad

% Inicializar variables:
EDP = nan(nStripes,A.nGM);
IM  = nan(nStripes,A.nGM);
IMc = zeros(1,A.nGM);

% Cargar datos:
for i = 1:A.nGM
    nStripes_i = length(A.IDA.(['Sa' num2str(i)]))-1;
    EDP(1:nStripes_i,i) = A.IDA.(['U' num2str(i)])(1:nStripes_i);
    IM(1:nStripes_i,i) =  A.IDA.(['Sa' num2str(i)])(1:nStripes_i)./g;
    IMc(i) = A.IDA.(['Sa' num2str(i)])(end)./g;
end

% Gurdar Backbone
Backbone_Data = A.Backbone;


