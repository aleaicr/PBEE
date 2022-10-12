function [EDP,IM] = getIdaCurves2(ResultsDir, ResultsName)

% Version    : 1.0
% Creado por : Alexis Contreras R. 
% Fecha      : 12/10/2022 (m/d/y)
%-------------------------------------------------------------------
% DESCRIPCION
%-------------------------------------------------------------------
% [EDP,IM] = getIdaCurves2(ResultsDir, ResultsName)
%
% Obtiene las curvas de EDP vs IM para todos los registros a partir de un
% analisis IDA realizado en II-DAP v1.3 o superior.
%-------------------------------------------------------------------
% FUNCIONES ADICIONALES LLAMADAS
%-------------------------------------------------------------------
% < ninguna >
%-------------------------------------------------------------------
% INPUTS
%-------------------------------------------------------------------
% ResultsDir:               String con el nombre (dirección) de la carpeta
%
% ResultsName :             String con el nombre del archivo
%-------------------------------------------------------------------
% OUTPUTS
%-------------------------------------------------------------------
% EDP(nFranjas,nRegistros): Vector con EDP
% 
% IM(nFranjas,nRegistros):  Vector con el IM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cargar resultados
IDA = readmatrix(ResultsDir + "\" + ResultsName);

%% Calcular numero de franjas:
% nStripes = length(IDA); 

%% Cargar EDP y IM en matrices:
% Inicializar variables:
EDP = IDA(:,1);
IM = IDA(:,2);
end

