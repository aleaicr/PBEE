function [IM_values,EDP_median] = nanEDP_median(EDP,IM,IM_interp1)
% Contreras - Sanguinetti
% Determinar la mediana de EDP cuando tiene valores NaN

% Inputs
% EDP:  Data(i).EDP: Matriz de valores de EDP devuelta por IIDAP
% IM:   Data(i).IM: Matriz de valores de IM devuelta por IIDAP
% IM_interp1: Vector con los valores de IM a interpolar (tarea3: 0.1:0.1:20)

% Outputs
% IM_values: Vector con IM de la media (los de IM_interp1 que sea capaz de interpolar
% para EDP)
% EDP_values: Vector con EDP de la media (los que alcance a interpolar)

%%
[~,n] = size(EDP);
EDP_vect = nan(size(EDP));
IM_vect = nan(size(EDP));
for i=1:n
    IM_rmm = rmmissing(IM(:,i));
    IM_rmm_length = length(IM_rmm);
    EDP_rmm = rmmissing(EDP(:,i));        
    
    EDP_vect(1:IM_rmm_length,i) = interp1(IM_rmm,EDP_rmm,IM_interp1(1:IM_rmm_length)','nearest');
    EDP_vals_len = length(rmmissing(EDP_vect(:,i)));
    IM_vect(1:EDP_vals_len,i) = IM_interp1(1:EDP_vals_len)';
 end
 
EDP_median = nanmedian(EDP_vect,2);
IM_values = IM_interp1(1:length(EDP_median))';



 


 
 
