TMD OPTIMIZATION COLLAPSE
CONTRERAS - SANGUINETTI
INGENIERIA SISMICA AVANZADA - USM 2022

En esta carpeta se genera el archivo .csv de los GroundMotions para ingresar a THAMDOF:
Todos los registros son escalados a las franjas que uno quiera, por lo que se tendrán la misma cantidad de registros para cada franja

En carpeta "Registros":
- Luego de seleccionar los registros, añadirlos como archivos NombreRegistro.txt como un vector con todas las aceleraciones en [g]
- Generar un archivo "GM Data.txt" (con ese nombre) con todos los nombres de los registros y los pasos temporales del sampling (muestreo, ej: "Maule2010.txt	0.01 [s]")


Ya agregados los registros y GM Data a la carpeta "Registros", solo ejecutar el script "txtGMtoCSVTHAMDOF.m" modificando los inputs a los correspondientes

GMFolder: Carpeta donde se encuentran los Registros
Cant_registros: Cantidad de registros a utilizar
GMDataName: Nombre del archivo "GM Data.txt", dejar como está
IM_SaT1: Cantidad de franjas para realizar el "IDA", mientras más mejor pero demora mucho más, por el momento, solo se tiene para poder seleccionar IMs como aceleraciones espectrales del primer modo
T: Periodo fundamental de la estructura, el de mayor masa traslacional
beta_newmark: beta del método de NewmarkBeta, dejar como está
xi: Fracción del amortiguamiento para realizar el espectro (cuidado al modificar, saber consecuencias)
t_extra: Tiempo extra luego de cada registro para que quede oscilando libremente la estructura
GMTHAMDOFName: Nombre del archivo .csv que se va a generar