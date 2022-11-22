% No terminado

Para realizar la optimización basada en la disminución del riesgo de colapso, seguir las siguientes instrucciones:

- Diseñar estructura, obtener los perfiles o utilizar un benchmark de una estructura existente
- Obtener la curva de capacidad (carga cíclica según un patrón de carga, pero se puede obtener de forma simplificada (e imprecisa) desde pushover en SAP2000)
- Dejar todos los datos de la estructura como están en el formato de THAMDOF (.csv)
- Calcular el periodo del primer modo (periodo fundamental) de la estructura
- Realizar selección de registro con script "Sa_avg_GMS.m" (Sa_avg Ground Motion Selection)
- Con cada registro, dejarlo en formato .csv como se muestra en THAMDOF, de tener los registros en Reg.txt y el archivo GM Data.txt se puede correr el script "GMThamdofFormat.m"

- Analisis Estructura SIN TMD
- Con los registros, escalarlos en franjas de 0.1g (IM = Sa(T1))
- 
