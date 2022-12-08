Optimización de TMD para la mejora del desempeño sísmico de edificios por medio de Evaluación del Riesgo de Colapso

- Lo que se planea realizar, es mejorar (disminuir) el lambda_c de una estructura utilizando un TMD, y encontrar el mejor TMD posible, es decir, encontrar aquel TMD que me disminuye lo máximo posible el lambda_c

Procedimiento:
1. Diseñar estructura
2. Obtener propiedades (rigidez, masa, amortiguamiento) de cada piso (para ingresar a THAMDOF)
3. Ajustar curvas de capacidad (monotónicas o cíclicas según un patrón de carga)
4. Seleccionar registros según un IM (el propuesto es Sa_avg ya que correlaciona bien con colapso)
5. Definir rango de IMs a utilizar (ej: 0.1:0.1:4 [g] para Sa_avg)
6. Dejar todo en formato THAMDOF (.csv)
    6.1 ConvertTxtToThamdof convierte los archivos desde formato .txt que otorga la carpeta GroundMotionSelection, a un archivo .csv con el formato que pide THAMDOF, con todos los IMs que uno quiera
7. Correr todos los registros escalados (para la estructura sin TMD, ojalá un IDA para ver qué rango de IMs son útiles)
8. Con respuestas de THAMDOF, realizar una evaluación del riesgo de colapso (obtener lambda_c)
9. 