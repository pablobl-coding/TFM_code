# Introducción al repositorio
Dentro de este repositorio, está disponible todo el código utilizado en la realización de mi Trabajo Fin de Máster para el Máster de Bioinformática y Biología Computacional de la Universidad Autónoma de Madrid.

Dicho esto, procedo a realizar la explicación del mismo, la cual es bastante sencilla. Encontramos cuatro carpetas llamadas "1. Pancreas", "2. Cerebro anterior", "3. Giro dentado del hipotálamo" y "4. Neutrofilos", y dentro de todas ellas encontraremos lo mismo: un fichero ".ipynb" con el código de Python necesario para llevar a cabo el análisis con Velocyto y scVelo; y un fichero ".R" en el que aparecerá el código utilizado para llevar a cabo el análisis usando los métodos de inferencia de trayectorias (*PAGA*, *PAGA_Tree*, *SCORPIUS* y *Slingshot*).

Por otro lado, además de estas cuatro carpetas, también vemos que nada más entrar en el repositorio aparece un archivo con el nombre "velocyto_script.sh". Este es el *script* en el que encontramos el código que permite el uso de datos obtenidos con BD Rhapsody en Velocyto, y consta de dos partes: una primera de procesamiento, la cual es necesaria para ajustar los datos a los parametros establecidos por la herramienta, y una segunda en la que se ejecuta el comando de Velocyto propiamente dicho y que permite obtener las matrices de transcritos maduros e inmaduros.
