# Introducción al repositorio
Dentro de este repositorio está disponible todo el código utilizado en la realización de mi Trabajo Fin de Máster para el Máster de Bioinformática y Biología Computacional de la Universidad Autónoma de Madrid.

Dicho esto, procedo a realizar la explicación del mismo, la cual es bastante sencilla. Encontramos cuatro carpetas llamadas "1. Pancreas", "2. Cerebro anterior", "3. Giro dentado del hipotálamo" y "4. Neutrofilos", y dentro de todas ellas encontraremos lo mismo: un fichero ".ipynb" con el código de Python necesario para llevar a cabo el análisis con Velocyto y scVelo; y un fichero ".R" en el que aparecerá el código utilizado para llevar a cabo el análisis usando los métodos de inferencia de trayectorias (*PAGA*, *PAGA_Tree*, *SCORPIUS* y *Slingshot*). 

**Nota**: *en el caso de la carpeta del giro dentado del hipotálamo, encontramos un fichero ".zip" porque el fichero ".ipynb" era demasiado grande si no. Además, en la carpeta de neutrofilos encontramos dos ficheros de ".R", uno con el análisis usando todos los genes y otro con el análisis usando solo los 2000 genes más importantes.*

Por otro lado, además de estas cuatro carpetas, también vemos que nada más entrar en el repositorio aparece un archivo con el nombre "velocyto_script.sh". Este es el *script* en el que encontramos el código que permite el uso de datos obtenidos con BD Rhapsody en Velocyto, y consta de dos partes: una primera de procesamiento, la cual es necesaria para ajustar los datos a los parametros establecidos por la herramienta, y una segunda en la que se ejecuta el comando de Velocyto propiamente dicho y que permite obtener las matrices de transcritos maduros e inmaduros.


Lista requisitos y requerimientos para poder llevar a cabo los análisis descritos al completo
* R packages
  + Seurat Disk--> https://github.com/mojaveazure/seurat-disk
  + Seurat --> https://satijalab.org/seurat/articles/install.html
  + dplyr --> https://www.r-project.org/nosvn/pandoc/dplyr.html
  + loomR --> https://github.com/mojaveazure/loomR
  + dyno --> https://dynverse.org/users/1-installation/
  + tidyverse --> https://www.tidyverse.org/packages/
  + Syngularity --> https://sylabs.io/docs/

* Python packages
  
    + velocyto --> http://velocyto.org/velocyto.py/install/index.html
    
	    * dependencies: numpy, scipy, cython, numba, matplotlib, scikit-learn, h5py, click, and pysam and loompy, which will be installed by "pip install velocyto" as dependencies.
      
    + scvelo --> https://scvelo.readthedocs.io/en/stable/installation/
    
	    * dependencies: anndata, scanpy, numpy, scipy, pandas, scikit-learn, igraph, matplotlib, louvain 

* Samtools --> http://www.htslib.org/
