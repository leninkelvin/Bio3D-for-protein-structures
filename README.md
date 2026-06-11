# Bio3D-for-multiple-structures
![A squirrel researcher looking stumped at a computer screen displaying a complex structural alignment.](GraphicalAbstract.png)

Este repositorio tiene varios scripts generales para ser usados en R/Rstudio para el análisis y clasificación de estructuras de proteínas.

El Abstract Gráfico desplegado en arriba corresponde al resumen que se sometió al 3er Congreso Nexus: de la bioinformática a la experimental.

El poster esta en [formato PDF](posterNEXUS3.pdf)

El material del poster se generó con el [script ProteinPCA-localsv2.0.R](ProteinPCA-localsv2.0.R). Este script usa archivos previamente descargados localmente y genera gráficas de PCA de 2 y 3 dimensiones interactivas. 

El [script ProteinPCA-specificv2.0.R](ProteinPCA-specificv2.0.R) realiza operaciones semejantes pero no genera gráficas de PCA interactivas. Generá un reporte de las características experimentales de las estructuras descargadas automaticamente. 

Contacto: lenin.dominguez.ramirez@proton.me

Notas:
La estructura 2oz7 tiene la mutación T877A (debería ser T878A)
Las estructuras 8fgy tiene las mutaciones L702H/H875Y/F877L/T878A y la 8fgz tiene las mutaciones L702H/H875Y/F877L/T878A.
Las estructuras 3v49 y 3v4a tienen peptido activador y in inhibidor.
