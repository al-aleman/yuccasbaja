# Análisis de datos genómicos de *Yucca* por secuenciación NextRad
## Para **María Clara Arteaga** | Por **Alberto Alemán**

*El presente documento está escrito en formato **Markdown**, se recomienda visualizar en un navegador web o un programa específico para este tipo de archivos*

**Notas previas:**

```
Cuando veamos una caja como esta, estamos refiriéndonos a comandos o resultados dentro de una terminal
```

`Y cuando veamos texto de ésta manera, será para resaltar funciones o características de los comandos`

Utilizaremos `$` para diferenciar los comandos que introducimos en la terminal de los resultados de estos comandos. Los comandos se copian **SIN** `$`. Si una línea no comienza con `$`, probablemente es que ese sea un output (resultado) de la terminal.

* Nuestro mejor y más importante amigo, para comenzar, será el editor de texto. Por si sí o por si no, tengamos siempre abierto uno.
* Nuestro  segundo mejor amigo será la internet. **Saber pedir ayuda** (y cómo pedirla) es una virtud de sabios.
* El presente taller considera que todos los programas, funciones, librerías, paquetes (*etc*) estén previamente instalados en nuestra computadora. En algún momento podemos explorar la instalación de programas.

**LA TERMINAL**

![](https://www.islabit.com/wp-content/uploads/2019/06/macOS-5.png)

Al abrir un terminal, ¿en dónde estamos y a dónde vamos?
```
$ pwd
/home/precisiontower5810
$ cd Documents/
$ pwd
/home/precisiontower5810/Documents
$ cd ..
/home/precisiontower5810
$ cd Downloads/bioinf/
$ pwd
/home/precisiontower5810/Downloads/bioinf
$ cd ~
$ pwd
/home/precisiontower5810
```
`pwd` imprime el directorio de trabajo
`cd` cambia de directorio
`cd ~` siempre nos llevará "a casa"

**Hint: **Una opción *user friendly* es buscar la carpeta en dónde deseamos trabajar, dar click derecho y abrir una terminal desde ahí. Hagamos eso y vayamos al directorio de trabajo.

¿Qué hay dentro de una carpeta con datos crudos?

```
$ ls
1295_AAGAGGCA-AAGGAGTA_S592_L008_R1_001.fastq.gz
1295_AAGAGGCA-ACTGCATA_S580_L008_R1_001.fastq.gz
1295_AAGAGGCA-AGAGTAGA_S556_L008_R1_001.fastq.gz
1295_AAGAGGCA-CTAAGCCT_S604_L008_R1_001.fastq.gz
...
```

`ls` sirve para **enlistar** lo que hay dentro de una carpeta

**Nota técnica:**La muestras provenientes de **SNPsaurus** vienen *demultiplexeadas*, lo que significa cada archivo corresponde a cada muestra (aka *individuo*) y sus respectivas secuencias. Las secuencias dentro de cada archivo **no tienen índices**. Los índices sólo están en el nombre de cada archivo y sirven para identificar quién es quién.

¿Qué tiene cada archivo?

```
zcat 1295_AGGCAGAA-AGAGTAGA_S548_L008_R1_001.fastq.gz | head -4
@K00337:84:HKVK7BBXX:8:1101:7953:1244 1:N:0:NGGCAGAA+NGAGTAGA
GTAAGCGCGGTGGGGGGTGACAGAGGACGTGCTCGTACGGTTCATAGTGGGGTATGATATGATATTCTCGTGGATCTGTCTCTTATACACATCTCCGAGCCCACGAGACAGGCAGAAATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAA
+
AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJF-AFFJJJJJJJJFJJJJAJJAJFFJJJJJJJJJ
```

`zcat` se usa para ver el contenido de un archivo comprimido sin necesidad de descomprimirlo.
`|` se conoce como *pipe*, sirve para dar una instrucción a realizar al resultado de un comando previo
`head` es el comando para imprimir el encabezado con una N cantidad de líneas (en el ejemplo anterior, 4) de un determinado archivo

**Hint:** Algunos **comandos** tienen **argumentos**, mismos que se "activan" o cuyos parámetros se establecen con el símbolo `-`
**Hint:** `tail` es la contraparte de `head`
**Hint:** La **mayoría del tiempo** los comandos se escriben con minúsculas
**Hint:** El **tabulador** del teclado sirve para autocompletar texto, y así parecer unos escritores muy rápidos
**Hint:** Las flechas **arriba** y **abajo** nos permiten ver los comandos más recientes

¿Qué pasa si utilizamos `more` en vez de `head`?

```
zcat 1295_AGGCAGAA-AGAGTAGA_S548_L008_R1_001.fastq.gz | more
@K00337:84:HKVK7BBXX:8:1101:7953:1244 1:N:0:NGGCAGAA+NGAGTAGA
GTAAGCGCGGTGGGGGGTGACAGAGGACGTGCTCGTACGGTTCATAGTGGGGTATGATATGATATTCTCGTGGATCTGTCTCTTATACACATCTCCGAGCCCACGAGACAGGCAGAAATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAA
+
AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJF-AFFJJJJJJJJFJJJJAJJAJFFJJJJJJJJJ
--More--
```

`more` sirve para imprimir un determinado archivo, permitiéndonos explorarlo de principio a fin. Existe un comando llamado `less`, que en términos prácticos, funciona igual
**Nota:** "salimos" de *more* con la letra `q`.

Después de la breve introducción, vamos a analizar la calidad de nuestras secuencias.

## FastQC

Comenzamos verificando la calidad en conjunto con **FastQC** y **MultiQC**.

Lo siguiente se realiza dentro de la carpeta en donde se encuentran las secuencias crudas:

```
$ mkdir calidadcrudos
```

`mkdir` crea un directorio nuevo, con el nombre de la palabra que pongamos después de un espacio (en el ejemplo, *calidadcrudos*)
Si ponemos varias palabras después de un espacio, creamos otros directorios en el mismo sitio

```
$ fastqc -t 6 *.fastq.gz -o ./calidadcrudos/
```

El comando anterior llama al programa **fastqc**, el argumento `t -6` establece con cuantos threads (núcleos) se realizará la acción (ojo, que hay que saber cuántos núcleos tiene la computadora para que funcione bien), `*.fastq.gz` llama a **todos** los archivos con extensión *fastq.gz*, el argumento `-o ./calidadcrudos/` envía el output a la carpeta *calidadcrudos*

**Hint:** el símbolo de `*` es un comodín para remplazar cualquier texto delante o detrás de su ubicación.

¿Qué hay ahora en la carpeta *calidadcrudos*?

```
ls ./calidadcrudos/
```

El archivo resultante, podemos visualizarlo en **multiqc**

```
$ cd calidadcrudos/
$ mkdir multiqc
```

¿Qué hicimos con los comandos anteriores?

```
$ conda activate
$ source activate multiqc_py2.7
$ multiqc ./*.zip -o ./multiqc --data-format json --export
```

`conda activate` activó un ambiente llamado **conda** (los ambientes los podemos ver a profundidad después)
`source activate multiqc_py2.7` activó un paquete, con nombre **multiqc_py2.7**
`multiqc ./*.zip -o ./multiqc` llamó al programa multiqc a abrir **todos** los archivos con  extensión *.zip* (que fueron resultado del análisis con **fastqc**), y dio instrucción de generar el o los *outputs* resultantes en la carpeta *multiqc*, mientas que los argumentos `--data-format json --export` dieron la instrucción de **exportar** el resultado en formato *json* (JavaScript Object Notation)

Antes de continuar, desactivamos el paquete **multiqc_py2.7** y el ambiente **conda**

```
$ conda deactivate
$ conda deactivate
```

Ahora sí, observemos los resultados en nuestro navegador:

```
$ cd multiQC
$ chromium multiqc_report.html
```

Una vez observada la calidad, es necesario hacer cortes para dejar todo en una calidad ideal y libarnos de adaptadores.

## BBduk (herramienta de BBTools)

En la herramienta **bbduk** llevamos a cabo cortes, filtrado y remoción de adaptadores (SNPsaurus usan *Nextera*). Los argumentos de bbduk se basan en las instrucciones de SNPsaurus y en que todas las secuencias midan lo mismo (**100** pb), que es un requerimiento del programa **Stacks**.

**Ojo:** hay un parámetro cuyo valor va a depender de cómo veamos las secuencias (*qtrim*)

```
$ for i in $(ls *.fastq.gz | sed 's/.fastq.gz//')
$ do
$ bbduk.sh in=$i\.fastq.gz out=$i\_cleaned.fastq.gz ktrim=r k=17 hdist=1 mink=8 ref=/home/precisiontower5810/.programs/bbmap/resources/nextera.fa.gz minlen=100 ow=t qtrim=rl trimq=10 forcetrimright=99
$ done
```
Lo anterior se conoce como **loop**.
`for i in $` "para el elemento *i* de la lista *$*"
(`ls *.fastq.gz | sed 's/.fastq.gz//'`) enlista `ls` los archivos con terminación `*.fastq.gz` y después `sed 's/.fastq.gz//'` elimina la terminación *fastq.gz*
`do` "apliquemos el siguiente comando"
`bbduk.sh` llama al programa **bbduk**
En `in=$i\.fastq.gz out=$i\_cleaned.fastq.gz` el *input* es el elemento *i* de la lista, (es decir, la instrucción se hará a cada uno en el listado) y el *output* es el elemento *i* de la lista (conservarán su nombre original) más la terminación "_cleaned.fastq.gz"

**Hint:** `sed` es un manipulador de texto que funciona desde el terminal. Tiene muchísimas funciones que iremos viendo poco a poco.

Los siguientes son argumentos (funciones específicas) del programa **bbduk:**
`ktrim=r k=17 hdist=1 mink=8` son parámetros establecidos por *SNPsaurus*, por lo que sus valores dependen del tipo de secuenciación. *ktrim* remueve bases coincidentes con las referencias a partir de la coincidencia y hacia la derecha. *k* es el tamaño de la longitud usada para buscar contaminantes (o adaptadores), dado que el tamaño mínimo es 1, nuestras 16 pb se representan como 16+1. *hdist* y *mink* se encargan de buscar y remover subsecuencias.
`ref=/home/precisiontower5810/.programs/bbmap/resources/nextera.fa.gz` es la ruta en la que están los adaptadores *Nextera*.
En `minlen=100 ow=t qtrim=rl trimq=10 forcetrimright=99` *minlength* hace que las lecturas más cortas que este tamaño (posterior) a la limpieza se descarten, *ow* permite que los archivos puedan ser sobreescritos (tener cuidado con ésta opción), *trimq* refiere a que las regiones con una calidad promedio debajo de éste número serán removidas, y *forcetrimright* forza las secuencias a medir 100pb (necesario para **Stacks**)

Idealmente, movemos los archivos *limpios* a una nueva carpeta

```
$ mkdir cleaned
$ mv *cleaned.fastq.gz ./cleaned
$ cd cleaned/
```

Y repetimos el análisis con **fastqc** y su revisión con **multiqc**

```
$ mkdir calidadclean
$ fastqc -t 6 *.fastq.gz -o ./calidadclean/
$ cd calidadclean/
$ mkdir multiQC
$ source activate multiqc_py2.7
$ multiqc ./*zip -o ./multiQC --data-format json --export
$ conda deactivate
$ cd multiQC
$ chromium multiqc_report.html
```

¿Cómo se ve ahora el análisis de calidad?

Pasamos al ensamble *de novo*

## Stacks (versión 2.4)

En Stacks se hace el ensamblado del mapa *de novo*. Para ello, necesitamos a) saber la dirección absoluta de nuestras muestras y de nuestro output (es decir, crear un directorio para poner el resultado), b) hacer un listado de nuestras muestras en un formato de texto plano y copiarlo al directorio en donde se encuentra instalado **Stacks 2.4** y, llamar al programa desde el directorio en el que se encuentra.

Los parámetros de inteŕes se enlistan a continuación:

- T: specify the number of threads to execute.
- m: specify a minimum number of identical, raw reads required to create a stack. (profundidad.)
- M: specify the number of mismatches allowed between loci when processing a single individual (default 2). ("diferencias en un mismo individuo")
- n: specify the number of mismatches allowed between loci when building the catalog (default 1). ("diferencias en una población")

**Nota:** El formato del listado de nuestras muestras para Stacks debe ser muy específico, según el manual "(format is *name* **TAB** *pop*, one sample per line)." De tal manera que, en nuestro ejemplo, el "popmap" es un archivo .txt que se encuentra dentro de la carpeta de Stacks y tiene la siguiente información:

```
1295_AAGAGGCA-AAGGAGTA_S592_L008_R1_001_cleaned	1
1295_AAGAGGCA-ACTGCATA_S580_L008_R1_001_cleaned	1
1295_AAGAGGCA-AGAGTAGA_S556_L008_R1_001_cleaned	1
1295_AAGAGGCA-CTAAGCCT_S604_L008_R1_001_cleaned	1
1295_AAGAGGCA-CTCTCTAT_S532_L008_R1_001_cleaned	1
1295_AAGAGGCA-GTAAGGAG_S568_L008_R1_001_cleaned	1
1295_AAGAGGCA-TAGATCGC_S520_L008_R1_001_cleaned	1
1295_AAGAGGCA-TATCCTCT_S544_L008_R1_001_cleaned	1
1295_AGGCAGAA-AGAGTAGA_S548_L008_R1_001_cleaned	1
1295_AGGCAGAA-GTAAGGAG_S560_L008_R1_001_cleaned	1
1295_ATTCGTTG-AAGGAGTA_S627_L008_R1_001_cleaned	1
1295_ATTCGTTG-ACTGCATA_S624_L008_R1_001_cleaned	1
1295_ATTCGTTG-AGAGTAGA_S618_L008_R1_001_cleaned	1
1295_ATTCGTTG-CTAAGCCT_S630_L008_R1_001_cleaned	1
1295_ATTCGTTG-CTCTCTAT_S612_L008_R1_001_cleaned	1
1295_ATTCGTTG-GTAAGGAG_S621_L008_R1_001_cleaned	1
1295_ATTCGTTG-TAGATCGC_S608_L008_R1_001_cleaned	1
1295_ATTCGTTG-TATCCTCT_S615_L008_R1_001_cleaned	1
1295_CAGAGAGG-CTAAGCCT_S601_L008_R1_001_cleaned	1
1295_CCATATGT-TAGATCGC_S609_L008_R1_001_cleaned	1
1295_CGAGGCTG-AAGGAGTA_S591_L008_R1_001_cleaned	1
1295_CGAGGCTG-ACTGCATA_S579_L008_R1_001_cleaned	1
1295_CGAGGCTG-AGAGTAGA_S555_L008_R1_001_cleaned	1
1295_CGAGGCTG-CTAAGCCT_S603_L008_R1_001_cleaned	1
1295_CGAGGCTG-GTAAGGAG_S567_L008_R1_001_cleaned	1
1295_CGAGGCTG-TATCCTCT_S543_L008_R1_001_cleaned	1
1295_GCTACGCT-TAGATCGC_S518_L008_R1_001_cleaned	1
1295_GTAGAGGA-AAGGAGTA_S593_L008_R1_001_cleaned	1
1295_GTAGAGGA-ACTGCATA_S581_L008_R1_001_cleaned	1
1295_GTAGAGGA-AGAGTAGA_S557_L008_R1_001_cleaned	1
1295_GTAGAGGA-CTCTCTAT_S533_L008_R1_001_cleaned	1
1295_GTAGAGGA-GTAAGGAG_S569_L008_R1_001_cleaned	1
1295_GTAGAGGA-TAGATCGC_S521_L008_R1_001_cleaned	1
1295_GTAGAGGA-TATCCTCT_S545_L008_R1_001_cleaned	1
1295_GTGCACAA-AAGGAGTA_S626_L008_R1_001_cleaned	1
1295_GTGCACAA-CTAAGCCT_S629_L008_R1_001_cleaned	1
1295_GTTGTGGC-AGAGTAGA_S616_L008_R1_001_cleaned	1
1295_GTTGTGGC-CTCTCTAT_S610_L008_R1_001_cleaned	1
1295_GTTGTGGC-TAGATCGC_S606_L008_R1_001_cleaned	1
1295_GTTGTGGC-TATCCTCT_S613_L008_R1_001_cleaned	1
```

El siguiente es el comando de **Stacks** que hemos utilizado para nuestro ejercicio. Debemos tomar en cuenta que lo ejecutamos desde la carpeta en donde se localiza Stacks:
```
$ ./denovo_map.pl -T 7 --samples /media/precisiontower5810/096ff30c-b17d-432d-bac8-7bd3839a4eb3/Jaime/Yucca/NextRad/SNPsaurus/Muestras/Analisis/Postcalibration/hibrida/cleaned -o /media/precisiontower5810/096ff30c-b17d-432d-bac8-7bd3839a4eb3/Jaime/Yucca/NextRad/SNPsaurus/Muestras/Analisis/Postcalibration/hibrida/cleaned/stacks -m 6 -M 3 -n 4 --popmap hibrida.txt
```

Si nos colocamos en la carpeta de Stacks y escribimos `./denovo_map.pl -h`, imprimirá el manual del programa, y además de entender lo que acabamos de hacer, podremos ver qué otras características permite modificar por medio de `--` *flags*

**REPASO HASTA AHORA:**

- Los comandos pueden tener argumentos que se activan/desactivan/personalizan con el uso de *flags* `-` o `--`
- `pwd` **imprime** el **directorio** de trabajo
- `cd` **cambia** de **directorio**
- `cd ~` siempre nos llevará "a **casa**"
- Existen rutas absolutas (la referencia completa a una carpeta o archivo) y relativas (la referencia a partir del sitio en el que nos encontramos)
- `ls` sirve para **enlistar** lo que hay dentro de una carpeta, `ls -a` permite enlistar las carpetas ocultas
- `zcat` se usa para ver el **contenido** de un archivo **comprimido** sin necesidad de descomprimirlo
- `|` se conoce como ***pipe***, sirve para dar una instrucción a realizar al resultado de un comando previo
- `head` es el comando para imprimir una N cantidad de líneas del **encabezado** de un determinado archivo
- `more` sirve para **imprimir un determinado archivo**, permitiéndonos explorarlo de principio a fin. Existe un comando llamado `less`, que en términos prácticos, funciona igual
- `mkdir` crea un **directorio nuevo**, con el nombre de la palabra que pongamos después de un espacio
- Para **copiar** y **pegar** en la terminal, tenemos que utilizar `Ctrl+Shift+C` y `Ctrl+Shift+V`
- En la terminal, la combinación de teclas `Ctrl+C` sirve para **cancelar** el trabajo en curso. Es una herramienta útil para **detener** lo que se esté llevando a cabo (principalmente, si hay un error)
- El símbolo `>` **redirige** el output de un comando a un archivo con el nombre y la extensión de la palabra que pongamos después de un espacio, eg. `> results.txt`
- En la mayoría de los casos, el flag `- h` permite visualizar el manual del programa desde la terminal

## Populations (herramienta de Stacks) y Plink(1.07)

```
$ mkdir /media/precisiontower5810/096ff30c-b17d-432d-bac8-7bd3839a4eb3/Jaime/Yucca/NextRad/SNPsaurus/Muestras/Analisis/Postcalibration/hibrida/cleaned/stacks/pops
```
**Pregunta:**¿Qué hicimos con el comando anterior?
Un paso preliminar a populations es verificar el archivo denovo_map.log. Ahí encontraremos información útil que nos ayudará a decidir, por ejemplo, a qué individuos mantener de aquí en adelante. El archivo se encuentra en la carpeta de resultados y, podemos abrirlo con cualquier editor de texto.

**Pregunta:** ¿Cuántos individuos superaron el umbral de 15x?
**Pregunta:** ¿Hay alguna asociación entre x y el número de lecturas limpias por individuo?
**Por lo tanto,** ¿tenemos que hacer alguna modificación al popmap?

**Hint:** De aquí en adelante, diversos programas, además del *output* deseado, generarán archivos con extensión `.log`. Este tipo de archivos, tiene muchas utilidades, tanto para saber qué ocurrió con nuestros datos, como buscar errores o en este ejemplo, darnos información sobre el ensamblado general.

En la carpeta donde se localiza Stacks:

```
$ sudo gedit hibrida.txt
$ ./populations -h
$ ./populations -P /media/precisiontower5810/096ff30c-b17d-432d-bac8-7bd3839a4eb3/Jaime/Yucca/NextRad/SNPsaurus/Muestras/Analisis/Postcalibration/hibrida/cleaned/stacks -M <nuevo_popmap>.txt --min_maf 0.05 -r 0.80 -O /media/precisiontower5810/096ff30c-b17d-432d-bac8-7bd3839a4eb3/Jaime/Yucca/NextRad/SNPsaurus/Muestras/Analisis/Postcalibration/hibrida/cleaned/stacks/pops/ -t 6 --plink
```

##### no. de loci que salen de populations
`sudo` es el usuario con derechos de administrador, `gedit` es es nombre de nuestro editor de texto.
¿Cuál fue la función de los `--` aquí?

Utilizando **Plink**, vamos a hacer la evaluación de los sitios ensamblados en la carpeta donde se localizan los outputs de **Populations**
**Nota:** Plink es una herramienta muy completa, con muchas funciones y cuyos formatos de archivos son utilizados en otros programas. En éste paso lo utilizamos para hacer una evaluación rápida de los sitios ensamblados. Es bueno tomar en cuenta leer **el manual físico** que tenemos en el laboratorio. Nos ayudará a entender mucho sobre nuestros datos.
Podemos (recomendado en posteriores análisis) cambiar los nombres de `<outfiles>` por alguno que nos agrade.

```
$ plink1 --ped populations.plink.ped --map populations.plink.map --make-bed --out outfile1 --noweb
28362 SNPs
$ plink1 --bfile outfile1 --hwe --make-bed --out outfile2 --noweb
 SNPs
$ plink1 --bfile outfile_2 --missing --out missingfile --noweb
$ plink1 --bfile outfile_2 --mind 0.2 --make-bed --out outfile_3 --noweb
 samples,  SNPs
$ plink1 --bfile outfile_3
 samples,  SNPs, r
# LISTA NEGRA
$ plink1 --bfile outfile_1 --write-snplist --out rawdata --noweb
$ plink1 --bfile outfile_3 --write-snplist --out cleandata
$ sort rawdata.snplist cleandata.snplist | uniq -u > preblacklist
$ sed 's/\_.*//' preblacklist | uniq > blacklist
```

Podemos revisar la lista negra para saber cuántos sitios quedaron fuera
`--ped populations.plink.ped` y `--map populations.plink.map` son flags para abrir los outputs de Stacks
**Pregunta:**  ¿Qué hacen los flags que utilizamos en Plink?

Regresamos a **Populations**, para hacer nuestra base de datos de SNP's

```
$ ./populations -P /media/precisiontower5810/096ff30c-b17d-432d-bac8-7bd3839a4eb3/Jaime/Yucca/NextRad/SNPsaurus/Muestras/Analisis/Postcalibration/hibrida/cleaned/stacks/pops -M <popmap> --min_maf 0.05 -r 0.80 -O /media/precisiontower5810/096ff30c-b17d-432d-bac8-7bd3839a4eb3/Jaime/Yucca/NextRad/SNPsaurus/Muestras/Analisis/Postcalibration/hibrida/cleaned/stacks/populations -B /media/precisiontower5810/096ff30c-b17d-432d-bac8-7bd3839a4eb3/Jaime/Yucca/NextRad/SNPsaurus/Muestras/Analisis/Postcalibration/hibrida/cleaned/stacks/pops/blacklist --write_random_snp -t 6 --vcf
```

El archivo .vcf tiene los nombres asignados para las muestras con base en los adaptadores, y en éste punto ya necesitamos saber quién es quién. Para ello, haremos uso de `sed`, nuestro manipulador de texto desde la terminal

```
$ sed -i 's/1295_AAGAGGCA-AAGGAGTA_S592_L008_R1_001_cleaned/H78_9/g' hibrida.vcf
sed -i 's/1295_AAGAGGCA-ACTGCATA_S580_L008_R1_001_cleaned/H77_17/g' hibrida.vcf
sed -i 's/1295_AAGAGGCA-AGAGTAGA_S556_L008_R1_001_cleaned/H76_9/g' hibrida.vcf
sed -i 's/1295_AAGAGGCA-CTAAGCCT_S604_L008_R1_001_cleaned/H79_12/g' hibrida.vcf
sed -i 's/1295_AAGAGGCA-CTCTCTAT_S532_L008_R1_001_cleaned/H76_2/g' hibrida.vcf
sed -i 's/1295_AAGAGGCA-GTAAGGAG_S568_L008_R1_001_cleaned/H77_13/g' hibrida.vcf
sed -i 's/1295_AAGAGGCA-TAGATCGC_S520_L008_R1_001_cleaned/H74_4/g' hibrida.vcf
sed -i 's/1295_AAGAGGCA-TATCCTCT_S544_L008_R1_001_cleaned/H76_6/g' hibrida.vcf
sed -i 's/1295_AGGCAGAA-AGAGTAGA_S548_L008_R1_001_cleaned/H11_6/g' hibrida.vcf
sed -i 's/1295_AGGCAGAA-GTAAGGAG_S560_L008_R1_001_cleaned/H11_12/g' hibrida.vcf
sed -i 's/1295_ATTCGTTG-AAGGAGTA_S627_L008_R1_001_cleaned/H94_7/g' hibrida.vcf
sed -i 's/1295_ATTCGTTG-ACTGCATA_S624_L008_R1_001_cleaned/H94_6/g' hibrida.vcf
sed -i 's/1295_ATTCGTTG-AGAGTAGA_S618_L008_R1_001_cleaned/H93_2/g' hibrida.vcf
sed -i 's/1295_ATTCGTTG-CTAAGCCT_S630_L008_R1_001_cleaned/H94_8/g' hibrida.vcf
sed -i 's/1295_ATTCGTTG-CTCTCTAT_S612_L008_R1_001_cleaned/H92_10/g' hibrida.vcf
sed -i 's/1295_ATTCGTTG-GTAAGGAG_S621_L008_R1_001_cleaned/H94_4/g' hibrida.vcf
sed -i 's/1295_ATTCGTTG-TAGATCGC_S608_L008_R1_001_cleaned/H92_8/g' hibrida.vcf
sed -i 's/1295_ATTCGTTG-TATCCTCT_S615_L008_R1_001_cleaned/H93_1/g' hibrida.vcf
sed -i 's/1295_CAGAGAGG-CTAAGCCT_S601_L008_R1_001_cleaned/H36_17/g' hibrida.vcf
sed -i 's/1295_CCATATGT-TAGATCGC_S609_L008_R1_001_cleaned/H94_9/g' hibrida.vcf
sed -i 's/1295_CGAGGCTG-AAGGAGTA_S591_L008_R1_001_cleaned/H72_2/g' hibrida.vcf
sed -i 's/1295_CGAGGCTG-ACTGCATA_S579_L008_R1_001_cleaned/H70_4/g' hibrida.vcf
sed -i 's/1295_CGAGGCTG-AGAGTAGA_S555_L008_R1_001_cleaned/H68_6/g' hibrida.vcf
sed -i 's/1295_CGAGGCTG-CTAAGCCT_S603_L008_R1_001_cleaned/H73_7/g' hibrida.vcf
sed -i 's/1295_CGAGGCTG-GTAAGGAG_S567_L008_R1_001_cleaned/H68_7/g' hibrida.vcf
sed -i 's/1295_CGAGGCTG-TATCCTCT_S543_L008_R1_001_cleaned/H67_1/g' hibrida.vcf
sed -i 's/1295_GCTACGCT-TAGATCGC_S518_L008_R1_001_cleaned/H36_20/g' hibrida.vcf
sed -i 's/1295_GTAGAGGA-AAGGAGTA_S593_L008_R1_001_cleaned/H88_8/g' hibrida.vcf
sed -i 's/1295_GTAGAGGA-ACTGCATA_S581_L008_R1_001_cleaned/H88_5/g' hibrida.vcf
sed -i 's/1295_GTAGAGGA-AGAGTAGA_S557_L008_R1_001_cleaned/H87_7/g' hibrida.vcf
sed -i 's/1295_GTAGAGGA-CTCTCTAT_S533_L008_R1_001_cleaned/H87_3/g' hibrida.vcf
sed -i 's/1295_GTAGAGGA-GTAAGGAG_S569_L008_R1_001_cleaned/H88_2/g' hibrida.vcf
sed -i 's/1295_GTAGAGGA-TAGATCGC_S521_L008_R1_001_cleaned/H92_3/g' hibrida.vcf
sed -i 's/1295_GTAGAGGA-TATCCTCT_S545_L008_R1_001_cleaned/H87_4/g' hibrida.vcf
sed -i 's/1295_GTGCACAA-AAGGAGTA_S626_L008_R1_001_cleaned/H92_6/g' hibrida.vcf
sed -i 's/1295_GTGCACAA-CTAAGCCT_S629_L008_R1_001_cleaned/H92_7/g' hibrida.vcf
sed -i 's/1295_GTTGTGGC-AGAGTAGA_S616_L008_R1_001_cleaned/H89_6/g' hibrida.vcf
sed -i 's/1295_GTTGTGGC-CTCTCTAT_S610_L008_R1_001_cleaned/H89_5/g' hibrida.vcf
sed -i 's/1295_GTTGTGGC-TAGATCGC_S606_L008_R1_001_cleaned/H88_3/g' hibrida.vcf
sed -i 's/1295_GTTGTGGC-TATCCTCT_S613_L008_R1_001_cleaned/H89_8/g' hibrida.vcf
```
## PGDSpider (conversión de formatos)

**PGDSpider** se volverá un gran aliado. Si bien es cierto que existen muchos conversores de formatos, éste es el más *user friendly*, de fácil instalación y además, tiene entorno gráfico (hecho que se agradece mucho).
La siguiente lista la necesitaremos en **PGDSpider**, pues necesitamos incluir un archivo con la localidad a la que pertenece cada individuo.
**Pregunta:** ¿Cómo obtuvimos ambas listas?

```
H36_17 H36
H36_20 H36
H67_1 H67
H68_6 H68
H68_7 H68
H70_4 H70
H72_2 H72
H73_7 H73
H74_4 H74
H76_2 H76
H76_6 H76
H76_9 H76
H77_13 H77
H77_17 H77
H78_9 H78
H79_12 H79
H87_3 H87
H87_4 H87
H87_7 H87
H88_2 H88
H88_3 H88
H88_5 H88
H88_8 H88
H89_5 H89
H89_6 H89
H89_8 H89
H92_10 H92
H92_3 H92
H92_6 H92
H92_7 H92
H92_8 H92
H93_1 H93
H93_2 H93
H94_4 H94
H94_6 H94
H94_7 H94
H94_8 H94
H94_9 H94
```

**Generemos** un archivo con formato **genepop** a partir de la conversión de nuestro archivo **.vcf** en **PGDSpider**.

* Se recomienda poner el VCF en una carpeta accesible, para facilitar las conversiones

Abriremos la terminal en la carpeta donde se encuentra **PGDSpider** y escribiremos:

```
$ ./PGDSpider2.sh
```

Elegiremos un destino y un nombre para nuestro archivo **genepop**, y por ahora solamente entraremos a R para ver que todo salió bien y es leído adecuadamente.

Llegar a este punto es el **PRIMER OBJETIVO**. Tuvimos la fortuna de iniciar con parámetros de ensamble preestablecidos funcionales y que nuestro set de individuos es corto. Para completar estos pasos a partir de secuenciaciones **nuevas**, es necesario hacer **optimizaciones** de Stacks con sets de una muestra (pequeña pero representativa) de todas las muestras.

A partir de éste momento tenemos el **primer set** de datos. Una vez realizadas todas las optimizaciones necesarias y que hayamos tomado un criterio final para el ensamble de Stacks, **se recomienda no volver a Stacks o Populations** (a menos que sea para cambiar parámetros de ensamblado o si vamos a hacer algo con los sitios monomórficos). De aquí en adelante, todo lo trabajaremos a partir de nuestro **VCF**.

El archivo **VCF** tiene todo lo que necesitamos para hacer los análisis de **genética de poblaciones** *per se*. De aquí en adelante solamente haremos conversiones de formato y le extraeremos la información que vayamos requiriendo.
