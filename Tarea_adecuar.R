setwd("C:/Users/52442/Desktop/")
# Con esto me muevo a mi escritorio antes de hacer cualquier cosa
# Ahi es donde se encuentra mi carpeta con todas las muestras sobre las que voy a trabajar

# Cargo las librerias, si hay de mas pero por si a caso 
library(sleuth)
library(rhdf5)
library(devtools)
library(edgeR)
library(ensembldb)

#Esta función permite mapear, a partir de la base de datos de
tx2gene <- function(){
  
  #     Dataset you want to use. To see the different datasets available within a biomaRt yo$
  #     host
  
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")  
  
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = "ensembl_transcript_id",
                       ens_gene = "ensembl_gene_id", ext_gene = "external_gene_name")
  return(t2g)
}

# Con esta funcion que acabamos de crear nos conectamos a la base de datos de Ensemble y obtenemos los id de los transcritos,
# el id de los genes, el nombre de los genes 

t2g <- tx2gene()

# Conectarse a la B.D de Ensemble
# Tarda por dos razones: esta saturado el servidor o algo anda mal
# Cambiamos celegans por hsapiens




base_dir <- "Practica_n/"
# Asignamos a un objeto la direccion en donde se encuentran nuestras muestras
# Mi carpeta se llama Practica_n y se encuentra en el escritorio




## Seleccionar las muestras
samples <- paste0("sample", c("1", "2", "3", "4"))
  # Asignamos a un objeto el nombre de las muestras asociadas a las carpetas que contienen las tablas de abundancias
  # De preferencia pocas porque es computacionalmente intensivo y tardaria mucho
  # Por lo que eligire solamente cuatro



kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
 # Con esto generamos la ruta de acceso a nuestras muestras





  # Selected samples
## Aqui decimos cuales son las condiciones experimentales
s2c <- data.frame(path = kal_dirs, sample = samples, muestras = c("sample5","sample5","sample6", "sample6"), stringsAsFactors=FALSE)
  # Aqui armamos una data frame que contenga el camino para acceder a las muestras, el nombre de las muestras y las condiciones experimentales
  # Las condiciones experimentales seran sample5 y 6
  # El script no corre si tienes solo sample5...6...7...8
  # Tienes que repetir mas de una condicion experimental para hacer el analisis diferencial de dos muestras en la misma condicion
  # De no hacerlo, al correr sleut_fit saldra un error en la terminal: valor ausente donde T/F es necesario


s2c # Esta es la base de datos generada


so <- sleuth_prep(s2c, ~ muestras, target_mapping = t2g,extra_bootstrap_summary = TRUE)
  # Con sleuth_prep operamos sobre la base de datos generada



so <- sleuth_fit(so)
  # Con sleuth_fit ajustamos modelos de error, y se estima la arianza de los remuestreos, la varianza biologica y el shrinkage...
  # No entiendo muy bien que es lo ultimo



so <- sleuth_wt(so, which_beta = "muestrassample6")
  # Con sleuth_wt se computa el test de Wald en el coeficiente beta de cada transcrito

sleuth_live(so)
  # Con el paquete shiny y corriendo sleuth_live podemos hacerlo interactivo
  # Aparecen varias pestañas:
    # Analisis: Aqui se observa algo semejante a un volcano plot
                # y algo acerca de el remuestreo que se realizo
    # Maps: graficas que muestran los componentes principales y 
                # una grafica de barras que nos dice el porcentaje de la variacion que explica 
                # cada componente
                # Y tambien un heatmap...
    # Summaries: Aqui podemos tener varias cosas
                # Distribucion de las abundancias
                # Distribucion de la longitud del fragmento... no se a que se refiera
                # Los datos que fueron procesados, semejante a lo que nos mostraste en clase
                # Y por ultimo la tabla de calisto, donde tenemos los id, est_counts.. etc
    # Diagnosticos: no se muy bien como describir esta seccion

# En conclusion es muy interactivo para visualizar los resultados de nuestro analisis


resultados <- read.table("Practica_n/test_table.csv",sep=",",
                       header=TRUE)
    # Test_table se supone que es una tabla
    # La tabla se obtiene de la seccion de analisis y se llama "test_table"
    # Se va a descargar un archivo csv que va a tener tooooodos los resultados hermosos de nuestro analisis
# De no haber leido que la tabla se llamaba como tal test_table no hubiera sabido que hacer despues

significativos <- which(resultados$qval < 0.1)  # Cual de los resultados tiene un valor menor a 0.1 es decir, son significativos
significativos <- resultados[significativos,]   # Agrega los resultados significativos en un renglon
upregulated <- which(significativos$b > 0)      # Todos los que tengan un valor mayor a 0 significa que estan sobreregulados 
upregulated <- significativos[upregulated,]     # Agrega los genes que son sobreregulados a un renglon
downregulated <- which(significativos$b < 0)    # Todos los que tengan un valor mayor a 0 significa que estan subregulados
downregulated <- significativos[downregulated,] # Agrega los genes que son subregulados a un renglon



write.table(upregulated,file = "Practica_n/Upregulated_S5vsS6.txt",sep="\t")
write.table(downregulated,file = "Practica_n/Downregulated_S5vsS6.txt",sep="\t")
# Aqui hacemos archivos txt que contienen los que estan sobreregulados y subregulados en las condiciones experimentales S5 y S6
# S5vsS6 
# Espero no equivocarme

# Y........ listo....
