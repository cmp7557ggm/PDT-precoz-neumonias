library(survival)
library(knitr)
library(dplyr)
library(gridExtra)
library(openintro)
library(plotrix)
library(formattable)
library(KMsurv)
library(MASS)
library(ggplot2)
library(survminer)
library(circlize)
library(readxl)
library(pillar)
library(kableExtra)
library(survival)
library(VSURF)
library(KMsurv)
library(randomForestSRC)
library(lubridate)
library(nlme)
library(fastDummies)
library(caret)
library(ggRandomForests)

neumonia<-read_excel("/Users/carlosmartinperez/Desktop/NEUMONIAS/BUENA/novedades/plantilla buena.xlsx",sheet=1)

emp<-as.duration(as.period(interval(neumonia$Fecha_h_neumo , neumonia$Fecha_h_empirico)))
emp<-ifelse(emp<0,0,emp)
emp<-round(emp/3600,1)
neumonia$emp<-emp

recepcion<-as.duration(as.period(interval(neumonia$Fecha_h_neumo , neumonia$Fecha_h_recepcion)))
recepcion<-ifelse(recepcion<0,0,recepcion)
recepcion<-round(recepcion/3600,1)
neumonia$recepcion<-recepcion
neumonia$recepcion

rapida_cruda<-as.duration(as.period(interval(neumonia$Fecha_h_neumo , neumonia$Fecha_h_rapido)))
rapida_cruda<-ifelse(rapida<0,0,rapida_cruda)
rapida_cruda<-round(rapida_cruda/3600,1)
neumonia$rapida_cruda<-rapida_cruda
neumonia$rapida_cruda

diagetiol_cruda<-as.duration(as.period(interval(neumonia$Fecha_h_neumo , neumonia$fecha_diag_etiol_verdad)))
diagetiol_cruda<-ifelse(diagetiol_cruda<0,0,diagetiol_cruda)
diagetiol_cruda<-round(diagetiol_cruda/3600,1)
neumonia$diagetiol_cruda<-diagetiol_cruda
neumonia$diagetiol_cruda

diagetiol2<-as.duration(as.period(interval(neumonia$Fecha_h_rapido , neumonia$fecha_diag_etiol_verdad)))
diagetiol2<-ifelse(diagetiol2<0,0,diagetiol2)
diagetiol2<-round(diagetiol2/3600,1)
neumonia$diagetiol2<-diagetiol2
neumonia$diagetiol2

final_cruda<-as.duration(as.period(interval(neumonia$Fecha_h_neumo , neumonia$Fecha_h_resultado_final)))
final_cruda<-ifelse(final_cruda<0,0,final_cruda)
final_cruda<-round(final_cruda/3600,1)
neumonia$final_cruda<-final_cruda
neumonia$final_cruda

alta_cruda<-as.duration(as.period(interval(neumonia$Fecha_h_neumo , neumonia$Fecha_h_alta)))
alta_cruda<-ifelse(alta_cruda<0,0,alta_cruda)
alta_cruda<-round(alta_cruda/3600,1)
neumonia$alta_cruda<-alta_cruda

provi_cruda<-as.duration(as.period(interval(neumonia$Fecha_h_neumo , neumonia$Fecha_h_provi)))
provi_cruda<-ifelse(provi_cruda<0,0,provi_cruda)
provi_cruda<-round(provi_cruda/3600,1)
neumonia$provi_cruda<-provi_cruda

cambiodat_cruda<-as.duration(as.period(interval(neumonia$Fecha_h_neumo , neumonia$Fecha_cambio)))
cambiodat_cruda<-ifelse(cambiodat_cruda<0,0,cambiodat_cruda)
cambiodat_cruda<-round(cambiodat_cruda/3600,1)
neumonia$cambiodat_cruda<-cambiodat_cruda
#ED0 es el 41
neumoRF<-dplyr::select(neumonia,c(4:8,10:27,34,36,38,39,41,42,59,60,64,65,67,68,71,73,75,76,78:80,
                                  86:89,91,92,95,96,100,107,109,111,113,114,115,116:124))
glimpse(neumoRF)
neumoRFsf<-dplyr::select(neumoRF,c(1:25,30,31,36,43,44,45,51:57,59,60,63,64))
glimpse(neumoRFsf)
neumoRFsf<-dplyr::select(neumoRFsf,-c(Microorganismo))
glimpse(neumoRFsf)

micro<-read_excel("/Users/carlosmartinperez/Desktop/NEUMONIAS/BUENA/novedades/plantilla buena.xlsx",sheet=2)
glimpse(micro)

neumoRFsf_dum<-cbind(neumoRFsf,micro)
glimpse(neumoRFsf_dum)
#BASE DE DATOS SIN FALTANTES

##################################################
#Localizar y eliminar variables con varianza cero o proxima a cero
##################################################

zv1<-neumoRFsf_dum %>%  nearZeroVar(saveMetrics = TRUE)

tablezv <- data.frame(zv1)
tablezv<-dplyr::select(tablezv,c(nzv))
tablezv$n<-seq(1:dim(tablezv)[1])
tablezv<-dplyr::filter(tablezv,nzv=="TRUE")
tablezv

datarf_nzv<-dplyr::select(neumoRFsf_dum,-c(11,15,17,29,42,43))
glimpse(datarf_nzv)

##################################################
#la base de datos pasa a llamarse df
##################################################

df<-datarf_nzv

##################################################
#AÑADIMOS VARIABLES PROPIAS
##################################################

df$ED1<-neumonia$ED1
df$provi<-neumonia$provi_cruda

#df<-dplyr::select(df,-c(Microorganismo))
glimpse(df)
##################################################
#Se quitan vay¡riables relacionadas con la intervención
##################################################
RF_sf<-dplyr::select(df,-c(Periodo,Pre_post,DER,PCR,DER_hecho,final_cruda))

#### Modelo RF #########################
set.seed(123)
mod1<-rfsrc(Surv(alta_cruda,Muertos )~.,data=RF_sf,ntree = 500, 
            importance = T)
mod1 

o <- tune(Surv(alta_cruda,Muertos)~.,data=RF_sf)
o$optimal

set.seed(123)
mod1<-rfsrc(Surv(alta_cruda,Muertos )~.,data=RF_sf,ntree = 500, nodesize = 3,mtry = 7,
            importance = T)
mod1

plot(gg_vimp(mod1))  

plot.variable(mod1, partial=T,plots.per.page = 1, xvar.names = c("ED1","Cancer","Hongo","Charlson_Index","NFR"))
plot.variable(mod1, partial=T,plots.per.page = 1, xvar.names = c("Edad","recepcion","Chronic_Kidney_disease","Provisional","rapida_cruda"))
plot.variable(mod1, partial=T,plots.per.page = 1, xvar.names = c("COVID","GRIPE"))


########################################
#Algoritmos para obtener todos los modelos posibles con su correspondiente error OOB
# y así seleccionar el modelo que menos error de predicción cometa

#_____funcion primera
data_select<-function(df,Yvar,Xvar=NULL){
  df<-df[,c(Yvar,Xvar)]
  Xvar2<-2:(length(Xvar)+1)
  return(df)
}

#_____funcion segunda (combinaciones tomandas de n en n)
combinar<-function(df,Xvar=NULL,numXVar=1){
  
  comb_Var<-t(combn(Xvar,numXVar));combinaciones<-comb_Var
  
  colnames(combinaciones)<-NULL
  
  options(warn=-1)
  
  get_names<-function(df,x){
    outcome<-matrix(data=NA,nrow=nrow(x),ncol=ncol(x))
    for(i in 1:nrow(x)){
      outcome[i,1:ncol(x)]<-colnames(df)[x[i,1:ncol(x)]]
    }
    return(outcome)
  }
  combn_names<-get_names(df,combinaciones)
  return(combn_names)
}

#_____funcion tercera (obtención de los modelos y su error)
Mconf<-function(df, x)
{
  filas<-nrow(ccc)
  materr2<-matrix(nrow=filas,ncol=2)
  colnames(materr2)<-c( "Variables", "Error")
  rownames(materr2)<-paste("Modelo",c(1:filas))
  for (i in 1:nrow(ccc)){
    xnam <- as.matrix(ccc[i,])
    fmla <- as.formula(paste("Surv(",colnames(b_datos[1]),",",colnames(b_datos[2]),") ~", paste(xnam, collapse= "+")))
    
    mod <- rfsrc(fmla, data = b_datos,ntree = 500)
    df <- mod$data
    prob = predict(mod,type = c("response"))
    df$prob = prob
    
    fmla<-as.character(fmla)
    #coeficiente<-summary(mod)$coefficients[2]
    #pv<-summary(mod)$coefficients[8]
    
    #materr2[i,1]<-nrow(xnam)
    materr2[i,1]<-fmla[3]
    materr2[i,2]<-get.cindex(mod$yvar[,1], mod$yvar[,2], mod$predicted.oob)
    #materr2[i,3]<-round(pv,4)
    
  }
  return(materr2)
}
b_datos<-dplyr::select(RF_sf,alta_cruda,Muertos,Cancer,Hongo,ED1,Charlson_Index,Edad,GRIPE,COVID)


#_________________________
ccc<-combinar(b_datos,c(3:9),1)
d1<-Mconf(b_datos,ccc)

ccc<-combinar(b_datos,c(3:9),2)
d2<-Mconf(b_datos,ccc)

ccc<-combinar(b_datos,c(3:9),3)
d3<-Mconf(b_datos,ccc)

ccc<-combinar(b_datos,c(3:9),4)
d4<-Mconf(b_datos,ccc)
d4
ccc<-combinar(b_datos,c(3:9),5)
d5<-Mconf(b_datos,ccc)

ccc<-combinar(b_datos,c(3:9),6)
d6<-Mconf(b_datos,ccc)

ccc<-combinar(b_datos,c(3:9),7)
d7<-Mconf(b_datos,ccc)

dfinal<-rbind(d1,d2,d3,d4,d5,d6,d7)
dfinal<-as.data.frame(dfinal)
dfinal<-arrange(dfinal,Error)
dfinalr<-slice(dfinal,1:10)
dfinalr%>%
  kbl(caption=" Todos los modelos posibles")%>%
  kable_styling(full_width = F) 

mod2<-rfsrc(Surv(alta_cruda,Muertos )~Cancer + Edad + ED1 +  GRIPE+COVID,data = RF_sf)
mod2
plot.variable(mod1, partial=T,plots.per.page = 1, xvar.names = c("ED1","Cancer","GRIPE","COVID","Edad"))


