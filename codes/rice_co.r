########################################################################
library(readxl)
library(XLConnect)
library(ggplot2)
library(dplyr)
library(scales)
library(cartography)

# library(ggpubr)
########################################################################

setwd("/home/alf/Scrivania/aa_lavori/lav_toscano_matrix")

########################################################################
# function

smooth_ecd_gray = function(dat, titleg,subtitleg="RICE",adj=0.1,e=0,labelx="Concentration_op") {
                      if ( labelx =="Concentration_op") { dat$x=dat$Concentration_op};
                      if ( labelx =="meanTot_op"){ dat$x=dat$meanTot_op};
                      dens = density(c(0,1), adjust=0.1, from=min(0,na.rm=T), to=max(1,na.rm=T),na.rm =T)
                       if (length(na.omit(dat$x)) > 1 ) {dens = density(dat$x, adjust=adj, from=min(dat$x,na.rm=T)-e, to=max(dat$x,na.rm=T) +e,na.rm =T)}
                      dens = data.frame(x=dens$x, y=dens$y)
                      ggplot(dat, aes(x)) +
                      stat_density(geom="line",colour="blue", adjust=adj, alpha=0.5) +
                      stat_smooth(data=dens, aes(x=x, y=cumsum(y)/sum(y)), method = "loess",colour="black",size=0.3, alpha=0.3) +
                      stat_ecdf(colour="red", size=0.6, alpha=0.6) +
                      scale_y_continuous(limits=c(0,1),oob = rescale_none)+
                      theme_gray()+
                      xlab(labelx)+
                      ylab("Density")+
                      ggtitle(titleg,subtitle = subtitleg)
}



info_extract=function(mat,param,type,plant) {
  mat$Concentration_op=mat$Concentration
  mat$ToTresValUncertSD_op=mat$ToTresValUncertSD
  mat$meanTot_op=mat$meanTot
  mat$POSresValUncertSD_op=mat$POSresValUncertSD
  mat$meanPos_op=mat$meanPos
  mat$min_op=mat$min
  mat$max_op=mat$max
  mat$LOD_op=mat$LOD
  names_mat=names(mat)
  mat=mat[which(mat$paramType==param),]
  mat=mat[which(mat$sampMatType==type),]
  if ( !is.na(plant)) {mat=mat[which(mat$sampMatbased==plant),]}
  mat=as.data.frame(mat)
  nrowmat=nrow(mat)
  if ( nrowmat==0) { mat[1,]=rep(0,length(names_mat))}
  vars_new=c("Concentration_op","ToTresValUncertSD_op","meanTot_op" ,"POSresValUncertSD_op","meanPos_op","min_op","max_op","LOD_op")
  mat_op=mat[vars_new]
  mat_op[mat_op==-1]=NA
  mat_op[mat_op==-2]=NA
  res=list()
  res$mat_op=mat_op
  res$mat=mat
  res$record=nrowmat
  res$N_conc=length(which(!is.na(mat_op$Concentration_op)))
  res$N_mtot=length(which(!is.na(mat_op$meanTot_op)))
  if ( nrowmat==0) { res$N_conc=0;res$N_mtot=0}
  res$sumsmapsize=sum(mat$sampSize,na.rm=T)
  res$tab=apply(mat_op,2,mean,na.rm=T)
  res$ls_tab=data.frame(data.frame(t(res$tab)),sumsmapsize=res$sumsmapsize,record=res$record,N_conc=res$N_conc,N_mtot=res$N_mtot)
  return(res)
}

info_extract_co_occur=function(mat,ref,type) {
  mat$Concentration_op=mat$Concentration
  mat$ToTresValUncertSD_op=mat$ToTresValUncertSD
  mat$meanTot_op=mat$meanTot
  mat$POSresValUncertSD_op=mat$POSresValUncertSD
  mat$meanPos_op=mat$meanPos
  mat$min_op=mat$min
  mat$max_op=mat$max
  mat$LOD_op=mat$LOD
  names_mat=names(mat)
  mat=mat[which(mat$sampMatType==type),]
  mat=mat[which(mat$Ref==ref),]
  mat=as.data.frame(mat)
  nrowmat=nrow(mat)
  if ( nrowmat==0) { mat[1,]=rep(0,length(names_mat))}
  
  vars_new=c("Concentration_op","ToTresValUncertSD_op","meanTot_op" ,"POSresValUncertSD_op","meanPos_op","min_op","max_op","LOD_op")
  mat_op=mat[vars_new]
  mat_op[mat_op==-1]=NA
  mat_op[mat_op==-2]=NA
  res=list()
  res$mat_op=mat_op
  res$mat=mat
  res$record=nrowmat
  res$N_conc=length(which(!is.na(mat_op$Concentration_op)))
  res$N_mtot=length(which(!is.na(mat_op$meanTot_op)))
  if ( nrowmat==0) { res$N_conc=0;res$N_mtot=0}
  res$sumsmapsize=sum(mat$sampSize,na.rm=T)
  res$tab=apply(mat_op,2,mean,na.rm=T)
  res$ls_tab=data.frame(data.frame(t(res$tab)),sumsmapsize=res$sumsmapsize,record=res$record,N_conc=res$N_conc,N_mtot=res$N_mtot)
  return(res)
}
############################################################################################################################################################

# RICE=read_xlsx("RICE.xlsx")
# RICE[RICE==-999]=NA
# saveRDS(RICE,"RICE.rds")



#########################################################################################################

RICE=readRDS("RICE.rds")
query="rice"
RICE_co=RICE[which(RICE$sampMatbased==query),]
RICE_co=as.data.frame(RICE_co[which(RICE_co$Co_occurrence==1),])

RICE_co$paramType=tolower(RICE_co$paramType)
RICE_co$sampMatType=tolower(RICE_co$sampMatType)
RICE_co$sampMatbased=tolower(RICE_co$sampMatbased)
RICE_co$sampMatCode=tolower(RICE_co$sampMatCode)




U_ref_co=unique(RICE_co$Ref)
U_type_co=unique(RICE_co$sampMatType)

#############################################################################################
res_co_data=list()
res_co_tab=list()

z=1
for ( i in seq_along(U_ref_co)) {
  for ( j in seq_along(U_type_co)){
      
    tab=info_extract_co_occur(RICE_co,U_ref_co[i],U_type_co[j])
    U_sampsize=unique(tab$mat$sampSize)
    
    for ( ind in seq_along(U_sampsize)) {
    res_co_tab[[z]]=tab$ls_mat
    matsize=tab$mat[which(tab$mat$sampSize==U_sampsize[ind]),]  
    res_co_data[[z]]=data.frame(sampMatbased=query,
                                sampSize=U_sampsize[ind],
                                paramType=paste(unique(matsize$paramType),collapse = "+"),
                                sampCountry=paste(unique(matsize$sampCountry),collapse=";"),
                                sampCountryorigin=paste(unique(matsize$sampCountry),collapse=";"),
                                Records=nrow(matsize),
                                Ref=U_ref_co[i])
                                
     z=z+1                           
    
    }

  }
  }
  
#########################################################################################################


# tabella RICE+Feed
# Paramtype	Sampsize	sampCountry	sampCountryOrigin	Records number
# DON+HT2+ADON	8	ES	ES

#########################################################################################################


data_tab=data.frame(do.call("rbind",res_co_data))
data_tab=data_tab[which(data_tab$paramType != 0),]
file.remove(paste0("Tabella_co_",query,".xls"))
XLConnect::writeWorksheetToFile(paste0("Tabella_co_",query,".xls"),data_tab,query)



#########################################################################################################
# References

# https://stats.stackexchange.com/questions/153725/plotting-a-ecdf-in-r-and-overlay-cdf
# https://stats.stackexchange.com/questions/197607/how-to-test-difference-between-times-series-does-time-series-anova-exist
# http://rstudio-pubs-static.s3.amazonaws.com/5554_25ed8319163a4df6bd644e68c6fd4b21.html
# https://cran.r-project.org/web/packages/codyn/vignettes/Temporal_Diversity_Indices.html
# https://www.r-bloggers.com/analysing-longitudinal-data-multilevel-growth-models-i/
# https://rstudio-pubs-static.s3.amazonaws.com/228019_f0c39e05758a4a51b435b19dbd321c23.html
