########################################################################

library(readxl)
library(XLConnect)
library(ggplot2)
library(dplyr)
library(scales)

# library(ggpubr)
########################################################################

setwd("/home/alf/Scrivania/lav_toscano_matrix")

########################################################################
# function

smooth_ecd_gray = function(dat, titleg,subtitleg="maize",adj=0.1,e=0,labelx="Concentration_op") {
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

############################################################################################################################################################

#MAIZE=read_xlsx("MaizeDB.xlsx")
#MAIZE[MAIZE==-999]=NA
#saveRDS(MAIZE,"MAIZE.rds")



#########################################################################################################

MAIZE=readRDS("MAIZE.rds")
query="glutenfree"
MAIZE_glutenfree=as.data.frame(MAIZE[grep(query,MAIZE$sampMatInfo),])
MAIZE_glutenfree$paramType=tolower(MAIZE_glutenfree$paramType)
MAIZE_glutenfree$sampMatType=tolower(MAIZE_glutenfree$sampMatType)
MAIZE_glutenfree$sampMatbased=tolower(MAIZE_glutenfree$sampMatbased)


U_param=unique(MAIZE_glutenfree$paramType)
U_type=unique(MAIZE_glutenfree$sampMatType)
U_plant=unique(MAIZE_glutenfree$sampMatbased)[1]

#########################################################################################################

res_plant_tab=list()
res_plant_data=list()
res_plant_names=list()

z=1
for ( i in seq_along(U_param)) {
      for ( j in seq_along(U_type)){
    
        tab=info_extract(MAIZE_glutenfree,U_param[i],U_type[j],U_plant)
        res_plant_data[[z]]=tab
        res_plant_tab[[z]]=tab$ls_tab
        res_plant_names[[z]]=paste0(U_param[i],";",U_type[j],";",query)
      
  #       ggplot(tab$mat_op, aes(meanTot_op))+
  #                  stat_ecdf(geom = "point")+
  #                  ylab("Density")+
  #                  ggtitle(paste(U_param[i],"on",U_type[j]),subtitle = query)
  # 
  #       ggsave(filename = paste0(U_param[i],"_on_",U_type[j],"_mtot_simpl_",query,".jpg"),dpi=300)
  # 
  # 
  # 
  #       ggplot(tab$mat_op, aes(Concentration_op))+
  #                  stat_ecdf(geom = "point")+
  #         ylab("Density")+
  #         ggtitle(paste(U_param[i],"on",U_type[j]),subtitle = query)
  #         ggsave(filename = paste0(U_param[i],"_on_",U_type[j],"_conc_simple_",query,".jpg"),dpi=300)
  # 
  # if (tab$N_conc > 0) {
  #                   smooth_ecd_gray(as.data.frame(tab$mat_op),titleg = paste(U_param[i],"on",U_type[j]),labelx = "Concentration_op")
  #                   ggsave(filename = paste0(U_param[i],"_on_",U_type[j],"_conc_full_",query,".jpg"),dpi=300)
  # 
  #                   }
  # if (tab$N_mtot > 0) {
  #                   smooth_ecd_gray(as.data.frame(tab$mat_op),titleg = paste(U_param[i],"on",U_type[j]),labelx = "meanTot_op")
  #                   ggsave(filename = paste0(U_param[i],"_on_",U_type[j],"_mtot_full_",query,".jpg"),dpi=300)
  # 
  #                   }
    z=z+1;
    }
}

data_tab=data.frame(names=do.call("rbind",res_plant_names),do.call("rbind",res_plant_tab))
file.remove(paste0("Tabella_",query,".xls"))
XLConnect::writeWorksheetToFile(paste0("Tabella_",query,".xls"),data_tab,query)
XLConnect::writeWorksheetToFile(paste0("Matrice_",query,".xls"),MAIZE_glutenfree,query)



#########################################################################################################
# References

# https://stats.stackexchange.com/questions/153725/plotting-a-ecdf-in-r-and-overlay-cdf
# https://stats.stackexchange.com/questions/197607/how-to-test-difference-between-times-series-does-time-series-anova-exist
# http://rstudio-pubs-static.s3.amazonaws.com/5554_25ed8319163a4df6bd644e68c6fd4b21.html
# https://cran.r-project.org/web/packages/codyn/vignettes/Temporal_Diversity_Indices.html
# https://www.r-bloggers.com/analysing-longitudinal-data-multilevel-growth-models-i/
# https://rstudio-pubs-static.s3.amazonaws.com/228019_f0c39e05758a4a51b435b19dbd321c23.html
