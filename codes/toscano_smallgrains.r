########################################################################
#
#
#
########################################################################

library(readxl)
library(XLConnect)
library(ggplot2)
library(dplyr)
library(scales)

# library(ggpubr)
########################################################################

setwd("/home/alf/Scrivania/aa_lavori/lav_toscano_matrix")

########################################################################
# function

smooth_ecd_gray = function(dat, titleg,subtitleg="SMALLGRAINS",adj=0.1,e=0,labelx="Concentration_op") {
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
  mat_op=as.data.frame(mat[vars_new])
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

# SMALLGRAINS=read_xlsx("small_grains.xlsx")
# SMALLGRAINS[SMALLGRAINS==-999]=NA
# saveRDS(SMALLGRAINS,"SMALLGRAINS.rds")



#########################################################################################################


SMALLGRAINS=readRDS("SMALLGRAINS.rds")

query_vec=unlist(strsplit("almond;barley;buckwheat;caraway;cereals;durum wheat;millet;oat;peanuts;rye;soft wheat;soy;spelt;triticale;walnuts;wheat",";"))



for ( query in query_vec) {

dir.create(paste0("ecdf_simple_SMALLGRAINS","_",query))
dir.create(paste0("ecdf_full_SMALLGRAINS","_",query))

SMALLGRAINSq=SMALLGRAINS[which(SMALLGRAINS$sampMatbased==query),]

SMALLGRAINSq$paramType=gsub("/-","_",tolower(SMALLGRAINSq$paramType))
SMALLGRAINSq$sampMatType=tolower(SMALLGRAINSq$sampMatType)
SMALLGRAINSq$sampMatbased=tolower(SMALLGRAINSq$sampMatbased)


U_param=unique(SMALLGRAINSq$paramType)
U_type=c("food","feed","raw materials")
U_plant=unique(SMALLGRAINSq$sampMatbased)[1]

#########################################################################################################

res_SMALLGRAINS_tab=list()
res_SMALLGRAINS_data=list()
res_SMALLGRAINS_names=list()


#########################################################################################################


z=1
for ( i in seq_along(U_param)) {
      for ( j in seq_along(U_type)){
    
        tab=info_extract(SMALLGRAINSq,U_param[i],U_type[j],query)
        res_SMALLGRAINS_data[[z]]=tab
        res_SMALLGRAINS_tab[[z]]=tab$ls_tab
        res_SMALLGRAINS_names[[z]]=paste0(U_param[i],";",U_type[j],";",paste(query,"on",U_type[j]))
      
        ggplot(tab$mat_op, aes(meanTot_op))+
                   stat_ecdf(geom = "point")+
                   ylab("Density")+
                   ggtitle(paste(U_param[i],"on",U_type[j]),subtitle = query)

        ggsave(filename = paste0(U_param[i],"_on_",U_type[j],"_mtot_simple_SMALLGRAINS_",query,".jpg"),dpi=300)



        ggplot(tab$mat_op, aes(Concentration_op))+
                   stat_ecdf(geom = "point")+
          ylab("Density")+
          ggtitle(paste(U_param[i],"on",U_type[j]),subtitle = query)
          ggsave(filename = paste0(U_param[i],"_on_",U_type[j],"_conc_simple_SMALLGRAINS_",query,".jpg"),dpi=300)

  if (tab$N_conc > 0) {
                    smooth_ecd_gray(as.data.frame(tab$mat_op),titleg = paste(U_param[i],"on",U_type[j]),subtitleg=query,labelx = "Concentration_op")
                    ggsave(filename = paste0(U_param[i],"_on_",U_type[j],"_conc_full_SMALLGRAINS_",query,".jpg"),dpi=300)

                    }
  if (tab$N_mtot > 0) {
                    smooth_ecd_gray(as.data.frame(tab$mat_op),titleg = paste(U_param[i],"on",U_type[j]),subtitleg=query,labelx = "meanTot_op")
                    ggsave(filename = paste0(U_param[i],"_on_",U_type[j],"_mtot_full_SMALLGRAINS_",query,".jpg"),dpi=300)

                    }
    z=z+1;
    }
}

SMALLGRAINSq_tab=data.frame(names=do.call("rbind",res_SMALLGRAINS_names),do.call("rbind",res_SMALLGRAINS_tab))
file.remove(paste0("Tabella_SMALLGRAINS_",query,".xls"))
XLConnect::writeWorksheetToFile(paste0("Tabella_SMALLGRAINS_",query,".xls"),SMALLGRAINSq_tab,query)

}

for ( query in query_vec) {
  system(paste0("mv *simple*",query,"*.jpg ecdf_simple_SMALLGRAINS_",query))
}

for ( query in query_vec) {
  system(paste0("mv *full*",query,"*.jpg ecdf_full_SMALLGRAINS_",query))
}
#########################################################################################################
# References

# https://stats.stackexchange.com/questions/153725/plotting-a-ecdf-in-r-and-overlay-cdf
# https://stats.stackexchange.com/questions/197607/how-to-test-difference-between-times-series-does-time-series-anova-exist
# http://rstudio-pubs-static.s3.amazonaws.com/5554_25ed8319163a4df6bd644e68c6fd4b21.html
# https://cran.r-project.org/web/packages/codyn/vignettes/Temporal_Diversity_Indices.html
# https://www.r-bloggers.com/analysing-longitudinal-data-multilevel-growth-models-i/
# https://rstudio-pubs-static.s3.amazonaws.com/228019_f0c39e05758a4a51b435b19dbd321c23.html