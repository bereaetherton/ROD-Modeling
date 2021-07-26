
library(readxl)
library(rgdal)
library(terra)
library(geosphere)

#------------------------------------------------------------#

#reading in the data

waipun<-read_excel("Waipunalei_2021_05_Shareable.xlsx",col_types = c("numeric"))
waipun<-as.matrix(waipun)
waipun<-waipun[c(1:304,306:8548),c(1,4:35)] #filter out broken data
colnames(waipun)<-c("Node","10/5/2009",	"1/6/2011","10/6/2014",	"4/18/2015",	"8/29/2015",	"5/4/2017",	"11/7/2017",	"12/28/2017",
                    "2/2/2018",	"3/21/2018",	"4/19/2018",	"5/15/2018",	"6/20/2018",	"7/17/2018",	"8/21/2018",	"10/5/2018",	"11/2/2018",	"12/11/2018",
                    "12/17/2018",	"1/2/2019",	"1/24/2019",	"2/28/2019",	"3/23/2019",	"3/27/2019",	"5/7/2019",	"6/7/2019",	"7/26/2019",	"8/30/2019",	"9/27/2019",	
                    "10/31/2019",	"12/10/2019",	"1/17/2020")


wp.utm<-read_excel("Waipunalei_2021_05_Shareable.xlsx")
wp.utm<-data.matrix(wp.utm)
wp.utm<-wp.utm[c(1:304,306:8548),2:3]
wp.size<-vect(wp.utm,crs="+proj=utm +zone=5n +datum=WGS84  +units=m") #convert UTM to lat/long
wp<-project(wp.size, "+proj=longlat +datum=WGS84")
wp<-geom(wp)[,c("x","y")]
colnames(wp)<-c("Long","Lat")

rm(wp.utm)
rm(wp.size)

size<-dim(wp)[1] #8547
time<-dim(waipun)[2]-1 #32

waipun[which(waipun!=4)]<-1 #making the data easier to work with 
waipun[which(waipun==4)]<-0 #by converting it to 0s and 1s
waipun[,1]<-1:size
dates<-list("Node","10/5/2009",	"1/6/2011","10/6/2014",	"4/18/2015",	"8/29/2015",	"5/4/2017",	"11/7/2017",	"12/28/2017",
            "2/2/2018",	"3/21/2018",	"4/19/2018",	"5/15/2018",	"6/20/2018",	"7/17/2018",	"8/21/2018",	"10/5/2018",	"11/2/2018",	"12/11/2018",
            "12/17/2018",	"1/2/2019",	"1/24/2019",	"2/28/2019",	"3/23/2019",	"3/27/2019",	"5/7/2019",	"6/7/2019",	"7/26/2019",	"8/30/2019",	"9/27/2019",	
            "10/31/2019",	"12/10/2019",	"1/17/2020")

waipun<-waipun[1:200,]
wp<-wp[1:200,]
size<-200

#------------------------------------------------------------#

#read in the arguments for beta and alpha

args<-commandArgs(TRUE)
ARG<-as.numeric(strsplit(args[1],",")[[1]])

ARG=c(2.6, 1.0)

print(ARG)

beta<-ARG[1]
alpha<-ARG[2]

#-------------------------------------------------------------#

#finding the days between the dates

days.list<-matrix(0,nrow=1,ncol=time)
days.list[1]<-0
for(i in 2:time){
  days.list[i]<-as.numeric(difftime(as.Date(dates[[i+1]],"%m/%d/%Y"),as.Date(dates[[i]],"%m/%d/%Y"),units ="days"))
  days.list[i]<-days.list[i-1]+days.list[i]
}


#----------------------------------------------------------#

waipun<-waipun[,7:33] #filtering out the timeline 
days<-988 #988 days between 5/4/2017 and 1/17/2020

#creating a distance matrix

wp.dist.mat<-matrix(0,ncol=size,nrow=size)
for(i in 1:size){
  for(j in 1:size){
    if(wp.dist.mat[j,i]==0){
      wp.dist.mat[i,j]<-distGeo(wp[(i),c(1,2)],wp[(j),c(1,2)])
    }
  }
}

ivp.wp.mat<-matrix(0,nrow=size,ncol=size)

#using distance matrix to creating the inverse power law probabilities

for(i in 1:size){
  for(j in 1:size){
    if(wp.dist.mat[i,j]==0){
      ivp.wp.mat[i,j]<-0
    }
    else{
      ivp.wp.mat[i,j]<-alpha*wp.dist.mat[i,j]^(-beta)
    }
  }
}

min<-min(ivp.wp.mat[which(ivp.wp.mat>0)])
max<-max(ivp.wp.mat)
mm<-max-min

#rescaling between 0 and 1

rs.ivp.wp<-matrix(0,nrow=size,ncol=size)

for(i in 1:size){
  for(j in 1:size){
    if(ivp.wp.mat[i,j]==0){
      rs.ivp.wp[i,j]<-0
    }
    else{
      rs.ivp.wp[i,j]<-((ivp.wp.mat[i,j]-min)/mm)
    }
  }
}

rm(ivp.wp.mat)
rm(max)
rm(min)
rm(mm)

rs.ivp.wp<-alpha*rs.ivp.wp
print("done with distance mat!")

#----------------------------------------------------------#

#1838 days between 06/01/2016 and 01/20/2021

wp.spread.mat<-matrix(0,nrow=size,ncol=days)
#wp.spread.mat[,1]<-init.dat.wp
wp.spread.mat[,1]<-waipun[,1]

bp.spread<-function(k,w){
  s<-sample(c(0,1),1,replace = TRUE,prob=c(1-rs.ivp.wp[k,w],rs.ivp.wp[k,w]))
  return(s)
}

for(t in 2:days){
  for(i in 1:size){
    if(wp.spread.mat[i,t-1]==1){
      wp.spread.mat[i,t]<-1
    }
    else{
      for(j in 1:size){  
        if(rs.ivp.wp[i,j]==0){
          if(bp.spread(j,i)==1 && wp.spread.mat[j,t-1]==1){
            wp.spread.mat[i,t]<-1
            break
          }
        }
        else{
          if(bp.spread(i,j)==1 && wp.spread.mat[j,t-1]==1){
            wp.spread.mat[i,t]<-1
            break
          }
        }
      }
    }
  }
  print(t)
  save.image("test.RData")
}

#-----------------------------------------------------------#

dl<-days.list[6:32]-2768
dl[1]<-1

#-----------------------------------------------------------#

NAME <- paste(ARG[1], ARG[2], sep="_")  
save.image(paste0("/blue/garrett/betherton/ROD/","WPIVP.test", NAME, ".RData"))

