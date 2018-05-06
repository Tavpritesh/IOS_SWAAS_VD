
rm(list=ls())

detach(dat)
source("./read.erg.R")


require(zoo)
require(RSEIS)


#dirs <- list.files()
#for(j in seq(along=dirs)) 
#{

files <- list.files(pattern=".erg$")


rname1 <- paste("VD_Indices.csv")
rname2 <- paste("VD_Spectra.csv")

#setwd(dirs[j])
#files <- list.files(pattern=".erg$")
labels=c("Asthmatics")
resultFile_Indices=c("Results_Indices")
resultFile_Spectra=c("Results_Spectra")
pdf("PSD_Plots.pdf")

for(i in seq(along=files))
{	
o_dat<-read.erg(files[i])
dat<-data.frame(o_dat[,c(5,8,11,17:28)])
attach(dat)


#a<-seq(1:length(R5),by=1)
#R5
#dat<-dat[,-c(1:7,29)]


#out<-array(0,c(1,5,15))		############ UnComment later ##############

#for(i in seq(1,15,by=1))
#{
		#R5<-dat[,5]			############ UnComment later ##############

#R5[which(R5=="**********")]<-NA
#R5<-as.numeric(R5)
#R5

#plot(R5,type="l")			############ UnComment later ##############

#R5<-R5-mean(R5)
even<-seq(2,length(R5),by=2)
odd<-seq(1,length(R5),by=2)

R5_even<-R5[even]

#plot(R5_even,type="l")		############ UnComment later ##############

R5_odd<-R5[odd]

#plot(R5_odd,type="l")		############ UnComment later ##############

### Outlier detection for Even Component

e_ts <- ts(R5_even,frequency=2.94)
e_ts<-ts(na.spline(e_ts),frequency=2.94)

#plot(e_ts,type="l")		############ UnComment later ##############

e_resid <- stl(e_ts,s.window="periodic",s.degree=1,t.window=1001)$time.series[,3]
e_resid.q <- quantile(e_resid,prob=c(0.25,0.75))
e_iqr <- diff(e_resid.q)
e_limits <- e_resid.q + 1.5*e_iqr*c(-1,1)
e_score <- abs(pmin((e_resid-e_limits[1])/e_iqr,0) + pmax((e_resid - e_limits[2])/e_iqr,0))
R5_even[e_score>0]<-NA
Even_Imputed<-na.spline(R5_even)
Even_Imputed<-Even_Imputed-(mean(Even_Imputed))

#plot(Even_Imputed,type="l")		############ UnComment later ##############

Even_Imputed[1:2]<-0
Even_Imputed[(length(Even_Imputed)-2):length(Even_Imputed)]<-0

#plot(Even_Imputed,type="l")

### Outlier detection for Odd Component

o_ts <- ts(R5_odd,frequency=2.94)
o_ts<-ts(na.spline(o_ts),frequency=2.94)

#plot(o_ts,type="l")			############ UnComment later ##############

o_resid <- stl(o_ts,s.window="periodic",s.degree=1,t.window=1001)$time.series[,3]
o_resid.q <- quantile(o_resid,prob=c(0.25,0.75))
o_iqr <- diff(o_resid.q)
o_limits <- o_resid.q + 1.5*o_iqr*c(-1,1)
o_score <- abs(pmin((o_resid-o_limits[1])/o_iqr,0) + pmax((o_resid - o_limits[2])/o_iqr,0))
R5_odd[o_score>0]<-NA
Odd_Imputed<-na.spline(R5_odd)
Odd_Imputed<-Odd_Imputed-(mean(Odd_Imputed))

#plot(Odd_Imputed,type="l")		############ UnComment later ##############

Odd_Imputed[1:2]<-0
Odd_Imputed[(length(Odd_Imputed)-2):length(Odd_Imputed)]<-0

#plot(Odd_Imputed,type="l")

### Reconstruct the merged R5 series

R5_Reconst<-c(rbind(Odd_Imputed,Even_Imputed))

#par(mfrow=c(2,1))						############ UnComment later ##############
#plot(R5,type="l",col="blue")			############ UnComment later ##############
#plot(R5_Reconst,type="l",col="red")	############ UnComment later ##############


##########			Calculation of PSD of R5		####################

y<-R5_Reconst
coef=2048
#win=256	#### Doesnt work for smaller data lengths   #########
win=64
inc=win/2
w<-setwelch(y, win=win, inc=inc, coef=coef, wintaper=0.2)
KK = apply(w$values, 2, FUN="mean")
dt<-0.17
fw=seq(from=0, to=0.5, length=coef)/(dt)
Wyy = (KK^2)/w$windowsize



############		Defining VD on the basis of a range			#########

VD_range = sum(Wyy[71:1393])			## Integral of PSD using 0-2 hertz
ID_range = sum(Wyy[1394:(length(Wyy))])	## Integral of PSD using 2-3 hertz
#legend("topleft", legend= "VD(range)=" VD_range)
#legend("topright", legend= "ID(range)=" ID_range)


############		Calculation of PSD of Volume				#########

yv<-Vol-mean(Vol)
wv<-setwelch(yv, win=win, inc=inc, coef=coef, wintaper=0.2)
KKv = apply(wv$values, 2, FUN="mean")
dt<-0.17
fw=seq(from=0, to=0.5, length=coef)/(dt)
Wvyy = (KKv^2)/w$windowsize
#plot(fw[71:length(fw)],Wvyy[71:length(Wvyy)] , col='red',type="l",main="PSD of Volume",xlab="Power (Watts/Hz)",ylab="frequency (Hz)")
#points(fw[Index_maxim],Wvyy[Index_maxim],col="red",pch=15)
VVD_range = sum(Wvyy[71:1393])

#VID = sum(Wvyy[(length(Wyy)/2):length(Wvyy)])
#legend("topleft", legend= VVD)
#legend("topright", legend= VID)

Norm_VD_range<-log(VD_range/VVD_range,base=10)

#Norm_ID<-log(ID/VID,base=10)



##########   Defining VD on the basis of a sharp peak    #############

maxim<-(max(Wvyy[71:length(fw)]))
Index_maxim<-which(Wvyy==maxim)
VD_sharp = sum(Wyy[(Index_maxim-70):(Index_maxim+70)])
VVD_sharp = sum(Wvyy[(Index_maxim-70):(Index_maxim+70)])
Norm_VD_sharp<-log(VD_sharp/VVD_sharp,base=10)

maxim_VD<-(max(Wyy[71:(length(fw)/2)]))

###########			Plotting PSDs			##############

par(mfrow=c(2,1))		############ UnComment later ##############
par(xpd=TRUE)
plot(fw[71:length(fw)],Wvyy[71:length(Wvyy)] , col='purple',type="l",main="PSD of Volume",ylab="Power(Watts/Hz)",xlab="frequency (Hz)",lwd=2)
points(fw[Index_maxim],Wvyy[Index_maxim],col="red",pch=15)

plot(fw[71:length(fw)],Wyy[71:length(Wyy)] , col='purple',type="l",main="PSD of R5",ylab="Power(Watts/Hz)",xlab="frequency (Hz)",lwd=2,ylim=c(0,(maxim_VD+0.08)))
#plot(fw[71:length(fw)],Wyy[71:length(Wyy)] , col='purple',type="l",main="PSD of R5",ylab="Power(Watts/Hz)",xlab="frequency (Hz)",lwd=2,ylim=c(0,(maxim_VD+(0.8*(maxim_VD)))))

points(fw[Index_maxim],Wyy[Index_maxim],col="red",pch=15)		############ UnComment later ##############



#x_line<-seq((Index_maxim-70),(Index_maxim+70),by=1)
#y_line<-rep(maxim_VD,141)
#lindat<-cbind(x_line,y_line)
#lines(lindat,lwd=2,col="red")

#lines(cumsum(Wyy),col="blue")

#legend("top", legend= c("Normalized VD(range)=", Norm_VD_range))		############ UnComment later ##############
#legend("bottom", legend= c("Normalized VD(sharp)=", Norm_VD_sharp))	############ UnComment later ##############


#maxim_VD<-(max(Wyy[71:(length(fw)/2)]))

Index_maxim_VD<-which(Wyy==maxim_VD)
#for(i in 1:141)
#{
	#points(fw[(Index_maxim-71)+i],(maxim_VD+0.005),col="blue",pch=15,pwd=0.8)
	##points(fw[Index_maxim+70],maxim_VD,col="blue",pch=15)
#}

##########  Uncommented_For Plots till here  ###########

#lines((fw[Index_maxim-70]:fw[Index_maxim+70],maxim_VD,col="blue",lwd=2)

#x_line<-seq((Index_maxim-70),(Index_maxim+70),by=1)
#y_line<-rep(maxim_VD,141)
#lindat<-cbind(x_line,y_line)
#lines(lindat,lwd=2,col="red")

#out[1,1,i]<-VD_range
#out[1,2,i]<-VD_sharp
#out[1,3,i]<-Norm_VD_range
#out[1,4,i]<-Norm_VD_sharp
#out[1,5,i]<-ID_range
result_Indices <- paste(files[i],VD_range,VD_sharp,Norm_VD_range,Norm_VD_sharp,ID_range)
result_Spectra<-paste(files[i],Wyy)
#out<-c(VD_range,VD_sharp,Norm_VD_range,Norm_VD_sharp,ID_range)

resultFile_Indices = c(resultFile_Indices,result_Indices)
resultFile_Spectra = c(resultFile_Spectra,result_Spectra)

detach(dat)
}

write.table(resultFile_Indices,rname1,row.names=F,col.names=F)
write.table(resultFile_Spectra,rname2,row.names=F,col.names=F)
dev.off()
#}