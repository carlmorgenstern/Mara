## libraries ----
require(tidyverse)
require(zoo)
require(WaveletComp)
require(xts)
require(gridExtra)
require(patchwork)
require(ggpmisc)
require(ggpubr)
require(plotly)

## runtime Options ----
options("max.contour.segments" = 500000)
options(digits=3)
Sys.setenv(TZ='UTC')
setwd("/Users/carlmorgenstern/Mara/Wavelets/R")

PDF <- FALSE
Plotly <- TRUE
CrossSpectrum <- FALSE

## Date Period ----
StartDate <- "2011-10-01" 
EndDate <- "2012-10-01"     ### MUST be at lease 1 year later!

Interval <- paste(StartDate,"/",EndDate,sep="")
FirstWaterYear <- year(StartDate)
LastWaterYear <- year(EndDate)
NumberofYears <- LastWaterYear -FirstWaterYear
WaterYearInterval <- paste0("Water Years Oct 1 ",FirstWaterYear," - Oct 1 ",LastWaterYear)
BodiesOfWater <- c("Lagoon")                    # Just the Lagoon
#BodiesOfWater <-c("Tides","Lagoon","Eta")      #datasets to analyze

## Constants and Frequencies ----
dT <- 1/2                      # 30 minutes between samples
dJ <- 1/128                    # octave resolution
upperPeriod <- 33        
lowerPeriod <- 4
spectrumWindowSize <- 5    
Nsim <- 10                     # Number of Wavelet simulations
trim <- 10                    # hack to truncate spikes of  tails of graph
#PlotlyFactor <- 2             # hack to thin out wavelet time axis for plot_ly rendering

Frequencies <- c(
0.0402306,
0.0410103,
0.0416667,
0.08333333,
0.16666667,
0.33333333
)  

Periods <- 1/(Frequencies)
PeriodIndex <- 0

if(PDF){pdf(paste0("Y",StartDate,"-",EndDate,"-dj" ,1/dJ,".pdf"))}

## Read Data and fix timestamps ----
Tides <- read_csv("~/Mara/csv/CabrilloTides_2000to2020.csv")
Lagoon <- read_csv("~/Mara/csv/laWY11_21.csv")
Tides$time <-as.POSIXct(Tides$time, format = " %d-%b-%Y %H:%M:%S", tz="UTC")  
Lagoon$time <- as.POSIXct(Lagoon$time, format = " %H:%M:%S %m/%d/%Y" , tz="UTC") 

##  create and merge time series ----
Tides <- xts(Tides$Inst,Tides$time)
Lagoon <- xts(Lagoon$Inst,Lagoon$time)
Tides <- Tides[Interval]
Lagoon <- Lagoon[Interval]
z <- merge(Tides,Lagoon)
z$Lagoon <- na.approx(z$Lagoon,rule=2)
Z <- z[index(Lagoon),]
Z$Eta <- Z$Lagoon - Z$Tides
Z <- na.omit(Z)
seriesLength <- length(Z)
# 
# ## create data frame W to analyze ----
# W <- data.frame(date=as.POSIXct(index(Z)), tides=Z$Tides,lagoon=Z$Lagoon,level=Z$Eta)
# row.names(W) <- NULL

WY <- function(date){ ifelse(month(date) %in% c(10,11,12),
                             year(date),year(date)-1)}

W <- data.frame(date=as.POSIXct(index(Z)), tides=Z$Tides,lagoon=Z$Lagoon,level=Z$Eta)
W <- transform(W,WaterYear=WY(date))
W <- W[,c(5,1,2,3,4)]
row.names(W) <- NULL

head(W)

## Plot the time series. ----
plot_data_column <- function (Data, column) {
    ggplot(Data, aes_string(x = W[,2],y=column)) +
        geom_line() +
        ggtitle(paste(WaterYearInterval, column)) +
            theme(plot.title = element_text(size=10),
                  axis.text.x=element_text(size=5.5),
                  plot.subtitle = element_text(size=3.5))+
        xlab("date") 
  #        scale_x_date(date_breaks = "months" , date_labels = "%b-%y")
}
# plot_lagoon_time_slice <- function(Date){
#                   X <- W[W$WaterYear == Date,]
#                   Y[ii] <<- ggplot(X,aes(x=date,y=Lagoon)) + geom_line() +
#                        ggtitle(paste0("Lagoon -- Water Year ",Date))
# }

g <- list()

for(i in FirstWaterYear:(LastWaterYear-1)){
     g[[i-FirstWaterYear+1]] <- ggplot(W[W$WaterYear==i,],aes(x=date,y=Lagoon)) + 
               geom_line() +
           ggtitle(paste0("Lagoon -- Water Year ",i))+
          theme(plot.title=element_text(size=10))
     }

 wrap_plots(lapply(colnames(W)[c(-1:-3,-5)], plot_data_column, Data = W),ncol=1)
#wrap_plots(lapply(colnames(W)[-1:-2], plot_data_column, Data = W[2:5]),ncol=1)

## Spectra ----
STides <- spec.pgram(W$Lagoon,spans=spectrumWindowSize,log="yes",
                     demean=TRUE,detrend=TRUE,plot=FALSE)
SLagoon <- spec.pgram(W$Lagoon,spans=spectrumWindowSize,log="yes",
                     dmean=TRUE,detrend=TRUE,plot=FALSE)
SEta <- spec.pgram(W$Lagoon,spans=spectrumWindowSize,log="yes",
                     demean=TRUE,detrend=TRUE,plot=FALSE)

S <- data.frame(freq=STides$freq,
                Period=1/(STides$freq/dT),
                TidesPower=log10(STides$spec),
                LagoonPower=log10(SLagoon$spec),
                EtaPower=log10(SEta$spec)
               )

plot_power_column = function (Data, column) {
    ggplot(Data, aes_string(x = S[,1],y=column)) +
        geom_line() +
        xlim(0.01,0.05) +
          ylim(-2.5,2.5) +
        xlab("Frequency (Hz)")   +
        ggtitle(label=str_replace(column,"Power", " -- Frequency"),
                subtitle = paste(WaterYearInterval) ) +
        theme(plot.title = element_text(size=12))+
        stat_peaks(color="red",span=100,geom="text",
                   ignore_threshold=0.70, show.legend=TRUE,
      #   position = position_nudge(x = 0.003, y = 0.2),
                   x.label.fmt = "%2.3f" )
}
plot_period_column = function (Data, column) {
        ggplot(Data, aes_string(x = S[,2],y=column)) +
        geom_line() +
        xlim(10,30) +
        xlab("Period (hrs)")   +
        ggtitle(label=str_replace(column,"Power", " -- Period"),
                subtitle = paste(WaterYearInterval) ) +
        theme(plot.title = element_text(size=12)) +
        stat_peaks(color="red",span=30,geom="text",
                   ignore_threshold=0.80, show.legend=TRUE,
                   position = position_nudge(x = 0.1, y = 0.5),
                   x.label.fmt = "%2.2f" )
}

H1 <- wrap_plots(lapply(colnames(S)[4], plot_power_column,Data=S),ncol=1)
H2 <- wrap_plots(lapply(colnames(S)[4], plot_period_column,Data=S),ncol=1)
H1+H2

## Wavelet Loop ----

for(i in BodiesOfWater){
eval(parse(text= paste0("Wavelet",i," <-","analyze.wavelet(W,",'"',i,'",',
                               "loess.span=0,
                                dt=dT,
                                dj=dJ,
                                make.pval = FALSE,
                                n.sim=Nsim,
                                date.format =", "'%F %T',
                                lowerPeriod=lowerPeriod,
                                upperPeriod = upperPeriod)")))
   
eval(parse(text=paste0("output",i,"Wavelet <- wt.image(Wavelet",i,
                             ",label.time.axis=TRUE,
                              periodlab  = 'period (hours)',
                              color.key = 'interval',
                              main=paste0(' ",i," wavelet') ,
                              date.format =", "'%F %T',
                              timelab=paste0(StartDate,'  ',EndDate),
                              show.date=TRUE,
                              verbose=TRUE) ")))
} # Wavelet loop

AA <- list()

for(i in BodiesOfWater){
  for(k in FirstWaterYear:(LastWaterYear-1)){
    B <- list()
    C <- list()
    C[[1]] <- g[[1]]  # Lagoon time series
    D <- data.frame()
    
    for(j in 1:length(Periods)){                                  # create vector of indices to plot over time
      AA[[j]] <- list()
      eval(parse(text= paste0("PeriodIndex[j] <- which.min(
                            abs(Wavelet",i,"$Period - Periods[j]))")))
                 
      eval(parse(text= paste0( "A <- data.frame(date=Wavelet",i,"$series$date,period=Periods[j],power=Wavelet",
                                  i,"$Power[PeriodIndex[j],],WaterYear = year(x=Wavelet",i,"$series$date)"," 
                                  )"
                             )))
          
       A <- A[(A$WaterYear == k | A$WaterYear == (k+1)),]
       D <- rbind(D,A)
          
       AA[[k-FirstWaterYear+1]][[j]] <- ggplot(A,aes(x=date,y=power)) + geom_line()  +
                xlab("Date") +
                ylab("Power") + 
                labs(title=paste0(i,"  Power at ", 
                            round(Periods[j],3)," hours.  "),
                            subtitle= paste0( WaterYearInterval, "  (index =",PeriodIndex[j],")")) + 
                theme(plot.title = element_text(size=10),
                      plot.subtitle = element_text(size=6.0)) # +

       B[[j]] <- AA[[k-FirstWaterYear+1]][[j]]
       C[[j+1]] <- B[[j]]
              
   }  # loop on j
       P24 <- filter(D,period == Periods[3])
       P244 <- filter(D,period==Periods[2])
       P249 <- filter(D,period==Periods[1])
       P24$totalPower24 <- P24$power+P244$power+P249$power
       
      G <- ggplot(P24,aes(x=date,y=totalPower24)) + geom_line() + 
                 ylab("total power") + 
                 labs(title = paste0("total power at ", round(Periods[[1]],3),
                                     " + ",round(Periods[[2]],3)," + ",round(Periods[[3]],3)," hours")) +
                 theme(plot.title = element_text(size=10),
                      plot.subtitle = element_text(size=6.0))
       
  # MG <-   marrangeGrob(C,nrow=5,ncol=1)
      MG <- C[[1]] / G + C[[5]] + C[[6]] + C[[7]]
   print(MG)
      
      
   #      reduce(B,`+`)
   #      grid.arrange(grobs=B,ncol=1,nrow=3,newpage=TRUE)
         
#        plot.pages <- ggarrange(plotlist=B,nrow=4,ncol=2)
#ggexport(plot.pages) #,filename='plots.pdf')

 # print(plot.pages)
         
 }  # loop on k
} # loop on i

### Reconstruct Time Series from Wavelets ----
# 
# colors <- c( "original" = "red", "reconstructed" = "green")
# 
# for(i in BodiesOfWater){
#      
# E <- paste0("RR",i,    "<- reconstruct(Wavelet",i,
#                        ", only.sig = FALSE,plot.rec=FALSE)")
# eval(parse(text=E))
# 
# E <-paste0( "R",i,    " <- data.frame(
#                         date=RR",i,"$series$date,
#                         original=RR",i,"$series$",i,
#                       " ,reconstructed=RR",i,"$series$",i,".r,
#                         error = RR",i,"$series$",i,"-RR",i,"$series$",i,".r)")
# eval(parse(text=E))
# 
# E <- paste0("ggplot(R",i,",aes(x=date)) +
#                         geom_line(aes(y=original,color='original'), linewidth=1.0,
#                         linetype='solid')  +
#                         geom_line(aes(y=reconstructed,color='reconstructed'),linewidth=0.4,
#                         linetype='dotted') +
#                         ggtitle('",i," -- Original and Reconstructed Time Series')")
# print(eval(parse(text=E)))
# 
#             
# E <- paste0("X <- R",i,"[trim:(nrow(R",i,")-trim),]" )
# eval(parse(text=E))
# 
# E <- paste0("ggplot(X,aes(x=date,y=error)) + 
#                     geom_line()  +
#                     ggtitle('",i, " Error  (note yaxis range)')")
# print(eval(parse(text=E)))
# 
# }   #for loop


if(CrossSpectrum){
#### compute cross wavelet transform, plot power and coherence ----
wc <-  analyze.coherency(W, my.pair = c("Tides", "Lagoon"), loess.span = 0,
                         dt = dT,
                         dj = dJ,
                         lowerPeriod = lowerPeriod,
                         upperPeriod = upperPeriod,
                         window.type.t = 1, window.type.s = 1,
                         window.size.t = 10, window.size.s = 1/32,
                         make.pval = FALSE,
                         n.sim = Nsim,
                         date.format = "%F %T", date.tz = NULL,
                         verbose = TRUE)

wc.image(wc, which.image="wp",                            #POWER
         main = " Tides and Lagoon: cross-wavelet power",
         plot.arrow = TRUE,
         n.levels = 1000,
         plot.ridge=TRUE,
    #     clear.area = TRUE,
    #     siglvl.area= 1.0,
         color.key="interval",
         legend.params = list(lab = "power levels"),
         label.time.axis = TRUE,
         show.date=TRUE,
         date.format = "%F %T",
         timelab=paste0(StartDate,"  ",EndDate),
         periodlab = "period (hours)")

 wc.image(wc, which.image="wc",                              #COHERENCY
         main = "Tides and Lagoon: cross-wavelet coherency",
         plot.arrow = TRUE,
         n.levels = 1000,
         plot.ridge = TRUE,
         clear.area = TRUE,
         siglvl.area=0.6,
         color.key="interval",
         legend.params = list(lab = "coherency levels"),
         label.time.axis = TRUE,
         show.date=TRUE,
         date.format = "%F %T",
         timelab=paste0(StartDate,"  ",EndDate),
         periodlab = "period (hours)")
 
####  Average cross series power and coherence ----

wc.avg(wc, which.avg = "wp", exponent = 1,
   show.siglvl = TRUE,
   siglvl = c(0.05, 0.1),
   sigcol = c("red", "blue"), sigpch = 20, sigcex = 1,
   minimum.level = NULL, maximum.level = NULL,
   label.avg.axis = TRUE,
   averagelab = NULL, averagetck = 0.02, averagetcl = 0.5,
   spec.avg.axis = list(at = NULL, labels = TRUE,
                        las = 1, hadj = NA, padj = NA),
   label.period.axis = TRUE,
   periodlab = NULL, periodtck = 0.02, periodtcl = 0.5,
   spec.period.axis = list(at = NULL, labels = TRUE,
                           las = 1, hadj = NA, padj = NA),
   show.legend = FALSE, legend.coords = "topright",
   main = "Tides x Lagoon average power",
   lwd = 1, col = 1,
   lwd.axis = 1,
   verbose = TRUE)

wc.avg(wc, which.avg = "wc", exponent = 1,
   show.siglvl = TRUE,
   siglvl = c(0.05, 0.1),
   sigcol = c("red", "blue"), sigpch = 20, sigcex = 1,
   minimum.level = NULL, maximum.level = NULL,
   label.avg.axis = TRUE,
   averagelab = NULL, averagetck = 0.02, averagetcl = 0.5,
   spec.avg.axis = list(at = NULL, labels = TRUE,
                        las = 1, hadj = NA, padj = NA),
   label.period.axis = TRUE,
   periodlab = NULL, periodtck = 0.02, periodtcl = 0.5,
   spec.period.axis = list(at = NULL, labels = TRUE,
                           las = 1, hadj = NA, padj = NA),
   show.legend = FALSE, legend.coords = "topright",
   main = "Tides x Lagoon -- average coherency",
   lwd = 1, col = 1,
   lwd.axis = 1,
   verbose = TRUE)

} # Cross Spectrum


### plotly ----
if(Plotly){
     
     C1 <- ggplotly(C[[1]])
     CG <- ggplotly(G)
     C5 <- ggplotly(C[[5]])
     C6 <- ggplotly(C[[6]])
     C7 <- ggplotly(C[[7]])
     
     t1 <- list(
  family = "Times New Roman",
  color = "red",
  size=20)
  
     CPlotly <- subplot(C1,CG,C5,C6,C7,nrows = 5,shareX = TRUE ) %>% 
                layout(title = list(text = paste0("Water Year ",FirstWaterYear),font = t1,x=10))
     CPlotly

for(i in BodiesOfWater){
     
# eval(parse(text= paste0("Wavelet",i," <-","analyze.wavelet(W,",'"',i,'",',
#                                "loess.span=0,
#                                 dt=dT,
#                                 dj=dJ,
#                                 make.pval = FALSE,
#                                 n.sim=Nsim,
#                                 date.format =", "'%F %T',
#                                 lowerPeriod=lowerPeriod,
#                                 upperPeriod = upperPeriod)")))
#    
# eval(parse(text=paste0("output",i,"Wavelet <- wt.image(Wavelet",i,",
#                               label.time.axis=TRUE,
#                               periodlab  = 'period (hours)',
#                               color.key = 'interval',
#                               main=paste0(' ",i," wavelet') ,
#                               date.format =", "'%F %T',
#                               timelab=paste0(StartDate,'  ',EndDate),
#                               show.date=TRUE,
#                               verbose=TRUE) ")))
     
P <- eval(parse(text= paste0("plot_ly(z = ~Wavelet",i,"$Power ,
        x=Wavelet",i,"$series$date,
        y=Wavelet",i,"$Period,
        type='contour') %>%
  #      add_surface() %>%
        layout(title = '",i," ',
             xaxis = list(title = 'Date'),
             yaxis = list(title = 'Period (Hrs)') ,
           showlegend=FALSE)  %>%
        hide_colorbar() "

            )))
}
     P
}


if(PDF){dev.off()}


