## libraries ----
require(tidyverse)
require(plotly)
require(shiny)
require(lubridate)
require(zoo)
require(WaveletComp)
require(xts)
require(readr)
require(gridExtra)
require(patchwork)
require(ggpmisc)

## Date Period ----
StartDate <- "2019-01-01" 
EndDate <- "2019-12-31"

## constants ----
options("max.contour.segments" = 500000)
options(digits=3)
Sys.setenv(TZ='UTC')
setwd("/Users/carlmorgenstern/Mara/Wavelets/R")
outputFileName <- "testFile.pdf"

dT <- 0.5                      # 30 minutes
dJ <- 1/128                    # octave resolution
upperPeriod <- 60          
lowerPeriod <- 6
spectrumWindowSize <-3    
Nsim <- 3                     # Number of Wavelet simulations
trim <- 10                    # hack to truncate spikes of  tails of graph


Interval <- paste(StartDate,"/",EndDate,sep="")
WaterYear <- year(StartDate)

Frequencies <- c(
# 0.0096826)
# 0.0104408,
# 0.0197724,
# 0.0201314,
# 0.0208280, 
# 0.0298059,
# 0.0305775,
# 0.0402306,
# 0.0410103)
  0.041666667,
  0.08333333,
  0.16666667
)  

Periods <- 1/(Frequencies)

## ----
pdf("testoutput.pdf")

## Read Data and fix timestamps ----
Tides <- read_csv("~/Mara/csv/CabrilloTides_2000to2020.csv")
Lagoon <- read_csv("~/Mara/csv/laWY11_21.csv")
Tides$time <-as.POSIXct(Tides$time, format = " %d-%b-%Y %H:%M:%S", tz="UTC")  
Lagoon$time <- as.POSIXct(Lagoon$time, format = " %H:%M:%S %m/%d/%Y" , tz="UTC") 

##  create and merge time series ---
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

## create data frame W to analyze ----
W <- data.frame(date=as.POSIXct(index(Z)), tides=Z$Tides,lagoon=Z$Lagoon,level=Z$Eta)
row.names(W) <- NULL

head(W)

## Plot the time series. ----
plot_data_column = function (Data, column) {
    ggplot(Data, aes_string(x = W[,1],y=column)) +
        geom_line() +
        ggtitle(paste(WaterYear, column)) +
            theme(plot.title = element_text(size=10),axis.text.x=element_text(size=5.5))+
        xlab("date") 
  #        scale_x_date(date_breaks = "months" , date_labels = "%b-%y")
}
wrap_plots(lapply(colnames(W)[-1], plot_data_column, Data = W),nrow=1)
wrap_plots(lapply(colnames(W)[-1], plot_data_column, Data = W),ncol=1)

## Spectra ----
STides <- spec.pgram(W$Tides,spans=spectrumWindowSize,log="yes",
                     demean=TRUE,detrend=TRUE,plot=FALSE)
SLagoon <- spec.pgram(W$Lagoon,spans=spectrumWindowSize,log="yes",
                     dmean=TRUE,detrend=TRUE,plot=FALSE)
SEta <- spec.pgram(W$Eta,spans=spectrumWindowSize,log="yes",
                     demean=TRUE,detrend=TRUE,plot=FALSE)

S <- data.frame(freq=STides$freq,
                Period=1/(STides$freq/dT),
                TidesPower=log10(STides$spec),
                LagoonPower=log10(SLagoon$spec),
                EtaPower=log10(SEta$spec)
               )

column_name <- function(column){
     str_replace(column,"Power"," -- Power")
}

plot_power_column = function (Data, column) {
    ggplot(Data, aes_string(x = S[,1],y=column)) +
        geom_line() +
        xlim(0,0.05) +
        xlab("Frequency (Hz)")   +
        ggtitle(paste(WaterYear, str_replace(column,"Power", " -- frequency"))) +
        theme(plot.title = element_text(size=10))+
        stat_peaks(color="red",span=5,geom="text",
                   ignore_threshold=0.75, show.legend=TRUE,
                   position = position_nudge(x = 0.01, y = 0.1),
                   x.label.fmt = "%2.3f" )
}

wrap_plots(lapply(colnames(S)[3:5], plot_power_column,Data=S),nrow=1)
wrap_plots(lapply(colnames(S)[3:5], plot_power_column,Data=S),ncol=1)


plot_period_column = function (Data, column) {
        ggplot(Data, aes_string(x = S[,2],y=column)) +
        geom_line() +
        xlim(10,30) +
        xlab("Period (hrs)")   +
        ggtitle(paste(WaterYear, str_replace(column,"Power", " -- period"))) +
        theme(plot.title = element_text(size=10)) +
        stat_peaks(color="red",span=30,geom="text",
                   ignore_threshold=0.80, show.legend=TRUE,
                   position = position_nudge(x = 0.1, y = 0.5),
                   x.label.fmt = "%2.2f" )
}

wrap_plots(lapply(colnames(S)[3:5], plot_period_column,Data=S),nrow=1)
wrap_plots(lapply(colnames(S)[3:5], plot_period_column,Data=S),ncol=1)


##  Tides Wavelet  ----
WaveletTides <- analyze.wavelet(W,"Tides",
                                loess.span=0,
                                dt=dT, 
                                dj=dJ, 
                                make.pval = FALSE,
                                n.sim=Nsim,
                                date.format = "%F %T",
                                lowerPeriod=lowerPeriod,
                                upperPeriod = upperPeriod)

WaveletTidesImage <- wt.image(WaveletTides,
                              label.time.axis=TRUE,
                              periodlab  = "period (hours)",
                              color.key = "interval",
                              main=paste0(WaterYear," Tides wavelet") ,
                              date.format = "%F %T",
                              timelab=paste0(StartDate,"  ",EndDate),
                              show.date=TRUE,
                              verbose=TRUE)

## Lagoon Wavelet  ----
WaveletLagoon <- analyze.wavelet(W,"Lagoon",
                                 loess.span=0,
                                 dt=dT,
                                 dj=dJ,
                                 make.pval = FALSE,
                                 n.sim=Nsim,
                                 date.format = "%F %T",
                                 lowerPeriod=lowerPeriod,
                                 upperPeriod = upperPeriod)

WaveletLagoonImage <- wt.image(WaveletLagoon,
                               label.time.axis=TRUE,
              #                 useRaster = FALSE,
                              color.key = "interval",
                               periodlab  = "period (hours)",
                               main=paste0(WaterYear,"  Lagoon wavelet") ,
                               date.format = "%F %T",
                               timelab=paste0(StartDate,"  ",EndDate),
                               show.date=TRUE,
                               verbose=TRUE)
 
### Eta Wavelet ----
WaveletEta <- analyze.wavelet(W,"Eta",
                                 loess.span=0,
                                 dt=dT,
                                 dj=dJ,
                                 make.pval = FALSE,
                                 n.sim=Nsim,
                                 date.format = "%F %T",
                                 lowerPeriod=lowerPeriod,
                                 upperPeriod = upperPeriod)

WaveletEtaImage <- wt.image(WaveletEta,
                               label.time.axis=TRUE,
              #                 useRaster = FALSE,
                              color.key = "interval",
                               periodlab  = "period (hours)",
                               main=paste0(WaterYear,"  Eta wavelet") ,
                               date.format = "%F %T",
                               timelab=paste0(StartDate,"  ",EndDate),
                               show.date=TRUE,
                               verbose=TRUE)
### Reconstruct Tides ----
RReconstructedTides  <- reconstruct(WaveletTides, only.sig = FALSE,plot.rec=FALSE)
ReconstructedTides <- data.frame(date=RReconstructedTides$series$date,
                                  original=RReconstructedTides$series$Tides,
                                  reconstructed=RReconstructedTides$series$Tides.r,
                                  error = RReconstructedTides$series$Tides-RReconstructedTides$series$Tides.r)

                            
colors <- c( "original" = "red", "reconstructed" = "green")
ggplot(ReconstructedTides,aes(x=date)) + 
               geom_line(aes(y=original,color="original"), linewidth=1.0, 
                         linetype="solid")  +
               geom_line(aes(y=reconstructed,color="reconstructed"),linewidth=0.4,
                         linetype="dotted") +
               ggtitle("Tides -- Original and Reconstructed Time Series")
X <-ReconstructedTides[trim:(length(ReconstructedTides$date)-trim),]
                ggplot(X,aes(x=date,y=error)) + geom_line() + 
                ggtitle("Tides Error  (note yaxis range)")

     
### Reconstruct Lagoon. ----
RReconstructedLagoon <- reconstruct(WaveletLagoon, only.sig = FALSE,plot.rec=FALSE)
ReconstructedLagoon <- data.frame(date=RReconstructedLagoon$series$date,
                                  original=RReconstructedLagoon$series$Lagoon,
                                  reconstructed=RReconstructedLagoon$series$Lagoon.r,
                                  error = RReconstructedLagoon$series$Lagoon-RReconstructedLagoon$series$Lagoon.r)
                            
ggplot(ReconstructedLagoon,aes(x=date)) + 
      geom_line(aes(y=original,color="original"), linewidth=1.0, 
                         linetype="solid")  +
               geom_line(aes(y=reconstructed,color="reconstructed"),linewidth=0.3,
                         linetype="dotted") +
      ggtitle("Lagoon -- Original and Reconstructed Time Series")

X <- ReconstructedLagoon[trim:(length(ReconstructedLagoon$date)-trim),]
      ggplot(X,aes(x=date,y=error)) + geom_line() + 
      ggtitle("Lagoon Error (note Yaxis range)")

### Reconstruct Eta ----
RReconstructedEta <- reconstruct(WaveletEta, only.sig = FALSE,plot.rec=FALSE)
ReconstructedEta <- data.frame(date=RReconstructedEta$series$date,
                                  original=RReconstructedEta$series$Eta,
                                  reconstructed=RReconstructedEta$series$Eta.r,
                                  Eta = RReconstructedEta$series$Eta-RReconstructedEta$series$Eta.r)
                            
ggplot(ReconstructedEta,aes(x=date)) + 
      geom_line(aes(y=original,color="original"), linewidth=1.0, 
                         linetype="solid")  +
      geom_line(aes(y=reconstructed,color="reconstructed"),linewidth=0.9,
                         linetype="dotted") + 
      ggtitle("Eta -- Original and Reconstructed Time Series")

X <- ReconstructedEta[trim:(length(ReconstructedEta$date)-trim),]
      ggplot(X,aes(x=date,y=Eta)) + geom_line() + 
      ggtitle("Eta Error (note Yaxis range)")



# ### Tides Phases  ---- 
#  
# TidesPhases <- wt.sel.phases(WaveletTides, my.series = 1,
#              sel.period = 24, sel.lower = 23, sel.upper = 25,
#               only.coi = TRUE,
#               only.sig = FALSE, siglvl = 0.05,
#               show.avg.phase = TRUE, phase.avg.col = "green",
#               label.time.axis = TRUE,
#               show.date = TRUE, date.format = "%F %T", date.tz = NULL,
#               timelab = NULL, timetck = 0.02, timetcl = 0.5,
#               spec.time.axis = list(at = NULL, labels = TRUE,
#                                     las = 1, hadj = NA, padj = NA),
#               label.phase.axis = TRUE,
#               phaselab = NULL, phasetck = 0.02, phasetcl = 0.5,
#               spec.phase.axis = list(at = NULL, labels = TRUE,
#                                      las = 1, hadj = NA, padj = NA),
#               main = "Tide Phase", sub = NULL,
#               lwd = 0.1, lwd.axis = 1,
#               verbose = TRUE) 
#       # G <- data.frame(date = TidesPhases$date, Phase = TidesPhases$Phase)
#       # ggplot(G[1:500,],aes(x=date,y=Phase)) + geom_line() 
# ### Lagoon Phases  ---- 
#  
# wt.sel.phases(WaveletLagoon, my.series = 1,
#               sel.period = 24, sel.lower = 23, sel.upper = 25,
#               only.coi = FALSE,
#               only.sig = FALSE, siglvl = 0.05,
#               show.avg.phase = FALSE, phase.avg.col = "black",
#               label.time.axis = TRUE,
#               show.date = TRUE, date.format = "%F %T", date.tz = NULL,
#               timelab = NULL, timetck = 0.02, timetcl = 0.5,
#               spec.time.axis = list(at = NULL, labels = TRUE,
#                                     las = 1, hadj = NA, padj = NA),
#               label.phase.axis = TRUE,
#               phaselab = NULL, phasetck = 0.02, phasetcl = 0.5,
#               spec.phase.axis = list(at = NULL, labels = TRUE,
#                                      las = 1, hadj = NA, padj = NA),
#               main = "Lagoon Phase", sub = NULL,
#               lwd = 1, lwd.axis = 1,
#               verbose = TRUE)
# ### Eta Phases  ---- 
#  
# wt.sel.phases(WaveletEta, my.series = 1,
#               sel.period = 24, sel.lower = 23, sel.upper = 25,
#               only.coi = FALSE,
#               only.sig = FALSE, siglvl = 0.05,
#               show.avg.phase = FALSE, phase.avg.col = "black",
#               label.time.axis = TRUE,
#               show.date = TRUE, date.format = "%F %T", date.tz = NULL,
#               timelab = NULL, timetck = 0.02, timetcl = 0.5,
#               spec.time.axis = list(at = NULL, labels = TRUE,
#                                     las = 1, hadj = NA, padj = NA),
#               label.phase.axis = TRUE,
#               phaselab = NULL, phasetck = 0.02, phasetcl = 0.5,
#               spec.phase.axis = list(at = NULL, labels = TRUE,
#                                      las = 1, hadj = NA, padj = NA),
#               main = "Eta Phase", sub = NULL,
#               lwd = 1, lwd.axis = 1,
#               verbose = TRUE)
# 
# 
# 

#####. compute cross wavelet transform, plot power and coherence ----
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
####  Average coss series power and coherence ----

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
   main = NULL, 
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
   main = NULL, 
   lwd = 1, col = 1, 
   lwd.axis = 1,
   verbose = TRUE)

## Cross series Coherency : Variable Periods----
# C <- analyze.coherency(W, my.pair = c("Tides", "Lagoon"), loess.span = 0,
#                   dt=dT,
#                   dj=dJ,
#                   upperPeriod = upperPeriod,
#                   lowerPeriod = lowerPeriod,
#                   window.type.t = 1, window.type.s = 1,       #Hamming Windows
#                   window.size.t = 10, window.size.s = 1/32,   # chamged from 1/4
#                   make.pval = FALSE, method = "white.noise", params = NULL,
#                   n.sim = Nsim,
#                   date.format = NULL, date.tz = NULL,
#                   verbose = TRUE)
# 
# 
# 
# wt.phase.image(C, my.series = 1,
#                plot.coi = TRUE,
#                plot.contour = TRUE,
#                siglvl = 0.1, col.contour = "white",
#                plot.ridge = TRUE, col.ridge = "black",
#                n.levels = 200,
#                color.palette = "rainbow(n.levels, start = 0, end = .7)",
#                useRaster = TRUE, max.contour.segments = 250000,
#                plot.legend = TRUE,
#                legend.params = list(width = 1.2, shrink = 0.9, mar = 5.1,
#                                     n.ticks = 6,
#                                     pi.style = TRUE,
#                                     label.digits = 1, label.format = "f",
#                                     lab = NULL, lab.line = 3),
#                label.time.axis = TRUE,
#                show.date = TRUE, date.format = "%F %T", date.tz = NULL,
#                timelab = paste0(StartDate,"  ",EndDate), 
#                                    timetck = 0.02, timetcl = 0.5,
#                spec.time.axis = list(at = NULL, labels = TRUE,
#                                      las = 1, hadj = NA, padj = NA),
#                label.period.axis = TRUE,
#                periodlab = NULL, periodtck = 0.02, periodtcl = 0.5,
#                spec.period.axis = list(at = NULL, labels = TRUE,
#                                        las = 1, hadj = NA, padj = NA),
#                main = paste0(WaterYear,"  Lagoon and Tide Coherence") ,
#                lwd = 2, lwd.axis = 1,
#                graphics.reset = TRUE,
#                verbose = TRUE)

 dev.off()

### plotly ----

plot_ly(z = ~WaveletTides$Power , 
        x=WaveletTides$series$date, 
        y=WaveletTides$Period,
        type="contour") %>% 
        layout(title = 'Tides ',
             xaxis = list(title = 'Date'),
             yaxis = list(title = 'Period (Hrs)',type='log' , tick=0.3010299957)) %>%
             colorbar( title =' Power ') 


plot_ly(z = ~WaveletLagoon$Power , 
        x=WaveletLagoon$series$date, 
        y=WaveletLagoon$Period,
        type="contour") %>% 
        layout(title = 'Lagoon ', 
             xaxis = list(title = 'Date'), 
             yaxis = list(title = 'Period (Hrs)',type='log',dtick=0.3010299957)) %>%
             colorbar( title =' Power ')

plot_ly(z = ~WaveletEta$Power , 
        x=WaveletEta$series$date, 
        y=WaveletEta$Period,
        type="contour") %>% 
        layout(title = 'Eta ', 
             xaxis = list(title = 'Date'), 
             yaxis = list(title = 'Period (Hrs)',type='log',dtick=0.3010299957)) %>%
             colorbar( title =' Power ')


#