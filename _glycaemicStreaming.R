library(data.table)

returnUnixDateTime<-function(date) {
  returnVal<-as.numeric(as.POSIXct(date, format="%Y-%m-%d", tz="GMT"))
  return(returnVal)
}

# diagnosisDataset<-read.csv("../GlCoSy/SDsource/diagnosisDateDeathDate.txt")
diagnosisDataset<-read.csv("~/R/GlCoSy/SDsource/demogALL.txt", quote = "", 
                           row.names = NULL, 
                           stringsAsFactors = FALSE)

diagnosisDatasetDT <- data.table(diagnosisDataset)

## hba1c import
## load in hba1c data
a1c_filename <- paste("~/R/GlCoSy/SD_workingSource/hba1cDTclean2.csv",sep="")
a1c_data<-read.csv(a1c_filename)
a1c_dataDT<-data.table(a1c_data)
a1c_dataDT<-a1c_dataDT[order(a1c_dataDT$LinkId,a1c_dataDT$dateplustime1),]

###

# DM subset for extracting hba1c data
T1DM_sub <- diagnosisDatasetDT[DiabetesMellitusType_Mapped == "Type 1 Diabetes Mellitus"]
T2DM_sub <- diagnosisDatasetDT[DiabetesMellitusType_Mapped == "Type 2 Diabetes Mellitus"]

T1DM_sub <- T2DM_sub

T1DM_sub_test <- data.frame(T1DM_sub$LinkId, T1DM_sub$DateOfDiagnosisDiabetes_Date)
colnames(T1DM_sub_test) <- c("LinkId", "diagnosisDate")
T1DM_sub_test$diagnosisDateUnix <- returnUnixDateTime(T1DM_sub_test$diagnosisDate)

T1DM_sub_test_a1c_merge <- merge(T1DM_sub_test, a1c_dataDT, by.x="LinkId", by.y="LinkId")
T1DM_sub_test_a1c_merge$interval_months <- (T1DM_sub_test_a1c_merge$dateplustime1 - T1DM_sub_test_a1c_merge$diagnosisDateUnix) / (60*60*24*(365.25/12))
T1DM_sub_test_a1c_merge <- data.table(T1DM_sub_test_a1c_merge)

findHbA1cValues <- function(dateplustime1, newNumeric, diagnosisDateUnix, timePoint_1, timePoint_2, window_1, window_2) {
  
  # dateplustime1 <- T1DM_sub_test_a1c_merge[LinkId=="2147484450"]$dateplustime1; newNumeric = T1DM_sub_test_a1c_merge[LinkId=="2147484450"]$newNumeric; diagnosisDateUnix <- T1DM_sub_test_a1c_merge[LinkId=="2147484450"]$diagnosisDateUnix;timePoint_1 = 12; timePoint_2 = 60; window_1 = 4; window_2 = 12
  
  diagnosisDateUnix <- diagnosisDateUnix[1]
  tableALL <- data.table(dateplustime1, newNumeric)
  
  timePoint_1_unix <- diagnosisDateUnix + (timePoint_1 * (60*60*24*(365.25/12)))
  timePoint_2_unix <- diagnosisDateUnix + (timePoint_2 * (60*60*24*(365.25/12)))
  
      value_1_sub <- tableALL[dateplustime1 > (timePoint_1_unix - (window_1 * (60*60*24*(365.25/12)))) & dateplustime1 < (timePoint_1_unix + (window_1 * (60*60*24*(365.25/12))))]
      
      if (nrow(value_1_sub)>0) {
        value_1_sub$flag <- ifelse(
          (sqrt((value_1_sub$dateplustime1 - timePoint_1_unix)^2)) == min((sqrt((value_1_sub$dateplustime1 - timePoint_1_unix)^2))), 1, 0
        )
        value_1 = value_1_sub[flag==1]$newNumeric
      }
      
      if (nrow(value_1_sub)==0) { value_1 = 0 }
      
      value_2_sub <- tableALL[dateplustime1 > (timePoint_2_unix - (window_2 * (60*60*24*(365.25/12)))) & dateplustime1 < (timePoint_2_unix + (window_2 * (60*60*24*(365.25/12))))]
      
      if (nrow(value_2_sub)>0) {
        value_2_sub$flag <- ifelse(
          (sqrt((value_2_sub$dateplustime1 - timePoint_2_unix)^2)) == min((sqrt((value_2_sub$dateplustime1 - timePoint_2_unix)^2))), 1, 0
        )
        value_2 = value_2_sub[flag==1]$newNumeric
      }
      
      if (nrow(value_2_sub)==0) { value_2 = 0 }
      
  return(list(value_1, value_2))
}

timePoint_1 = 24 # months
timePoint_2 = 60 # months
window_1 = 4     # window (half)
window_2 = 6     # window (half)

T1DM_sub_test_a1c_merge[, c("a1c_point1", "a1c_point2") := findHbA1cValues(dateplustime1, newNumeric, diagnosisDateUnix, timePoint_1, timePoint_2, window_1, window_2) , by=.(LinkId)]

T1DM_sub_test_a1c_merge <- T1DM_sub_test_a1c_merge[order(T1DM_sub_test_a1c_merge$LinkId),]
T1DM_sub_test_a1c_merge$idDiff <- 1
T1DM_sub_test_a1c_merge$idDiff[2:nrow(T1DM_sub_test_a1c_merge)] <- diff(T1DM_sub_test_a1c_merge$LinkId)

# bring down to 1 row per ID
T1DM_sub_test_a1c_merge_plot <- T1DM_sub_test_a1c_merge[idDiff==1]

# sub with values at both time points:
a1c_2point_plot <- T1DM_sub_test_a1c_merge_plot[a1c_point1>0 & a1c_point2>0]

# 
plot(a1c_2point_plot$a1c_point1, a1c_2point_plot$a1c_point2)
  fit <- lm(a1c_2point_plot$a1c_point2 ~ a1c_2point_plot$a1c_point1)
  abline(fit, col="red")
boxplot(a1c_2point_plot$a1c_point2 ~ cut(a1c_2point_plot$a1c_point1, breaks=seq(0,200,5)), varwidth=T, xlab=paste("a1c at ",timePoint_1," months", sep=""), ylab=paste("a1c at ",timePoint_2," months", sep=""), main=paste("glycemic streaming. T1DM. n=",nrow(a1c_2point_plot), sep=""))












timepoints <- c(0,3,6,12,24,36,60,120)
window1 <- 4
window2 <- 6

reportingFrame <- as.data.frame(matrix(0,nrow=1,ncol=(length(timepoints)+0)))
colnames(reportingFrame) <- as.character(timepoints)
reportingFrame$LinkId = 0

findHba1c <- function(linkid, diagnosisDateUnix) {
  
  a1c_sub <- a1c_dataDT[LinkId == linkid]
  
  for (jj in timepoints) {
    window <- ifelse(jj < 24, window1, window2)
    windowSec <- window*(60*60*24*(365.25/12))
    
    timePoint <- diagnosisDateUnix + (jj * (60*60*24*(365.25/12)))
    
    hba1c_window_set <- a1c_sub[dateplustime1 > (timePoint - windowSec) & dateplustime1 < (timePoint + windowSec)]
    
    if (nrow(hba1c_window_set)>0) {
      hba1c_window_set$flag <- ifelse(
      (sqrt((hba1c_window_set$dateplustime1 - timePoint)^2)) == min((sqrt((hba1c_window_set$dateplustime1 - timePoint)^2))), 1, 0
    )
      returnVal = hba1c_window_set[flag==1]$newNumeric
    }
    
    if (nrow(hba1c_window_set)==0) {returnVal=0}

    
    
  }
}