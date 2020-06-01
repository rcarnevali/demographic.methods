
# Contains functions related to the Brass Growth Balance method

##########BGB Columns##########################################################

#' @title Estimate death registration coverage for a single year/sex/region using the Brass Growth Balance method
#'
#' @description Given one census and an average annual number of deaths in each age class, we can use stable population assumptions to estimate the completeness of the reporting of deaths relative to an estimate of the population. The method makes use of the observation that in a stable population that is closed to migration and has accurately reported data, the growth rate, r, is equal to the birth rate, b, less the death rate, d.
#' 
#' @details \code{bgbcolumns()}
#' 
#' @param Data \code{data.frame} with columns, \code{$Age}, \code{$Pop}, \code{$Deaths}.
#' @param censysyear Year of the census
#' @param minA The minimum of the age range searched. Default 5
#' @param maxA The maximum of the age range searched. Default 69
#' @param deaths.summed Logical. Is the deaths column given as the total per age in the intercensal period (\code{TRUE}). By default we assume \code{FALSE}.
#' @param deathsyear1 The first reference year for the deaths' distribution. Instead of transforming the exact reference date for the deaths' distribution, we used only the year and considered the middle point for this given year.
#' @param deathsyear2 The final reference year for the deaths' distribution. If the deaths refer to just one year, this parameter may be left blank 
#' 
#'  
#' @export
#' 
#' @example
#' ## El Salvador 1961 Females Census data:
#' 
#' elsalvador <- data.frame("Age" = seq(0, 75, 5), 
#'                          "population" = c(214089, 190234, 149538, 125040, 113490, 91663, 77711, 72936, 56942, 46205, 38616, 26154, 29273, 14964, 11205, 16193),
#'                          "death" = c(6909, 610, 214, 266, 291, 271, 315, 349, 338, 357, 385, 387, 647, 449, 504, 1360))
#' 
#' bgbcolumns(Data = elsalvador, censusyear = 1961, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1 = 1961, deathsyear2 = 1961)

 bgbcolumns <- function (Data = Data, censusyear = year, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1, deathsyear2){
   require(dplyr)
   
   if (abs(censusyear-deathsyear1) > 5) {
     
     warning("The census' reference year and the deaths' first reference year are more than 5 years apart!")
   }
   
   ## 1) Change the name of deaths column to "Deaths" to fit the code
   tempcols             <- tolower(colnames(Data))
   dcol                 <- grepl("deaths", tempcols)
   stopifnot(any(dcol))
   colnames(Data)[dcol] <- "Deaths"
   
   ## 2) Change the name of population column to "Pop" to fit the code
   tempcols             <- tolower(colnames(Data))
   pcol                 <- grepl("population", tempcols)
   stopifnot(any(pcol))
   colnames(Data)[pcol] <- "Pop"
   
   ## 3) If the deaths should be considered an average of the period or they refer to just a year.
   avgDeaths.alt <- function(Data, deaths.summed, deathsyear1, deathsyear2){
     if (!"deathsAvg" %in% colnames(Data)){
       if (deaths.summed){
         deathsmid      <- ((deathsyear1 + 0.5) + (deathsyear2 + 0.5)) / (deathsyear2 - deathsyear1)
         Data$deathsAvg <- Data$Deaths / deathmid
       } else {
         Data$Deaths <- Data$Deaths
       }
     }
     Data
   }
   
   Data <- avgDeaths.alt(Data = Data, deaths.summed = deaths.summed, deathsyear1 = deathsyear1, deathsyear2 = deathsyear2)
   
   
   ## 4) Detects the interval between deaths years
   detectDeathsPeriod <- function (Data, year1, year2) {
     if (year1 > year2){
       warning("The start year of deaths registration is larger than final year of deaths registration")
       return(NA)
     } else {
       if (year1 == year2) {
         deathsperiod <- 1
       } else {
         deathsperiod <- year2 - year1
       }
     }
     deathsperiod
   }
   
   deathsperiod <- detectDeathsPeriod (Data = Data, year1 = deathsyear1, year2 = deathsyear2)
   
   ## 5) Detects the midpoint of deaths
   deathsmid <- ifelse (deathsyear1 == deathsyear2,
                        deathsmid <- deathsyear1 + 0.5,
                        deathsmid <- ((deathsyear2 + deathsyear1) / 2) + 0.5)
   
   
   ## 6) This detect the age interval within the data
   detectAgeInterval <- function(Data, MinAge, MaxAge, ageColumn){
     
     stopifnot(tolower(ageColumn) %in% tolower(colnames(Data)))
     
     colnames(Data)[grepl(tolower(ageColumn),tolower(colnames(Data)))] <- "Age"
     
     Ages.int <- with(Data, unique(Age[Age >= MinAge & Age <= MaxAge]))
     Interval <- unique(diff(sort(Ages.int)))
     
     if (length(Interval) > 1){
       warning("You have more than one interval!")
       return(NA)
     }
     Interval
   }
   
   AgeInt <- detectAgeInterval(Data = Data, MinAge = minA, MaxAge = maxA, ageColumn = "Age")
   
   ## 7) If there is a 0-1 and 1-4 age groups, this part turns them into 0-4 group
   Ages  <- Data$Age
   
   if (1 %in% Ages){
     ind0           <- Ages == 0
     ind1           <- Ages == 1
     cnames         <- c("Pop","Deaths")
     
     if ("deathsAvg" %in% colnames(X)){
       cnames       <- c(cnames,"deathsAvg")
     }
     
     X[ind0,cnames] <- colSums(X[ ind0 | ind1, cnames])
     X              <- X[!ind1, ]
   }
   
   N                <- nrow(Data)
   Data$Pop         <- as.double(Data$Pop)
   Data$Deaths      <- as.double(Data$Deaths)
   Data$PopCum      <- rev ( cumsum( rev ( Data$Pop ) ) )
   Data$DeathCum    <- rev ( cumsum( rev ( Data$Deaths ) ) )
   Data$PYL         <- Data$PopCum * deathsperiod
   Data$Nx          <- (c ( 0, rep(sqrt(Data$Pop[ -N  ] * Data$Pop[ -1 ] ), length.out = N-2 ), 0) * deathsperiod ) / AgeInt
   Data$bx          <- Data$Nx / Data$PYL
   Data$X           <- ifelse(Data$Nx == 0, Data$X <- 0, Data$X <- Data$DeathCum/Data$PYL)
   Data$Y           <- ifelse(Data$Age >= minA,
                              ifelse((Data$Age - 1) < maxA, Data$Y <- Data$bx, Data$Y <- NA),
                              Data$Y <- NA)
   
   temp             <- filter(Data, !(Age < minA), !(Age >= maxA), !is.na(Y))
   b                <- as.double( sd(temp$Y) / sd(temp$X))
   r                <- as.double( mean(temp$Y) - mean(temp$X) * b)
   c                <- (1 / b) * (exp (r * (censusyear - deathsmid)))
   
   Data$abx         <- ifelse((Data$Age - 1) < maxA, Data$abx <- r + b * Data$X, Data$abx <- NA)
   Data <- Data %>%
     group_by(Age) %>%
     mutate(residual = ifelse(Age >= minA,
                              ifelse(Age <= maxA, residual <- (Y-abx), residual <- NA),
                              residual <- NA)) %>%
     ungroup()
   
   Data             <- naniar::replace_with_na_at(Data, .vars = c("Nx","bx", "Y", "abx", "residuals"), condition = ~.x == 0)
   
   Data
 }

##########BGB Coverage#########################################################

#' @details \code{bgbcoverage()}
#' 
#' @param Data a chunk of data (single sex, year, etc) with all columns required by \code{bgbcolumns()}
#' @param censusyear Year of the census
#' @param minA The minimum of the age range searched. Default 5
#' @param maxA The maximum of the age range searched. Default 69
#' @param deaths.summed Logical. is the deaths column given as the total per age in the intercensal period (\code{TRUE}). By default we assume \code{FALSE}.
#' @param deathsyear1 The first reference year for the deaths' distribution. Instead of transforming the exact reference date for the deaths' distribution, we used only the year and considered the middle point for this given year.
#' @param deathsyear2 The final reference year for the deaths' distribution. If the deaths refer to just one year, this parameter may be left blank 
#' 
#'@return A table with two informations: 1) The completeness relative to population at midpoint, and 2) The annual growth rate of stable population
#'
#'@export
#'
#'@example
#'El Salvador 1961 Females Census data:
#' 
#' elsalvador <- data.frame("Age" = seq(0, 75, 5), 
#'                          "population" = c(214089, 190234, 149538, 125040, 113490, 91663, 77711, 72936, 56942, 46205, 38616, 26154, 29273, 14964, 11205, 16193),
#'                          "death" = c(6909, 610, 214, 266, 291, 271, 315, 349, 338, 357, 385, 387, 647, 449, 504, 1360))
#' 
#' bgbcoverage(Data = elsalvador, censusyear = 1961, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1 = 1961, deathsyear2 = 1961)
#' 
 bgbcoverage <- function(Data = Data, censusyear = year, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1, deathsyear2) {
   require(dplyr)
   
   if (abs(censusyear-deathsyear1) > 5) {
     
     warning("The census' reference year and the deaths' first reference year are more than 5 years apart!")
   }
   
   ## 1) Change the name of deaths column to "Deaths" to fit the code
   tempcols             <- tolower(colnames(Data))
   dcol                 <- grepl("deaths", tempcols)
   stopifnot(any(dcol))
   colnames(Data)[dcol] <- "Deaths"
   
   ## 2) Change the name of population column to "Pop" to fit the code
   tempcols             <- tolower(colnames(Data))
   pcol                 <- grepl("population", tempcols)
   stopifnot(any(pcol))
   colnames(Data)[pcol] <- "Pop"
   
   ## 3) If the deaths should be considered an average of the period or they refer to just a year.
   
   avgDeaths.alt <- function(Data, deaths.summed, deathsyear1, deathsyear2){
     if (!"deathsAvg" %in% colnames(Data)){
       if (deaths.summed){
         deathsmid      <- ((deathsyear1 + 0.5) + (deathsyear2 + 0.5)) / (deathsyear2 - deathsyear1)
         Data$deathsAvg <- Data$Deaths / deathmid
       } else {
         Data$Deaths <- Data$Deaths
       }
     }
     Data
   }
   
   Data <- avgDeaths.alt(Data = Data, deaths.summed = deaths.summed, deathsyear1 = deathsyear1, deathsyear2 = deathsyear2)
   
   
   ## 4) Detects the interval between deaths years
   detectDeathsPeriod <- function (Data, year1, year2) {
     if (year1 > year2){
       warning("The start year of deaths registration is larger than final year of deaths registration")
       return(NA)
     } else {
       if (year1 == year2) {
         deathsperiod <- 1
       } else {
         deathsperiod <- year2 - year1
       }
     }
     deathsperiod
   }
   
   deathsperiod <- detectDeathsPeriod (Data = Data, year1 = deathsyear1, year2 = deathsyear2)
   
   ## 5) Detects the midpoint of deaths
   deathsmid <- ifelse (deathsyear1 == deathsyear2,
                        deathsmid <- deathsyear1 + 0.5,
                        deathsmid <- ((deathsyear2 + deathsyear1) / 2) + 0.5)
   
   
   ## 6) This detect the age interval within the data
   detectAgeInterval <- function(Data, MinAge, MaxAge, ageColumn){
     
     stopifnot(tolower(ageColumn) %in% tolower(colnames(Data)))
     
     colnames(Data)[grepl(tolower(ageColumn),tolower(colnames(Data)))] <- "Age"
     
     Ages.int <- with(Data, unique(Age[Age >= MinAge & Age <= MaxAge]))
     Interval <- unique(diff(sort(Ages.int)))
     
     if (length(Interval) > 1){
       warning("You have more than one interval!")
       return(NA)
     }
     Interval
   }
   
   AgeInt <- detectAgeInterval(Data = Data, MinAge = minA, MaxAge = maxA, ageColumn = "Age")
   
   ## 7) If there is a 0-1 and 1-4 age groups, this part turns them into 0-4 group
   Ages  <- Data$Age
   
   if (1 %in% Ages){
     ind0           <- Ages == 0
     ind1           <- Ages == 1
     cnames         <- c("Pop","Deaths")
     
     if ("deathsAvg" %in% colnames(X)){
       cnames       <- c(cnames,"deathsAvg")
     }
     
     X[ind0,cnames] <- colSums(X[ ind0 | ind1, cnames])
     X              <- X[!ind1, ]
   }
   
   N                <- nrow(Data)
   Data$Pop         <- as.double(Data$Pop)
   Data$Deaths      <- as.double(Data$Deaths)
   Data$PopCum      <- rev ( cumsum( rev ( Data$Pop ) ) )
   Data$DeathCum    <- rev ( cumsum( rev ( Data$Deaths ) ) )
   Data$PYL         <- Data$PopCum * deathsperiod
   Data$Nx          <- (c ( 0, rep(sqrt(Data$Pop[ -N  ] * Data$Pop[ -1 ] ), length.out = N-2 ), 0) * deathsperiod ) / AgeInt
   Data$bx          <- Data$Nx / Data$PYL
   Data$X           <- ifelse(Data$Nx == 0, Data$X <- 0, Data$X <- Data$DeathCum/Data$PYL)
   Data$Y           <- ifelse(Data$Age >= minA,
                              ifelse((Data$Age - 1) < maxA, Data$Y <- Data$bx, Data$Y <- NA),
                              Data$Y <- NA)
   
   temp             <- filter(Data, !(Age < minA), !(Age >= maxA), !is.na(Y))
   b                <- as.double( sd(temp$Y) / sd(temp$X))
   r                <- as.double( mean(temp$Y) - mean(temp$X) * b)
   c                <- (1 / b) * (exp (r * (censusyear - deathsmid)))
   
   result <- rbind("Completeness relative to population at midpoint (%)" = round(c, 4)*100,
                   "Annual growth rate of stable population (%)" = round(r,4)*100)
   
   result
 }
 
##########BGB Fitted Plot####################################################

#' @details \code{bgb.fitted.plot()}
#' 
#' @param Data a chunk of data (single sex, year, region, etc) with all columns required by \code{bgbcolumns()}
#' @param censysyear Year of the census
#' @param minA The minimum of the age range searched. Default 5
#' @param maxA The maximum of the age range searched. Default 69
#' @param deaths.summed Logical. is the deaths column given as the total per age in the intercensal period (\code{TRUE}). By default we assume \code{FALSE}, i.e. that the average annual was given.
#' @param deathsyear1 The first reference year for the deaths' distribution. Instead of transforming the exact reference date for the deaths' distribution, we used only the year and considered the middle point for this given year.
#' @param deathsyear2 The final reference year for the deaths' distribution. If the deaths refer to just one year, this parameter may be left blank 
#' @param countryname The country that the estimations are being made. A character parameter
#' 
#'@return A plot with the observed information as dots and the fitted information as lines estimated from the Brass Growth Balance Method
#'@export
#'
#'@example
#'El Salvador 1961 Females Census data:
#' 
#' elsalvador <- data.frame("Age" = seq(0, 75, 5), 
#'                          "population" = c(214089, 190234, 149538, 125040, 113490, 91663, 77711, 72936, 56942, 46205, 38616, 26154, 29273, 14964, 11205, 16193),
#'                          "death" = c(6909, 610, 214, 266, 291, 271, 315, 349, 338, 357, 385, 387, 647, 449, 504, 1360))
#'   
#' bgb.fitted.plot(Data = elsalvador, censusyear = 1961, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1 = 1961, deathsyear2 = 1961, countryname = 'El Salvador')

bgb.fitted.plot <- function(Data = Data, censusyear = year, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1, deathsyear2, countryname = country) {
  require(dplyr)
  require(ggplot2)
  
  if (abs(censusyear-deathsyear1) > 5) {
    
    warning("The census' reference year and the deaths' first reference year are more than 5 years apart!")
  }
  
  ## 1) Change the name of deaths column to "Deaths" to fit the code
  tempcols             <- tolower(colnames(Data))
  dcol                 <- grepl("deaths", tempcols)
  stopifnot(any(dcol))
  colnames(Data)[dcol] <- "Deaths"
  
  ## 2) Change the name of population column to "Pop" to fit the code
  tempcols             <- tolower(colnames(Data))
  pcol                 <- grepl("population", tempcols)
  stopifnot(any(pcol))
  colnames(Data)[pcol] <- "Pop"
  
  ## 3) If the deaths should be considered an average of the period or they refer to just a year.
  
  avgDeaths.alt <- function(Data, deaths.summed, deathsyear1, deathsyear2){
    if (!"deathsAvg" %in% colnames(Data)){
      if (deaths.summed){
        deathsmid      <- ((deathsyear1 + 0.5) + (deathsyear2 + 0.5)) / (deathsyear2 - deathsyear1)
        Data$deathsAvg <- Data$Deaths / deathmid
      } else {
        Data$Deaths <- Data$Deaths
      }
    }
    Data
  }
  
  Data <- avgDeaths.alt(Data = Data, deaths.summed = deaths.summed, deathsyear1 = deathsyear1, deathsyear2 = deathsyear2)
  
  
  ## 4) Detects the interval between deaths years
  detectDeathsPeriod <- function (Data, year1, year2) {
    if (year1 > year2){
      warning("The start year of deaths registration is larger than final year of deaths registration")
      return(NA)
    } else {
      if (year1 == year2) {
        deathsperiod <- 1
      } else {
        deathsperiod <- year2 - year1
      }
    }
    deathsperiod
  }
  
  deathsperiod <- detectDeathsPeriod (Data = Data, year1 = deathsyear1, year2 = deathsyear2)
  
  ## 5) Detects the midpoint of deaths
  deathsmid <- ifelse (deathsyear1 == deathsyear2,
                       deathsmid <- deathsyear1 + 0.5,
                       deathsmid <- ((deathsyear2 + deathsyear1) / 2) + 0.5)
  
  
  ## 6) This detect the age interval within the data
  detectAgeInterval <- function(Data, MinAge, MaxAge, ageColumn){
    
    stopifnot(tolower(ageColumn) %in% tolower(colnames(Data)))
    
    colnames(Data)[grepl(tolower(ageColumn),tolower(colnames(Data)))] <- "Age"
    
    Ages.int <- with(Data, unique(Age[Age >= MinAge & Age <= MaxAge]))
    Interval <- unique(diff(sort(Ages.int)))
    
    if (length(Interval) > 1){
      warning("You have more than one interval!")
      return(NA)
    }
    Interval
  }
  
  AgeInt <- detectAgeInterval(Data = Data, MinAge = minA, MaxAge = maxA, ageColumn = "Age")
  
  ## 7) If there is a 0-1 and 1-4 age groups, this part turns them into 0-4 group
  Ages  <- Data$Age
  
  if (1 %in% Ages){
    ind0           <- Ages == 0
    ind1           <- Ages == 1
    cnames         <- c("Pop","Deaths")
    
    if ("deathsAvg" %in% colnames(X)){
      cnames       <- c(cnames,"deathsAvg")
    }
    
    X[ind0,cnames] <- colSums(X[ ind0 | ind1, cnames])
    X              <- X[!ind1, ]
  }
  
  N                <- nrow(Data)
  Data$Pop         <- as.double(Data$Pop)
  Data$Deaths      <- as.double(Data$Deaths)
  
  ## Add if there is a ife expectation vector to use
  ## 8) This method demands to make a prior estimation of ex. Here this estimation is made through Brass General Growth Balance method
  Data$PopCum      <- rev ( cumsum( rev ( Data$Pop ) ) )
  Data$DeathCum    <- rev ( cumsum( rev ( Data$Deaths ) ) )
  Data$PYL         <- Data$PopCum * deathsperiod
  Data$Nx          <- (c ( 0, rep(sqrt(Data$Pop[ -N  ] * Data$Pop[ -1 ] ), length.out = N-2 ), 0) * deathsperiod ) / AgeInt
  Data$bx          <- Data$Nx / Data$PYL
  Data$X           <- ifelse(Data$Nx == 0, Data$X <- 0, Data$X <- Data$DeathCum/Data$PYL)
  Data$Y           <- ifelse(Data$Age >= minA,
                             ifelse((Data$Age - 1) < maxA, Data$Y <- Data$bx, Data$Y <- NA),
                             Data$Y <- NA)
  
  temp             <- filter(Data, !(Age < minA), !(Age >= maxA), !is.na(Y))
  b                <- as.double( sd(temp$Y) / sd(temp$X))
  r                <- as.double( mean(temp$Y) - mean(temp$X) * b)
  c                <- (1 / b) * (exp (r * (censusyear - deathsmid)))
  
  Data$abx         <- ifelse((Data$Age - 1) < maxA, Data$abx <- r + b * Data$X, Data$abx <- NA)
  Data <- Data %>%
    group_by(Age) %>%
    mutate(residual = ifelse(Age >= minA,
                             ifelse(Age <= maxA, residual <- (Y-abx), residual <- NA),
                             residual <- NA)) %>%
    ungroup()
  
  Data             <- naniar::replace_with_na_at(Data, .vars = c("Nx","bx", "Y", "abx", "residuals"), condition = ~.x == 0)

  title2 <- paste('(', countryname, ': ' , as.character(censusyear), ')', sep = " ")
  plot   <- Data %>%
    filter(!is.na(abx)) %>%
    ggplot(mapping = aes(x = X)) + 
    geom_path(aes(y = abx), color =  "#0072B2", size = 0.8) + 
    geom_point(aes(y = Y), color =  "#D55E00", size = 2.5) +                        
    labs(x = "Partial Death Rates - d(x+)",
         y = "Partial Birth Rates - b(x+)") +
    scale_color_manual(name = 'Estimation', labels = c("Fitted", "Observed"), values = c("#0072B2", "#D55E00")) +
    theme_bw() +
    ggtitle(paste("Brass Growth Balance Diagnostic Plot", title2, sep = "\n ")) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "right")
  
  plot
}

#########Residuals Plot########################################################

#' @details \code{bgb.residuals.plot()}
#' 
#' @param Data a chunk of data (single sex, year, region, etc) with all columns required by \code{bgbcolumns()}
#' @param censysyear Year of the census
#' @param minA The minimum of the age range searched. Default 5
#' @param maxA The maximum of the age range searched. Default 69
#' @param deaths.summed Logical. is the deaths column given as the total per age in the intercensal period (\code{TRUE}). By default we assume \code{FALSE}.
#' @param deathsyear1 The first reference year for the deaths' distribution. Instead of transforming the exact reference date for the deaths' distribution, we used only the year and considered the middle point for this given year.
#' @param deathsyear2 The final reference year for the deaths' distribution. If the deaths refer to just one year, this parameter may be left blank 
#' @param countryname The country that the estimations are being made. A character parameter
#' 
#'@return A plot with residuals estimated from the Brass Growth Balance Method
#'@export
#'
#'@example
#'El Salvador 1961 Females Census data:
#' 
#' elsalvador <- data.frame("Age" = seq(0, 75, 5), 
#'                          "population" = c(214089, 190234, 149538, 125040, 113490, 91663, 77711, 72936, 56942, 46205, 38616, 26154, 29273, 14964, 11205, 16193),
#'                          "death" = c(6909, 610, 214, 266, 291, 271, 315, 349, 338, 357, 385, 387, 647, 449, 504, 1360))
#'
#' bgb.residuals.plot(Data = elsalvador, censusyear = 1961, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1 = 1961, deathsyear2 = 1961, countryname = 'El Salvador')

bgb.residuals.plot <- function(Data = Data, censusyear = year, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1, deathsyear2, countryname = country) {
  require(dplyr)
  require(ggplot2)
  
  if (abs(censusyear-deathsyear1) > 5) {
    
    warning("The census' reference year and the deaths' first reference year are more than 5 years apart!")
  }
  
  ## 1) Change the name of deaths column to "Deaths" to fit the code
  tempcols             <- tolower(colnames(Data))
  dcol                 <- grepl("deaths", tempcols)
  stopifnot(any(dcol))
  colnames(Data)[dcol] <- "Deaths"
  
  ## 2) Change the name of population column to "Pop" to fit the code
  tempcols             <- tolower(colnames(Data))
  pcol                 <- grepl("population", tempcols)
  stopifnot(any(pcol))
  colnames(Data)[pcol] <- "Pop"
  
  ## 3) If the deaths should be considered an average of the period or they refer to just a year.
  
  avgDeaths.alt <- function(Data, deaths.summed, deathsyear1, deathsyear2){
    if (!"deathsAvg" %in% colnames(Data)){
      if (deaths.summed){
        deathsmid      <- ((deathsyear1 + 0.5) + (deathsyear2 + 0.5)) / (deathsyear2 - deathsyear1)
        Data$deathsAvg <- Data$Deaths / deathmid
      } else {
        Data$Deaths <- Data$Deaths
      }
    }
    Data
  }
  
  Data <- avgDeaths.alt(Data = Data, deaths.summed = deaths.summed, deathsyear1 = deathsyear1, deathsyear2 = deathsyear2)
  
  
  ## 4) Detects the interval between deaths years
  detectDeathsPeriod <- function (Data, year1, year2) {
    if (year1 > year2){
      warning("The start year of deaths registration is larger than final year of deaths registration")
      return(NA)
    } else {
      if (year1 == year2) {
        deathsperiod <- 1
      } else {
        deathsperiod <- year2 - year1
      }
    }
    deathsperiod
  }
  
  deathsperiod <- detectDeathsPeriod (Data = Data, year1 = deathsyear1, year2 = deathsyear2)
  
  ## 5) Detects the midpoint of deaths
  deathsmid <- ifelse (deathsyear1 == deathsyear2,
                       deathsmid <- deathsyear1 + 0.5,
                       deathsmid <- ((deathsyear2 + deathsyear1) / 2) + 0.5)
  
  
  ## 6) This detect the age interval within the data
  detectAgeInterval <- function(Data, MinAge, MaxAge, ageColumn){
    
    stopifnot(tolower(ageColumn) %in% tolower(colnames(Data)))
    
    colnames(Data)[grepl(tolower(ageColumn),tolower(colnames(Data)))] <- "Age"
    
    Ages.int <- with(Data, unique(Age[Age >= MinAge & Age <= MaxAge]))
    Interval <- unique(diff(sort(Ages.int)))
    
    if (length(Interval) > 1){
      warning("You have more than one interval!")
      return(NA)
    }
    Interval
  }
  
  AgeInt <- detectAgeInterval(Data = Data, MinAge = minA, MaxAge = maxA, ageColumn = "Age")
  
  ## 7) If there is a 0-1 and 1-4 age groups, this part turns them into 0-4 group
  Ages  <- Data$Age
  
  if (1 %in% Ages){
    ind0           <- Ages == 0
    ind1           <- Ages == 1
    cnames         <- c("Pop","Deaths")
    
    if ("deathsAvg" %in% colnames(X)){
      cnames       <- c(cnames,"deathsAvg")
    }
    
    X[ind0,cnames] <- colSums(X[ ind0 | ind1, cnames])
    X              <- X[!ind1, ]
  }
  
  N                <- nrow(Data)
  Data$Pop         <- as.double(Data$Pop)
  Data$Deaths      <- as.double(Data$Deaths)
  Data$PopCum      <- rev ( cumsum( rev ( Data$Pop ) ) )
  Data$DeathCum    <- rev ( cumsum( rev ( Data$Deaths ) ) )
  Data$PYL         <- Data$PopCum * deathsperiod
  Data$Nx          <- (c ( 0, rep(sqrt(Data$Pop[ -N  ] * Data$Pop[ -1 ] ), length.out = N-2 ), 0) * deathsperiod ) / AgeInt
  Data$bx          <- Data$Nx / Data$PYL
  Data$X           <- ifelse(Data$Nx == 0, Data$X <- 0, Data$X <- Data$DeathCum/Data$PYL)
  Data$Y           <- ifelse(Data$Age >= minA,
                             ifelse((Data$Age - 1) < maxA, Data$Y <- Data$bx, Data$Y <- NA),
                             Data$Y <- NA)
  
  temp             <- filter(Data, !(Age < minA), !(Age >= maxA), !is.na(Y))
  b                <- as.double( sd(temp$Y) / sd(temp$X))
  r                <- as.double( mean(temp$Y) - mean(temp$X) * b)
  c                <- (1 / b) * (exp (r * (censusyear - deathsmid)))
  
  Data$abx         <- ifelse((Data$Age - 1) < maxA, Data$abx <- r + b * Data$X, Data$abx <- NA)
  Data <- Data %>%
    group_by(Age) %>%
    mutate(residual = ifelse(Age >= minA,
                             ifelse(Age <= maxA, residual <- (Y-abx), residual <- NA),
                             residual <- NA)) %>%
    ungroup()
  
  Data             <- naniar::replace_with_na_at(Data, .vars = c("Nx","bx", "Y", "abx", "residuals"), condition = ~.x == 0)

  title1 <- paste('Residuals of the Census Completeness Estimation', 'Using Brass Growth Balance Method', sep ='\n')
  title2 <- paste('(', countryname, ': ' , as.character(censusyear), ')', sep = "")

  Plot   <- Data %>%
    filter(Age >= minA,
           !is.na(residual)) %>%
    ggplot(mapping = aes(x = X,
                         y = residual)) +
    geom_path(aes(), colour = "#0072B2", size = 0.5) +
    geom_point(aes(), colour = "#0072B2", size = 2)+
    labs(x="d(x+)",
         y = "residual") +
    theme_bw()+
    geom_hline(yintercept = 0, 
               show.legend = NA, 
               linetype = "dashed",
               color = "#D55E00", 
               size = 1) +
    ggtitle(paste(title1, title2, sep = "\n ")) +
    theme(plot.title = element_text(hjust = 0.5,face="bold"), legend.position = "right")
  
  Plot
}

###############################################################################
# Run all results
  #Sample data
  elsalvador <- data.frame("Age" = seq(0, 75, 5), 
                          "population" = c(214089, 190234, 149538, 125040, 113490, 91663, 77711, 72936, 56942, 46205, 38616, 26154, 29273, 14964, 11205, 16193),
                          "deaths" = c(6909, 610, 214, 266, 291, 271, 315, 349, 338, 357, 385, 387, 647, 449, 504, 1360))

    # Census Year: 1961
    # Deaths Year: 1961
 
  #Brass estimations   
  bgbcolumns (Data = elsalvador, censusyear = 1961, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1 = 1961, deathsyear2 = 1961)
  bgbcoverage (Data = elsalvador, censusyear = 1961, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1 = 1961, deathsyear2 = 1961)
  bgb.fitted.plot (Data = elsalvador, censusyear = 1961, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1 = 1961, deathsyear2 = 1961, countryname = 'El Salvador')
  bgb.residuals.plot (Data = elsalvador, censusyear = 1961, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1 = 1961, deathsyear2 = 1961, countryname = 'El Salvador')
