# Contains functions related to the Preston and Coale method

## Refer to: http://demographicestimation.iussp.org/content/preston-coale-method

#############Preston and Coale Columns##########################

#' @title Estimate death registration coverage for a single year/sex/region using the Preston and Coale method
#'
#' @description Given one census and an average annual number of deaths in each age class, we can use stable population and closed to migration assumptions to estimate the completeness of the reporting of deaths relative to an estimate of the population. This is a particular case of the more general Synthetic Extinct Generations method, which requires estimates of the population at two points in time but does not require that the population is stable.
#'  
#' @details \code{prestoncoalecolumns()}
#' 
#' @param Data \code{data.frame} with columns, \code{$Age}, \code{$Pop}, \code{$Deaths}.
#' @param censysyear Year of the census
#' @param lowerage The minimum of the age range searched. Default 15
#' @param upperage The maximum of the age range searched. Default 64
#' @param deaths.summed Logical. Is the deaths column given as the total per age in the intercensal period (\code{TRUE}). By default we assume \code{FALSE}.
#' @param deathsyear1 The first reference year for the deaths' distribution. Instead of transforming the exact reference date for the deaths' distribution, we used only the year and considered the middle point for this given year.
#' @param deathsyear2 The final reference year for the deaths' distribution. If the deaths refer to just one year, this parameter may be left blank 
#' @param minA The minimum of the age range used for estimating the ex from the Brass Growth Balance. Default 5
#' @param maxA The maximum of the age range used for estimating the ex from the Brass Growth Balance. Default 74
#' @param bgbmaxage The max age to estimate the life expectancy with Brass Growth Balance Method
#' @param model.LT \code{data.frame} in the wide format with the relational logit Model Life Table values
#' @param sex Information on Model Life Table sex 
#' @param family Model Life Table Families such as the UN's 'General' family,	The Priceton's 'East',	'North',	'South',	'West' families, and	'AIDS'. Defaut 'West'
#' @param model.agefrom Lower Bound Age for the Model Life Table fitting. Default 45
#' @param model.ageto Upper Bound Age for the Model Life Table fitting. Default 75
#'  
#' @export
#' 
#' @example
#' ## El Salvador 1961 Females Census data:
#' 
#' prestoncoalecolumns(Data = elsalvador, censusyear = 1961, lowerage = 15, upperage = 60,  deaths.summed = FALSE, deathsyear1 = 1961, deathsyear2 = 1961,
#'                      minA = 5, maxA = 74, bgbmaxage = 69, model.LT = Model.LT.logits, model.sex = 'female', model.family = 'West', model.agefrom = 45, model.ageto = 75)

prestoncoalecolumns <- function (Data = Data, censusyear, lowerage = 15, upperage = 60,  deaths.summed = FALSE, 
                                 deathsyear1, deathsyear2, minA = 5, maxA = 74, bgbmaxage = 69,
                                 model.LT = Model.LT.logits, model.sex = sex, model.family = 'West', model.agefrom = 45, model.ageto = 75) {
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
  
  ## Add if there is a ife expectation vector to use
  ## 8) This method demands to make a prior estimation of ex. Here this estimation is made through Brass General Growth Balance method
  Data$PopCum      <- rev ( cumsum( rev ( Data$Pop ) ) )
  Data$DeathCum    <- rev ( cumsum( rev ( Data$Deaths ) ) )
  Data$PYL         <- Data$PopCum * deathsperiod
  Data$Nx          <- (c ( 0, rep(sqrt(Data$Pop[ -N  ] * Data$Pop[ -1 ] ), length.out = N-2 ), 0) * deathsperiod ) / AgeInt
  Data$bx          <- Data$Nx / Data$PYL
  Data$X           <- ifelse(Data$Nx == 0, Data$X <- 0, Data$X <- Data$DeathCum/Data$PYL)
  Data$Y           <- ifelse(Data$Age >= minA,
                             ifelse((Data$Age - 1) < bgbmaxage, Data$Y <- Data$bx, Data$Y <- NA),
                             Data$Y <- NA)
  
  temp             <- subset(Data, !(Age < minA) & !(Age >= maxA) & !is.na(Y))
  b                <- as.double( sd(temp$Y) / sd(temp$X))
  r                <- as.double( mean(temp$Y) - mean(temp$X) * b)
  c                <- (1 / b) * (exp (r * (censusyear - deathsmid)))
  
  Data$abx         <- ifelse((Data$Age - 1) < bgbmaxage, Data$abx <- r + b * Data$X, Data$abx <- NA)
  age_by           <- by(Data, 
                         INDICES = Data$Age, 
                         FUN = function(x){
                           x$residual <- ifelse(x$Age >= minA,
                                                ifelse(x$Age <= maxA, x$residual <- (x$Y - x$abx), x$residual <- NA),
                                                x$residual <- NA)
                           return(x)
                         })
  
  Data             <- do.call(rbind,age_by)
  Data             <- naniar::replace_with_na_at(Data, .vars = c("Nx","bx", "Y", "abx", "residual"), condition = ~.x == 0)
  Data$Pop.Ajust   <- ifelse(Data$Age >= minA, Data$Pop * exp( -r * (censusyear - deathsmid)), NA)
  Data$Death.Ajust <- ifelse(Data$Age >= minA, Data$Deaths / c, NA)
  Data$PYL.Ajust   <- ifelse(Data$Age >= minA, Data$Pop.Ajust * deathsperiod, NA)
  Data$mx.Ajust    <- ifelse(Data$Age >= minA, Data$Death.Ajust / Data$PYL.Ajust, NA)
  Data$qx          <- ifelse(Data$Age <= maxA+5,
                             ifelse(Data$Age >= minA, AgeInt * Data$mx.Ajust / ( 1 + 2.5 * Data$mx.Ajust), NA),
                             NA)
  
  Data$lx <- NA
  for(i in 1:nrow(Data)){
    ifelse(Data$Age[i] < minA,
           Data$lx[i] <- NA,
           ifelse(Data$Age[i] == minA,
                  Data$lx[i] <- 1.0, 
                  Data$lx[i] <- Data$lx[i-1] * (1-Data$qx[i-1])))
  }
  
  Data$Obs.Y       <- ifelse(!is.na(Data$lx),
                             ifelse(Data$lx != 0, 0.5 * (log((1 - Data$lx) / Data$lx)), NA),
                             NA)
  Data             <- mutate_if(Data, is.numeric, list(~na_if(., -Inf)))
  
  LT               <- gather(Model.LT.logits, fam, Ys, -c(sex, e0, Age))
  
  temp             <- subset(LT, sex == model.sex & fam == model.family)
  temp$Std         <- 1/ (1 + (exp(2 * temp$Ys)))
  temp$lx.Std      <- temp$Std / temp$Std [1]
  
  x                <- data.frame('Age' = seq(5, 100, by = 5))
  temp             <- merge(temp, x, by = "Age", all.y = TRUE)
  
  for(i in 1:nrow(temp)){
    if (temp$Age[i] > 85){
      temp$lx.Std[i] <- temp$lx.Std[i-1] * exp((log(temp$lx.Std[i-1] / temp$lx.Std[i-2])^2) / log(temp$lx.Std[i-2] / temp$lx.Std[i-3]))
    } else {
      next
    }
  }
  
  temp$Cdn.Ys      <- 0.5 * log ((temp$lx.Std [1] - temp$lx.Std) / temp$lx.Std)
  temp             <- mutate_if (temp, is.numeric, list(~na_if(., -Inf)))
  
  Data             <- full_join (Data, temp , by = "Age")
  
  model            <- subset(Data, Age >= model.agefrom & Age <= model.ageto)
  model            <- lm (formula = Obs.Y ~ Cdn.Ys, data = model)
  alp              <- as.double(model$coefficients [1])
  bet              <- as.double(model$coefficients [2])
  
  Data$Fit.Ys      <- alp + bet * Data$Cdn.Ys
  Data$Fit.lx      <- ifelse(Data$Age == minA, Data$Fit.lx <- 1, Data$Fit.lx <- 1 / (1 + exp(2 * Data$Fit.Ys)))
  
  for(i in nrow(Data):1){
    ifelse(Data$Age[i] == max(Data$Age),
           Data$Tx[i] <- 2.5 * Data$Fit.lx[i],
           ifelse(Data$Age[i] >= minA,
                  Data$Tx[i] <- Data$Tx[i + 1] + 2.5*(Data$Fit.lx[i + 1] + Data$Fit.lx[i]),
                  Data$Tx[i] <- NA))
  }
  
  Data$ex          <- ifelse(Data$Age >= minA, Data$Tx/Data$Fit.lx, Data$Tx <- NA)
  
  # Estimating Est.Nx 
  for(i in nrow(Data):1){
    
    if (Data$Age[i] < minA){
      next
    }
    
    ifelse (Data$Age[i] > 64, 
            ifelse(maxA < Data$Age[i],
                   Data$Est.Nx[i] <- (exp(r * Data$ex[i]) - (r * Data$ex[i] ) ^ 2 / 6) * Data$DeathCum[i],
                   Data$Est.Nx[i] <- Data$Est.Nx[i + 1] * exp(5 * r) + Data$Deaths[i] * exp(2.5 * r)),
            Data$Est.Nx[i] <- Data$Est.Nx[i + 1] * exp(5 * r) + Data$Deaths[i] * exp(2.5 * r))
  }
  
  # Estimating Est.5Nx
  for(i in nrow(Data):1){
    
    if (Data$Age[i] < minA){
      next
    }
    
    ifelse (Data$Age[i] < maxA,
            Data$Est.5Nx[i] <- round(2.5 * (Data$Est.Nx[i] + Data$Est.Nx[i+1]),1),
            Data$Est.5Nx[i] <- NA)
  }
  
  Data$Obs.5Nx   <- ifelse(maxA > Data$Age, 
                           Data$Pop*deathsperiod, 
                           NA)
  Data$cNx       <- ifelse(maxA > Data$Age,
                           Data$Est.5Nx/Data$Obs.5Nx,
                           NA)
  
  # Estimating a.xNx
  Data$a.xNx  <-  NA
  for(i in 1:nrow(Data)){
    if (Data$Age[i] < minA){
      next
    }
    ifelse (Data$Age[i] < maxA,
            Data$a.xNx[i] <- sum(Data$Est.5Nx[i:N],na.rm = T)/sum(Data$Obs.5Nx[i:N], na.rm = T),
            Data$a.xNx[i] <- NA)
  }  
  
  Method <- Data[(Data$Age <= model.ageto), (names(Data) %in% c("Age","Pop","Deaths", "Est.Nx", "Est.5Nx", "Obs.5Nx", "cNx", "a.xNx"))]
  
  Method
}

#############Preston and Coale Coverage####################

#' @details \code{prestoncoalecov()}
#' 
#' @param Data \code{data.frame} with columns, \code{$Age}, \code{$Pop}, \code{$Deaths}.
#' @param censysyear Year of the census
#' @param lowerage The minimum of the age range searched. Default 15
#' @param upperage The maximum of the age range searched. Default 64
#' @param deaths.summed Logical. Is the deaths column given as the total per age in the intercensal period (\code{TRUE}). By default we assume \code{FALSE}.
#' @param deathsyear1 The first reference year for the deaths' distribution. Instead of transforming the exact reference date for the deaths' distribution, we used only the year and considered the middle point for this given year.
#' @param deathsyear2 The final reference year for the deaths' distribution. If the deaths refer to just one year, this parameter may be left blank 
#' @param minA The minimum of the age range used for estimating the ex from the Brass Growth Balance. Default 5
#' @param maxA The maximum of the age range used for estimating the ex from the Brass Growth Balance. Default 74
#' @param bgbmaxage The max age to estimate the life expectancy with Brass Growth Balance Method
#' @param model.LT \code{data.frame} in the wide format with the relational logit Model Life Table values
#' @param sex Information on Model Life Table sex 
#' @param family Model Life Table Families such as the UN's 'General' family,	The Priceton's 'East',	'North',	'South',	'West' families, and	'AIDS'. Defaut 'West'
#' @param model.agefrom Lower Bound Age for the Model Life Table fitting. Default 45
#' @param model.ageto Upper Bound Age for the Model Life Table fitting. Default 75
#'  
#' @export
#' 
#' @example
#' ## El Salvador 1961 Females Census data:
#' 
#' prestoncoalecov (Data = elsalvador, censusyear = 1961, lowerage = 15, upperage = 60,  deaths.summed = FALSE, deathsyear1 = 1961, deathsyear2 = 1961,
#'                      minA = 5, maxA = 74, bgbmaxage = 69, model.LT = Model.LT.logits, model.sex = 'female', model.family = 'West', model.agefrom = 45, model.ageto = 75)

prestoncoalecov <- function (Data = Data, censusyear, lowerage = 15, upperage = 60,  deaths.summed = FALSE, 
                             deathsyear1, deathsyear2, minA = 5, maxA = 74, bgbmaxage = 69,
                             model.LT = Model.LT.logits, model.sex = sex, model.family = 'West', model.agefrom = 45, model.ageto = 75) {
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
  
  ## Add if there is a ife expectation vector to use
  ## 8) This method demands to make a prior estimation of ex. Here this estimation is made through Brass General Growth Balance method
  Data$PopCum      <- rev ( cumsum( rev ( Data$Pop ) ) )
  Data$DeathCum    <- rev ( cumsum( rev ( Data$Deaths ) ) )
  Data$PYL         <- Data$PopCum * deathsperiod
  Data$Nx          <- (c ( 0, rep(sqrt(Data$Pop[ -N  ] * Data$Pop[ -1 ] ), length.out = N-2 ), 0) * deathsperiod ) / AgeInt
  Data$bx          <- Data$Nx / Data$PYL
  Data$X           <- ifelse(Data$Nx == 0, Data$X <- 0, Data$X <- Data$DeathCum/Data$PYL)
  Data$Y           <- ifelse(Data$Age >= minA,
                             ifelse((Data$Age - 1) < bgbmaxage, Data$Y <- Data$bx, Data$Y <- NA),
                             Data$Y <- NA)
  
  temp             <- subset(Data, !(Age < minA) & !(Age >= maxA) & !is.na(Y))
  b                <- as.double( sd(temp$Y) / sd(temp$X))
  r                <- as.double( mean(temp$Y) - mean(temp$X) * b)
  c                <- (1 / b) * (exp (r * (censusyear - deathsmid)))
  
  Data$abx         <- ifelse((Data$Age - 1) < bgbmaxage, Data$abx <- r + b * Data$X, Data$abx <- NA)
  age_by           <- by(Data, 
                         INDICES = Data$Age, 
                         FUN = function(x){
                           x$residual <- ifelse(x$Age >= minA,
                                                ifelse(x$Age <= maxA, x$residual <- (x$Y - x$abx), x$residual <- NA),
                                                x$residual <- NA)
                           return(x)
                         })
  
  Data             <- do.call(rbind,age_by)
  Data             <- naniar::replace_with_na_at(Data, .vars = c("Nx","bx", "Y", "abx", "residual"), condition = ~.x == 0)
  Data$Pop.Ajust   <- ifelse(Data$Age >= minA, Data$Pop * exp( -r * (censusyear - deathsmid)), NA)
  Data$Death.Ajust <- ifelse(Data$Age >= minA, Data$Deaths / c, NA)
  Data$PYL.Ajust   <- ifelse(Data$Age >= minA, Data$Pop.Ajust * deathsperiod, NA)
  Data$mx.Ajust    <- ifelse(Data$Age >= minA, Data$Death.Ajust / Data$PYL.Ajust, NA)
  Data$qx          <- ifelse(Data$Age <= maxA+5,
                             ifelse(Data$Age >= minA, AgeInt * Data$mx.Ajust / ( 1 + 2.5 * Data$mx.Ajust), NA),
                             NA)
  Data$lx <- NA
  for(i in 1:nrow(Data)){
    ifelse(Data$Age[i] < minA,
           Data$lx[i] <- NA,
           ifelse(Data$Age[i] == minA,
                  Data$lx[i] <- 1.0, 
                  Data$lx[i] <- Data$lx[i-1] * (1-Data$qx[i-1])))
  }
  
  Data             <- as.data.frame(Data)
  
  Data$Obs.Y       <- ifelse(!is.na(Data$lx),
                             ifelse(Data$lx != 0, 0.5 * (log((1 - Data$lx) / Data$lx)), NA),
                             NA)
  Data             <- mutate_if(Data, is.numeric, list(~na_if(., -Inf)))
  
  LT               <- gather(Model.LT.logits, fam, Ys, -c(sex, e0, Age))
  
  temp             <- subset(LT, sex == model.sex & fam == model.family)
  temp$Std         <- 1/ (1 + (exp(2 * temp$Ys)))
  temp$lx.Std      <- temp$Std / temp$Std [1]
  
  x                <- data.frame('Age' = seq(5, 100, by = 5))
  temp             <- merge(temp, x, by = "Age", all.y = TRUE)
  
  for(i in 1:nrow(temp)){
    if (temp$Age[i] > 85){
      temp$lx.Std[i] <- temp$lx.Std[i-1] * exp((log(temp$lx.Std[i-1] / temp$lx.Std[i-2])^2) / log(temp$lx.Std[i-2] / temp$lx.Std[i-3]))
    } else {
      next
    }
  }
  
  temp$Cdn.Ys      <- 0.5 * log ((temp$lx.Std [1] - temp$lx.Std) / temp$lx.Std)
  temp             <- mutate_if (temp, is.numeric, list(~na_if(., -Inf)))
  
  Data             <- full_join (Data, temp , by = "Age")
  
  model            <- subset(Data, Age >= model.agefrom & Age <= model.ageto)
  model            <- lm (formula = Obs.Y ~ Cdn.Ys, data = model)
  alp              <- as.double(model$coefficients [1])
  bet              <- as.double(model$coefficients [2])
  
  Data$Fit.Ys      <- alp + bet * Data$Cdn.Ys
  Data$Fit.lx      <- ifelse(Data$Age == minA, Data$Fit.lx <- 1, Data$Fit.lx <- 1 / (1 + exp(2 * Data$Fit.Ys)))
  
  for(i in nrow(Data):1){
    ifelse(Data$Age[i] == max(Data$Age),
           Data$Tx[i] <- 2.5 * Data$Fit.lx[i],
           ifelse(Data$Age[i] >= minA,
                  Data$Tx[i] <- Data$Tx[i + 1] + 2.5*(Data$Fit.lx[i + 1] + Data$Fit.lx[i]),
                  Data$Tx[i] <- NA))
  }
  
  Data$ex          <- ifelse(Data$Age >= minA, Data$Tx/Data$Fit.lx, Data$Tx <- NA)
  
  # Estimating Est.Nx 
  for(i in nrow(Data):1){
    
    if (Data$Age[i] < minA){
      next
    }
    
    ifelse (Data$Age[i] > 64, 
            ifelse(maxA < Data$Age[i],
                   Data$Est.Nx[i] <- (exp(r * Data$ex[i]) - (r * Data$ex[i] ) ^ 2 / 6) * Data$DeathCum[i],
                   Data$Est.Nx[i] <- Data$Est.Nx[i + 1] * exp(5 * r) + Data$Deaths[i] * exp(2.5 * r)),
            Data$Est.Nx[i] <- Data$Est.Nx[i + 1] * exp(5 * r) + Data$Deaths[i] * exp(2.5 * r))
  }
  
  # Estimating Est.5Nx
  for(i in nrow(Data):1){
    
    if (Data$Age[i] < minA){
      next
    }
    
    ifelse (Data$Age[i] < maxA,
            Data$Est.5Nx[i] <- round(2.5 * (Data$Est.Nx[i] + Data$Est.Nx[i+1]),1),
            Data$Est.5Nx[i] <- NA)
  }
  
  Data$Obs.5Nx   <- ifelse(maxA > Data$Age, 
                           Data$Pop*deathsperiod, 
                           NA)
  Data$cNx       <- ifelse(maxA > Data$Age,
                           Data$Est.5Nx/Data$Obs.5Nx,
                           NA)
  
  # Estimating a.xNx
  Data$a.xNx  <-  NA
  for(i in 1:nrow(Data)){
    if (Data$Age[i] < minA){
      next
    }
    ifelse (Data$Age[i] < maxA,
            Data$a.xNx[i] <- sum(Data$Est.5Nx[i:N],na.rm = T)/sum(Data$Obs.5Nx[i:N], na.rm = T),
            Data$a.xNx[i] <- NA)
  }  
  
  temp        <- subset(Method, Method$Age >= lowerage & Method$Age <= upperage) 
  comp        <- ((median(temp$cNx)*0.5) + as.double((quantile(temp$cNx, 0.25) + quantile(temp$cNx, 0.75)))*0.25) *
    exp(r * (deathsmid - censusyear))
  comp.lower  <- as.double(comp - (quantile(temp$cNx, 0.75) - quantile(temp$cNx, 0.25)))
  comp.upper  <- as.double(comp + (quantile(temp$cNx, 0.75) - quantile(temp$cNx, 0.25)))
  
  result <- rbind("Lower 25% (%)" = round(comp.lower,4)*100,
                  "Completeness (%)" = round(comp, 4)*100,
                  "Upper 25% (%)" = round(comp.upper, 4)*100)
  result
}


#############Ratios Plot###########################################

#' @details \code{prestoncoale.ratiosplot()}
#' 
#' @param Data \code{data.frame} with columns, \code{$Age}, \code{$Pop}, \code{$Deaths}.
#' @param censysyear Year of the census
#' @param lowerage The minimum of the age range searched. Default 15
#' @param upperage The maximum of the age range searched. Default 64
#' @param deaths.summed Logical. Is the deaths column given as the total per age in the intercensal period (\code{TRUE}). By default we assume \code{FALSE}.
#' @param deathsyear1 The first reference year for the deaths' distribution. Instead of transforming the exact reference date for the deaths' distribution, we used only the year and considered the middle point for this given year.
#' @param deathsyear2 The final reference year for the deaths' distribution. If the deaths refer to just one year, this parameter may be left blank 
#' @param minA The minimum of the age range used for estimating the ex from the Brass Growth Balance. Default 5
#' @param maxA The maximum of the age range used for estimating the ex from the Brass Growth Balance. Default 74
#' @param countryname The country that the estimations are being made. A character parameter
#' @param bgbmaxage The max age to estimate the life expectancy with Brass Growth Balance Method
#' @param model.LT \code{data.frame} in the wide format with the relational logit Model Life Table values
#' @param sex Information on Model Life Table sex 
#' @param family Model Life Table Families such as the UN's 'General' family,	The Priceton's 'East',	'North',	'South',	'West' families, and	'AIDS'. Defaut 'West'
#' @param model.agefrom Lower Bound Age for the Model Life Table fitting. Default 45
#' @param model.ageto Upper Bound Age for the Model Life Table fitting. Default 75
#'  
#' @export
#' 
#' @example
#' ## El Salvador 1961 Females Census data:
#' 
#' prestoncoale.ratiosplot (Data = elsalvador, censusyear = 1961, lowerage = 15, upperage = 60,  deaths.summed = FALSE, deathsyear1 = 1961, deathsyear2 = 1961,
#'                      minA = 5, maxA = 74, bgbmaxage = 69, countryname ='El Salvador', model.LT = Model.LT.logits, model.sex = 'female', model.family = 'West', model.agefrom = 45, model.ageto = 75)

prestoncoale.ratiosplot <- function (Data = Data, censusyear, lowerage = 15, upperage = 60,  deaths.summed = FALSE, 
                                     deathsyear1, deathsyear2, minA = 5, maxA = 74, bgbmaxage = 69, countryname = country,
                                     model.LT = Model.LT.logits, model.sex = sex, model.family = 'West', model.agefrom = 45, model.ageto = 75) {
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
                             ifelse((Data$Age - 1) < bgbmaxage, Data$Y <- Data$bx, Data$Y <- NA),
                             Data$Y <- NA)
  
  temp             <- subset(Data, !(Age < minA) & !(Age >= maxA) & !is.na(Y))
  b                <- as.double( sd(temp$Y) / sd(temp$X))
  r                <- as.double( mean(temp$Y) - mean(temp$X) * b)
  c                <- (1 / b) * (exp (r * (censusyear - deathsmid)))
  
  Data$abx         <- ifelse((Data$Age - 1) < bgbmaxage, Data$abx <- r + b * Data$X, Data$abx <- NA)
  age_by           <- by(Data, 
                         INDICES = Data$Age, 
                         FUN = function(x){
                           x$residual <- ifelse(x$Age >= minA,
                                                ifelse(x$Age <= maxA, x$residual <- (x$Y - x$abx), x$residual <- NA),
                                                x$residual <- NA)
                           return(x)
                         })
  
  Data             <- do.call(rbind,age_by)
  Data             <- naniar::replace_with_na_at(Data, .vars = c("Nx","bx", "Y", "abx", "residual"), condition = ~.x == 0)
  Data$Pop.Ajust   <- ifelse(Data$Age >= minA, Data$Pop * exp( -r * (censusyear - deathsmid)), NA)
  Data$Death.Ajust <- ifelse(Data$Age >= minA, Data$Deaths / c, NA)
  Data$PYL.Ajust   <- ifelse(Data$Age >= minA, Data$Pop.Ajust * deathsperiod, NA)
  Data$mx.Ajust    <- ifelse(Data$Age >= minA, Data$Death.Ajust / Data$PYL.Ajust, NA)
  Data$qx          <- ifelse(Data$Age <= maxA+5,
                             ifelse(Data$Age >= minA, AgeInt * Data$mx.Ajust / ( 1 + 2.5 * Data$mx.Ajust), NA),
                             NA)
  Data$lx <- NA
  for(i in 1:nrow(Data)){
    ifelse(Data$Age[i] < minA,
           Data$lx[i] <- NA,
           ifelse(Data$Age[i] == minA,
                  Data$lx[i] <- 1.0, 
                  Data$lx[i] <- Data$lx[i-1] * (1-Data$qx[i-1])))
  }

  Data$Obs.Y       <- ifelse(!is.na(Data$lx),
                             ifelse(Data$lx != 0, 0.5 * (log((1 - Data$lx) / Data$lx)), NA),
                             NA)
  Data             <- mutate_if(Data, is.numeric, list(~na_if(., -Inf)))
  
  LT               <- gather(Model.LT.logits, fam, Ys, -c(sex, e0, Age))
  
  temp             <- subset(LT, sex == model.sex & fam == model.family)
  temp$Std         <- 1/ (1 + (exp(2 * temp$Ys)))
  temp$lx.Std      <- temp$Std / temp$Std [1]
  
  x                <- data.frame('Age' = seq(5, 100, by = 5))
  temp             <- merge(temp, x, by = "Age", all.y = TRUE)
  
  for(i in 1:nrow(temp)){
    if (temp$Age[i] > 85){
      temp$lx.Std[i] <- temp$lx.Std[i-1] * exp((log(temp$lx.Std[i-1] / temp$lx.Std[i-2])^2) / log(temp$lx.Std[i-2] / temp$lx.Std[i-3]))
    } else {
      next
    }
  }
  
  temp$Cdn.Ys      <- 0.5 * log ((temp$lx.Std [1] - temp$lx.Std) / temp$lx.Std)
  temp             <- mutate_if (temp, is.numeric, list(~na_if(., -Inf)))
  
  Data             <- full_join (Data, temp , by = "Age")
  
  model            <- subset(Data, Age >= model.agefrom & Age <= model.ageto)
  model            <- lm (formula = Obs.Y ~ Cdn.Ys, data = model)
  alp              <- as.double(model$coefficients [1])
  bet              <- as.double(model$coefficients [2])
  
  Data$Fit.Ys      <- alp + bet * Data$Cdn.Ys
  Data$Fit.lx      <- ifelse(Data$Age == minA, Data$Fit.lx <- 1, Data$Fit.lx <- 1 / (1 + exp(2 * Data$Fit.Ys)))
  
  for(i in nrow(Data):1){
    ifelse(Data$Age[i] == max(Data$Age),
           Data$Tx[i] <- 2.5 * Data$Fit.lx[i],
           ifelse(Data$Age[i] >= minA,
                  Data$Tx[i] <- Data$Tx[i + 1] + 2.5*(Data$Fit.lx[i + 1] + Data$Fit.lx[i]),
                  Data$Tx[i] <- NA))
  }
  
  Data$ex          <- ifelse(Data$Age >= minA, Data$Tx/Data$Fit.lx, Data$Tx <- NA)
  
  # Estimating Est.Nx 
  for(i in nrow(Data):1){
    
    if (Data$Age[i] < minA){
      next
    }
    
    ifelse (Data$Age[i] > 64, 
            ifelse(maxA < Data$Age[i],
                   Data$Est.Nx[i] <- (exp(r * Data$ex[i]) - (r * Data$ex[i] ) ^ 2 / 6) * Data$DeathCum[i],
                   Data$Est.Nx[i] <- Data$Est.Nx[i + 1] * exp(5 * r) + Data$Deaths[i] * exp(2.5 * r)),
            Data$Est.Nx[i] <- Data$Est.Nx[i + 1] * exp(5 * r) + Data$Deaths[i] * exp(2.5 * r))
  }
  
  # Estimating Est.5Nx
  for(i in nrow(Data):1){
    
    if (Data$Age[i] < minA){
      next
    }
    
    ifelse (Data$Age[i] < maxA,
            Data$Est.5Nx[i] <- round(2.5 * (Data$Est.Nx[i] + Data$Est.Nx[i+1]),1),
            Data$Est.5Nx[i] <- NA)
  }
  
  Data$Obs.5Nx   <- ifelse(maxA > Data$Age, 
                           Data$Pop*deathsperiod, 
                           NA)
  Data$cNx       <- ifelse(maxA > Data$Age,
                           Data$Est.5Nx/Data$Obs.5Nx,
                           NA)
  
  # Estimating a.xNx
  Data$a.xNx  <-  NA
  for(i in 1:nrow(Data)){
    if (Data$Age[i] < minA){
      next
    }
    ifelse (Data$Age[i] < maxA,
            Data$a.xNx[i] <- sum(Data$Est.5Nx[i:N],na.rm = T)/sum(Data$Obs.5Nx[i:N], na.rm = T),
            Data$a.xNx[i] <- NA)
  }  
  
  Method <- Data[(Data$Age <= model.ageto), (names(Data) %in% c("Age","Pop","Deaths", "Est.Nx", "Est.5Nx", "Obs.5Nx", "cNx", "a.xNx"))]
  
  title  <- paste('Ratios of Census Completeness by Age', 'Using Preston-Coale Method', sep ='\n')
  
  title2 <- paste('(', countryname, ': ' , as.character(censusyear), ')', sep = "")
  
  ratios.plot <- ggplot(subset(Method, !is.na(a.xNx & cNx)), aes(Age)) + 
    geom_path(aes(y = cNx,  color = "darkred"), size = 1.5,  linetype="twodash") + 
    geom_path(aes(y = a.xNx, color = "steelblue"), size = 1.5) +
    labs(x = 'Ages',
         y = 'Ratio of Completeness') +
    scale_color_manual(name = "Ratios", labels = c("c: N(x to x+5)", "c: N(x to A)"), values = c("darkred", "steelblue")) +
    theme_bw() +
    ggtitle(paste(title, title2, sep = "\n ")) + 
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
  
  ratios.plot
}

#############Run all results#################################################
# Run all results

## Estimations
prestoncoalecolumns(Data = elsalvador, censusyear = 1961, lowerage = 15, upperage = 60,  deaths.summed = FALSE, 
                    deathsyear1 = 1961, deathsyear2 = 1961, minA = 5, maxA = 74, bgbmaxage = 69, 
                    model.LT = Model.LT.logits, model.sex = 'female', model.family = 'West', model.agefrom = 45, model.ageto = 75)

prestoncoalecov (Data = elsalvador, censusyear = 1961, lowerage = 15, upperage = 60,  deaths.summed = FALSE, 
                 deathsyear1 = 1961, deathsyear2 = 1961, minA = 5, maxA = 74, bgbmaxage = 69, 
                 model.LT = Model.LT.logits, model.sex = 'female', model.family = 'West', model.agefrom = 45, model.ageto = 75)

prestoncoale.ratiosplot (Data = elsalvador, censusyear = 1961, lowerage = 15, upperage = 60,  deaths.summed = FALSE, 
                         deathsyear1 = 1961, deathsyear2 = 1961, minA = 5, maxA = 74, bgbmaxage = 69, countryname ='El Salvador', 
                         model.LT = Model.LT.logits, model.sex = 'female', model.family = 'West', model.agefrom = 45, model.ageto = 75)   
