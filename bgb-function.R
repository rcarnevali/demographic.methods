# Contains functions related to the Brass Growth Balance method

###############################################################################

#' @title Estimate death registration coverage for a single year/sex/region using the Brass Growth Balance method
#'
#' @description Given one census and an average annual number of deaths in each age class, we can use stable population assumptions to estimate the completeness of the reporting of deaths relative to an estimate of the population. The method makes use of the observation that in a stable population that is closed to migration and has accurately reported data, the growth rate, r, is equal to the birth rate, b, less the death rate, d.
#' 
#' @details \code{bgbcolumns()}
#' 
#' @param Data \code{data.frame} with columns, \code{$Age}, \code{$Population}, \code{$Deaths}.
#' @param censysyear Year of the census
#' @param minA The minimum of the age range searched. Default 5
#' @param maxA The maximum of the age range searched. Default 69
#' @param deaths.summed Logical. is the deaths column given as the total per age in the intercensal period (\code{TRUE}). By default we assume \code{FALSE}, i.e. that the average annual was given.
#' @param deathsyear1 The first reference year for the deaths' distribution
#' @param deathsyear2 The final reference year for the deaths' distribution. If the deaths refer to just one year, this parameter may be left blank 
#' 
#'  
#' @export
#' 
#' @example
#' ## El Salvador 1961 Females Census data:
#' 
#' elsalvador <- data.frame("Age" = seq(0, 75, 5), 
#' "Pop" = c(214089, 190234, 149538, 125040, 113490, 91663, 77711, 72936, 56942, 46205, 38616, 26154, 29273, 14964, 11205, 16193),
#' "Deaths" = c(6909, 610, 214, 266, 291, 271, 315, 349, 338, 357, 385, 387, 647, 449, 504, 1360))
#' 
#' bgbcolumns(Data = elsalvador, censusyear = 1961, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1 = 1961)
#' 

bgbcolumns <- function (Data = Data, censusyear = year, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1, deathsyear2){
  require(dplyr)
  require(DDM)
  
  if (abs(censusyear-deathsyear1)>5)
    
    warning("The census' reference year and the deaths' first reference year are more than 5 years apart")
  
  #Took the codes from Tim's DDM functions.
  #Obs: Must check how to properly  give credit where credit is due
    
  #1) Change the name of deaths column to "Deaths" to fit the code
  tempcols <- tolower(colnames(Data))
  dcol     <- grepl("death", tempcols)
  stopifnot(any(dcol))
  colnames(Data)[dcol] <- "Deaths"
  
  #2) Change the name of population column to "Pop" to fit the code
  tempcols <- tolower(colnames(Data))
  pcol     <- grepl("population", tempcols)
  stopifnot(any(pcol))
  colnames(Data)[pcol] <- "Pop"

  #3) If the deaths should be considered an average of the period or they refer to just a year. Variation from Tim Riffe's code
    #Hence 'alt' refer to alteration. 
    #Not the best solution, but Tim has very complex code and this is an easy way out
    
  avgDeaths.alt <- function(Data, deaths.summed = FALSE){
    if (!"deathsAvg" %in% colnames(Data)){
      if (deaths.summed){
        dif <- ((deathsyear1+0.5)+(deathsyear2+0.5))/(deathsyear2-deathsyear1)
        Data$deathsAvg <- Data$deaths / dif
        
      } else {
        
        Data$deathsAvg <- Data$death
        deathsyear2 <- deathsyear1
        deathmid <- deathsyear1 + 0.5
      }
    }
    Data
  }
  
    
  Data <- avgDeaths.alt(Data = Data, deaths.summed = deaths.summed)
  
  #4) This detect the age interval within the data   
  AgeInt <- DDM::detectAgeInterval(Dat = Data, MinAge =  minA, MaxAge = maxA, ageColumn = "age")

    
  #5) If there is a 0-1 and 1-4 age groups, this part turns them into 0-4 group
  ages  <- Data$age
  
    if (1 %in% ages){
      ind0 <- ages == 0
      ind1 <- ages == 1
      cnames <- c("Pop","Deaths")
      if ("deathsAvg" %in% colnames(X)){
        cnames <- c(cnames,"deathsAvg")
      }
      
      # TR: look for age interval column and do that too?
      X[ind0,cnames] <- colSums(X[ ind0 | ind1, cnames])
      
      X <- X[!ind1, ]
    }
  
  N <- nrow(Data)
  
  Data$Pop <- as.double(Data$Pop)
  
  Data$Deaths <- as.double(Data$Deaths)

  Data <- Data %>%
    mutate(PopCum = rev(cumsum(rev(Pop))),
           DeathCum = rev(cumsum(rev(Deaths))),
           PYL = PopCum*deathsperiod,
           Nx = (c(0, rep(sqrt(Pop[ -N  ] * Pop[ -1 ]), length.out = N-2), 0)*deathsperiod)/AgeInt, #tem alguma coisa dando errado no ultimo grupo
           bx = Nx/PYL,
           X = ifelse(Nx == 0, #Change to X
                      X <- 0,
                      X <- DeathCum/PYL),
           Y = ifelse(age-1 < maxA,
                      Y <- bx,
                      Y <- NA))
  
  b <- Data %>%
    select(age, Y, X) %>%
    filter(!(minA  > age) & !(age >= maxA)) %>%
    select(X, Y) %>%
    summarise(sd = sd(Y)/sd(X))%>%
    as.double()
  
  a <- Data %>%
    select(age, Y, X) %>%
    filter(!(minA  > age) & !(age >= maxA)) %>%
    select(X, Y) %>%
    summarise(a = mean(Y)-mean(X)*b) %>%
    as.double()
  
  c <- (1/b)*(exp(a*(censusyear-deathsmid)))
  
  Data <- Data %>%
    mutate(abx = ifelse(age-1 < maxA,
                        abx <- a+b*X,
                        abx <- NA)) %>%
    group_by(age) %>%
    mutate(residual = ifelse(age >= minA,
                             ifelse(age <= maxA,
                                    residual <- (Y-abx),
                                    residual <- NA),
                             residual <- NA)) %>%
    naniar::replace_with_na_at(.vars = c("Nx","bx", "Y", "abx", "residuals"),
                               condition = ~.x == 0) %>%
    ungroup() %>%
    as.data.table()

  Data
}


###############################################################################

#' @details \code{bgbcoverage()}
#' 
#' @param Data a chunk of data (single sex, year, etc) with all columns required by \code{bgbcolumns()}
#' @param censysyear Year of the census
#' @param minA The minimum of the age range searched. Default 5
#' @param maxA The maximum of the age range searched. Default 69
#' @param deaths.summed Logical. is the deaths column given as the total per age in the intercensal period (\code{TRUE}). By default we assume \code{FALSE}, i.e. that the average annual was given.
#' @param deathsyear1 The first reference year for the deaths' distribution
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
#' "Pop" = c(214089, 190234, 149538, 125040, 113490, 91663, 77711, 72936, 56942, 46205, 38616, 26154, 29273, 14964, 11205, 16193),
#' "Deaths" = c(6909, 610, 214, 266, 291, 271, 315, 349, 338, 357, 385, 387, 647, 449, 504, 1360))
#' 
#' bgbcoverage(Data = elsalvador, censusyear = 1961, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1 = 1961)
#' 

bgbcoverage <- function(Data = Data, censusyear = year, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1, deathsyear2) {
  require(dplyr)
  require(DDM)
  
  values <- bgbcolumns (Data = Data, censusyear = year, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1, deathsyear2)
  
  results <- formattable::percent(t(data.table("Annual growth rate of stable population" = a, 
                                               "Completeness relative to population at midpoint"= c)), 
                                  digits = 2)
  results
}

###############################################################################

#' @details \code{bgb.fitted.graph()}
#' 
#' @param Data a chunk of data (single sex, year, region, etc) with all columns required by \code{bgbcolumns()}
#' @param censysyear Year of the census
#' @param minA The minimum of the age range searched. Default 5
#' @param maxA The maximum of the age range searched. Default 69
#' @param deaths.summed Logical. is the deaths column given as the total per age in the intercensal period (\code{TRUE}). By default we assume \code{FALSE}, i.e. that the average annual was given.
#' @param deathsyear1 The first reference year for the deaths' distribution
#' @param deathsyear2 The final reference year for the deaths' distribution. If the deaths refer to just one year, this parameter may be left blank 
#' 
#'@return A graph with the observed information as dots and the fitted information as lines estimated from the Brass Growth Balance Method
#'@export
#'
#'@example
#'El Salvador 1961 Females Census data:
#' 
#' elsalvador <- data.frame("Age" = seq(0, 75, 5), 
#' "Pop" = c(214089, 190234, 149538, 125040, 113490, 91663, 77711, 72936, 56942, 46205, 38616, 26154, 29273, 14964, 11205, 16193),
#' "Deaths" = c(6909, 610, 214, 266, 291, 271, 315, 349, 338, 357, 385, 387, 647, 449, 504, 1360))
#'     
#' bgb.fitted.graph(Data = elsalvador, censusyear = 1961, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1 = 1961)

bgb.fitted.graph <- function(Data = Data, censusyear = year, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1, deathsyear2) {
  
  values <- bgbcolumns (Data = Data, censusyear = year, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1, deathsyear2)
  
  graph <- values %>%
    select(X, Y, abx) %>%
    ggplot(mapping = aes(x = X)) + 
    geom_path(aes(y = abx), color =  "#0072B2", size = 0.8) + 
    geom_point(aes(y = Y), color =  "#D55E00", size = 2.5) +                        
    labs(x="d(x+)",
         y = "b(x+)",
         title ="TEST: Brass Growth Balance - El Salvador (1961)") +
    #scale_color_hue(labels = c("Fitted", "Observed")) + #does not work!!
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5, face="bold"))
  
  graph
}

###############################################################################

#' @details \code{bgb.residuals.graph()}
#' 
#' @param Data a chunk of data (single sex, year, region, etc) with all columns required by \code{bgbcolumns()}
#' @param censysyear Year of the census
#' @param minA The minimum of the age range searched. Default 5
#' @param maxA The maximum of the age range searched. Default 69
#' @param deaths.summed Logical. is the deaths column given as the total per age in the intercensal period (\code{TRUE}). By default we assume \code{FALSE}, i.e. that the average annual was given.
#' @param deathsyear1 The first reference year for the deaths' distribution
#' @param deathsyear2 The final reference year for the deaths' distribution. If the deaths refer to just one year, this parameter may be left blank 
#' 
#'@return A graph with residuals estimated from the Brass Growth Balance Method
#'@export
#'
#'@example
#'El Salvador 1961 Females Census data:
#' 
#' elsalvador <- data.frame("Age" = seq(0, 75, 5), 
#' "Pop" = c(214089, 190234, 149538, 125040, 113490, 91663, 77711, 72936, 56942, 46205, 38616, 26154, 29273, 14964, 11205, 16193),
#' "Deaths" = c(6909, 610, 214, 266, 291, 271, 315, 349, 338, 357, 385, 387, 647, 449, 504, 1360))
#' 
#' bgb.residuals.graph(Data = elsalvador, censusyear = 1961, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1 = 1961)

bgb.residuals.graph <- function(Data = Data, censusyear = year, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1, deathsyear2) {
  
  values <- bgbcolumns (Data = Data, censusyear = year, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1, deathsyear2)
  
  values %>%
    filter(Age >= agemin) %>%
    select(X, residual) %>%
    ggplot(mapping = aes(x = X, 
                         y = residual)) + 
    geom_line(aes(), colour = "#0072B2", size = 0.5) +
    geom_point(aes(), colour = "#0072B2", size = 2)+
    labs(x="d(x+)",
         y = "residual",
         title ="TEST: Residual estimation for Brass Growth Balance - \n El Salvador (1961)") +
    theme_bw()+
    geom_hline(yintercept = 0, show.legend = NA, linetype="dashed", 
               color = "#D55E00", size = 1)+
    theme(plot.title = element_text(hjust = 0.5, face="bold"))
}

###############################################################################
# Run all results
  #Sample data
    elsalvador <- data.frame("Age" = seq(0, 75, 5), 
     "Pop" = c(214089, 190234, 149538, 125040, 113490, 91663, 77711, 72936, 56942, 46205, 38616, 26154, 29273, 14964, 11205, 16193),
     "Deaths" = c(6909, 610, 214, 266, 291, 271, 315, 349, 338, 357, 385, 387, 647, 449, 504, 1360))
    # Census Year: 1961
    # Deaths Year: 1961

 
  #Brass estimations   
  bgbcolumns (Data = Data, censusyear = year, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1, deathsyear2)
  bgbcoverage (Data = Data, censusyear = year, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1, deathsyear2)
  bgb.fitted.graph (Data = Data, censusyear = year, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1, deathsyear2)
  bgb.residuals.graph (Data = Data, censusyear = year, minA = 5, maxA = 69, deaths.summed = FALSE, deathsyear1, deathsyear2)
