#Estimation of adult mortality using the Growth Balance method

#Please refer to: [BGB](http://demographicestimation.iussp.org/content/brass-growth-balance-method/)
#this code works as a base for the functionalized method. The data used here is just an example and is not used in the functionalized version

################
## Packages ##
# List of packages for session
.packages = c("devtools", "data.table", "psych", "tidyverse", "formattable", "lubridate", "naniar", "Hmisc", "janitor", "dplyr")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0)
  install.packages(.packages[!.inst],dependencies = T)

# Load packages into session
lapply(.packages, require, character.only = T)

#############################################################################################################
## Example Data ##
Data <- data.frame("age" = seq(0, 75, 5), 
                   "Pop" = c(214089, 190234, 149538, 125040, 113490, 91663, 77711, 72936, 56942, 46205, 
                             38616, 26154, 29273, 14964, 11205, 16193),
                   "Deaths" = c(6909, 610, 214, 266, 291, 271, 315, 349, 338, 357, 385, 387, 647, 449, 504, 1360))


 # Converting integers into doubles
 Data$Pop <- as.double(Data$Pop)
 Data$Deaths <- as.double(Data$Deaths

# Read census date as fraction
#This was later changed in the function. Intead of using the census' full reference date, the function uses only the reference year to create a 'mid point'.
censusdate <- as.Date("05/05/1961", format = "%d/%m/%Y") #defining the date of the census as a 'date' format
                          
fraccensusdate <- yrfraction(censusdate) + as.numeric(lubridate::year(censusdate)) 

# Define period of deaths
#This was later changed in the function. Intead of using the death's full reference date, the function uses only the difference between two reference years to create a 'mid point'.
deathsperiod <- round((as.numeric((as.Date("31/12/1961", format = "%d/%m/%Y")) - (as.Date("01/01/1961", format = "%d/%m/%Y"))) / 365), 
                        digits = 2) #Estimating the death's reference period

deathsmid <- (deathsperiod / 2) + year #Finding the mid point

# Estimating the age interval
AgeInt <- DDM::detectAgeInterval(Dat = Data, ageColumn = "age") 
                          
N <- nrow(Data)

agemin <-	as.double(5)
agemax <- as.double(69)

#Age Limits used for the mortality estimation
minA <- 45
maxA <- 75

#################################################################################################
## Estimating Completeness and Growth Rate ##
Data <- Data %>%
  mutate(PopCum = rev (cumsum (rev (Pop))), #Cumulate population
         DeathCum = rev (cumsum (rev (Deaths))), #Cumulate deaths
         PYL = PopCum * deathsperiod, #person-years of life lived
         Nx = (c (0, rep(sqrt(Pop[ -N  ] * Pop[ -1 ]), length.out = N - 2), 0) * deathsperiod) / AgeInt, #number of people who turned x
         bx = Nx / PYL, #partial ‘birth’ rate
         X = ifelse(Nx == 0, #partial death rate
                    X <- 0,
                    X <- DeathCum/PYL),
         Y = ifelse(age - 1 < agemax,
                    Y <- bx,
                    Y <- NA)
        )

b <- Data %>% 
  select(age, Y, X) %>%
  filter(!(agemin  > age) & 
         !(age >= agemax)
        ) %>%
  select(X, Y) %>%
  summarise(sd = sd(Y) / sd(X))%>%
  as.double()

a <- Data %>% #annual growth rate of the stable population
  select(age, Y, X) %>%
  filter(!(agemin  > age) & 
         !(age >= agemax)
        ) %>%
  select(X, Y) %>%
  summarise(a = mean(Y) - mean(X) * b) %>%
  as.double()

c <- (1 / b) * (exp (a * (fraccensusdate - deathsmid))) #completeness of reporting of deaths

Data <- Data %>%
  mutate(abx = ifelse(age - 1 < agemax,
                      abx <- a + b * X,
                      abx <- NA)
        ) %>%
  group_by(age) %>%
  mutate(residual = ifelse(age >= agemin, #residuals = y-(a+bx)
                           ifelse(age <= agemax,
                                  residual <- (Y-abx),
                                  residual <- NA),
                           residual <- NA)
        ) %>%
  ungroup() %>%
  naniar::replace_with_na_at(.vars = c("Nx","bx", "Y", "abx", "residuals"),
                             condition = ~.x == 0) %>%
  mutate(Pop.Ajust = ifelse(age >= agemin,  #population at the mid-point of the period over which the deaths were recorded
                            Pop * exp(- a * (fraccensusdate - deathsmid)),
                            NA),
         Death.Ajust = ifelse(age >= agemin, #adjusted number of deaths for incompleteness
                              as.integer(Deaths/c),
                              NA),
         PYL.Ajust = ifelse(age >= agemin, #person-years of exposure
                            Pop.Ajust * deathsperiod,
                            NA),
         mx.Ajust = ifelse(age >= agemin,  #mortality rates adjusted for incompleteness of reporting of deaths
                           Death.Ajust / PYL.Ajust,
                           NA),
         qx = ifelse(age <= agemax + 5, #probabilities of death from the adjusted rates of mortality
                     ifelse(age >= agemin,
                            AgeInt * mx.Ajust / (1 + 2.5 * mx.Ajust),
                            0),
                     1),
         px = 1 - qx, #survivourship probabilities
         lx = cumprod(px), #survivours at the exact age x
         lx = 1 * c(1, rep(lx, length.out = N - 1)), #put one in front of everything
         lx = lx [1:length(px)],
         Obs.Y = ifelse(!is.na(lx),
                        ifelse(lx != 0,
                               0.5 * (log ((1 - lx) / lx)),
                               NA),
                        NA)
        ) %>%
  mutate_if(is.numeric, list(~na_if(., -Inf)))

################
# Getting the Model life table (LT) information
LT <- melt(Model.LT.logits, na.rm = FALSE,
               value.name = "Ys", 
               measure.vars = c("UN.General",	"Princeton.East",	"Princeton.North", "Princeton.South",	"Princeton.West", "AIDS"))

# Smooth using relational logit model life table
              
  #Obs: From here onwards is the estimation using the logit of model Life Tables

  #Note: Must find a way to handle the model LT choice. When estimating the completeness level it isn't necessary, however, when estimating the preston-coale method we need the life expectancy estimated through bgb (or in another way, but we need the life expactancy anyway)
                          
Data <- LT %>%
  filter(sex == "female", #selecting the model LT sex based on the census information
         variable == "Princeton.West") %>% #selecting the LT model > How should we do to select the e0 level?
  mutate(Std = 1 / (1 + exp (2 * Ys)),
         lx.Std = Std / Std[1],
         Cdn.Ys = 0.5 * log((lx.Std[1] - lx.Std) / lx.Std) #logit transformations of the proportions surviving
        ) %>%
  mutate_if(is.numeric, list(~na_if(., -Inf))) %>%
  right_join(Data, by = "age") %>%
  select(age, Obs.Y, Cdn.Ys)

model <- Data %>%
  filter(age >= minA & 
           age <= maxA) %>%
  lm(formula = Obs.Y ~ Cdn.Ys, data = .) #fitting the relational logit model

#intercept and slope, respectively, of the straight line fitted to the logit transformations over the range of ages chosen by the user (45 and 75 in this case)
alpha <- as.double(model$coefficients[1])
beta <- as.double(model$coefficients[2])

Data <- Data %>%
  mutate(Fit.Ys = alpha + beta * Cdn.Ys, 
         Fit.lx = 1 / (1 + exp (2 * Fit.Ys)),
         Fit.lx = replace_na(Fit.lx, 1),
         Tx = (5 * Fit.lx), #conditional years of life lived
         Tx = c(rep (rev (cumsum (rev (Tx))), length.out = N - 1), (5 * Fit.lx[N]) + 1), #This estimation is very problematic
         ex = Tx / Fit.lx #Life expectancy
        )  

Data[1, c("Fit.lx","Tx", "ex")] <- NA
