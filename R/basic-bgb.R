#Estimation of adult mortality using the Growth Balance method

#Refer to: http://demographicestimation.iussp.org/content/brass-growth-balance-method

#############################################################################################################
## Packages ##
# List of packages for session
.packages = c("devtools", "data.table", "psych", "tidyverse", "formattable", "rvest","lubridate", "naniar",
              "cwhmisc","Hmisc", "janitor", "dplyr")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0)
  install.packages(.packages[!.inst],dependencies = T)

# Load packages into session
lapply(.packages, require, character.only = T)

#############################################################################################################
Data <- data.frame("age" = seq(0, 75, 5), 
                         "Pop" = c(214089, 190234, 149538, 125040, 113490, 91663, 77711, 72936, 56942, 46205, 
                                   38616, 26154, 29273, 14964, 11205, 16193),
                         "Deaths" = c(6909, 610, 214, 266, 291, 271, 315, 349, 338, 357, 385, 387, 647, 449, 504, 1360))

fraccensusdate <- 1961.3397260274
deathsmid <- 1961.5
deathsperiod <- 1
AgeInt <- 5
N <- nrow(Data)

agemin <-	as.double(5)
agemax <- as.double(69)

#Mortality Rates Estimation Data

#Model Life Tables
Model.LT.logits <- data.frame("sex" = rep(c("female", "male"), each = 17),
                              "e0" = 60,
                              "age" = rep(seq(5, 85, 5), times =2),
                              "UN General" = c(-1.0558, -1.0060, -0.9779, -0.9381, -0.8876, -0.8322, -0.7716, -0.7054, -0.6323, -0.5472,
                                               -0.4432, -0.3141, -0.1535, 0.0472, 0.2983, 0.6132, 1.0182, -1.1183, -1.0744, -1.0466, 
                                               -1.0067, -0.9539, -0.8987, -0.8389, -0.7690, -0.6845, -0.5809, -0.4533, -0.2991, 
                                               -0.1120, 0.1135, 0.3839, 0.7102, 1.1176),
                              "Princeton East" = c(-1.0035, -0.9621, -0.9369, -0.9014, -0.8565, -0.8085, -0.7572, -0.7014, -0.6406, 
                                                   -0.5697, -0.4806, -0.3653, -0.2118, -0.0065, 0.2668, 0.6303, 1.1332, -1.0839, 
                                                   -1.0473, -1.0232, -0.9821, -0.9268, -0.8758, -0.8247, -0.7668, -0.6954, -0.6029, 
                                                   -0.4816, -0.3307, -0.1485, 0.0744, 0.3569, 0.7245, 1.2272),
                              "Princeton North" = c(-1.0613, -0.9787, -0.9346, -0.8858, -0.8306, -0.7720, -0.7117, -0.6490, -0.5801, 
                                                    -0.5086, -0.4228, -0.3210, -0.1876, -0.0105, 0.2241, 0.5305, 0.9446, -1.1186, 
                                                    -1.0394, -0.9981, -0.9416, -0.8692, -0.8021, -0.7373, -0.6722, -0.6004, -0.5206, 
                                                    -0.4210, -0.3068, -0.1606, 0.0250, 0.2639, 0.5773, 1.0014),
                              "Princeton South" = c(-0.9088, -0.8743, -0.8527, -0.8221, -0.7841, -0.7431, -0.7021, -0.6578, -0.6081, 
                                                    -0.5521, -0.4801, -0.3886, -0.2591, -0.0798, 0.1754, 0.5294, 1.0331, -0.9780, 
                                                    -0.9478, -0.9272, -0.8981, -0.8562, -0.8170, -0.7719, -0.7223, -0.6604, -0.5838, 
                                                    -0.4846, -0.3597, -0.2018, -0.0036, 0.2626, 0.6309, 1.1510),
                              "Princeton West" = c(-1.0842, -1.0329, -0.9956, -0.9449, -0.8841, -0.8214, -0.7564, -0.6886, -0.6163, 
                                                   -0.5354, -0.4371, -0.3187, -0.1664, 0.0260, 0.2779, 0.6087, 1.0567, -1.1444, 
                                                   -1.0968, -1.0631, -1.0130, -0.9485, -0.8877, -0.8259, -0.7570, -0.6756, -0.5772, 
                                                   -0.4558, -0.3074, -0.1262, 0.0928, 0.3652, 0.7161, 1.1877),
                              "AIDS" = c(-0.8518, -0.7992, -0.7789, -0.7400, -0.6675, -0.5550, -0.4395, -0.3475, -0.2797, -0.2268, 
                                         -0.1732, -0.1105, -0.0246, 0.1003, 0.2644, 0.4789, 0.7683, -0.7983, -0.7369, -0.6961, 
                                         -0.6654, -0.6193, -0.5383, -0.4240, -0.3054, -0.1910, -0.0913, -0.0058, 0.0805, 0.1841, 
                                         0.3261, 0.5119, 0.7559, 1.0995))

#Age Limits
minA <- 45
maxA <- 75

#################################################################################################
## Essa parte calcula a aba "Method" da planilha da IUSSP ##
#Ela roda sem problemas e pode rodar ela toda tambÃ©m

Data <- Data %>%
  mutate(PopCum = rev(cumsum(rev(Pop))),
         DeathCum = rev(cumsum(rev(Deaths))),
         PYL = PopCum*deathsperiod,
         Nx = (c(0, rep(sqrt(Pop[ -N  ] * Pop[ -1 ]), length.out = N-2), 0)*deathsperiod)/AgeInt, 
         bx = Nx/PYL,
         X = ifelse(Nx == 0, #Change to X
                    X <- 0,
                    X <- DeathCum/PYL),
         Y = ifelse(age-1 < agemax,
                    Y <- bx,
                    Y <- NA))

b <- Data %>%
  select(age, Y, X) %>%
  filter(!(agemin  > age) & !(age >= agemax)) %>%
  select(X, Y) %>%
  summarise(sd = sd(Y)/sd(X))%>%
  as.double()

a <- Data %>%
  select(age, Y, X) %>%
  filter(!(agemin  > age) & !(age >= agemax)) %>%
  select(X, Y) %>%
  summarise(a = mean(Y)-mean(X)*b) %>%
  as.double()


c <- (1/b)*(exp(a*(fraccensusdate-deathsmid)))

Data <- Data %>%
  mutate(abx = ifelse(age-1 < agemax,
                      abx <- a+b*X,
                      abx <- NA)) %>%
  group_by(age) %>%
  mutate(residual = ifelse(age >= agemin,
                           ifelse(age <= agemax,
                                  residual <- (Y-abx),
                                  residual <- NA),
                           residual <- NA)) %>%
  ungroup() %>%
  naniar::replace_with_na_at(.vars = c("Nx","bx", "Y", "abx", "residuals"),
                             condition = ~.x == 0) %>%
  mutate(Pop.Ajust = ifelse(age >= agemin,
                            Pop * exp(-a*(fraccensusdate-deathsmid)),
                            NA),
         Death.Ajust = ifelse(age >= agemin,
                              as.integer(Deaths/c),
                              NA),
         PYL.Ajust = ifelse(age >= agemin,
                            Pop.Ajust*deathsperiod,
                            NA),
         mx.Ajust = ifelse(age >= agemin,
                           Death.Ajust/PYL.Ajust,
                           NA),
         qx = ifelse(age <= agemax+5,
                     ifelse(age >= agemin,
                            AgeInt*mx.Ajust/(1+2.5*mx.Ajust),
                            0),
                     1),
         px = 1-qx,
         lx = cumprod(px),
         lx = 1 * c(1, rep(lx, length.out = N-1)), #Put one in front of everything
         lx = lx[1:length(px)],
         Obs.Y = ifelse(!is.na(lx),
                        ifelse(lx != 0,
                               0.5*(log((1-lx)/lx)),
                               NA),
                        NA)) %>%
  mutate_if(is.numeric, list(~na_if(., -Inf)))

LT <- melt(Model.LT.logits, na.rm = FALSE,
               value.name = "Ys", 
               measure.vars = c("UN.General",	"Princeton.East",	"Princeton.North", "Princeton.South",	"Princeton.West", "AIDS"))


Data <- LT %>%
  filter(sex == "female",
         variable == "Princeton.West") %>%
  mutate(Std = 1/(1+exp(2*Ys)),
         lx.Std = Std/Std[1],
         Cdn.Ys = 0.5*log((lx.Std[1]-lx.Std)/lx.Std)) %>%
  mutate_if(is.numeric, list(~na_if(., -Inf))) %>%
  right_join(Data, by = "age") %>%
  select(age, Obs.Y, Cdn.Ys)


model <- Data %>%
  filter(age >= minA & 
           age <= maxA) %>%
  lm(formula = Obs.Y ~ Cdn.Ys, data = .)

alpha <- as.double(model$coefficients[1])
beta <- as.double(model$coefficients[2])


Data <- Data %>%
  mutate(Fit.Ys = alpha+beta*Cdn.Ys, 
         Fit.lx = 1/(1 + exp(2*Fit.Ys)),
         Fit.lx = replace_na(Fit.lx, 1),
         Tx = (5*Fit.lx),
         Tx = c(rep(rev(cumsum(rev(Tx))), length.out = N-1), (5*Fit.lx[N])+1),
         ex = Tx/Fit.lx) 

Data[1,c("Fit.lx","Tx", "ex")] <- NA
