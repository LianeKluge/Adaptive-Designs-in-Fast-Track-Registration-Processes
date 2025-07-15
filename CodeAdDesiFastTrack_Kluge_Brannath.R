# Code for plots of Paper 
# "Adaptive Designs in Fast-Track Registration Processes"

# authors: Liane Kluge and Werner Brannath

# PLEASE NOTE: To ensure that the paper contains accurate diagrams,
# the step size was often set to a small value. This means that the numerical
# calculations take some time. To obtain results more quickly,
# please reduce the step sizes! This applies in particular to the vectors
# relativeI1vec_1 to relativeI1vec_6.

# PLEASE NOTE: When running the code, images are saved in the current working 
# directory!

###################### Code Structure #########################################
# For clarity, the code is divided into five parts.

# 1) The first part provides the functions needed to determine the minimum needed 
# (relative) informationI_2,min for the second stage

# 2) The second part provides the functions needed to  determine the constant 
# information I_2,const 

# 3) The third part gives the functions to compare the efficiency of the 
# different designs 

# 4) The fourth part is needed to process the input parameters for the functions
# of part 1-3

# 5) The last part runs the functions and provides the plots of the paper


################### Explanations/ overview of input parameters and assumptions 
# of the main-functions #######################################################

# A: conditional error function, either A_IN, A_F, A_Z or the separate study
# design (Function A_C), see Part 5

# c_A: vector of level constants (c_A[i] corresponds to relativeI1vec[i])
# for the conditional error functions A_IN, A_F, A_Z or A_C
# (see part 5 for definition of the functions)

# iNMethodChangingWeights: Indicates whether A_Z is considered or one of the
# other functions A_IN. A_F or A_C. A_Z is considered iff 
# iNMethodChangingWeights == TRUE

# relativeI1vec: vector of all considered relative first stage informations
# I_1/I_delta where I_delta is the information needed to achieve a 
# power of 1-beta to reject H_0 with a single fixed size sample test

# relativeI_2Const: a numeric vector of the same length as relativeI1vec,
# giving for each considered (relative) information relativeI1vec[i]
# the (relative) constant information used for the second stage
# when no conditional registration is achieved. Determination of 
# relativeI_2Const via the function determination_I_2Const.

# xi: xi=delta/delta_rel is the relative a priori assumed or true healthcare 
# effect 

# z_er: vector of early rejection boundaries on the z-scale. z_er=Inf holds for 
# the examples in the paper. z_er[i] corresponds to relativeI1vec[i]. z_er can also
# have length 1 if for all considered first stage informations z_er is the same.
# z_er >= pmax(z_f,z_s) is assumed. 

# z_f: vector of boundaries on the z-scale to achieve early rejection 
# (iff Z_1>= z_f). z_f[i] corresponds to relativeI1vec[i]. z_f < z_er is assumed

# z_s: vector of early futility boundaries on the z-scale. The two-staged design is
# stopped at interim iff Z_1< z_s. z_s < z_er is assumed. z_s[i] corresponds to 
# relativeI1vec[i]. z_s can also have length 1 if for all considered first stage 
# informations z_s is the same. 

# p_1: for the examples in the paper with required early rejection, 
# p_1 = 1-beta is set and gives the aimed overall success probability.
# For the examples without required early rejection, 
# p_1[i] = (1-beta)*P_delta(Z_1 >= z_f[i]) is set to achieve a conditional 
# success probability of 1-beta. p_1[i] corresponds to relativeI1vec[i].

# p_2 : only needed for the examples in the paper without required conditional
# registration. It holds p_2[i]=P_delta(Z_1<z_f[i]) and p_2[i] corresponds to
# relativeI1vec[i].

# uppBoundSearch_I_2Min: 0 and uppBoundSearch_I_2Min give the limits for the 
# bisection search to determine I_2,min. 

# uppBoundSearch_I_2Const: 0 and uppBoundSearch_I_2Const give the limits for the 
# bisection search to determine I_2,Const. uppBoundSearch_I_2Const should be 
# greater 1. 

# tolBisecDeterm_I_2Min: error tolerance of bisection method to determine 
# I_2,Min

# tolBisecDeterm_I_2Const: error tolerance of bisection method to determine
# I_2,Const

# stabConstIntegral: required to stabilize integration with normal densities
# for more accurate results.
# stabConstIntegral should be the lowest positive z-value such that 
# dnorm(stabConstIntegral) = 0 numerically. When dnorm returns the value 0 
# depends on the floating point representation. We determined 38.56802  for 
# the IEEE 754 double-precision format. Please note that the stability constant 
# may need to be adjusted if you are using a different floating point 
# representation!
# Numerical calculations show that very similar results can be achieved even 
# without stabilization. Therefore, it is not impossible to set
# stabConstIntegral <- Inf

# searchIntervalDeterLevelConst: vector of length 2 where the first entry must
# be lower than the second one. It gives the bounds for the determination
# of the level constants for the inverse normal method with fixed weights
# and Fisher's product test. This is necessary because the conditional error
# functions are truncated with 0.5 and may also have an early futility stop.
# Determining the level constant is important for achieving the 
# Type I error rate alpha. searchIntervalDeterLevelConst[1] must be chosen such
# the the Type I error rate is controlled. The vector is needed for the 
# determineLevelConst-function.

# tolBisecDetermLevelConst: error tolerance of bisection method to determine
# the level constants with the determineLevelConst-function or the
# determineLevelConstA_Z- function.


###############################################################################
###############################################################################
###############################################################################
# PART 1: functions to determine the minimum needed  (relative) information 
# I_2,min for the second stage
###############################################################################

############################# helper functions of PART 1 #####################


# The function calculates the power to receive both a conditional
# and a final registration (for designs without early rejection as
# in the paper).  The function is  needed to determine the minimum needed 
# (relative) Information I_2,min for  the second stage (with the function
# determinationI_2Min) to achieve a certain overall power  (if conditional
# registration is required) or to achieve a certain probability of success
# condition on Z_1>=z_f (if conditional registration is not required).
# The function can also be used  to calculate for designs with early rejection
# (not explicitly considered in the paper), the probability to achieve an
# registration at interim (conditional or final) and staying listed. 
# In this case, z_er_i must be set less than infinity. 
powerTwiceSucc <- function(relativeI1, relativeI_2Min, xi, Ared, z_er_i, z_f_i,
                           z_s_i, stabConstIntegral) {
  noncenP <- sqrt(relativeI1)*eta_f
  lower <- max(z_f_i, z_s_i, -stabConstIntegral+noncenP) # shifted stability
  # constants due to normal distribution with non centrality parameter noncenP
  upper <- min(z_er_i, stabConstIntegral+noncenP)
  if(lower < upper){
    powerTwiceSuccPart1 <-
      integrate(Vectorize(integrandPower(xi, relativeI1, relativeI_2Min, Ared)),
                lower = lower,
                upper = upper)
    powerTwiceSuccValuePart1 <- powerTwiceSuccPart1$value  
  }else{
    powerTwiceSuccValuePart1 <- 0    
  }
  
  if (z_er_i < Inf) { # z_er_i == Inf without early rejection as in the paper
    powerTwiceSuccValuePart2 <- 1-pnorm(z_er_i-sqrt(relativeI1)*eta_f)

  }else{
    powerTwiceSuccValuePart2 <- 0
  }
  powerTwiceSuccValue <- powerTwiceSuccValuePart1 + powerTwiceSuccValuePart2
  return(powerTwiceSuccValue)
}

# function gives the integrand for the calculation of the integral in
# the powerTwiceSucc-function:
integrandPower <- function(xi, relativeI1, relativeI_2Min, A) {
  
  relativeI_2A <- function(z1) { # calculation of second stage information
    relativeSampSizeFormZTest <-
      relativeI1 * (qnorm(1-beta) + qnorm(1 - A(z1))) ^ 2 /
      (max(z1, (sqrt(relativeI1) * eta_f) / xi)) ^ 2
    # (sqrt(relativeI1) * eta_f) / xi = sqrt(I_1)delta_rel
    return(max(relativeI_2Min, relativeSampSizeFormZTest))
  }
  
  integrand <- function(z1) {
    cp_delta <- # conditional power at delta
      (1 - pnorm(qnorm(1 - A(z1)) - eta_f * sqrt(relativeI_2A(z1))))
    return(cp_delta * dnorm(z1 - sqrt(relativeI1) * eta_f))
    # sqrt(relativeI1) * eta_f = sqrt(I_1)*delta
  }
  return(integrand)
}

############################ main function of PART 1 ##########################

# the function returns a vector relativeI_2Min of relative informations 
# I_2,min (relative to I_delta) where relativeI_2Min[i] corresponds 
# to the relative first-stage information relativeI1vec[i]. 
# If relativeI_2Min[i] == NaN, the search interval, i.e. the variable 
# uppBoundSearch_I_2Min, should be (significantly) increased. If this 
# does not help, z_f[i] or z_s[i] are chosen too large.
determinationI_2Min <- function(relativeI1vec,
                              p_1,
                              xi,
                              A,
                              c_A,
                              z_er,
                              z_f,
                              z_s,
                              uppBoundSearch_I_2Min,
                              tolBisecDeterm_I_2Min,
                              stabConstIntegral,
                              relativeI_2Const=NULL,
                              iNMethodChangingWeights=FALSE) {
  numI1 <- length(relativeI1vec)
  relativeI_2Min <- numeric(numI1) # Initializing output vector

# process the parameters:
  if(length(p_1)== 1){
    p_1 <- rep(p_1, numI1)
  }
  if (length(c_A) == 1) {
    c_A <- rep(c_A, numI1)
  }
  if (length(z_er) == 1) {
    z_er <- rep(z_er, numI1)
  }
  if (length(z_s) == 1) {
    z_s <- rep(z_s, numI1)
  }


  for (i in 1:numI1) {
    searchInterval <- c(0, uppBoundSearch_I_2Min)
    relativeI1 <- relativeI1vec[i]

    # set the conditional error function corresponding to i or relativeI1vec[i]:
    if(iNMethodChangingWeights){ # conditional error function depends three
      # parameters, this is only the case for A_Z 
      # in the second part of the paper
      Ared <- function(z1) {
        return(A(z1, c_A[i], relativeI1, relativeI_2Const[i]))
      }
    }else{ # conditional error function only depends on the level constant, this 
      # is the case for all considered conditional error functions except A_Z
      Ared <- function(z1) {
        return(A(z1, c_A[i]))
      }
    }
    
    # check whether without a minimum information the desired power
    # is at least achieved:
    powerTwiceSuccNoMinInform<-
      powerTwiceSucc(relativeI1, searchInterval[1], xi, Ared,
                            z_er[i], z_f[i], z_s[i], stabConstIntegral)

    if (powerTwiceSuccNoMinInform >= p_1[i]) { # no minimum information for
      # the second stage is needed
      relativeI_2Min[i] <- 0
    } else{
      # check whether search interval is large enough:
      powerTwiceSuccGreatestConsideredInform <-
        powerTwiceSucc(relativeI1, searchInterval[2], xi, Ared, z_er[i],
                       z_f[i], z_s[i], stabConstIntegral)

      if (powerTwiceSuccGreatestConsideredInform < p_1[i]) {
        relativeI_2Min[i] <- NaN # If this appears in the output,
        # the search interval, i.e. the variable uppBoundSearch_I_2Min should be
        # increased. 
      }else{
        diff <- searchInterval[2] - searchInterval[1]

        # a bisection search is used to find/ approximate I_2,min:
        while (diff > tolBisecDeterm_I_2Min) {
          midSearchInterv <- sum(searchInterval)/2
          powermMidSearchInterv <- powerTwiceSucc(relativeI1,
                                                       midSearchInterv,
                                                       xi,
                                                       Ared,
                                                       z_er[i],
                                                       z_f[i],
                                                       z_s[i],
                                                       stabConstIntegral)

          if (powermMidSearchInterv >= p_1[i]) {
            searchInterval[2] <- midSearchInterv
          } else{
            searchInterval[1] <- midSearchInterv
          }
          diff <- searchInterval[2] - searchInterval[1]
        }
        relativeI_2Min[i] <- searchInterval[2] # An upper approximation of
        # the minimum needed information for the second stage is returned
      }
    }
  }
  return(relativeI_2Min)
}

###############################################################################
###############################################################################
###############################################################################
# PART 2: functions to determine the constant information I_2,const for the 
# second stage for the case that no conditional registration is achieved
# (i.e. Z_1<z_f). This is only relevant for the considered scenarios in the
# paper without required registration.
###############################################################################

############################ helper function of PART 2: ######################

# function gives the integrand for the calculation of the integrals in
# the determination_I_2Const-function:

integrandPowerFixedNoncenP <- function(relativeI1, A, noncenP) {
  integrand <- function(z1) {
    cp_delta <-
      (1 - pnorm(qnorm(1 - A(z1)) - noncenP))
    return(cp_delta * dnorm(z1 - sqrt(relativeI1) * eta_f))
  }
  return(integrand)
}

######################### main function of PART 2: ###########################

# the function returns a vector relativeI_2Const of relative informations 
# I_2,const (relative to I_delta) where relativeI_2Const[i] corresponds 
# to the relative first-stage information relativeI1vec[i]. 
# If relativeI_2Const[i] == NaN, the search interval, i.e. the variable 
# uppBoundSearch_I_2Const, should be increased. 
# The function is used to find I_2,const/I_delta such that the equation
# P_delta(Z_2 >= qnorm(1-A(z_1)) | z_s<= Z_1 < z_f)>= 1-beta is fulfilled.
# This is equivalent to 
# P_delta(Z_2 >= qnorm(1-A(z_1)), z_s<= Z_1 < z_f)
# >= (1-beta)* P_delta(z_s<= Z_1 < z_f).
# The determination_I_2Const-function determines I_2,const such that
# this equation is fulfilled.
# In the paper it holds z_s=-Inf.

determination_I_2Const <- function(p_2, A, c_A, 
                                z_er, z_f, z_s, relativeI1vec, xi,
                                uppBoundSearch_I_2Const,
                                tolBisecDeterm_I_2Const, stabConstIntegral,
                                iNMethodChangingWeights=FALSE) {
  numI1 <- length(relativeI1vec)

  # process the parameters:
  if(length(p_2)== 1){
    p_1 <- rep(p_2, numI1)
  }
  if (length(c_A) == 1) {
    c_A <- rep(c_A, numI1)
  }
  if (length(z_er) == 1) {
    z_er <- rep(z_er, numI1)
  }
  if(length(z_s)==1){
    z_s <- rep(z_s, numI1)
  }
  uppBoundSearch_I_2Const <- ceiling(uppBoundSearch_I_2Const)
  relativeI_2Const <- numeric(numI1) # Initializing output vector

  # determination of I_2Const[i] for each relativeI1vec[i] considered:
  for(i in 1:numI1){
    relativeI1 <- relativeI1vec[i]
    
    if (z_s[i] >= z_f[i]) {
      if(p_2[i]>0){
        relativeI_2Const[i] <- "Condition not fulfilable"
      }else{
        relativeI_2Const[i] <- 0
      }
    }else{
      # set the current conditional error function:
      if (iNMethodChangingWeights){
        Ared <- function(z1) {
          return(A(z1, c_A[i], relativeI1, uppBoundSearch_I_2Const)) # for the inverse
          # normal method with weights chosen depending on I_2,const, the
          # conditional error function will change during the (bisection) search
          # for I_2,const.
        }
      }else{
        Ared <- function(z1) {
          return(A(z1, c_A[i]))
        }
      }
      # check whether numerical calculations are possible:
      noncenP <- sqrt(relativeI1)*eta_f
      lower <- max(z_s[i], -stabConstIntegral+noncenP)
      upper <- min(z_f[i], stabConstIntegral+noncenP)
      if(lower >= upper){
        relativeI_2Const[i] <- 0 # if lower >= upper, there are two possible
        # reasons. The first one is that z_s[i] >= z_f[i] (not the case in 
        # the paper). In this case, the fast track is not continued
        # if Z_1 < z_f[i]. Thus, relativeI_2Const[i] <- 0
        # The second possible reason is that z_s[i]>= stabConstIntegral+noncenP
        # or -stabConstIntegral+noncenP >= z_f[i]. In both cases this implies
        # that # P_delta(Z_2 >= qnorm(1-A(z_1)), z_s<= Z_1 < z_f) == 0 numerically
        # and P_delta(z_s<= Z_1 < z_f) == 0 numerically (if stabConstIntegral is 
        # chosen as described) regardless of
        # relativeI_2Const[i]. Thus, relativeI_2Const[i] <- 0 can be chosen.
      }else{
        # check whether the search interval (0,uppBoundSearch_I_2Const) is 
        # large enough: 
        noncenP <- sqrt(uppBoundSearch_I_2Const)*eta_f
        checkSearchInterval <- integrate(
          Vectorize(integrandPowerFixedNoncenP(relativeI1, Ared, noncenP)),
          lower = max(z_s[i], -stabConstIntegral+noncenP),
          upper = min(z_f[i], stabConstIntegral+noncenP)
        )$value
        
        if(checkSearchInterval < p_2[i]){
          relativeI_2Const[i] <- NaN # If this appears in the output,
          # the search interval, i.e. the variable uppBoundSearch_I_2Const, should 
          # be increased. 
        }else{
          searchInterval <- numeric(2)
          
          # first narrow the search range
          # (to have less iterations when uppBoundSearch_I_2Const is large):
          for (j in 1:uppBoundSearch_I_2Const) {
            # If the weights of the inverse normal method depend on the constant
            # information I_2,const, set the current conditional error function
            # as follows :
            if(iNMethodChangingWeights){
              Ared <- function(z1) {
                return(A(z1, c_A[i], relativeI1, j)) # j is current
                # considered I_2,const/I_delta, i.e. t_xi(I_2,const)
              }
            }
            noncenP <- sqrt(j)*eta_f
            PowerNoCondRegistrAndRejectH_0 <- integrate(
              Vectorize(integrandPowerFixedNoncenP(relativeI1, Ared, noncenP)),
              lower = max(z_s[i],-stabConstIntegral+noncenP),
              upper = min(z_f[i], stabConstIntegral+noncenP)
            )$value
            if(PowerNoCondRegistrAndRejectH_0 >= p_2[i]) {
              searchInterval <- c(j - 1, j)
              diff <- 1
              # search within the specified interval for 
              # the relative information I_2const/I_delta:
              while (diff > tolBisecDeterm_I_2Const) {
                midInterval <- sum(searchInterval) / 2
                # If the weights of the inverse normal method depend on the constant
                # information I_2,const, set the current conditional error function
                # as follows :
                if(iNMethodChangingWeights){
                  Ared <- function(z1) {
                    return(A(z1, c_A[i], relativeI1, midInterval))
                  }
                }
                noncenP <- sqrt(midInterval)*eta_f
                PowerNoCondRegistrAndRejectH_0 <- integrate(
                  Vectorize(integrandPowerFixedNoncenP(relativeI1, Ared, noncenP)),
                  lower = max(z_s[i],-stabConstIntegral+noncenP),
                  upper = min(z_f[i], stabConstIntegral+noncenP)
                )$value
                if(PowerNoCondRegistrAndRejectH_0 >= p_2[i]) {
                  searchInterval[2] <- midInterval
                } else{
                  searchInterval[1] <- midInterval
                }
                diff <- searchInterval[2]-searchInterval[1]
              }
              break
            }
          }
          relativeI_2Const[i] <- searchInterval[2] 
        }  
      }
    }
  }
  return(relativeI_2Const)
}
###############################################################################
###############################################################################
###############################################################################
# PART 3: functions to compare the efficiency of the different designs 
###############################################################################


############################ helper function of PART 3: ######################

# the function gives the integrand for the calculation of the mean information
# with the function investigateEfficiency 
integrandMeanInformSuccCondRegis <- function(xi, relativeI1, relativeI_2Min_i, A) {
  relativeInfo <- function(z1) {
    sampSizeFormZTest <- relativeI1 * (qnorm(1-beta) + qnorm(1 - A(z1))) ^ 2 /
      (max(z1, eta_f * sqrt(relativeI1) / xi)) ^ 2
    # (eta_f * sqrt(relativeI1) ) / xi = sqrt(I_1)*delta_rel
    return(max(relativeI_2Min_i, sampSizeFormZTest))
  }
  integrand <- function(z1) {
    return(relativeInfo(z1) * dnorm(z1 - sqrt(relativeI1) * eta_f))
  }
  return(integrand)
}

############################ main function of PART 3: ######################

# function returns a list containing the maximum information 
# (vector relativeMaxInform) and mean information (vector relativeMeanInform) 
# of the fast track procedure for each considered relativeI1vec[i].
# Additionally, t_xi(I_2(z_f)) is returned by the vector maxSamSizeFormulaZTest.

investigateEfficiency <- function(relativeI_2Min,
                      xi,
                      A,
                      c_A,
                      z_er,
                      z_f,
                      z_s,
                      relativeI1vec,
                      stabConstIntegral,
                      relativeI_2Const=NULL,
                      iNMethodChangingWeights=FALSE) {
  numI1 <- length(relativeI1vec)

  # for the output:
  relativeMaxInform <- numeric(numI1)
  maxSamSizeFormulaZTest <- numeric(numI1)
  relativeMeanInform <- numeric(numI1)

  # process parameters:
  if (length(c_A) == 1) {
    c_A <- rep(c_A, numI1)
  }
  if (length(z_er) == 1) {
    z_er <- rep(z_er, numI1)
  }
  if (length(z_s) == 1) {
    z_s <- rep(z_s, numI1)
  }

  for (i in 1:numI1) {
    relativeI1 <- relativeI1vec[i]

    # set the current conditional error function:
    if(iNMethodChangingWeights){
      Ared <- function(z1) {
        return(A(z1, c_A[i], relativeI1, relativeI_2Const[i]))
      }
    }else{
      Ared <- function(z1) {
        return(A(z1, c_A[i]))
      }
    }


    ############## determine the (relative) maximum information:  #############
      pointMaxSamSizeFormZTest <- max(z_f[i], z_s[i])
      maxSamSizeFormulaZTest[i] <- relativeI1 *
        (qnorm(1-beta) + qnorm(1 - Ared(pointMaxSamSizeFormZTest))) ^ 2 /
        (pointMaxSamSizeFormZTest) ^ 2
     if(is.null(relativeI_2Const)){ # means no continuation after non sucessful
       # conditional registration
       relativeMaxInform[i] <- max(relativeI_2Min[i], maxSamSizeFormulaZTest[i])  
     }else{
       relativeMaxInform[i] <- max(relativeI_2Min[i], maxSamSizeFormulaZTest[i],
                                   relativeI_2Const[i]) # The first two parts of
       # the maximum belong to the case of a successful conditional registration
       # and the third part belongs to the case of no successful conditional 
       # registration  
     }
      

    ############## determine the (relative) mean information:  ################

    if(is.null(relativeI_2Const)){ # no continuation of the 
      # fast-track procedure after unsuccessful conditional registration
      relMeanInformNoCondRegis <- 0
    }else{
      if(z_s[i] >= z_f[i]){ # no continuation of the 
        # fast-track procedure after unsuccessful conditional registration
        relMeanInformNoCondRegis <- 0 
      }else{
        relMeanInformNoCondRegis <- relativeI_2Const[i] * (pnorm(z_f[i] - sqrt(relativeI1) * eta_f) -
                                                             pnorm(z_s[i] - sqrt(relativeI1) * eta_f)) 
      }
    }
      noncenP <- sqrt(relativeI1)*eta_f
      lower <- max(z_f[i], z_s[i], -stabConstIntegral+noncenP)
      upper <- min(z_er[i], stabConstIntegral+noncenP)
      if(lower < upper){
        relMeanInformSuccCondRegis<- integrate(
          Vectorize(integrandMeanInformSuccCondRegis(xi, relativeI1, relativeI_2Min[i], Ared)),
          lower = lower,
          upper = upper
        )$value  
      }else{
        # if lower >= upper due to the assumptions on z_f, z_s, z_er, it must
        # hold max(z_f[i], z_s[i])>= stabConstIntegral+noncenP or
        # z_er[i], <= -stabConstIntegral+noncenP (not possible in the paper
        # because z_er[i]=Inf). In both cases the probability 
        # P_delta(max(z_f[i], z_s[i]) <= Z_1 < z_er[i])==0 numerically
        # (if stabConstIntegral is chosen as described). Thus,
        # this part of the mean information of the second stage can be set
        # equal to zero.
        relMeanInformSuccCondRegis <- 0  
      }
      
    relativeMeanInform[i] <- relMeanInformSuccCondRegis+ relMeanInformNoCondRegis
  }
  return(list("relativeMaxInform" = relativeMaxInform,
              "maxSamSizeFormulaZTest" = maxSamSizeFormulaZTest,
              "relativeMeanInform" = relativeMeanInform))
}

###############################################################################
###############################################################################
###############################################################################
# PART 4: functions to process input parameters for functions from PART 1 - 3

# the function is used to determine I_1,min/I_delta in case of required 
# early rejection
determine_I1min <- function(xi, beta, alpha_c){
  firstPart <- (qnorm(1-beta)/((xi-1)*eta_f))^2 
  secondPart <- ((qnorm(1-alpha_c)+qnorm(1-beta))/(xi*eta_f))^2
  return((max(firstPart, secondPart))*xi^2)
}

##############################################################################
# the function is used to determine I_2,max/I_delta
determine_I1max<- function(xi, alpha){
  return(qnorm(1-alpha)^2*xi^2/eta_f^2)
}

###############################################################################

# the function is used to calculate (1-beta)*P_{delta}(Z_1>=z_f) such that the 
# function determinationI_2Min can be used to determine I_2,min to achieve a 
# certain probability of success condition on Z_1>=z_f 
# (second part of the paper, if conditional registration is not required).

create_p <- function(relativeI1vec, z_f, z_s){
  numI1 <- length(relativeI1vec)
  if(length(z_s)==1){
    z_s <- rep(z_s, numI1)
  }
  pi_1 <- numeric(numI1)
  for(i in 1:numI1){
    probCondRegisSucessful <- 1 - pnorm(max(z_f[i], z_s[i])-sqrt(relativeI1vec[i]) *eta_f)
    pi_1[i] <- probCondRegisSucessful*(1-beta)
  }
  return(pi_1)
}
###############################################################################

# The function is used to determine the level constants for the inverse normal 
# method with fixed weights and Fisher's product test. This is necessary because
# the conditional error functions are truncated with 0.5 and may also have an 
# early futility stop.
# Determining the level constant is important for achieving the Type I error 
# rate alpha. The vector levelConst_vec is returned.
# If NaN is contained in the returned vector, try changing 
# searchIntervalDeterLevelConst.
# If NA==levelConst_vec[i], z_s[i] and z_er[i] are not sensibly chosen.

determineLevelConst <- function(A, relativeI1vec, z_s, z_er,
                               searchIntervalDeterLevelConst,
                               tolBisecDetermLevelConst,
                               stabConstIntegral){
  numI1 <- length(relativeI1vec)
  
  # process parameters:
  if(length(z_s)==1){
    z_s <- rep(z_s, numI1)
  }
  if(length(z_er)==1){
    z_er <- rep(z_er, numI1)
  }
  
  levelConst_vec <- numeric(numI1)
  for(i in 1:numI1){
    searchArea <- searchIntervalDeterLevelConst
    alpha_1 <- 1 - pnorm(z_er[i]-sqrt(relativeI1vec[i])*eta_f)
    
    # check whether numerical determination is possible:
    if(max(z_s[i], -stabConstIntegral) >= min(z_er[i], stabConstIntegral)){
      levelConst_vec[i] <- NA  
    }else{
      # check lower search bound:
      c <- searchArea[1]
      A_red <- function(z1){
        return(A(z1, c))
      }
      integrand <- function(z1){
        return(A_red(z1)*dnorm(z1))
      }
      checkLowerSearchB <- integrate(Vectorize(integrand),
                                     lower = max(z_s[i], -stabConstIntegral),
                                     upper = min(z_er[i], stabConstIntegral))$value
      if(checkLowerSearchB > alpha-alpha_1){
        levelConst_vec[i] <- NaN # if this appears in the output, 
        # the searchIntervalDeterLevelConst must be changed because the 
        # lower bound is too large
      }else{
        # check upper search bound:
        c <- searchArea[2]
        A_red <- function(z1){
          return(A(z1, c))
        }
        integrand <- function(z1){
          return(A_red(z1)*dnorm(z1))
        }
        checkUpperSearchB <- integrate(Vectorize(integrand),
                                       lower = max(z_s[i], -stabConstIntegral),
                                       upper = min(z_er[i], stabConstIntegral))$value
        if(checkUpperSearchB <= alpha-alpha_1){
          levelConst_vec[i] <- searchArea[2] # if searchArea[2] fulfills
          # the level condition, it should be returned.
        }else{ # bisection search to find level constant:
          diff <- searchArea[2]- searchArea[1]
          while(diff > tolBisecDetermLevelConst){
            c <- (searchArea[2]+ searchArea[1])/2
            A_red <- function(z1){
              return(A(z1, c))
            }
            integrand <- function(z1){
              return(A_red(z1)*dnorm(z1))
            }
            cer <- integrate(Vectorize(integrand),
                             lower = max(z_s[i], -stabConstIntegral),
                             upper = min(z_er[i], stabConstIntegral))$value
            if(cer <= alpha-alpha_1){
              searchArea[1] <- c
            }else{
              searchArea[2] <- c
            }
            diff <- searchArea[2]-searchArea[1]
          }
          levelConst_vec[i] <- searchArea[1]
        }  
      }
    }
  }
  return(levelConst_vec)
}




##############################################################################

# The function is used to determine the level constants (second part of the 
# paper) for the inverse normal method with changing weights depending on I_2,const. 
# This corresponds to the A_Z- conditional error function in the paper.
# The function gives the level constant for the upper branch. For the lower
# branch the level constant is alpha (or alpha-alpha_1 with early rejection).
# This determination of the level constant is necesarry because the conditional 
# error functions are truncated with 0.5 (and may also have an early futility
# stop, not considered explicitely in the paper).
# Determining the level constant is important for achieving the Type I error 
# rate alpha. 

determineLevelConstA_Z <- function(A_Z, relativeI1vec, relativeI2Const, xi, 
                                   z_s, z_f, z_er,
                                   tolBisecDetermLevelConst,
                                   stabConstIntegral){
  numI1 <- length(relativeI1vec)
  # process parameters:
  if(length(z_s)==1){
    z_s <- rep(z_s, numI1)
  }
  if(length(z_er)==1){
    z_er <- rep(z_er, numI1)
  }
  
  levelConst_vec <- numeric(numI1)
  
  for(i in 1:numI1){
    alpha_1 <- 1 - pnorm(z_er[i]-sqrt(relativeI1vec[i])*eta_f) # ==0 in the paper
    # because no early rejection
    A_Z_red_part1 <- function(z1){
      return(A_Z(z1, alpha-alpha_1, relativeI1vec[i], relativeI2Const[i]))
    }
    integrand <- function(z1){
      return(A_Z_red_part1(z1)* dnorm(z1))
    }
    lower <- max(z_s[i],-stabConstIntegral)
    upper <- min(z_f[i],stabConstIntegral)
    if(lower < upper){
      alpha_part1 <- integrate(Vectorize(integrand), lower=lower, 
                               upper=upper)$value  # Type I error rate of lower
      # branch
    }else{
      alpha_part1 <- 0  # if lower >= upper integral is numerically zero
    }
    ######### determine level constant for upper branch:##############
    searchArea <- c(alpha-alpha_1,1)
    
    # check upper limit of search interval:
    c <- searchArea[2]
    A_Z_red_part2 <- function(z1){
      return(A_Z(z1, c , relativeI1vec[i], relativeI2Const[i]))
    }
    integrand <- function(z1){
      return(A_Z_red_part2(z1)* dnorm(z1))
    }
    lower <- max(z_f[i], z_s[i], -stabConstIntegral)
    upper <- min(z_er[i], stabConstIntegral)
    if(lower < upper){
       cer <- integrate(Vectorize(integrand),
                     lower = max(z_f[i], z_s[i], -stabConstIntegral),
                     upper = min(z_er[i], stabConstIntegral))$value  
    }else{
      cer <- 0
    }
    if(cer <= alpha-alpha_part1-alpha_1){ # if greatest considered level 
      # constant fulfils level condition, it should be returned
      levelConst_vec[i] <- searchArea[2] 
    }else{ # bisection search for level constant upper branch:
      diff <- searchArea[2]-searchArea[1]
      while(diff > tolBisecDetermLevelConst){
        c <- (searchArea[2]+ searchArea[1])/2
        A_Z_red_part2 <- function(z1){
          return(A_Z(z1, c , relativeI1vec[i], relativeI2Const[i]))
        }
        integrand <- function(z1){
          return(A_Z_red_part2(z1)* dnorm(z1))
        }
        cer <- integrate(Vectorize(integrand),
                         lower = max(z_f[i], z_s[i], -stabConstIntegral),
                         upper = min(z_er[i], stabConstIntegral))$value
        if(cer <= alpha-alpha_part1-alpha_1){
          searchArea[1] <- c
        }else{
          searchArea[2] <- c
        }
        diff <- searchArea[2]-searchArea[1]
      }
      levelConst_vec[i] <- searchArea[1]  
    }
  }
  return(levelConst_vec)
}

###############################################################################
###############################################################################
###############################################################################
# PART 5: Creation of plots
###############################################################################

###################### common input parameters: #############################
alpha <- 0.025
beta <- 0.2
tol <- 1 / 10 ^ 3
eta_f <- qnorm(1 - beta) + qnorm(1 - alpha)
stabConstIntegral <- 38.56802# for greater values, dnorm becomes zero
uppBoundSearch_I_2Min <- 10
tolBisecDeterm_I_2Min <- 1/10^3
uppBoundSearch_I_2Const <- 10
tolBisecDeterm_I_2Const <- 1/10^3
tolBisecDetermLevelConst <- 1/10^6

#################### conditional error functions: #############################
# the level constants of the conditional error functions can depend on
# I_1/I_delta if z_f, z_s or z_er depend on I_1/I_delta.
# This is the case, for example, if the early futility limit depends on 
# I_1/I_delta as in the first part of the paper. The level constants
# are calculated and plugged into the following functions to
# specify the conditional error function.

w1 <- 1 / sqrt(2)
w2 <- 1 / sqrt(2)

A_IN <- function(z1, c) { # inverse normal method with fixed equal weights
  A_def <- 1 - pnorm((qnorm(1 - c) - w1 * z1) / w2)
  return(min(A_def, 0.5))
}

A_F <- function(z1, c) { # Fishers product test
  A_def <- c/(1-pnorm(z1))
  return(min(A_def, 0.5))
}

# conditional error function for design with separate studies:
A_C <- function(z1, c){
  return(alpha)
}


# Note that the conditional error function A_Z has two level constants,
# one for the 'lower branch' (i.e. Z_1 < z_f), which is given by alpha and one 
# for the 'upper branch' (i.e. Z_1 >= z_f), alpha', which must be calculated.
# In the following we use A_Z with c=alpha or c=alpha'.
# For the determinationI_2Min-function and investigateEfficiency-function
# c=alpha' must be set. To determine alpha', the  determineLevelConstA_Z-
# function can be used. This requires the specification of I_2,const/I_delta
# via the determination_I_2Const-function. For the 
# determination_I_2Const-function c=alpha must be set.

A_Z <- function(z1, c , relativeI1, relativeI2C){ # conditional error
  # function for second part of the paper: inverse normal method with weights
  # such that in case of a unsuccessful conditional registration,
  # a fixed size sample test is performed
  return(min(1-pnorm((qnorm(1-c)*sqrt(relativeI1+relativeI2C)-sqrt(relativeI1)*z1)/sqrt(relativeI2C)),0.5))
}


###############################################################################
##### Plots for first part of the paper: required conditional registration ####
###############################################################################

##############################################################################
#### xi = 2####################################################################
xi_1 = 2
alpha_c <- 0.15
rel_I1_min_xi_1<- determine_I1min(xi_1, beta, alpha_c)
rel_I1_max_xi_1 <- determine_I1max(xi_1, alpha)
relativeI1vec_1 <- seq(rel_I1_min_xi_1+1/10^6, rel_I1_max_xi_1 -1/10^6,
                       by = 0.01) 
numI1 <- length(relativeI1vec_1)

# set important boundaries: 
z_er1 <- Inf # no early rejection
z_f1 <- pmax(sqrt(relativeI1vec_1)*eta_f/xi_1, qnorm(1-0.15))
z_s1 <- z_f1 # due to the early futility bound

# determine level constants:
c_AIN1 <- determineLevelConst(A_IN, relativeI1vec_1, z_s1, z_er1,
                             searchIntervalDeterLevelConst= c(alpha,1),
                             tolBisecDetermLevelConst,
                             stabConstIntegral)
c_AF1 <- determineLevelConst(A_F, relativeI1vec_1, z_s1, z_er1,
                            searchIntervalDeterLevelConst= c(0.0038,1),
                            tolBisecDetermLevelConst,
                            stabConstIntegral)
c_AC1 <- rep(alpha,numI1)

# set further input parameter:
p_1_1 <- 0.8

#################### determination of I_2,min/I_delta: ########################

relativeI_2MinIN1 <- determinationI_2Min(relativeI1vec_1, p_1_1, xi_1, A_IN,
                                        c_AIN1, z_er1, z_f1, z_s1,
                                        uppBoundSearch_I_2Min, 
                                        tolBisecDeterm_I_2Min,
                                        stabConstIntegral)
relativeI_2MinClas1 <- determinationI_2Min(relativeI1vec_1, p_1_1, xi_1, A_C,
                                          c_AC1, z_er1, z_f1, z_s1,
                                          uppBoundSearch_I_2Min, 
                                          tolBisecDeterm_I_2Min,
                                          stabConstIntegral)
relativeI_2MinFisher1 <- determinationI_2Min(relativeI1vec_1, p_1_1, xi_1, A_F,
                                            c_AF1, z_er1, z_f1, z_s1,
                                            uppBoundSearch_I_2Min,
                                            tolBisecDeterm_I_2Min,
                                            stabConstIntegral)


########   Plot of I_2,min/I_delta, i.e. t_xi(I_2,min):   #####################

png("BindFutiliMinInform_2.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(3.5,1,0))
plot(relativeI1vec_1, relativeI_2MinClas1, type = 'l', col = "red", lty= 3, lwd=2,
     xlab = expression(t[xi](I[1])), ylab=expression(t[xi](I[2*","*min])),
     ylim=c(0,3), cex.lab=1.7, xaxt="n")
axis(1, at= seq(0.75,1.75, by=0.25))
axis(1, at=c(0.6704^2,4*0.4894), labels = c(expression(t[xi](I[1*","*min])),
                                            expression(t[xi](I[1*","*max]))),
     cex.axis=1.4)
points(relativeI1vec_1, relativeI_2MinFisher1, type = 'l', col = "blue", lty= 1,
       lwd=2)
points(relativeI1vec_1, relativeI_2MinIN1, type = 'l', col = "black", lty=2,
       lwd=2)
abline(h = 1, col = "black") 
legend(0.7, 3, legend=c("separate studies", "inverse normal method",
                        "Fisher's product test"),
       col=c("red", "black", "blue"), lty=c(3,2,1), cex=1.25, lwd=c(2,2,2),
       text.width = 0.75)
dev.off()

########## Comparison of efficiency between the considered designs designs #####

# efficiency of inverse normal method with fixed equal weights:
efficiencyResultIN1 <- investigateEfficiency(relativeI_2MinIN1, xi_1, A_IN, 
                                            c_AIN1, z_er1, z_f1, z_s1,
                                            relativeI1vec_1,
                                            stabConstIntegral)

relativeMeanInformIN1 <- efficiencyResultIN1$relativeMeanInform
relativeMaxInformIN1 <- efficiencyResultIN1$relativeMaxInform

# efficiency of Fisher's product test:
efficiencyResultFisher1 <- investigateEfficiency(relativeI_2MinFisher1, xi_1, A_F, 
                                                c_AF1, z_er1, z_f1, z_s1, 
                                                relativeI1vec_1, 
                                                stabConstIntegral)

relativeMeanInformFisher1 <- efficiencyResultFisher1$relativeMeanInform
relativeMaxInformFisher1 <- efficiencyResultFisher1$relativeMaxInform

# efficiency of design with two separate studies:
efficiencyResultClas1 <- investigateEfficiency(relativeI_2MinClas1, xi_1, A_C, 
                                              c_AC1, z_er1, z_f1, z_s1, 
                                              relativeI1vec_1,
                                              stabConstIntegral)
relativeMeanInformClas1 <- efficiencyResultClas1$relativeMeanInform
relativeMaxInformClas1 <- efficiencyResultClas1$relativeMaxInform

######### plot of relative mean information for the second stage: ###########

png("BindFutilMeanInformSeconStage_2.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(3.5,1,0))
plot(relativeI1vec_1, relativeMeanInformClas1, col = "red", type ='l', lty= 3,
     lwd=2, ylim=c(0,3), xlab = expression(t[xi](I[1])),
     ylab=expression(t[xi](E(I[2]))), cex.lab=1.7, cex.main=1.5, xaxt="n")
axis(1, at= seq(0.75,1.75, by=0.25))
axis(1, at=c(0.6704^2,4*0.4894), labels = c(expression(t[xi](I[1*","*min])),
                                            expression(t[xi](I[1*","*max]))),
     cex.axis=1.4)
points(relativeI1vec_1, relativeMeanInformFisher1, col = "blue", lty= 1, lwd=2,
       type ='l')
points(relativeI1vec_1, relativeMeanInformIN1, col = "black", type ='l', lty=2,
       lwd=2)
legend(0.7, 2.75, legend=c("separate studies", "inverse normal method",
                           "Fisher's product test"),
       col=c("red", "black", "blue"), lty=c(3,2,1), cex=1.25, lwd=c(2,2,2),
       text.width = 0.75)
abline(h = 1, col = "black") 
dev.off()


######### plot of relative maximum information for the second stage: ###########
png("BindFutilMaxInformSeconStage_2.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(3.5,1,0))
plot(relativeI1vec_1, relativeMaxInformClas1, col = "red", type ='l', lty= 3,
     lwd=2, ylim =c(0,5), xlab = expression(t[xi](I[1])), cex.lab=1.7,
     ylab=expression(t[xi](I[2*","*max])), xaxt="n")
axis(1, at= seq(0.75,1.75, by=0.25))
axis(1, at=c(0.6704^2,4*0.4894), labels = c(expression(t[xi](I[1*","*min])),
                                            expression(t[xi](I[1*","*max]))),
     cex.axis=1.4)
lines(relativeI1vec_1, relativeMaxInformFisher1, col = "blue", type='l', lty= 1,
       lwd=2)
points(relativeI1vec_1, relativeMaxInformIN1, col = "black", type='l',lty=2,
       lwd=2)
abline(h = 1, col = "black")
legend(1, 3.75, legend=c("separate studies", "inverse normal method",
                         "Fisher's product test"),
       col=c("red", "black", "blue"), lty=c(3,2,1), cex=1.2, lwd=c(2,2,2),
       text.width = 0.70)
dev.off()


######### plot of relative mean information over both stage: ###########
png("BindFutilMeanInformBothStages_2.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(3.5,1,0))
plot(relativeI1vec_1, relativeMeanInformClas1+relativeI1vec_1, col = "red",
     type ='l', lty= 3, lwd=2, ylim=c(0.75,3.75),
     xlab = expression(t[xi](I[1])), cex.lab=1.7,
     ylab=expression(t[xi](I[1])+t[xi](E(I[2]))), cex.main=1.5, xaxt="n")
axis(1, at= seq(0.75,1.75, by=0.25))
axis(1, at=c(0.6704^2,4*0.4894), labels = c(expression(t[xi](I[1*","*min])),
                                            expression(t[xi](I[1*","*max]))),
     cex.axis=1.4)
points(relativeI1vec_1, relativeMeanInformFisher1+relativeI1vec_1, col = "blue", 
       lty= 1, lwd=2, type ='l')
points(relativeI1vec_1, relativeMeanInformIN1+relativeI1vec_1, col = "black", 
       type ='l', lty=2, lwd=2)
legend(0.49, 3.75, legend=c("separate studies", "inverse normal method", 
                            "Fisher's product test"),
       col=c("red", "black", "blue"), lty=c(3,2,1), cex=1.2, lwd=c(2,2,2),
       text.width = 0.75)
abline(h = 1, col = "black") 
dev.off()


######### plot of relative maximum information over both stages: ###########
png("BindFutilMaxInformBothStages_2.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(3.5,1,0))
plot(relativeI1vec_1, relativeMaxInformClas1+relativeI1vec_1, col = "red",
     type ='l', lty= 3, lwd=2, ylim =c(1,6), xlab = expression(t[xi](I[1])),
     cex.lab=1.7, ylab=expression(t[xi](I[1])+t[xi](I[2*","*max])), xaxt="n")
axis(1, at= seq(0.75,1.75, by=0.25))
axis(1, at=c(0.6704^2,4*0.4894), labels = c(expression(t[xi](I[1*","*min])),
                                            expression(t[xi](I[1*","*max]))),
     cex.axis=1.4)
points(relativeI1vec_1, relativeMaxInformFisher1+relativeI1vec_1, col = "blue",
       lty= 1, lwd=2, type ='l')
points(relativeI1vec_1, relativeMaxInformIN1+relativeI1vec_1, col = "black",
       type ='l', lty=2, lwd=2)
legend(0.9, 4.75, legend=c("separate studies", "inverse normal method", 
                           "Fisher's product test"),
       col=c("red", "black", "blue"), lty=c(3,2,1), cex=1.25, lwd=c(2,2,2),
       text.width = 0.75)
dev.off()

################################################################################
########### xi_2=1.75 #########################################################
xi_2 = 1.75
alpha_c <- 0.15
rel_I1_min_xi_2<- determine_I1min(xi_2, beta, alpha_c)
rel_I1_max_xi_2 <- determine_I1max(xi_2, alpha)
relativeI1vec_2 <- seq(rel_I1_min_xi_2+1/10^6, rel_I1_max_xi_2-1/10^6,
                       by = 0.01) 
numI1_2 <- length(relativeI1vec_2)

# set important boundaries:
z_er2 <- Inf # no early rejection
z_f2 <- pmax(sqrt(relativeI1vec_2)*eta_f/xi_2, qnorm(1-0.15))
z_s2 <- z_f2

# determine level constants:
c_AIN2 <- determineLevelConst(A_IN, relativeI1vec_2, z_s2, z_er2,
                              searchIntervalDeterLevelConst= c(alpha,1),
                              tolBisecDetermLevelConst,
                              stabConstIntegral)
  
c_AF2 <- determineLevelConst(A_F, relativeI1vec_2, z_s2, z_er2,
                             searchIntervalDeterLevelConst= c(0.0038,1),
                             tolBisecDetermLevelConst,
                             stabConstIntegral)
c_AC2 <- rep(alpha,numI1_2)

# set further input parameter:
p_1_2 <- 0.8


#################### determination of I_2,min/I_delta: ########################
relativeI_2MinIN2 <- determinationI_2Min(relativeI1vec_2, p_1_2, xi_2, A_IN,
                                         c_AIN2, z_er2, z_f2, z_s2,
                                         uppBoundSearch_I_2Min,
                                         tolBisecDeterm_I_2Min,
                                         stabConstIntegral)
relativeI_2MinClas2 <- determinationI_2Min(relativeI1vec_2, p_1_2, xi_2, A_C,
                                           c_AC2, z_er2, z_f2, z_s2,
                                           uppBoundSearch_I_2Min,
                                           tolBisecDeterm_I_2Min,
                                           stabConstIntegral)
relativeI_2MinFisher2 <- determinationI_2Min(relativeI1vec_2, p_1_2, xi_2, A_F,
                                             c_AF2, z_er2, z_f2, z_s2,
                                             uppBoundSearch_I_2Min,
                                             tolBisecDeterm_I_2Min,
                                             stabConstIntegral)

########   Plot of I_2,min/I_delta, i.e. t_xi(I_2,min):   #####################

png("BindFutiliMinInform1.75.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(3.5,1,0))
plot(relativeI1vec_2, relativeI_2MinClas2, type = 'l', col = "red", lty= 3,
     lwd=2, xaxt="n", xlab = expression(t[xi](I[1])), 
     ylab=expression(t[xi](I[2*","*min])), ylim=c(0,3), cex.lab=1.7)
points(relativeI1vec_2, relativeI_2MinFisher2, type = 'l', col = "blue", lty= 1,
       lwd=2)
points(relativeI1vec_2, relativeI_2MinIN2, type = 'l', col = "black", lty=2,
       lwd=2)
axis(1, at= seq(0.75,1.25, by=0.25))
axis(1, at=c(0.7010^2,1.75^2*0.4894),labels = c(expression(t[xi](I[1*","*min])),
                                                expression(t[xi](I[1*","*max]))),
     cex.axis=1.4)
abline(h = 1, col = "black") 
legend(0.7, 3, legend=c("separate studies", "inverse normal method", 
                        "Fisher's product test"),
       col=c("red", "black", "blue"), lty=c(3,2,1), cex=1.25, lwd=c(2,2,2), 
       text.width = 0.5)
dev.off()


########## Comparison of efficiency between the considered designs designs #####

efficiencyResultIN2 <- investigateEfficiency(relativeI_2MinIN2, xi_2, A_IN,
                                             c_AIN2, z_er2, z_f2, z_s2,
                                             relativeI1vec_2, stabConstIntegral)
relativeMeanInformIN2 <- efficiencyResultIN2$relativeMeanInform
relativeMaxInformIN2 <- efficiencyResultIN2$relativeMaxInform

efficiencyResultFisher2 <- investigateEfficiency(relativeI_2MinFisher2, xi_2,
                                                 A_F, c_AF2, z_er2, z_f2, z_s2,
                                                 relativeI1vec_2,
                                                 stabConstIntegral)
relativeMeanInformFisher2 <- efficiencyResultFisher2$relativeMeanInform
relativeMaxInformFisher2 <- efficiencyResultFisher2$relativeMaxInform


efficiencyResultClas2 <- investigateEfficiency(relativeI_2MinClas2, xi_2, A_C,
                                               c_AC2, z_er2, z_f2, z_s2, 
                                               relativeI1vec_2,
                                               stabConstIntegral)
relativeMeanInformClas2 <- efficiencyResultClas2$relativeMeanInform
relativeMaxInformClas2 <- efficiencyResultClas2$relativeMaxInform

######### plot of relative mean information for the second stage: ###########

png("BindFutilMeanInformSeconStage_1.75.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(3.5,1,0))
plot(relativeI1vec_2, relativeMeanInformClas2, col = "red", type ='l', lty= 3, 
     lwd=2, ylim=c(0,3), xlab = expression(t[xi](I[1])), xaxt="n",
     ylab=expression(t[xi](E(I[2]))), cex.lab=1.7)
axis(1, at= seq(0.75,1.25, by=0.25))
axis(1, at=c(0.7010^2,1.75^2*0.4894),labels = c(expression(t[xi](I[1*","*min])),
                                                expression(t[xi](I[1*","*max]))),
     cex.axis=1.4)
points(relativeI1vec_2, relativeMeanInformFisher2, col = "blue", lty= 1,
       lwd=2, type ='l')
points(relativeI1vec_2, relativeMeanInformIN2, col = "black", type ='l',
       lty=2, lwd=2)
legend(0.7, 2.75, legend=c("separate studies", "inverse normal method",
                           "Fisher's product test"),
       col=c("red", "black", "blue"), lty=c(3,2,1), cex=1.25, lwd=c(2,2,2),
       text.width = 0.5)
abline(h = 1, col = "black")
dev.off()

######### plot of relative maximum information for the second stage: ###########

png("BindFutilMaxInformSeconStage_1.75.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(3.5,1,0))
plot(relativeI1vec_2, relativeMaxInformIN2, col = "black", type ='l', lty=2,
     lwd=2, ylim =c(0,5), xlab = expression(t[xi](I[1])), cex.lab=1.7, xaxt="n",
     ylab=expression(t[xi](I[2*","*max])))
axis(1, at= seq(0.75,1.25, by=0.25))
axis(1, at=c(0.7010^2,1.75^2*0.4894),labels = c(expression(t[xi](I[1*","*min])),
                                                expression(t[xi](I[1*","*max]))),
     cex.axis=1.4)
points(relativeI1vec_2, relativeMaxInformFisher2, col = "blue", lty= 1, lwd=2,
       type ='l')
points(relativeI1vec_2, relativeMaxInformClas2, col = "red", type ='l', lty= 3,
       lwd=2)
abline(h = 1, col = "black")
legend(0.7, 5, legend=c("separate studies", "inverse normal method",
                        "Fisher's product test"),
       col=c("red", "black", "blue"), lty=c(3,2,1), cex=1.25, lwd=c(2,2,2),
       text.width = 0.5)
dev.off()

############ plot of relative mean information over both stage: ###############

png("BindFutilMeanInformBothStages_1.75.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(3.5,1,0))
plot(relativeI1vec_2, relativeMeanInformClas2+relativeI1vec_2, col = "red", 
     type ='l', lty= 3, lwd=2,  ylim=c(0.75,3.75),
     xlab = expression(t[xi](I[1])), xaxt="n",
     ylab=expression(t[xi](I[1])+t[xi](E(I[2]))), cex.lab=1.7)
axis(1, at= seq(0.75,1.25, by=0.25))
axis(1, at=c(0.7010^2,1.75^2*0.4894),labels = c(expression(t[xi](I[1*","*min])),
                                                expression(t[xi](I[1*","*max]))),
     cex.axis=1.4)
points(relativeI1vec_2, relativeMeanInformFisher2+relativeI1vec_2, col = "blue",
       lty= 1, lwd=2, type ='l')
points(relativeI1vec_2, relativeMeanInformIN2+relativeI1vec_2, col = "black",
       type ='l', lty=2, lwd=2)
legend(0.6, 3.75, legend=c("separate studies", "inverse normal method",
                           "Fisher's product test"),
       col=c("red", "black", "blue"), lty=c(3,2,1), cex=1.25, lwd=c(2,2,2),
       text.width = 0.5)
abline(h = 1, col = "black") 
dev.off()

######### plot of relative maximum information over both stages: ###########

png("BindFutilMaxInformBothStages_1.75.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(3.5,1,0))
plot(relativeI1vec_2, relativeMaxInformIN2+relativeI1vec_2, col = "black",
     type ='l', lty=2, lwd=2, ylim =c(1.5,6),
     xlab = expression(t[xi](I[1])), cex.lab=1.7, xaxt="n",
     ylab=expression(t[xi](I[1])+t[xi](I[2*","*max])))
axis(1, at= seq(0.75,1.25, by=0.25))
axis(1, at=c(0.7010^2,1.75^2*0.4894), labels = c(expression(t[xi](I[1*","*min])),
                                                 expression(t[xi](I[1*","*max]))),
     cex.axis=1.4)
points(relativeI1vec_2, relativeMaxInformFisher2+relativeI1vec_2, col = "blue",
       lty= 1, lwd=2, type ='l')
points(relativeI1vec_2, relativeMaxInformClas2+relativeI1vec_2, col = "red",
       type ='l', lty= 3, lwd=2)
legend(0.55, 6, legend=c("separate studies", "inverse normal method",
                         "Fisher's product test"),
       col=c("red", "black", "blue"), lty=c(3,2,1), cex=1.25, lwd=c(2,2,2),
       text.width = 0.5)
dev.off()

###############################################################################
### Plots for second part of the paper:  no required conditional registration ##
###############################################################################


###############################################################################
###### xi = 1.25  #############################################################
xi_3 = 1.25
alpha_c <- 0.15
rel_I1_max_xi_3 <- determine_I1max(xi_3, alpha)
relativeI1vec_3 <- seq(1/10^6, rel_I1_max_xi_3 -1/10^6, by = 0.01)
numI1_3 <- length(relativeI1vec_3)

# set important boundaries:
z_er3 <- Inf # no early rejection
z_f3 <- pmax(sqrt(relativeI1vec_3)*eta_f/xi_3, qnorm(1-0.15))
z_s3 <- -Inf # no early futility stop

# determine level constants:
c_AIN3  <- determineLevelConst(A_IN, relativeI1vec_3[1], z_s3, z_er3,
                              searchIntervalDeterLevelConst= c(alpha,1),
                              tolBisecDetermLevelConst,
                              stabConstIntegral) # the level constants are the 
# same for all entries of relativeI1vec_3, since z_er and z_s do not depend 
# on I_1/I_delta

c_AF3 <- determineLevelConst(A_F, relativeI1vec_3[1], z_s3, z_er3,
                            searchIntervalDeterLevelConst= c(0.0038,1),
                            tolBisecDetermLevelConst, 
                            stabConstIntegral) # the level constants are the 
# same for all entries of relativeI1vec_3, since z_er and z_s do not depend 
# on I_1/I_delta

c_AC3 <- rep(alpha,numI1_3)

# Note that the conditional error function A_Z has two level constants,
# one for the 'lower branch' (z_1 < z_f), which is given by alpha and one 
# for the 'upper branch' (i.e. Z_1 >= z_f) alpha', which must be calculated. 
# To calculate alpha', I_2,const must already have been calculated. 
# This means that after calculating I_2,const/I_delta 
# using the determination_I_2Const-function,
# alpha' can be calculated using the determineLevelConstA_Z-function.

# set further input parameters:
p_1_3 <- create_p(relativeI1vec_3, z_f3, z_s3)
p_2_3 <- 0.8 - p_1_3

#################### determination of I_2,const/I_delta: ######################
relativeI_2ConstIN3 <- determination_I_2Const(p_2_3, A_IN, c_AIN3,
                                              z_er3, z_f3, z_s3,
                                              relativeI1vec_3, xi_3, 
                                              uppBoundSearch_I_2Const, 
                                              tolBisecDeterm_I_2Const,
                                              stabConstIntegral)
relativeI_2ConstFisher3 <- determination_I_2Const(p_2_3, A_F, c_AF3, 
                                                  z_er3, z_f3, z_s3,
                                                  relativeI1vec_3, xi_3,
                                                  uppBoundSearch_I_2Const,
                                                  tolBisecDeterm_I_2Const,
                                                  stabConstIntegral)

relativeI_2ConstClas3 <- rep(1, numI1_3) # follows from the definition
# of I_2,const/I_delta for the design with separate studies

relativeI_2ConstA_Z3 <- determination_I_2Const(p_2_3, A_Z, c_A=alpha,
                                               z_er3, z_f3, z_s3,
                                               relativeI1vec_3, xi_3,
                                               uppBoundSearch_I_2Const, 
                                               tolBisecDeterm_I_2Const, 
                                               stabConstIntegral,
                                               iNMethodChangingWeights=TRUE)
# note that c_A=alpha because the level constant of A_Z for the lower branch is
# alpha


#################### determination of I_2,min/I_delta: ########################
relativeI_2MinIN3 <- determinationI_2Min(relativeI1vec_3, p_1_3, xi_3, A_IN,
                                         c_AIN3, z_er3, z_f3, z_s3,
                                         uppBoundSearch_I_2Min,
                                         tolBisecDeterm_I_2Min,
                                         stabConstIntegral)
relativeI_2MinFisher3 <- determinationI_2Min(relativeI1vec_3, p_1_3, xi_3, A_F,
                                             c_AF3, z_er3, z_f3, z_s3,
                                             uppBoundSearch_I_2Min, 
                                             tolBisecDeterm_I_2Min,
                                             stabConstIntegral)
relativeI_2MinClas3 <- determinationI_2Min(relativeI1vec_3, p_1_3, xi_3, A_C, 
                                           c_AC3, z_er3, z_f3, z_s3,
                                           uppBoundSearch_I_2Min, 
                                           tolBisecDeterm_I_2Min, 
                                           stabConstIntegral)

c_AZ3 <- determineLevelConstA_Z(A_Z, relativeI1vec_3, relativeI_2ConstA_Z3, xi_3,
                                z_s3, z_f3, z_er3,
                                tolBisecDetermLevelConst,
                                stabConstIntegral)
relativeI_2MinA_Z3 <- determinationI_2Min(relativeI1vec_3, p_1_3, xi_3, A_Z, 
                                          c_AZ3, z_er3, z_f3, z_s3,
                                          uppBoundSearch_I_2Min,
                                          tolBisecDeterm_I_2Min,
                                          stabConstIntegral,
                                          relativeI_2ConstA_Z3,
                                          iNMethodChangingWeights=TRUE)
# Note that c_AZ3=alpha' because the upper branch (Z_1>=z_f) is considered

########   Plot of I_2,min/I_delta, i.e. t_xi(I_2,min):   #####################

png("NonBindingMinInformSecondStage_1.25.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(3.8,1,0))
plot(relativeI1vec_3, relativeI_2MinClas3, type = 'l', col = "red", lty=3,lwd=2,
     xlab = expression(t[xi](I[1])), ylab=expression(t[xi](I[2*","*min])),
     ylim=c(0,2.5), cex.lab=1.7, yaxt="n")
axis(1, at=c(rel_I1_max_xi_3 -1/10^6),
     labels = c(expression(t[xi](I[1*","*max]))), cex.axis=1.4)
axis(2, at = c(0, 0.5, 1, 2, 2.5), las=1)
axis(2, at=c(xi_3^2),
     labels = c(expression(t[xi](I[rel]))), cex.axis=1.4, las=1)
points(relativeI1vec_3, relativeI_2MinFisher3, type = 'l', col = "blue",
       lty=1, lwd=2)
points(relativeI1vec_3, relativeI_2MinIN3, type = 'l', col = "black",lty=2,
       lwd=2)
points(relativeI1vec_3, relativeI_2MinA_Z3, type = 'l', col = "purple",lty=5,
       lwd=2)
abline(v=0.1368, col="darkgrey", lwd=2, lty=1)
abline(h = 1, col = "black")
abline(h = 0.5, col = "black")
abline(h = 1.25^2, col = "black")
par(new=TRUE)
probCondRegistr <- 1-pnorm(pmax(sqrt(relativeI1vec_3)*eta_f/xi_3,
                                qnorm(1-0.15))-sqrt(relativeI1vec_3)*eta_f)
plot(relativeI1vec_3, probCondRegistr, type='l', lwd=2, lty=4, col="darkgreen",
     axes=FALSE, ylim=c(0,2.5), xlab="", ylab="")
axis(side=4, at = seq(0,1, by=0.1), col="darkgreen", col.axis="darkgreen",
     las=1.1)
mtext("Power conditional registration", side=4, line=3, col="darkgreen", 
      cex=1.3)
legend(0.2, 2.52,  legend=c("separate studies", "inverse normal method",
                            "Fisher's product test", expression(A[Z]),
                            expression(P[delta](Z[1]>=z[f])),
                            "gambling sep. studies"),
       col=c("red", "black", "blue", "purple", "darkgreen", "darkgrey"), 
       lty=c(3,2,1,5,4,1), cex=1.25, lwd=c(2,2,2,2,2,2), text.width = 0.4)
dev.off()
########   Plot of I_2,const/I_delta, i.e. t_xi(I_2,min):   ###################

png("NonBindingConstInform_1.25.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(4,1,0))
plot(relativeI1vec_3, relativeI_2ConstClas3, type = 'l', col = "red", lty=3,
     lwd=2, xlab = expression(t[xi](I[1])),
     ylab=expression(t[xi](I[2*","*const])),
     ylim=c(0,3), cex.lab=1.7, yaxt='n')
axis(1, at=c(rel_I1_max_xi_3 -1/10^6),
     labels = c(expression(t[xi](I[1*","*max]))), cex.axis=1.4)
axis(2, at= c(0,0.5,1,2,2.5,3), las=1)
axis(2, at=c(xi_3^2),
     labels = c(expression(t[xi](I[rel]))), cex.axis=1.4, las=1)
points(relativeI1vec_3, relativeI_2ConstFisher3, type = 'l', col = "blue",
       lty=1, lwd=2)
points(relativeI1vec_3, relativeI_2ConstIN3, type = 'l', col = "black",lty=2,
       lwd=2)
points(relativeI1vec_3, relativeI_2ConstA_Z3, type = 'l', col = "purple",lty=5,
       lwd=2)
abline(v=0.1368, col="darkgrey", lwd=2, lty=1)
abline(h = 1, col = "black")
abline(h = 0.5, col = "black")
abline(h = 1.25^2, col = "black")
par(new=TRUE)
probCondRegistr <- 1-pnorm(pmax(sqrt(relativeI1vec_3)*eta_f/xi_3,
                                qnorm(1-0.15))-sqrt(relativeI1vec_3)*eta_f)
plot(relativeI1vec_3, probCondRegistr, type='l', lwd=2, lty=4, col="darkgreen",
     axes=FALSE, ylim=c(0,3), xlab="", ylab="")
axis(side=4, at = seq(0,1, by=0.1), col="darkgreen", col.axis="darkgreen",
     las=1.1)
mtext("Power conditional registration", side=4, line=3, col="darkgreen", 
      cex=1.3)
legend(0.25, 3.1,  legend=c("separate studies", "inverse normal method",
                            "Fisher's product test", expression(A[Z]),
                            expression(P[delta](Z[1]>=z[f])),
                            "gambling sep. studies"),
       col=c("red", "black", "blue", "purple", "darkgreen", "darkgrey"), 
       lty=c(3,2,1,5,4,1), cex=1.25, lwd=c(2,2,2,2,2,2), text.width = 0.38)
dev.off()

########## Comparison of efficiency between the considered designs designs #####
efficiencyResultIN3 <- investigateEfficiency(relativeI_2MinIN3, xi_3, A_IN, 
                                             c_AIN3, z_er3, z_f3, z_s3,
                                             relativeI1vec_3, stabConstIntegral,
                                             relativeI_2ConstIN3)
relativeMeanInformIN3 <- efficiencyResultIN3$relativeMeanInform
relativeMaxInformIN3 <- efficiencyResultIN3$relativeMaxInform
maxSamSizeFormulaZTestIN3 <- efficiencyResultIN3$maxSamSizeFormulaZTest

efficiencyResultFisher3 <- investigateEfficiency(relativeI_2MinFisher3, xi_3, 
                                                 A_F, c_AF3, z_er3, z_f3, z_s3,
                                                 relativeI1vec_3,
                                                 stabConstIntegral,
                                                 relativeI_2ConstFisher3)
relativeMeanInformFisher3 <- efficiencyResultFisher3$relativeMeanInform
relativeMaxInformFisher3 <- efficiencyResultFisher3$relativeMaxInform
maxSamSizeFormulaZTestFisher3 <- efficiencyResultFisher3$maxSamSizeFormulaZTest


efficiencyResultClas3 <- investigateEfficiency(relativeI_2MinClas3, xi_3, A_C, 
                                               c_AC3, z_er3, z_f3, z_s3,
                                               relativeI1vec_3,
                                               stabConstIntegral,
                                               relativeI_2ConstClas3)
relativeMeanInformClas3 <- efficiencyResultClas3$relativeMeanInform
relativeMaxInformClas3 <- efficiencyResultClas3$relativeMaxInform
maxSamSizeFormulaZTestClas3 <- efficiencyResultClas3$maxSamSizeFormulaZTest


efficiencyResultA_Z3 <- investigateEfficiency(relativeI_2MinA_Z3, xi_3, A_Z, 
                                              c_AZ3, z_er3, z_f3, z_s3,
                                              relativeI1vec_3,
                                              stabConstIntegral,
                                              relativeI_2ConstA_Z3,
                                              iNMethodChangingWeights=TRUE)
relativeMeanInformA_Z3 <- efficiencyResultA_Z3$relativeMeanInform
relativeMaxInformA_Z3 <- efficiencyResultA_Z3$relativeMaxInform
maxSamSizeFormulaZTestA_Z3 <- efficiencyResultA_Z3$maxSamSizeFormulaZTest



######### plot of relative maximum information for the second stage: ###########
png("NonBindFutilMaxInformSeconStage_1.25.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(4.1,1,0))
plot(relativeI1vec_3, relativeMaxInformIN3, col = "black", type ='l', lty=2, 
     lwd=2, ylim=c(0,4.5), xlab = expression(t[xi](I[1])),  xaxt="n", yaxt="n",
     ylab=expression(t[xi](I[2*","*max])), cex.lab=1.7)
axis(1, at= c(0,0.2,0.4,0.6))
axis(1, at=c(0.76), labels = c(expression(t[xi](I[1*","*max]))), cex.axis=1.4)
axis(2, at= c(0,0.5,1,2,2.5,3.0,3.5,4,4.5), las=1)
axis(2, at=c(1.5625), labels = c(expression(t[xi](I[rel]))), las=1, cex.axis=1.4)
points(relativeI1vec_3, relativeMaxInformA_Z3, col = "purple", lty= 5, lwd=2, 
       type ='l')
points(relativeI1vec_3, relativeMaxInformFisher3, col = "blue", lty= 1, lwd=2, 
       type ='l')
abline(h = 1.25^2, col = "black")
abline(h = 1, col = "black") 
points(relativeI1vec_3, relativeMaxInformClas3, col = "red", type ='l', lty= 3, 
       lwd=2)
par(new=TRUE)
probCondRegistr <- 
  1-pnorm(pmax(sqrt(relativeI1vec_3)*eta_f/xi_3,qnorm(1-0.15))
          -sqrt(relativeI1vec_3)*eta_f)
plot(relativeI1vec_3, probCondRegistr, type='l', lwd=2, lty=4, col="darkgreen",
     axes=FALSE, ylim=c(0,2.25), xlab="", ylab="")
axis(side=4, at = seq(0,1, by=0.1), col="darkgreen", col.axis="darkgreen", 
     las=1)
mtext("Power conditional registration", side=4, line=3, col="darkgreen",
      cex=1.3)
abline(v=0.1368, col="darkgrey", lwd=2, lty=1)
legend(0.25, 2.3, legend=c("separate studies", "inverse normal method",
                           "Fisher's product test",  expression(A[Z]),
                           expression(P[delta](Z[1]>=z[f])),
                           "gambling sep. studies"),
       col=c("red", "black", "blue", "purple","darkgreen","darkgrey"),
       lty=c(3,2,1,5,4,1), cex=1.25, lwd=c(2,2,2,2,2,2), text.width = 0.38)
dev.off()

######### plot of relative mean information for the second stage: #############

png("NonBindingMeanInformSecondStage_1.25.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(4.1,1,0))
plot(relativeI1vec_3, relativeMeanInformClas3, col = "red", type ='l', lty= 3,
     lwd=2, ylim=c(0,3.6), xlab = expression(t[xi](I[1])),  xaxt="n", yaxt="n",
     ylab=expression(t[xi](E(I[2]))), cex.lab=1.7)
axis(1, at= c(0,0.2,0.4,0.6))
axis(1, at=c(0.76), labels = c(expression(t[xi](I[1*","*max]))), cex.axis=1.4)
axis(2, at= c(0,0.5,1.0,2,2.5,3.0, 3.5),
     las=1)
axis(2, at=c(1.5625), labels = c(expression(t[xi](I[rel]))), las=1, 
     cex.axis=1.4)
abline(h = 1.25^2, col = "black")
abline(h = 1, col = "black") 
points(relativeI1vec_3, relativeMeanInformFisher3, col = "blue", lty= 1, lwd=2,
       type ='l')
points(relativeI1vec_3, relativeMeanInformIN3, col = "black", type ='l', lty=2,
       lwd=2)
points(relativeI1vec_3, relativeMeanInformA_Z3, col = "purple", lty= 5, lwd=2,
       type ='l')
abline(v=0.1368, col="darkgrey", lwd=2, lty=1)
par(new=TRUE)
probCondRegistr <- 
  1-pnorm(pmax(sqrt(relativeI1vec_3)*eta_f/xi_3,qnorm(1-0.15))
          -sqrt(relativeI1vec_3)*eta_f)
plot(relativeI1vec_3, probCondRegistr, type='l', lwd=2, lty=4, col="darkgreen",
     axes=FALSE, ylim=c(0,1.8), xlab="", ylab="")
axis(side=4, at = seq(0,1, by=0.1), col="darkgreen", col.axis="darkgreen",
     las=1)
mtext("Power conditional registration", side=4, line=3, col="darkgreen",
      cex=1.3)
legend(0.25, 1.88, legend=c("separate studies", "inverse normal method",
                            "Fisher's product test", expression(A[Z]),
                            expression(P[delta](Z[1]>=z[f])),
                            "gambling sep. studies"),
       col=c("red", "black", "blue", "purple", "darkgreen", "darkgrey"),
       lty=c(3,2,1,5,4,1), cex=1.25, lwd=c(2,2,2,2,2,2), text.width = 0.38)
dev.off()

######### plot of relative maximum information over both stages: ##############

png("NonBindFutilMaxInformBothStages_1.25.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(4,1,0))
plot(relativeI1vec_3, relativeMaxInformIN3+relativeI1vec_3, col = "black",
     type ='l', lty=2, lwd=2, ylim=c(1,4.5), xlab = expression(t[xi](I[1])),
     xaxt="n", yaxt="n", ylab=expression(t[xi](I[1])+t[xi](I[2*","*max])),
     cex.lab=1.7)
axis(1, at= c(0,0.2,0.4,0.6))
axis(1, at=c(0.76), labels = c(expression(t[xi](I[1*","*max]))), cex.axis=1.4)
axis(2, at= c(1.0,2,2.5,3.0,3.5,4.0,4.5),las=1)
axis(2, at=c(1.5625), labels = c(expression(t[xi](I[rel]))), 
     cex.axis=1.4, las=1)
abline(h = 1.25^2, col = "black")
abline(h = 1, col = "black") 
points(relativeI1vec_3, relativeMaxInformA_Z3+relativeI1vec_3, col = "purple",
       lty= 5, lwd=2, type ='l')
points(relativeI1vec_3, relativeMaxInformFisher3+relativeI1vec_3, col = "blue",
       lty= 1, lwd=2, type ='l')
points(relativeI1vec_3, relativeMaxInformClas3+relativeI1vec_3, col = "red",
       type ='l', lty= 3, lwd=2)
abline(v=0.1368, col="darkgrey", lwd=2, lty=1)
par(new=TRUE)
probCondRegistr <- 
  1-pnorm(pmax(sqrt(relativeI1vec_3)*eta_f/xi_3,qnorm(1-0.15))
          -sqrt(relativeI1vec_3)*eta_f)
plot(relativeI1vec_3, probCondRegistr, type='l', lwd=2, lty=4, col="darkgreen",
     axes=FALSE, ylim=c(0,3.12), xlab="", ylab="")
axis(side=4, at = seq(0,1, by=0.1), col="darkgreen", col.axis="darkgreen", 
     las=1)
mtext("Power conditional registration", side=4, line=3, col="darkgreen",
      cex=1.3)
legend(0.25, 3.14, legend=c("separate studies", "inverse normal method",
                            "Fisher's product test", expression(A[Z]),
                            expression(P[delta](Z[1]>=z[f])),
                            "gambling sep. studies"),
       col=c("red", "black", "blue", "purple","darkgreen", "darkgrey"), 
       lty=c(3,2,1,5,4,1), cex=1.25, lwd=c(2,2,2,2,2,2), text.width = 0.38)
dev.off()


######### plot of relative mean information over both stages: #################
png("NonBindingMeanInformBothStages_1.25.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(4.1,1,0))
plot(relativeI1vec_3, relativeMeanInformClas3+relativeI1vec_3, col = "red",
     type ='l', xaxt="n", lty= 3, lwd=2, ylim=c(0.5,3.5),
     xlab = expression(t[xi](I[1])), cex.lab=1.7, yaxt="n",
     ylab=expression(t[xi](I[1])+t[xi](E(I[2]))))
axis(1, at= c(0,0.2,0.4,0.6))
axis(1, at=c(0.76), labels = c(expression(t[xi](I[1*","*max]))), cex.axis=1.4)
axis(2, at= c(0.5,1.0,2,2.5,3.0, 3.5), las=1)
axis(2, at=c(1.5625), labels = c(expression(t[xi](I[rel]))),  cex.axis=1.4, las=1)
abline(h = 1.25^2, col = "black")
abline(h = 1, col = "black") 
points(relativeI1vec_3, relativeMeanInformA_Z3+relativeI1vec_3, col = "purple",
       lty= 5, lwd=2, type ='l')
points(relativeI1vec_3, relativeMeanInformFisher3+relativeI1vec_3, col = "blue",
       lty= 1, lwd=2, type ='l')
points(relativeI1vec_3, relativeMeanInformIN3+relativeI1vec_3, col = "black",
       type ='l', lty=2, lwd=2)
abline(v=0.1368, col="darkgrey", lwd=2, lty=1)
probCondRegistr <- 
  1-pnorm(pmax(sqrt(relativeI1vec_3)*eta_f/xi_3,qnorm(1-0.15))
          -sqrt(relativeI1vec_3)*eta_f)
par(new=TRUE)
plot(relativeI1vec_3, probCondRegistr, type='l', lwd=2, lty=4, col="darkgreen",
     axes=FALSE, ylim=c(0,3), xlab="", ylab="")
axis(side=4, at = seq(0,1, by=0.1), col="darkgreen", col.axis="darkgreen", las=1)
mtext("Power conditional registration", side=4, line=3, col="darkgreen", cex=1.3)
legend(0.15, 3, legend=c("separate studies", "inverse normal method", 
                         "Fisher's product test", expression(A[Z]), 
                         expression(P[delta](Z[1]>=z[f])),"gambling sep. studies"),
       col=c("red", "black", "blue", "purple", "darkgreen", "darkgrey"),
       lty=c(3,2,1,5,4,1), cex=1.25, lwd=c(2,2,2,2,2), text.width = 0.38)
dev.off()

############# Detailed plot for each design: ##################################

###################### Inverse-Normal method: ##################################

png("NonBindingMaxInformationSplit_1.25InverseNormal.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(4,1,0))
plot(relativeI1vec_3, pmax(relativeI_2MinIN3,maxSamSizeFormulaZTestIN3),
     col = "blue", type ='l', lty=2, lwd=2, xaxt="n", yaxt="n",
     ylim =c(0,3.0), xlab = expression(t[xi](I[1])), cex.lab=1.7,
     ylab=expression("information on " *t[xi]*" -scale"))
axis(1, at= c(0,0.2,0.4,0.6))
axis(1, at=c(0.76), labels = c(expression(t[xi](I[1*","*max]))), cex.axis=1.4)
axis(2, at= c(0,0.5,1.0,2,2.5,3.0), las=1)
axis(2, at=c(1.5625), labels = c(expression(t[xi](I[rel]))), las=1,
     cex.axis=1.4)
points(relativeI1vec_3, relativeI_2MinIN3, col = "black", type ='l', lty=1,
       lwd=2)
points(relativeI1vec_3, relativeI_2ConstIN3, col = "red", type ='l',
     lty=3, lwd=2)
abline(h = 1.25^2, col = "black")
abline(h = 1, col = "black") 
probCondRegistr <- 
  1-pnorm(pmax(sqrt(relativeI1vec_3)*eta_f/xi_3,qnorm(1-0.15))
          -sqrt(relativeI1vec_3)*eta_f)
par(new=TRUE)
plot(relativeI1vec_3, probCondRegistr, type='l', lwd=2, lty=4, col="darkgreen",
     axes=FALSE, ylim=c(0,1.5), xlab="", ylab="")
axis(side=4, at = seq(0,1, by=0.1), col="darkgreen", col.axis="darkgreen",
     las=1)
mtext("Power conditional registration", side=4, line=3, col="darkgreen",
      cex=1.3)
legend(0.1, 1.5, legend=c(expression(t[xi](I[2*","*min])),
                          expression(t[xi](I[2](z[f]))),
                          expression(t[xi](I[2*", const"])),
                          expression(P[delta](Z[1]>=z[f]))),
       col=c("black", "blue", "red", "darkgreen"), lty=c(1,2,3,4), cex=1.25,
       lwd=c(2,2,2,2), text.width = 0.25)
dev.off()


#################### Fisher's product test: ###################################
png("NonBindingMaxInformationSplit_1.25Fisher.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(4,1,0))
plot(relativeI1vec_3, pmax(maxSamSizeFormulaZTestFisher3,relativeI_2MinFisher3),
     col = "blue", type ='l', lty=2, lwd=2, xaxt="n", yaxt="n",
      ylim =c(0,3.0), xlab = expression(t[xi](I[1])), cex.lab=1.7,
     ylab=expression("information on " * t[xi]*" -scale"))
axis(1, at= c(0,0.2,0.4,0.6))
axis(1, at=c(0.76), labels = c(expression(t[xi](I[1*","*max]))), cex.axis=1.4)
axis(2, at= c(0,0.5,1.0,2,2.5,3), las=1)
axis(2, at=c(1.5625), labels = c(expression(t[xi](I[rel]))), las=1, 
     cex.axis=1.4)
points(relativeI1vec_3, relativeI_2MinFisher3, col = "black", type ='l', lty=1,
       lwd=2)
points(relativeI1vec_3, relativeI_2ConstFisher3, col = "red", type ='l',
       lty=3, lwd=2)
abline(h = 1.25^2, col = "black")
abline(h = 1, col = "black") 
probCondRegistr <- 
  1-pnorm(pmax(sqrt(relativeI1vec_3)*eta_f/xi_3,qnorm(1-0.15))-
            sqrt(relativeI1vec_3)*eta_f)
par(new=TRUE)
plot(relativeI1vec_3, probCondRegistr, type='l', lwd=2, lty=4, col="darkgreen",
     axes=FALSE, ylim=c(0,1.5), xlab="", ylab="")
axis(side=4, at = seq(0,1, by=0.1), col="darkgreen", col.axis="darkgreen", 
     las=1)
mtext("Power conditional registration", side=4, line=3, col="darkgreen",
      cex=1.3)
legend(0.1, 1.5, legend=c(expression(t[xi](I[2*","*min])),
                          expression(t[xi](I[2](z[f]))),
                          expression(t[xi](I[2*", const"])),
                          expression(P[delta](Z[1]>=z[f]))),
       col=c("black", "blue", "red", "darkgreen"), lty=c(1,2,3,4), cex=1.25,
       lwd=c(2,2,2,2), text.width = 0.25)
dev.off()


################### design with two separate studies: #########################
png("NonBindingMaxInformationSplit_1.25NonAdaptive.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(4,1,0))
plot(relativeI1vec_3, pmax(maxSamSizeFormulaZTestClas3,relativeI_2MinClas3),
     col = "blue", type ='l',, lty=2, lwd=2, xaxt="n", yaxt="n",
      ylim =c(0,3.0), xlab = expression(t[xi](I[1])), cex.lab=1.7,
     ylab=expression("information on " * t[xi]*" -scale"))
axis(1, at= c(0,0.2,0.4,0.6))
axis(1, at=c(0.76), labels = c(expression(t[xi](I[1*","*max]))), cex.axis=1.4)
axis(2, at= c(0,0.5,1.0,2,2.5,3.0), las=1)
axis(2, at=c(1.5625), labels = c(expression(t[xi](I[rel]))), las=1,
     cex.axis=1.4)
points(relativeI1vec_3, relativeI_2MinClas3, col = "black", type ='l',
       lty=1, lwd=2)
abline(h = 1, col = "black") 
points(relativeI1vec_3, relativeI_2ConstClas3, col = "red", type ='l',
       lty=3, lwd=2)
abline(h = 1.25^2, col = "black")
abline(v=0.1368, col="darkgrey", lwd=2, lty=1)
probCondRegistr <- 
  1-pnorm(pmax(sqrt(relativeI1vec_3)*eta_f/xi_3,qnorm(1-0.15))
          -sqrt(relativeI1vec_3)*eta_f)
par(new=TRUE)
plot(relativeI1vec_3, probCondRegistr, type='l', lwd=2, lty=4, col="darkgreen",
     axes=FALSE, ylim=c(0,1.5), xlab="", ylab="")
axis(side=4, at = seq(0,1, by=0.1), col="darkgreen", col.axis="darkgreen",
     las=1)
mtext("Power conditional registration", side=4, line=3, col="darkgreen",
      cex=1.3)
legend(0.2, 1.6, legend=c(expression(t[xi](I[2*","*min])),
                          expression(t[xi](I[2](z[f]))),
                          expression(t[xi](I[2*", const"])),
                          expression(P[delta](Z[1]>=z[f])),
                          "gambling sep. studies"),
       col=c("black", "blue", "red", "darkgreen","darkgrey"), 
       lty=c(1,2,3,4,1), cex=1.25, lwd=c(2,2,2,2,2), text.width = 0.35)
dev.off()

############ inverse normal method with changing weights; A_Z: #################
png("NonBindingMaxInformationSplit_1.25INChangingWeights_A_Z.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(4,1,0))
plot(relativeI1vec_3, pmax(maxSamSizeFormulaZTestA_Z3,relativeI_2MinA_Z3),
     col = "blue", type ='l', lty=2, lwd=2, xaxt="n", yaxt="n",
     ylim =c(0,3.0), xlab = expression(t[xi](I[1])), cex.lab=1.7,
     ylab=expression("information on " *t[xi]*" -scale"))
axis(1, at= c(0,0.2,0.4,0.6))
axis(1, at=c(0.76), labels = c(expression(t[xi](I[1*","*max]))), cex.axis=1.4)
axis(2, at= c(0,0.5,1.0,2,2.5,3.0), las=1)
axis(2, at=c(1.5625), labels = c(expression(t[xi](I[rel]))), las=1,
     cex.axis=1.4)
points(relativeI1vec_3, relativeI_2MinA_Z3, col = "black", type ='l', lty=1,
       lwd=2)
points(relativeI1vec_3, relativeI_2ConstA_Z3, col = "red", type ='l',
       lty=3, lwd=2)
abline(h = 1.25^2, col = "black")
abline(h = 1, col = "black")
probCondRegistr <- 
  1-pnorm(pmax(sqrt(relativeI1vec_3)*eta_f/xi_3,qnorm(1-0.15))
          -sqrt(relativeI1vec_3)*eta_f)
par(new=TRUE)
plot(relativeI1vec_3, probCondRegistr, type='l', lwd=2, lty=4, col="darkgreen",
     axes=FALSE, ylim=c(0,1.5), xlab="", ylab="")
axis(side=4, at = seq(0,1, by=0.1), col="darkgreen", col.axis="darkgreen",
     las=1)
mtext("Power conditional registration", side=4, line=3, col="darkgreen", 
      cex=1.3)
legend(0.1, 1.5, legend=c(expression(t[xi](I[2*","*min])),
                          expression(t[xi](I[2](z[f]))) , 
                          expression(t[xi](I[2*", const"])),
                          expression(P[delta](Z[1]>=z[f]))),
       col=c("black", "blue", "red", "darkgreen"), lty=c(1,2,3,4), cex=1.25,
       lwd=c(2,2,2,2), text.width = 0.25)
dev.off()

###############################################################################
#### xi = xi_ min =1.43 ########################### ###########################
xi_4 = 1.43
alpha_c <- 0.15
rel_I1_max_xi_4 <- determine_I1max(xi_4, alpha)
relativeI1vec_4 <- seq(1/10^6, rel_I1_max_xi_4-1/10^6, by = 0.0001)
numI1_4 <- length(relativeI1vec_4)

# set important boundaries:
z_er4 <- Inf # no early rejection
z_f4 <- pmax(sqrt(relativeI1vec_4)*eta_f/xi_4, qnorm(1-0.15))
z_s4 <- -Inf # no early futility stop

# determine level constants:
c_AIN4  <- determineLevelConst(A_IN, relativeI1vec_4[1], z_s4, z_er4,
                               searchIntervalDeterLevelConst= c(alpha,1),
                               tolBisecDetermLevelConst,
                               stabConstIntegral) # the level constants are the 
# same for all entries of relativeI1vec_4, since z_er and z_s do not depend 
# on I_1/I_delta

c_AF4 <- determineLevelConst(A_F, relativeI1vec_4[1], z_s4, z_er4,
                             searchIntervalDeterLevelConst= c(0.0038,1),
                             tolBisecDetermLevelConst, 
                             stabConstIntegral) # the level constants are the 
# same for all entries of relativeI1vec_4, since z_er and z_s do not depend 
# on I_1/I_delta

c_AC4 <- rep(alpha,numI1_3)

# Note that the conditional error function A_Z has two level constants,
# one for the 'lower branch' (z_1 < z_f), which is given by alpha and one 
# for the 'upper branch' (i.e. Z_1 >= z_f), alpha', which must be calculated. 
# To calculate alpha', I_2,const must already have been calculated. 
# This means that after calculating I_2,const/I_delta using 
# the determination_I_2Const-function, alpha' can be calculated using the 
# determineLevelConstA_Z-function.

# set further input parameters:

p_1_4 <- create_p(relativeI1vec_4, z_f4, z_s4)
p_2_4 <- 0.8 - p_1_4

#################### determination of I_2,const/I_delta: ######################
relativeI_2ConstIN4 <- determination_I_2Const(p_2_4, A_IN, c_AIN4, 
                                              z_er4, z_f4, z_s4,
                                              relativeI1vec_4, xi_4,
                                              uppBoundSearch_I_2Const,
                                              tolBisecDeterm_I_2Const,
                                              stabConstIntegral)

relativeI_2ConstFisher4 <- determination_I_2Const(p_2_4, A_F, c_AF4, 
                                                  z_er4, z_f4, z_s4,
                                                  relativeI1vec_4, xi_4,
                                                  uppBoundSearch_I_2Const,
                                                  tolBisecDeterm_I_2Const,
                                                  stabConstIntegral)

relativeI_2ConstClas4 <- rep(1, numI1_3) # follows from the definition
# of I_2,const/I_delta for the design with separate studies

relativeI_2ConstA_Z4 <- determination_I_2Const(p_2_4, A_Z, c_A=alpha,
                                               z_er4, z_f4, z_s4,
                                               relativeI1vec_4, xi_4,
                                               uppBoundSearch_I_2Const, 
                                               tolBisecDeterm_I_2Const, 
                                               stabConstIntegral,
                                               iNMethodChangingWeights=TRUE)
# note that c_A=alpha because the level constant of A_Z for the lower branch is
# alpha


#################### determination of I_2,min/I_delta: ########################

relativeI_2MinIN4 <- determinationI_2Min(relativeI1vec_4, p_1_4, xi_4, A_IN,
                                         c_AIN4, z_er4, z_f4, z_s4,
                                         uppBoundSearch_I_2Min,
                                         tolBisecDeterm_I_2Min,
                                         stabConstIntegral)
relativeI_2MinFisher4 <- determinationI_2Min(relativeI1vec_4, p_1_4, xi_4, A_F,
                                             c_AF4, z_er4, z_f4, z_s4,
                                             uppBoundSearch_I_2Min, 
                                             tolBisecDeterm_I_2Min, 
                                             stabConstIntegral)
relativeI_2MinClas4 <- determinationI_2Min(relativeI1vec_4, p_1_4, xi_4, A_C,
                                           c_AC4, z_er4, z_f4, z_s4,
                                           uppBoundSearch_I_2Min,
                                           tolBisecDeterm_I_2Min, 
                                           stabConstIntegral)

c_AZ4 <- determineLevelConstA_Z(A_Z, relativeI1vec_4, relativeI_2ConstA_Z4,
                                xi_4, z_s4, z_f4, z_er4,
                                tolBisecDetermLevelConst,
                                stabConstIntegral)
relativeI_2MinA_Z4<- determinationI_2Min(relativeI1vec_4, p_1_4, xi_4, A_Z,
                                         c_AZ4, z_er4, z_f4, z_s4,
                                         uppBoundSearch_I_2Min, 
                                         tolBisecDeterm_I_2Min,
                                         stabConstIntegral,
                                         relativeI_2Const=relativeI_2ConstA_Z4, 
                                         iNMethodChangingWeights=TRUE)
# Note that c_AZ3=alpha' because the upper branch (Z_1>=z_f) is considered


########   Plot of I_2,min/I_delta, i.e. t_xi(I_2,min):   #####################

png("NonBindingMinInformSecondStage_1.43.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(3.5,1,0))
plot(relativeI1vec_4, relativeI_2MinClas4, type = 'l', col = "red", lty=3,lwd=2,
     xlab = expression(t[xi](I[1])), ylab=expression(t[xi](I[2*","*min])),
     ylim=c(0,2.5), cex.lab=1.7, yaxt="n", xaxt="n")
axis(1, at=c(0,0.2,0.4,0.6,0.8))
axis(1, at=c(rel_I1_max_xi_4 -1/10^6),
     labels = c(expression(t[xi](I[1*","*max]))), cex.axis=1.4)
axis(2, at = c(0, 0.5, 1, 1.5, 2.5), las=1)
axis(2, at=c(xi_4^2),
     labels = c(expression(t[xi](I[rel]))), cex.axis=1.4, las=1)
points(relativeI1vec_4, relativeI_2MinFisher4, type = 'l', col = "blue",
       lty=1, lwd=2)
points(relativeI1vec_4, relativeI_2MinIN4, type = 'l', col = "black",lty=2,
       lwd=2)
points(relativeI1vec_4, relativeI_2MinA_Z4, type = 'l', col = "purple",lty=5,
       lwd=2)
abline(v=0.1368, col="darkgrey", lwd=2, lty=1)
abline(h = 1, col = "black")
abline(h = 0.5, col = "black")
#abline(h = 1.43^2, col = "black")
par(new=TRUE)
probCondRegistr <- 1-pnorm(pmax(sqrt(relativeI1vec_4)*eta_f/xi_4,
                                qnorm(1-0.15))-sqrt(relativeI1vec_4)*eta_f)
plot(relativeI1vec_4, probCondRegistr, type='l', lwd=2, lty=4, col="darkgreen",
     axes=FALSE, ylim=c(0,2.5), xlab="", ylab="")
axis(side=4, at = seq(0,1, by=0.1), col="darkgreen", col.axis="darkgreen",
     las=1.1)
mtext("Power conditional registration", side=4, line=3, col="darkgreen", 
      cex=1.3)
legend(0.3, 2.52,  legend=c("separate studies", "inverse normal method",
                            "Fisher's product test", expression(A[Z]),
                            expression(P[delta](Z[1]>=z[f])),
                            "gambling sep. studies"),
       col=c("red", "black", "blue", "purple", "darkgreen", "darkgrey"), 
       lty=c(3,2,1,5,4,1), cex=1.25, lwd=c(2,2,2,2,2,2), text.width = 0.5)
dev.off()
########   Plot of I_2,const/I_delta, i.e. t_xi(I_2,min):   ###################

png("ConstInformNonBinding_1.43.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(4,1,0))
plot(relativeI1vec_4, relativeI_2ConstClas4, type = 'l', col = "red", lty=3,
     lwd=2, xlab = expression(t[xi](I[1])),
     ylab=expression(t[xi](I[2*","*const])),
     ylim=c(0,3), cex.lab=1.7, yaxt='n', xaxt='n')
axis(1, at=c(0,0.2,0.4,0.6,0.8))
axis(1, at=c(rel_I1_max_xi_4 -1/10^6),
     labels = c(expression(t[xi](I[1*","*max]))), cex.axis=1.4)
axis(2, at= c(0,0.5,1,1.5,2.5,3), las=1)
axis(2, at=c(xi_4^2),
     labels = c(expression(t[xi](I[rel]))), cex.axis=1.4, las=1)
points(relativeI1vec_4, relativeI_2ConstFisher4, type = 'l', col = "blue",
       lty=1, lwd=2)
points(relativeI1vec_4, relativeI_2ConstIN4, type = 'l', col = "black",lty=2,
       lwd=2)
points(relativeI1vec_4, relativeI_2ConstA_Z4, type = 'l', col = "purple",lty=5,
       lwd=2)
abline(v=0.1368, col="darkgrey", lwd=2, lty=1)
abline(h = 1, col = "black")
abline(h = 0.5, col = "black")
abline(h = xi_4^2, col = "black")
par(new=TRUE)
probCondRegistr <- 1-pnorm(pmax(sqrt(relativeI1vec_4)*eta_f/xi_4,
                                qnorm(1-0.15))-sqrt(relativeI1vec_4)*eta_f)
plot(relativeI1vec_4, probCondRegistr, type='l', lwd=2, lty=4, col="darkgreen",
     axes=FALSE, ylim=c(0,3), xlab="", ylab="")
axis(side=4, at = seq(0,1, by=0.1), col="darkgreen", col.axis="darkgreen",
     las=1.1)
mtext("Power conditional registration", side=4, line=3, col="darkgreen", 
      cex=1.3)
legend(0.35, 3.2,  legend=c("separate studies", "inverse normal method",
                            "Fisher's product test", expression(A[Z]),
                            expression(P[delta](Z[1]>=z[f])),
                            "gambling sep. studies"),
       col=c("red", "black", "blue", "purple", "darkgreen", "darkgrey"), 
       lty=c(3,2,1,5,4,1), cex=1.25, lwd=c(2,2,2,2,2,2), text.width = 0.5)
dev.off()


########## Comparison of efficiency between the considered designs designs #####
efficiencyResultIN4 <- investigateEfficiency(relativeI_2MinIN4, xi_4, A_IN, 
                                             c_AIN4, z_er4, z_f4, z_s4,
                                             relativeI1vec_4, stabConstIntegral,
                                             relativeI_2ConstIN4)
  
relativeMeanInformIN4 <- efficiencyResultIN4$relativeMeanInform
relativeMaxInformIN4 <- efficiencyResultIN4$relativeMaxInform
maxSamSizeFormulaZTestIN4 <- efficiencyResultIN4$maxSamSizeFormulaZTest

efficiencyResultFisher4 <- investigateEfficiency(relativeI_2MinFisher4, xi_4, 
                                                 A_F, c_AF4, z_er4, z_f4, z_s4,
                                                 relativeI1vec_4,
                                                 stabConstIntegral,
                                                 relativeI_2ConstFisher4)

relativeMeanInformFisher4 <- efficiencyResultFisher4$relativeMeanInform
relativeMaxInformFisher4 <- efficiencyResultFisher4$relativeMaxInform
maxSamSizeFormulaZTestFisher4<- efficiencyResultFisher4$maxSamSizeFormulaZTest


efficiencyResultClas4 <- investigateEfficiency(relativeI_2MinClas4, xi_4, A_C, 
                                               c_AC4, z_er4, z_f4, z_s4,
                                               relativeI1vec_4,
                                               stabConstIntegral,
                                               relativeI_2ConstClas4)

relativeMeanInformClas4 <- efficiencyResultClas4$relativeMeanInform
relativeMaxInformClas4 <- efficiencyResultClas4$relativeMaxInform
maxSamSizeFormulaZTestClas4 <- efficiencyResultClas4$maxSamSizeFormulaZTest

efficiencyResultA_Z4 <- investigateEfficiency(relativeI_2MinA_Z4, xi_4, A_Z, 
                                              c_AZ4, z_er4, z_f4, z_s4,
                                              relativeI1vec_4,
                                              stabConstIntegral,
                                              relativeI_2ConstA_Z4,
                                              iNMethodChangingWeights=TRUE)
relativeMeanInformA_Z4 <- efficiencyResultA_Z4$relativeMeanInform
relativeMaxInformA_Z4 <- efficiencyResultA_Z4$relativeMaxInform
maxSamSizeFormulaZTestA_Z4 <- efficiencyResultA_Z4$maxSamSizeFormulaZTest


######### plot of relative maximum information for the second stage: ##########
png("NonBindingMaxInformationSecondStage_1.43.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(4.1,1,0))
plot(relativeI1vec_4, relativeMaxInformIN4, col = "black", type ='l', lty=2, 
     lwd=2, ylim=c(0,4.5), xlab = expression(t[xi](I[1])), yaxt="n", xaxt="n",
     ylab=expression(t[xi](I[2*","*max])), cex.lab=1.7)
axis(1, at=c(1), labels = c(expression(t[xi](I[1*","*max]))), cex.axis=1.4)
axis(1, at=c(0,0.2,0.4,0.6,0.8))
axis(2, at= c(seq(0,1.75,by=0.5),c(0,0.5,1,1.5,2.5,3,3.5,4,4.5)), las=1)
axis(2, at=c(1.43^2), labels = c(expression(t[xi](I[rel]))), las=1, 
     cex.axis=1.4)
points(relativeI1vec_4, relativeMaxInformA_Z4, col = "purple", lty= 5, lwd=2, 
       type ='l')
abline(h = 1.43^2, col = "black")
abline(h = 1, col = "black") #
points(relativeI1vec_4, relativeMaxInformFisher4, col = "blue", lty= 1, lwd=2, 
       type ='l')
points(relativeI1vec_4, relativeMaxInformClas4, col = "red", type ='l', lty= 3,
       lwd=2)
abline(v=0.1428, col="darkgrey", lwd=2, lty=1)
par(new=TRUE)
probCondRegistr <-
  1-pnorm(pmax(sqrt(relativeI1vec_4)*eta_f/xi_4,qnorm(1-0.15))
           -sqrt(relativeI1vec_4)*eta_f)
plot(relativeI1vec_4, probCondRegistr, type='l', lwd=2, lty=4, col="darkgreen",
     axes=FALSE, ylim=c(0,2.25), xlab="", ylab="")
axis(side=4, at = seq(0,1, by=0.1), col="darkgreen", col.axis="darkgreen",
     las=1)
mtext("Power conditional registration", side=4, line=3, col="darkgreen",
      cex=1.3)
legend(0.25, 2.4, legend=c("separate studies", "inverse normal method",
                           "Fisher's product test", expression(A[Z]),
                           expression(P[delta](Z[1]>=z[f])),
                           "gambling sep. studies"),
       col=c("red", "black", "blue", "purple", "darkgreen","darkgrey"), 
       lty=c(3,2,1,5,4,1), cex=1.25, lwd=c(2,2,2,2,2,2), text.width = 0.5)
dev.off()


######### plot of relative mean information for the second stage: #############

png("NonBindingMeanInformSecondStage_1.43.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(4.1,1,0))
plot(relativeI1vec_4, relativeMeanInformClas4, col = "red", type ='l', lty= 3,
     lwd=2, ylim=c(0,3.5), xlab = expression(t[xi](I[1])),  yaxt="n", xaxt="n",
     ylab=expression(t[xi](E(I[2]))), cex.lab=1.7)
axis(1, at=c(0,0.2,0.4,0.6,0.8))
axis(1, at=c(1), labels = c(expression(t[xi](I[1*","*max]))), cex.axis=1.4)
axis(2, at= c(seq(0,1.5,by=0.5),c(2.5,3,3.5)), las=1)
axis(2, at=c(1.43^2), labels = c(expression(t[xi](I[rel]))), las=1,
     cex.axis=1.4)
abline(h = 1.43^2, col = "black")
abline(h = 1, col = "black") 
points(relativeI1vec_4, relativeMeanInformFisher4, col = "blue", lty= 1, lwd=2,
       type ='l')
points(relativeI1vec_4, relativeMeanInformIN4, col = "black", type ='l', lty=2,
       lwd=2)
points(relativeI1vec_4, relativeMeanInformA_Z4, col = "purple", lty= 5, lwd=2,
       type ='l')
abline(v=0.1428, col="darkgrey", lwd=2, lty=1)
probCondRegistr <- 
  1-pnorm(pmax(sqrt(relativeI1vec_4)*eta_f/xi_4,qnorm(1-0.15))
          -sqrt(relativeI1vec_4)*eta_f)
par(new=TRUE)
plot(relativeI1vec_4, probCondRegistr, type='l', lwd=2, lty=4, col="darkgreen",
     axes=FALSE, ylim=c(0,1.75), xlab="", ylab="")
axis(side=4, at = seq(0,1, by=0.1), col="darkgreen", col.axis="darkgreen",
     las=1)
mtext("Power conditional registration", side=4, line=3, col="darkgreen",
      cex=1.3)
legend(0.25, 1.88, legend=c("separate studies", "inverse normal method",
                            "Fisher's product test", expression(A[Z]),
                            expression(P[delta](Z[1]>=z[f])),
                            "gambling sep. studies"),
       col=c("red", "black", "blue", "purple", "darkgreen","darkgrey"),
       lty=c(3,2,1,5,4,1), cex=1.25, lwd=c(2,2,2,2,2,2), text.width = 0.5)
dev.off()

######### plot of relative maximum information over both stages: ##############

png("NonBindFutilMaxInformBothStages_1.43.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(4,1,0))
plot(relativeI1vec_4, relativeMaxInformIN4+relativeI1vec_4, col = "black", 
     type ='l', lty=2, lwd=2, ylim=c(1,4.5), xlab = expression(t[xi](I[1])),
     yaxt="n", xaxt="n", ylab=expression(t[xi](I[1])+t[xi](I[2*","*max])),
     cex.lab=1.7)
axis(1, at=c(0,0.2,0.4,0.6,0.8))
axis(1, at=c(1), labels = c(expression(t[xi](I[1*","*max]))), cex.axis=1.4)
axis(2, at= c(0,0.5,1.0,1.5,2.5,3.0,3.5,4,4.5), las=1)
axis(2, at=c(1.43^2), labels = c(expression(t[xi](I[rel]))), cex.axis=1.4, las=1)
abline(h = 1.43^2, col = "black")
abline(h = 1, col = "black") 
abline(v=0.1428, col="darkgrey", lwd=2)
points(relativeI1vec_4, relativeMaxInformFisher4+relativeI1vec_4, col = "blue",
       lty= 1, lwd=2, type ='l')
points(relativeI1vec_4, relativeMaxInformClas4+relativeI1vec_4, col = "red",
       type ='l', lty= 3, lwd=2)
points(relativeI1vec_4, relativeMaxInformA_Z4+relativeI1vec_4, col = "purple", 
       lty= 5, lwd=2,
       type ='l')
par(new=TRUE)
probCondRegistr <- 
  1-pnorm(pmax(sqrt(relativeI1vec_4)*eta_f/xi_4,qnorm(1-0.15))
          -sqrt(relativeI1vec_4)*eta_f)
plot(relativeI1vec_4, probCondRegistr, type='l', lwd=2, lty=4, col="darkgreen",
     axes=FALSE, ylim=c(0,1.68), xlab="", ylab="")
axis(side=4, at = seq(0,1, by=0.1), col="darkgreen", col.axis="darkgreen", las=1)
mtext("Power conditional registration", side=4, line=3, col="darkgreen", cex=1.3)
legend(0, 1.82, legend=c("separate studies", "inverse normal method",
                         "Fisher's product test", expression(A[Z]),
                         expression(P[delta](Z[1]>=z[f])),
                         "gambling sep. studies"),
       col=c("red", "black", "blue", "purple","darkgreen","darkgrey"),
       lty=c(3,2,1,5,4,1), cex=1.25, lwd=c(2,2,2,2,2,2), text.width = 0.5)
dev.off()


######### plot of relative mean information over both stages: ##############
png("NonBindingMeanInformBothStages_1.43.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(4,1,0))
plot(relativeI1vec_4, relativeMeanInformClas4+relativeI1vec_4, col = "red", 
     type ='l', xaxt="n", lty= 3, lwd=2, ylim=c(0.5,3.5), 
     xlab = expression(t[xi](I[1])), cex.lab=1.7, yaxt="n",
     ylab=expression(t[xi](I[1])+t[xi](E(I[2]))))
axis(1, at= c(0,0.2,0.4,0.6,0.8))
axis(1, at=c(1,1.43^2), labels = c(expression(t[xi](I[1*","*max])),
                                   expression(t[xi](I[rel]))),cex.axis=1.4)
axis(2, at= c(seq(0.5,1.5,by=0.5),c(2.5,3, 3.5)), las=1)
axis(2, at=c(1.43^2), labels = c(expression(t[xi](I[rel]))), las=1, cex.axis=1.4)
abline(h = 1.43^2, col = "black")
abline(h = 1, col = "black") 
abline(v=0.1428, col="darkgrey", lwd=2)
points(relativeI1vec_4, relativeMeanInformFisher4+relativeI1vec_4, col = "blue",
       lty= 1, lwd=2, type ='l')
points(relativeI1vec_4, relativeMeanInformIN4+relativeI1vec_4, col = "black",
       type ='l', lty=2, lwd=2)
points(relativeI1vec_4, relativeMeanInformA_Z4+relativeI1vec_4, col = "purple",
       lty= 5, lwd=2,
       type ='l')
probCondRegistr <- 
  1-pnorm(pmax(sqrt(relativeI1vec_4)*eta_f/xi_4,qnorm(1-0.15))
          -sqrt(relativeI1vec_4)*eta_f)
par(new=TRUE)
plot(relativeI1vec_4, probCondRegistr, type='l', lwd=2, lty=4, col="darkgreen",
     axes=FALSE, ylim=c(0,3), xlab="", ylab="")
axis(side=4, at = seq(0,1, by=0.1), col="darkgreen", col.axis="darkgreen", las=1)
mtext("Power conditional registration", side=4, line=3, col="darkgreen", cex=1.3)
legend(0, 3.2, legend=c("separate studies", "inverse normal method",
                        "Fisher's product test", expression(A[Z]),
                        expression(P[delta](Z[1]>=z[f])),
                        "gambling sep. studies"),
       col=c("red", "black", "blue", "purple","darkgreen","darkgrey"),
       lty=c(3,2,1,5,4,1), cex=1.25, lwd=c(2,2,2,2,2,2), text.width = 0.5)
dev.off()


############# Detailed plot for each design: ##################################

###################### Inverse-Normal method: ##################################
png("NonBindingMaxInformationSplit_1.43InverseNormal.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(4,1,0))
plot(relativeI1vec_4, pmax(maxSamSizeFormulaZTestIN4, relativeI_2MinIN4),
     col = "blue", type ='l', lty=2, lwd=2, ylim =c(0,4),
     xlab = expression(t[xi](I[1])), cex.lab=1.7,
     ylab=expression("information on " *t[xi]*" -scale"),xaxt="n", yaxt="n")
axis(1, at= c(0,0.2,0.4,0.6,0.8))
axis(1, at=c(1,1.43^2), labels = c(expression(t[xi](I[1*","*max])),
                                   expression(t[xi](I[rel]))), cex.axis=1.4)
axis(2, at= c(seq(0,1.5,by=0.5),c(2.5,3,3.5,4)), las=1)
axis(2, at=c(1.43^2), labels = c(expression(t[xi](I[rel]))), las=1, cex.axis=1.4)
points(relativeI1vec_4, relativeI_2MinIN4, col = "black", type ='l', 
       lty=1, lwd=2)
points(relativeI1vec_4, relativeI_2ConstIN4, col = "red", type ='l',
       lty=3, lwd=2)
abline(h = 1.43^2, col = "black")
abline(h = 1, col = "black") 
probCondRegistr <- 
  1-pnorm(pmax(sqrt(relativeI1vec_4)*eta_f/xi_4,qnorm(1-0.15))
          -sqrt(relativeI1vec_4)*eta_f)
par(new=TRUE)
plot(relativeI1vec_4, probCondRegistr, type='l', lwd=2, lty=4, col="darkgreen",
     axes=FALSE, ylim=c(0,2), xlab="", ylab="")
axis(side=4, at = seq(0,1, by=0.1), col="darkgreen", col.axis="darkgreen", las=1)
mtext("Power conditional registration", side=4, line=3, col="darkgreen", cex=1.3)
legend(0.02, 2, legend=c(expression(t[xi](I[2*","*min])),
                            expression(t[xi](I[2](z[f]))),
                            expression(t[xi](I[2*", const"])), 
                            expression(P[delta](Z[1]>=z[f]))),
       col=c("black", "blue", "red", "darkgreen"), lty=c(1,2,3,4), cex=1.25,
       lwd=c(2,2,2,2), text.width = 0.3)
dev.off()

##################### Fisher's product test: ##################################
png("NonBindingMaxInformationSplit_1.43Fisher.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(4,1,0))
plot(relativeI1vec_4, pmax(maxSamSizeFormulaZTestFisher4, relativeI_2MinFisher4), 
     col = "blue", type ='l', lty=2, lwd=2, ylim =c(0,4), 
     xlab = expression(t[xi](I[1])), cex.lab=1.7,
     ylab=expression("information on " * t[xi]*" -scale"),xaxt="n", yaxt="n")
axis(1, at= c(0,0.2,0.4,0.6,0.8))
axis(1, at= c(1,1.43^2), labels = c(expression(t[xi](I[1*","*max])),
                                    expression(t[xi](I[rel]))), cex.axis=1.4)
axis(2, at= c(seq(0,1.5,by=0.5),c(2.5,3, 3.5,4)), las=1)
axis(2, at=c(1.43^2), labels = c(expression(t[xi](I[rel]))), las=1, cex.axis=1.4)
points(relativeI1vec_4, relativeI_2MinFisher4, col = "black", type ='l', 
       lty=1, lwd=2)
points(relativeI1vec_4, relativeI_2ConstFisher4, col = "red", type ='l',
       lty=3, lwd=2)
abline(h = 1.43^2, col = "black")
abline(h = 1, col = "black") 
probCondRegistr <- 
  1-pnorm(pmax(sqrt(relativeI1vec_4)*eta_f/xi_4,qnorm(1-0.15))
          -sqrt(relativeI1vec_4)*eta_f)
par(new=TRUE)
plot(relativeI1vec_4, probCondRegistr, type='l', lwd=2, lty=4, col="darkgreen",
     axes=FALSE, ylim=c(0,2), xlab="", ylab="")
axis(side=4, at = seq(0,1, by=0.1), col="darkgreen", col.axis="darkgreen", las=1)
mtext("Power conditional registration", side=4, line=3, col="darkgreen", cex=1.3)
legend(0, 2, legend=c(expression(t[xi](I[2*","*min])),
                         expression(t[xi](I[2](z[f]))), 
                         expression(t[xi](I[2*", const"])), 
                         expression(P[delta](Z[1]>=z[f]))),
       col=c("black", "blue", "red", "darkgreen"), lty=c(1,2,3,4), cex=1.25,
       lwd=c(2,2,2,2), text.width = 0.3)
dev.off()


################### design with two separate studies: #########################
png("NonBindingMaxInformationSplit_1.43NonAdaptive.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(4,1,0))
plot(relativeI1vec_4, pmax(maxSamSizeFormulaZTestClas4,relativeI_2MinClas4),
     col = "blue", type ='l', lty=2, lwd=2, ylim =c(0,4),
     xlab = expression(t[xi](I[1])), cex.lab=1.7,
     ylab=expression("information on " * t[xi]*" -scale"), yaxt="n", xaxt="n")
axis(1, at= c(0,0.2,0.4,0.6,0.8))
axis(1, at= c(1,1.43^2), labels = c(expression(t[xi](I[1*","*max])),
                                    expression(t[xi](I[rel]))), cex.axis=1.4)
axis(2, at= c(seq(0,1.5,by=0.5),c(2.5,3, 3.5,4)), las=1)
axis(2, at=c(1.43^2), labels = c(expression(t[xi](I[rel]))), las=1, cex.axis=1.4)
points(relativeI1vec_4, relativeI_2MinClas4, col = "black", type ='l', xaxt="n",
       yaxt="n", lty=1, lwd=2)
abline(h = 1, col = "black") 
points(relativeI1vec_4, relativeI_2ConstClas4, col = "red", type ='l',
       lty=3, lwd=2)
abline(h = 1.43^2, col = "black")
abline(v=0.1428, col="darkgrey", lwd=2)
probCondRegistr <- 
  1-pnorm(pmax(sqrt(relativeI1vec_4)*eta_f/xi_4,qnorm(1-0.15))
          -sqrt(relativeI1vec_4)*eta_f)
par(new=TRUE)
plot(relativeI1vec_4, probCondRegistr, type='l', lwd=2, lty=4, col="darkgreen",
     axes=FALSE, ylim=c(0,2), xlab="", ylab="")
axis(side=4, at = seq(0,1, by=0.1), col="darkgreen", col.axis="darkgreen", las=1)
mtext("Power conditional registration", side=4, line=3, col="darkgreen", cex=1.3)
legend(0.25, 2.1, legend=c(expression(t[xi](I[2*","*min])), 
                           expression(t[xi](I[2](z[f]))), 
                           expression(t[xi](I[2*","* const])),
                           expression(P[delta](Z[1]>=z[f])),
                           "gambling sep. studies"),
       col=c("black", "blue", "red", "darkgreen", "darkgrey"), 
       lty=c(1,2,3,4,1), cex=1.25, lwd=c(2,2,2,2,2), text.width = 0.5)
dev.off()

############ inverse normal method with changing weights; A_Z: #################
png("NonBindingMaxInformationSplit_1.43INChangingWeights_A_Z.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(4,1,0))
plot(relativeI1vec_4, pmax(maxSamSizeFormulaZTestA_Z4,relativeI_2MinA_Z4),
     col = "blue", type ='l', lty=2, lwd=2, xaxt="n", yaxt="n",
     ylim =c(0,4), xlab = expression(t[xi](I[1])), cex.lab=1.7,
     ylab=expression("information on " *t[xi]*" -scale"))
axis(1, at= c(0,0.2,0.4,0.6,0.8))
axis(1, at=c(1), labels = c(expression(t[xi](I[1*","*max]))), cex.axis=1.4)
axis(2, at= c(0,0.5,1.0,1.5,2.5,3.0,3.5,4), las=1)
axis(2, at=c(1.43^2), labels = c(expression(t[xi](I[rel]))), las=1,
     cex.axis=1.4)
points(relativeI1vec_4, relativeI_2MinA_Z4, col = "black", type ='l', lty=1,
       lwd=2)
points(relativeI1vec_4, relativeI_2ConstA_Z4, col = "red", type ='l',
       lty=3, lwd=2)
abline(h = 1.43^2, col = "black")
abline(h = 1, col = "black") 
probCondRegistr <- 
  1-pnorm(pmax(sqrt(relativeI1vec_4)*eta_f/xi_4,qnorm(1-0.15))
          -sqrt(relativeI1vec_4)*eta_f)
par(new=TRUE)
plot(relativeI1vec_4, probCondRegistr, type='l', lwd=2, lty=4, col="darkgreen",
     axes=FALSE, ylim=c(0,2), xlab="", ylab="")
axis(side=4, at = seq(0,1, by=0.1), col="darkgreen", col.axis="darkgreen",
     las=1)
mtext("Power conditional registration", side=4, line=3, col="darkgreen", 
      cex=1.3)
legend(0.1, 2, legend=c(expression(t[xi](I[2*","*min])),
                          expression(t[xi](I[2](z[f]))) , 
                          expression(t[xi](I[2*", const"])),
                          expression(P[delta](Z[1]>=z[f]))),
       col=c("black", "blue", "red", "darkgreen"), lty=c(1,2,3,4), cex=1.25,
       lwd=c(2,2,2,2), text.width = 0.25)
dev.off()

###############################################################################
###############################################################################
# The following plots correspond to the case that at the time of the planning 
# of a fast-track procedure by the applicants, a successful conditional 
# registration is a requirement for permanent registration. 
# In contrast to the above plots for xi_1=2 and xi_2=1.75, the conditional error 
# functions are in the following not chosen in a way that there is a binding 
# futility stop after unsuccessful conditional registration.
# Instead, we deal with this case like an unbinding futility stop and 
# set z_s=-Inf. These are somewhat conservative designs with somewhat 
# larger sample sizes. However, such designs provide the opportunity to continue
# with a second stage for a permanent registration whenever there are good
# reasons for and sufficient resources have been achieved. 

########### xi_5=1.75 #########################################################
xi_5 = 1.75
alpha_c <- 0.15
rel_I1_min_xi_5<- determine_I1min(xi_5, beta, alpha_c)
rel_I1_max_xi_5 <- determine_I1max(xi_5, alpha)
relativeI1vec_5 <- seq(rel_I1_min_xi_5+1/10^6, rel_I1_max_xi_5-1/10^6,
                       by = 0.01) 
numI1_5 <- length(relativeI1vec_5)

# set important boundaries:
z_er5 <- Inf # no early rejection
z_f5 <- pmax(sqrt(relativeI1vec_5)*eta_f/xi_5, qnorm(1-0.15))
z_s5 <- -Inf

# determine level constants:
c_AIN5 <- determineLevelConst(A_IN, relativeI1vec_5[1], z_s5, z_er5,
                              searchIntervalDeterLevelConst= c(alpha,1),
                              tolBisecDetermLevelConst,
                              stabConstIntegral) # the level constants are the 
# same for all entries of relativeI1vec_5, since z_er and z_s do not depend 
# on I_1/I_delta


c_AF5 <- determineLevelConst(A_F, relativeI1vec_5[1], z_s5, z_er5,
                             searchIntervalDeterLevelConst= c(0.0038,1),
                             tolBisecDetermLevelConst,
                             stabConstIntegral) # the level constants are the 
# same for all entries of relativeI1vec_5, since z_er and z_s do not depend 
# on I_1/I_delta

c_AC5 <- rep(alpha,numI1_5)

# set further input parameter:
p_1_5 <- 0.8


#################### determination of I_2,min/I_delta: ########################
relativeI_2MinIN5 <- determinationI_2Min(relativeI1vec_5, p_1_5, xi_5, A_IN,
                                         c_AIN5, z_er5, z_f5, z_s5,
                                         uppBoundSearch_I_2Min,
                                         tolBisecDeterm_I_2Min,
                                         stabConstIntegral)
relativeI_2MinClas5 <- determinationI_2Min(relativeI1vec_5, p_1_5, xi_5, A_C,
                                           c_AC5, z_er5, z_f5, z_s5,
                                           uppBoundSearch_I_2Min,
                                           tolBisecDeterm_I_2Min,
                                           stabConstIntegral)
relativeI_2MinFisher5 <- determinationI_2Min(relativeI1vec_5, p_1_5, xi_5, A_F,
                                             c_AF5, z_er5, z_f5, z_s5,
                                             uppBoundSearch_I_2Min,
                                             tolBisecDeterm_I_2Min,
                                             stabConstIntegral)

########   Plot of I_2,min/I_delta, i.e. t_xi(I_2,min):   #####################

png("NonBindFutiliMinInform_1.75.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(3.5,1,0))
plot(relativeI1vec_5, relativeI_2MinClas5, type = 'l', col = "red", lty= 3,
     lwd=2, xaxt="n", xlab = expression(t[xi](I[1])), 
     ylab=expression(t[xi](I[2*","*min])), ylim=c(0,3), cex.lab=1.7)
points(relativeI1vec_5, relativeI_2MinFisher5, type = 'l', col = "blue", lty= 1,
       lwd=2)
points(relativeI1vec_5, relativeI_2MinIN5, type = 'l', col = "black", lty=2,
       lwd=2)
axis(1, at= seq(0.75,1.25, by=0.25))
axis(1, at=c(0.7010^2,1.75^2*0.4894),labels = c(expression(t[xi](I[1*","*min])),
                                                expression(t[xi](I[1*","*max]))),
     cex.axis=1.4)
abline(h = 1, col = "black") 
legend(0.7, 3, legend=c("separate studies", "inverse normal method", 
                        "Fisher's product test"),
       col=c("red", "black", "blue"), lty=c(3,2,1), cex=1.25, lwd=c(2,2,2), 
       text.width = 0.5)
dev.off()


########## Comparison of efficiency between the considered designs designs #####

efficiencyResultIN5 <- investigateEfficiency(relativeI_2MinIN5, xi_5, A_IN,
                                             c_AIN5, z_er5, z_f5, z_s5,
                                             relativeI1vec_5, stabConstIntegral)
relativeMeanInformIN5 <- efficiencyResultIN5$relativeMeanInform
relativeMaxInformIN5 <- efficiencyResultIN5$relativeMaxInform

efficiencyResultFisher5 <- investigateEfficiency(relativeI_2MinFisher5, xi_5,
                                                 A_F, c_AF5, z_er5, z_f5, z_s5,
                                                 relativeI1vec_5,
                                                 stabConstIntegral)
relativeMeanInformFisher5 <- efficiencyResultFisher5$relativeMeanInform
relativeMaxInformFisher5 <- efficiencyResultFisher5$relativeMaxInform


efficiencyResultClas5 <- investigateEfficiency(relativeI_2MinClas5, xi_5, A_C,
                                               c_AC5, z_er5, z_f5, z_s5, 
                                               relativeI1vec_5,
                                               stabConstIntegral)
relativeMeanInformClas5 <- efficiencyResultClas5$relativeMeanInform
relativeMaxInformClas5 <- efficiencyResultClas5$relativeMaxInform

######### plot of relative mean information for the second stage: ###########

png("NonBindFutilMeanInformSeconStage_1.75.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(3.5,1,0))
plot(relativeI1vec_5, relativeMeanInformClas5, col = "red", type ='l', lty= 3, 
     lwd=2, ylim=c(0,3), xlab = expression(t[xi](I[1])), xaxt="n",
     ylab=expression(t[xi](E(I[2]))), cex.lab=1.7)
axis(1, at= seq(0.75,1.25, by=0.25))
axis(1, at=c(0.7010^2,1.75^2*0.4894),labels = c(expression(t[xi](I[1*","*min])),
                                                expression(t[xi](I[1*","*max]))),
     cex.axis=1.4)
points(relativeI1vec_5, relativeMeanInformFisher5, col = "blue", lty= 1,
       lwd=2, type ='l')
points(relativeI1vec_5, relativeMeanInformIN5, col = "black", type ='l',
       lty=2, lwd=2)
legend(0.7, 2.75, legend=c("separate studies", "inverse normal method",
                           "Fisher's product test"),
       col=c("red", "black", "blue"), lty=c(3,2,1), cex=1.25, lwd=c(2,2,2),
       text.width = 0.5)
abline(h = 1, col = "black")
dev.off()

######### plot of relative maximum information for the second stage: ###########

png("NonBindFutilMaxInformSeconStage_1.75.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(3.5,1,0))
plot(relativeI1vec_5, relativeMaxInformIN5, col = "black", type ='l', lty=2,
     lwd=2, ylim =c(0,5), xlab = expression(t[xi](I[1])), cex.lab=1.7, xaxt="n",
     ylab=expression(t[xi](I[2*","*max])))
axis(1, at= seq(0.75,1.25, by=0.25))
axis(1, at=c(0.7010^2,1.75^2*0.4894),labels = c(expression(t[xi](I[1*","*min])),
                                                expression(t[xi](I[1*","*max]))),
     cex.axis=1.4)
points(relativeI1vec_5, relativeMaxInformFisher5, col = "blue", lty= 1, lwd=2,
       type ='l')
points(relativeI1vec_5, relativeMaxInformClas5, col = "red", type ='l', lty= 3,
       lwd=2)
abline(h = 1, col = "black")
legend(0.7, 5, legend=c("separate studies", "inverse normal method",
                        "Fisher's product test"),
       col=c("red", "black", "blue"), lty=c(3,2,1), cex=1.25, lwd=c(2,2,2),
       text.width = 0.5)
dev.off()

############ plot of relative mean information over both stage: ###############

png("NonBindFutilMeanInformBothStages_1.75.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(3.5,1,0))
plot(relativeI1vec_5, relativeMeanInformClas5+relativeI1vec_5, col = "red", 
     type ='l', lty= 3, lwd=2,  ylim=c(0.75,3.75),
     xlab = expression(t[xi](I[1])), xaxt="n",
     ylab=expression(t[xi](I[1])+t[xi](E(I[2]))), cex.lab=1.7)
axis(1, at= seq(0.75,1.25, by=0.25))
axis(1, at=c(0.7010^2,1.75^2*0.4894),labels = c(expression(t[xi](I[1*","*min])),
                                                expression(t[xi](I[1*","*max]))),
     cex.axis=1.4)
points(relativeI1vec_5, relativeMeanInformFisher5+relativeI1vec_5, col = "blue",
       lty= 1, lwd=2, type ='l')
points(relativeI1vec_5, relativeMeanInformIN5+relativeI1vec_5, col = "black",
       type ='l', lty=2, lwd=2)
legend(0.6, 3.75, legend=c("separate studies", "inverse normal method",
                           "Fisher's product test"),
       col=c("red", "black", "blue"), lty=c(3,2,1), cex=1.25, lwd=c(2,2,2),
       text.width = 0.5)
abline(h = 1, col = "black") 
dev.off()

######### plot of relative maximum information over both stages: ###########

png("NonBindFutilMaxInformBothStages_1.75.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(3.5,1,0))
plot(relativeI1vec_5, relativeMaxInformIN5+relativeI1vec_5, col = "black",
     type ='l', lty=2, lwd=2, ylim =c(1.5,6),
     xlab = expression(t[xi](I[1])), cex.lab=1.7, xaxt="n",
     ylab=expression(t[xi](I[1])+t[xi](I[2*","*max])))
axis(1, at= seq(0.75,1.25, by=0.25))
axis(1, at=c(0.7010^2,1.75^2*0.4894), labels = c(expression(t[xi](I[1*","*min])),
                                                 expression(t[xi](I[1*","*max]))),
     cex.axis=1.4)
points(relativeI1vec_5, relativeMaxInformFisher5+relativeI1vec_5, col = "blue",
       lty= 1, lwd=2, type ='l')
points(relativeI1vec_5, relativeMaxInformClas5+relativeI1vec_5, col = "red",
       type ='l', lty= 3, lwd=2)
legend(0.55, 6, legend=c("separate studies", "inverse normal method",
                         "Fisher's product test"),
       col=c("red", "black", "blue"), lty=c(3,2,1), cex=1.25, lwd=c(2,2,2),
       text.width = 0.5)
dev.off()

###############################################################################
#### xi = 2####################################################################
xi_6 = 2
alpha_c <- 0.15
rel_I1_min_xi_6<- determine_I1min(xi_6, beta, alpha_c)
rel_I1_max_xi_6 <- determine_I1max(xi_6, alpha)
relativeI1vec_6 <- seq(rel_I1_min_xi_6+1/10^6, rel_I1_max_xi_6 -1/10^6,
                       by = 0.01) 
numI6 <- length(relativeI1vec_6)

# set important boundaries: 
z_er6 <- Inf # no early rejection
z_f6 <- pmax(sqrt(relativeI1vec_6)*eta_f/xi_6, qnorm(1-0.15))
z_s6 <- -Inf

# determine level constants:
c_AIN6 <- determineLevelConst(A_IN, relativeI1vec_6[1], z_s6, z_er6,
                              searchIntervalDeterLevelConst= c(alpha,1),
                              tolBisecDetermLevelConst,
                              stabConstIntegral) # the level constants are the 
# same for all entries of relativeI1vec_6, since z_er and z_s do not depend 
# on I_1/I_delta


c_AF6 <- determineLevelConst(A_F, relativeI1vec_6[1], z_s6, z_er6,
                             searchIntervalDeterLevelConst= c(0.0038,1),
                             tolBisecDetermLevelConst,
                             stabConstIntegral) # the level constants are the 
# same for all entries of relativeI1vec_6, since z_er and z_s do not depend 
# on I_1/I_delta

c_AC6 <- rep(alpha,numI6)

# set further input parameter:
p_1_6 <- 0.8

#################### determination of I_2,min/I_delta: ########################

relativeI_2MinIN6 <- determinationI_2Min(relativeI1vec_6, p_1_6, xi_6, A_IN,
                                         c_AIN6, z_er6, z_f6, z_s6,
                                         uppBoundSearch_I_2Min, 
                                         tolBisecDeterm_I_2Min,
                                         stabConstIntegral)
relativeI_2MinClas6 <- determinationI_2Min(relativeI1vec_6, p_1_6, xi_6, A_C,
                                           c_AC6, z_er6, z_f6, z_s6,
                                           uppBoundSearch_I_2Min, 
                                           tolBisecDeterm_I_2Min,
                                           stabConstIntegral)
relativeI_2MinFisher6 <- determinationI_2Min(relativeI1vec_6, p_1_6, xi_6, A_F,
                                             c_AF6, z_er6, z_f6, z_s6,
                                             uppBoundSearch_I_2Min,
                                             tolBisecDeterm_I_2Min,
                                             stabConstIntegral)


########   Plot of I_2,min/I_delta, i.e. t_xi(I_2,min):   #####################

png("NonBindFutiliMinInform_2.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(3.5,1,0))
plot(relativeI1vec_6, relativeI_2MinClas6, type = 'l', col = "red", lty= 3, lwd=2,
     xlab = expression(t[xi](I[1])), ylab=expression(t[xi](I[2*","*min])),
     ylim=c(0,3), cex.lab=1.7, xaxt="n")
axis(1, at= seq(0.75,1.75, by=0.25))
axis(1, at=c(0.6704^2,4*0.4894), labels = c(expression(t[xi](I[1*","*min])),
                                            expression(t[xi](I[1*","*max]))),
     cex.axis=1.4)
points(relativeI1vec_6, relativeI_2MinFisher6, type = 'l', col = "blue", lty= 1,
       lwd=2)
points(relativeI1vec_6, relativeI_2MinIN6, type = 'l', col = "black", lty=2,
       lwd=2)
abline(h = 1, col = "black") 
legend(0.7, 3, legend=c("separate studies", "inverse normal method",
                        "Fisher's product test"),
       col=c("red", "black", "blue"), lty=c(3,2,1), cex=1.25, lwd=c(2,2,2),
       text.width = 0.75)
dev.off()

########## Comparison of efficiency between the considered designs designs #####

# efficiency of inverse normal method with fixed equal weights:
efficiencyResultIN6 <- investigateEfficiency(relativeI_2MinIN6, xi_6, A_IN, 
                                             c_AIN6, z_er6, z_f6, z_s6,
                                             relativeI1vec_6,
                                             stabConstIntegral)

relativeMeanInformIN6 <- efficiencyResultIN6$relativeMeanInform
relativeMaxInformIN6 <- efficiencyResultIN6$relativeMaxInform

# efficiency of Fisher's product test:
efficiencyResultFisher6 <- investigateEfficiency(relativeI_2MinFisher6, xi_6, A_F, 
                                                 c_AF6, z_er6, z_f6, z_s6, 
                                                 relativeI1vec_6, 
                                                 stabConstIntegral)

relativeMeanInformFisher6 <- efficiencyResultFisher6$relativeMeanInform
relativeMaxInformFisher6 <- efficiencyResultFisher6$relativeMaxInform

# efficiency of design with two separate studies:
efficiencyResultClas6 <- investigateEfficiency(relativeI_2MinClas6, xi_6, A_C, 
                                               c_AC6, z_er6, z_f6, z_s6, 
                                               relativeI1vec_6,
                                               stabConstIntegral)
relativeMeanInformClas6 <- efficiencyResultClas6$relativeMeanInform
relativeMaxInformClas6 <- efficiencyResultClas6$relativeMaxInform

######### plot of relative mean information for the second stage: ###########

png("NonBindFutilMeanInformSeconStage_2.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(3.5,1,0))
plot(relativeI1vec_6, relativeMeanInformClas6, col = "red", type ='l', lty= 3,
     lwd=2, ylim=c(0,3), xlab = expression(t[xi](I[1])),
     ylab=expression(t[xi](E(I[2]))), cex.lab=1.7, cex.main=1.5, xaxt="n")
axis(1, at= seq(0.75,1.75, by=0.25))
axis(1, at=c(0.6704^2,4*0.4894), labels = c(expression(t[xi](I[1*","*min])),
                                            expression(t[xi](I[1*","*max]))),
     cex.axis=1.4)
points(relativeI1vec_6, relativeMeanInformFisher6, col = "blue", lty= 1, lwd=2,
       type ='l')
points(relativeI1vec_6, relativeMeanInformIN6, col = "black", type ='l', lty=2,
       lwd=2)
legend(0.7, 2.75, legend=c("separate studies", "inverse normal method",
                           "Fisher's product test"),
       col=c("red", "black", "blue"), lty=c(3,2,1), cex=1.25, lwd=c(2,2,2),
       text.width = 0.75)
abline(h = 1, col = "black") 
dev.off()


######### plot of relative maximum information for the second stage: ###########
png("NonBindFutilMaxInformSeconStage_2.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(3.5,1,0))
plot(relativeI1vec_6, relativeMaxInformClas6, col = "red", type ='l', lty= 3,
     lwd=2, ylim =c(0,5), xlab = expression(t[xi](I[1])), cex.lab=1.7,
     ylab=expression(t[xi](I[2*","*max])), xaxt="n")
axis(1, at= seq(0.75,1.75, by=0.25))
axis(1, at=c(0.6704^2,4*0.4894), labels = c(expression(t[xi](I[1*","*min])),
                                            expression(t[xi](I[1*","*max]))),
     cex.axis=1.4)
lines(relativeI1vec_6, relativeMaxInformFisher6, col = "blue", type='l', lty= 1,
      lwd=2)
points(relativeI1vec_6, relativeMaxInformIN6, col = "black", type='l',lty=2,
       lwd=2)
abline(h = 1, col = "black")
legend(0.6, 1.5, legend=c("separate studies", "inverse normal method",
                         "Fisher's product test"),
       col=c("red", "black", "blue"), lty=c(3,2,1), cex=1.2, lwd=c(2,2,2),
       text.width = 0.70)
dev.off()


######### plot of relative mean information over both stage: ###########
png("NonBindFutilMeanInformBothStages_2.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(3.5,1,0))
plot(relativeI1vec_6, relativeMeanInformClas6+relativeI1vec_6, col = "red",
     type ='l', lty= 3, lwd=2, ylim=c(0.75,3.75),
     xlab = expression(t[xi](I[1])), cex.lab=1.7,
     ylab=expression(t[xi](I[1])+t[xi](E(I[2]))), cex.main=1.5, xaxt="n")
axis(1, at= seq(0.75,1.75, by=0.25))
axis(1, at=c(0.6704^2,4*0.4894), labels = c(expression(t[xi](I[1*","*min])),
                                            expression(t[xi](I[1*","*max]))),
     cex.axis=1.4)
points(relativeI1vec_6, relativeMeanInformFisher6+relativeI1vec_6, col = "blue", 
       lty= 1, lwd=2, type ='l')
points(relativeI1vec_6, relativeMeanInformIN6+relativeI1vec_6, col = "black", 
       type ='l', lty=2, lwd=2)
legend(0.49, 3.75, legend=c("separate studies", "inverse normal method", 
                            "Fisher's product test"),
       col=c("red", "black", "blue"), lty=c(3,2,1), cex=1.2, lwd=c(2,2,2),
       text.width = 0.75)
abline(h = 1, col = "black") 
dev.off()


######### plot of relative maximum information over both stages: ###########
png("NonBindFutilMaxInformBothStages_2.png", width = 2000, height = 1500, res = 300)
par(mar = c(5, 4+2, 4, 2+2), mgp=c(3.5,1,0))
plot(relativeI1vec_6, relativeMaxInformClas6+relativeI1vec_6, col = "red",
     type ='l', lty= 3, lwd=2, ylim =c(1,6), xlab = expression(t[xi](I[1])),
     cex.lab=1.7, ylab=expression(t[xi](I[1])+t[xi](I[2*","*max])), xaxt="n")
axis(1, at= seq(0.75,1.75, by=0.25))
axis(1, at=c(0.6704^2,4*0.4894), labels = c(expression(t[xi](I[1*","*min])),
                                            expression(t[xi](I[1*","*max]))),
     cex.axis=1.4)
points(relativeI1vec_6, relativeMaxInformFisher6+relativeI1vec_6, col = "blue",
       lty= 1, lwd=2, type ='l')
points(relativeI1vec_6, relativeMaxInformIN6+relativeI1vec_6, col = "black",
       type ='l', lty=2, lwd=2)
legend(0.75, 3, legend=c("separate studies", "inverse normal method", 
                           "Fisher's product test"),
       col=c("red", "black", "blue"), lty=c(3,2,1), cex=1.25, lwd=c(2,2,2),
       text.width = 0.75)
dev.off()

