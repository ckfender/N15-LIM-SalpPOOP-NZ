#Note that first you have to set the directory: Session>Set Working Directory
rm(list=ls())
library(R.matlab)
library(limSolve)
library(MASS)
setwd('C:/Users/Admin/Dropbox/SalpLIM/N15-LIM-SalpPOOP-main/')
source('xsampleN15outputsNZ.r')


#Choose which Lagrangian Cycle to run model for
Cycle <- 1
if(Cycle==1){
  source('ExternalFunctionsNZC1.r')
} else{
  source('ExternalFunctionsNZ.r')
}

#Read in Matlab data containing the data matrices
tmp <- paste('N15NZInverseCycle',Cycle,'.mat',sep="")
data <- readMat(tmp, fixNames=TRUE)
#Note these are all transposed relative to the excel matrix (since that's how Mike formatted his original CRD excel)
h <- data[['h']]  #The values for which Gx is greater than or equal to (Gx>==h)
b <- data[['b']]  #The values to which all exact and approximate equalities are equal
ba <- data[['ba']] #The values to which the approximate equalities are equal (rate measurements)
be <- data[['be']] #The values to which the exact equalities are equal (mass balances so all 0)
Ae <- data[['Ae']] #The matrix of compartments*flows (Ex) that codifies the mass balances (the first section of the excel spreadsheet)
Aa <- data[['Aa']] #The matrix of rate measurements*flows (Ax) that codifies the approximate equalities
Aa0 <- data[['Aa0']]   #Aa prior to the application of weights
G <- data[['G']]
sdba <- data[['sdba']]
InputCol <- data[['InputCol']]
A <- data[['A']]
Inputs <- data[['Inputs']]
d15NInputs <- data[['d15NInputs']]
wts <- data[['wts']]
#W <- rep(1,39)
#Set run length
BurninLength1 <- 1000
BurninLength2 <- 1000
IterLength <- 1000000
##Set what proportion of model solutions to keep
OutputLength <- IterLength/1000
#jmpLength determines how quickly the solution space is sampled and should be tuned to individual datasets
#If the acceptance rate is above 20%, make the jump length longer and vice versa
if(Cycle==1){
  #jmpLength <- 0.05
  jmpLength <- 0.03
  jmpLength15 <- 0.2
  #jmpLength15 <- 0.15
}else if(Cycle==2){
  #jmpLength <- 0.075
  jmpLength <- 0.035
  jmpLength15 <- 0.02
}else if(Cycle==3){
  #jmpLength=.07
  jmpLength=.02
  jmpLength15 <- 0.1
  #jmpLength15 <- 0.02
}else if(Cycle==4){
  #jmpLength <- 0.14
  jmpLength <- 0.055
  jmpLength15 <- 0.1
}else if(Cycle==5){
  jmpLength <- 0.09
  #jmpLength <- 0.006
  jmpLength15 <- 0.02
}

#Aa2 is approximate equations without 15N (used for initial non-15N model to determine starting location)
if(Cycle==1){
  Aa2 <- Aa[1:36,] #This seems to be just the rate measurements, no 15N for initial L2MN. Note that any columns missing measurements will have been removed in matlab and thus may not match Excel
  ba2 <- ba[1:36]
  sdba2 <- sdba[1:36]
}else{
  Aa2 <- Aa[1:39,] 
  ba2 <- ba[1:39]
  sdba2 <- sdba[1:39]
}


#Initializing unknown 15N values
if(Cycle==1){
  d15N0 <- c(0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0) #Setting unknown d15N values. Mike did 2 extra for his RSDetInputs and RSDetInputd in External Function, I need only 1 extra for unknown RupNO3
}else{
  d15N0 <- c(0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0)
}
#Here Mike started off at the 15N of upwelled nitrate but that is unknown for me so staying at 0
#Setting known 15N values
Meso15N=d15NInputs[1,InputCol]
Macro15N=d15NInputs[2,InputCol]
Salp15N=d15NInputs[3,InputCol]
del15Nknown = c(Meso15N,Macro15N,Salp15N)


#L2MN Solution
test <- lsei(A = Aa2, B = ba2, E = Ae, F = be, G = G, H = h, type=2)
X <- test[['X']]
solNorm <- test[['solutionNorm']]
lseisol <- as.matrix(X)
rm(test, X)
fileout <- paste('N15NZInverseLSEICycle',Cycle,'Routputs.mat',sep="")
writeMat(fileout, Aa = Aa2, ba = ba2, Ae = Ae, be = be, G = G, h = h, lseisol = lseisol, wts=wts)



#Calculating Central value to use as starting location for MCMC algorithm
center <- xranges(E = Ae, F = be, G = G, H = h,
                  ispos = FALSE, tol = 1e-8, central = TRUE, full=FALSE)
centralval <- center[,3]


#Initial Burn-in for MCMC without 15N (large standard deviation to ensure that we move away from starting value)
test2 <- xsample(A = Aa2, B = ba2, E = Ae, F = be, G = G,
                 H = h, sdB = sdba2*10, iter = BurninLength1, type="mirror", jmp=jmpLength+runif(1, 0, 1)*jmpLength/5, x0 = centralval, fulloutput='TRUE')
Burninmat <- test2[['X']]
Startpt <- Burninmat[BurninLength1,]
rm(Burninmat)

##Second Burn-in for MCMC without 15N (for a proper burn-in period)
# test2 <- xsample(A = Aa2, B = ba2, E = Ae, F = be, G = G,
#                  H = h, sdB = sdba2*10, iter = BurninLength2, type="mirror", jmp=jmpLength+runif(1, 0, 1)*jmpLength/5, x0 = Startpt, fulloutput='TRUE')
# Burninmat <- test2[['X']]
# Startpt <- Burninmat[BurninLength2,]
# rm(Burninmat)

#MCMC without 15N
test2 <- xsample(A = Aa2, B = ba2, E = Ae, F = be, G = G,
                 H = h, sdB = sdba2, iter = IterLength/10, outputlength = OutputLength, type="mirror", jmp=jmpLength+runif(1, 0, 1)*jmpLength/5, x0 = Startpt, fulloutput='TRUE')

#MCMC solution (without 15N)
MCMCmatplain <- test2[['X']]
MCMCmatmean <- colMeans(MCMCmatplain)


#Initializing values needed for 15N solution
Eps_Remin <- -1
Eps_Eg <- -2
Eps_NH4up <- -10
sdbafactor <- -Eps_Remin
if(Cycle==1){
  sdba[37:60] <- sdbafactor #setting sd of 15N approximate equalities
} else{
  sdba[40:64] <- sdbafactor
}


#Initial Burn-in for MCMC with 15N (large standard deviation to ensure that we move away from starting value)
test2 <- xsampleN15outputs(A = Aa0, B = ba, E = Ae, F = be, G = G,
                 H = h, sdB = sdba*10, wts = wts, iter = BurninLength2, burninlength = BurninLength1, type="mirror", jmp=jmpLength+runif(1, 0, 1)*jmpLength/5, jmp2=jmpLength15/5, x0 = MCMCmatmean, del15N1=d15N0, del15Nknown=del15Nknown, fulloutput='TRUE')
Burninmat <- test2[['X']]
del15Ntemp <- test2[['del15Ntrack']]
d15N0 <- del15Ntemp[BurninLength2,]
Startpt <- Burninmat[BurninLength2,]
rm(del15Ntemp, Burninmat)

#MCMC with 15N
test2 <- xsampleN15outputs(A = Aa0, B = ba, E = Ae, F = be, G = G,
                 H = h, sdB = sdba, wts=wts, iter = IterLength, outputlength = OutputLength, type="mirror", jmp=jmpLength+runif(1, 0, 1)*jmpLength/5, jmp2=jmpLength15, x0 = Startpt, del15N1=d15N0, del15Nknown=del15Nknown, fulloutput='TRUE')
# Xa1 <- test2[['Xa1']]   #Not sure what all this is
# Xa2 <- test2[['Xa2']]
# Aa1 <- test2[['Aa1']]
# Aa2 <- test2[['Aa2']]
# SSRq1 <- test2[['SSRq1']]
# SSRq2 <- test2[['SSRq2']]

#MCMC with 15N solution
MCMCmat <- test2[['X']]
MCMCmat_rejects <- test2[['x_rejects']]
del15N <- test2[['del15Ntrack']]
del15N_rejects <- test2[['del15Ntrack_rejects']]
randomnumber <- test2[['randomnumber']] #Possibly a RNG seed output that's no longer implemented?
p <- test2[['p']]
p_rejects <- test2[['p_rejects']]
err <- test2[['err']] #This also no longer exists in the code
acceptedratio <- test2[['acceptedratio']]
fileout <- paste('N15NZInverseCycle',Cycle,'Routputs.mat',sep="")
writeMat(fileout, MCMCmat = MCMCmat, MCMCmat_rejects = MCMCmat_rejects, MCMCmatplain = MCMCmatplain, del15N=del15N, del15N_rejects=del15N_rejects, Aa = Aa, ba = ba, Ae = Ae, be = be, G = G, h = h, sdba = sdba, Inputs = Inputs, InputCol = InputCol, jmpLength = jmpLength, jmpLength15 = jmpLength15, acceptedratio = acceptedratio, Startpt = Startpt, centralval = centralval, p=p, p_rejects=p_rejects, err=err, del15Nknown=del15Nknown, lseisol=lseisol, wts=wts, IterLength=IterLength)

# N <- 1:(IterLength/100)
# N <- N*100
# MCMCmat <- MCMCmattemp[N,]
# rm(MCMCmattemp,N)
# 
# 
#writeMat('N15InverseModel.NEMURO.Coastal.ROutputs.mat', MCMCmat = MCMCmat, del15N=del15N, Xa1=Xa1, Xa2=Xa2, Aa1=Aa1, Aa2=Aa2, SSRq1=SSRq1, SSRq2=SSRq2, ba=ba, sdba=sdba, randomnumber=randomnumber)
#plot(MCMCmat[,1], type="o", col="blue")

# Aa3 = ResetRN15(Aa,del15N[OutputLength,],del15Nknown[1],del15Nknown[2],del15Nknown[3],del15Nknown[4],del15Nknown[5])
# errors=(Aa3%*%MCMCmat[OutputLength,]-ba)/sdba
# sum(errors^2)