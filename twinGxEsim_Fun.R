require(OpenMx)
require(MASS)


##############################################################################################################
##############################################################################################################
##                                                                                                          ##
##  Function to simulate genome-wide and twin data                                                          ##
##                                                                                                          ##
##############################################################################################################
##############################################################################################################


gxeSim <- function(N, DZr, MZr, nSNPs, nModSNPs, binMod = F, b0, b2, snpSIGN = T, modSIGN = T, effDist = "fix") {

snpBETA <- as.matrix(rep(sqrt(1/nSNPs), nSNPs))
modBETA <- as.matrix(c(rep(sqrt(1/nSNPs), nModSNPs), rep(0, nSNPs-nModSNPs)))

if(snpSIGN == T & effDist == "fix") {
snpSign <- rbinom(nSNPs, 1, .5)
snpBETA <- ifelse(snpSign==1, snpBETA, -1*snpBETA)
}

if(modSIGN == T & effDist == "fix") {
modSign <- rbinom(nSNPs, 1, .5)
modBETA <- ifelse(modSign==1, modBETA, -1*modBETA)
}

if(effDist == "norm") {
snpBETA <- rnorm(nSNPs, 0, sqrt(1/nSNPs))
modBETA <- rnorm(nSNPs, 0, sqrt(1/nSNPs))
}


###################################
###################################
## Simulate data for DZ twins 

#Maternal Genotype

DZmomA1 <- matrix(NA, N, nSNPs)
DZmomA2 <- matrix(NA, N, nSNPs)

for(i in 1:nSNPs){
DZmomA1[,i] <- rbinom(N, 1, .5)
DZmomA2[,i] <- rbinom(N, 1, .5)
}
DZmomG <- DZmomA1 + DZmomA2

#Paternal Genotype
DZdadA1 <- matrix(NA, N, nSNPs)
DZdadA2 <- matrix(NA, N, nSNPs)

for(i in 1:nSNPs){
DZdadA1[,i] <- rbinom(N, 1, .5)
DZdadA2[,i] <- rbinom(N, 1, .5)
}
DZdadG <- DZdadA1 + DZdadA2

# Parental Unique Environments 
DZmomE <- rnorm(N)
DZdadE <- rnorm(N)

#Parental Phenotypes
DZmom <- DZmomG%*%snpBETA + DZmomE * sqrt(.5)
DZdad <- DZdadG%*%snpBETA + DZdadE * sqrt(.5)

# Selecting a random allele from mom and dad for each twin
mom_dz1 <- matrix(NA, N, nSNPs)
dad_dz1 <- matrix(NA, N, nSNPs)

mom_dz2 <- matrix(NA, N, nSNPs)
dad_dz2 <- matrix(NA, N, nSNPs)

for(i in 1:nSNPs){
mom_dz1[,i] <- rbinom(N, 1, .5)
dad_dz1[,i] <- rbinom(N, 1, .5)

mom_dz2[,i] <- rbinom(N, 1, .5)
dad_dz2[,i] <- rbinom(N, 1, .5)
}

dz1_A1 <- matrix(NA, N, nSNPs)
dz1_A2 <- matrix(NA, N, nSNPs)

dz2_A1 <- matrix(NA, N, nSNPs)
dz2_A2 <- matrix(NA, N, nSNPs)

for(i in 1:nSNPs){
dz1_A1[,i] <- ifelse(mom_dz1[,i] == 0, DZmomA1[,i], DZmomA2[,i])
dz1_A2[,i] <- ifelse(dad_dz1[,i] == 0, DZdadA1[,i], DZdadA2[,i])

dz2_A1[,i] <- ifelse(mom_dz2[,i] == 0, DZmomA1[,i], DZmomA2[,i])
dz2_A2[,i] <- ifelse(dad_dz2[,i] == 0, DZdadA1[,i], DZdadA2[,i])
}

# DZ twin genotypes
dz1G <- dz1_A1 + dz1_A2
dz2G <- dz2_A1 + dz2_A2

# DZ twin moderators
DZmods <- as.data.frame(mvrnorm(N, c(0,0), matrix(c(1, DZr, DZr, 1),2,2), empirical = T))
colnames(DZmods) <- c("dz1M", "dz2M")

if (binMod == T) DZmods <- as.data.frame(ifelse(DZmods>0, 1, 0))

dz1M <- DZmods$dz1M
dz2M <- DZmods$dz2M

# DZ twin unique environments 
dz1E <- rnorm(N)
dz2E <- rnorm(N)

# define the DZ phenotypes
dz1 <- b0 + dz1G%*%snpBETA + b2*dz1M + (dz1G* dz1M)%*%modBETA + dz1E* sqrt(.5) 
dz2 <- b0 + dz2G%*%snpBETA + b2*dz2M + (dz2G* dz2M)%*%modBETA + dz2E* sqrt(.5) 

###################################
###################################
## Simulate data for MZ twins 

#Maternal Genotype
MZmomA1 <- matrix(NA, N, nSNPs)
MZmomA2 <- matrix(NA, N, nSNPs)

for(i in 1:nSNPs){
MZmomA1[,i] <- rbinom(N, 1, .5)
MZmomA2[,i] <- rbinom(N, 1, .5)
}
MZmomG <- MZmomA1 + MZmomA2

#Paternal Genotype
MZdadA1 <- matrix(NA, N, nSNPs)
MZdadA2 <- matrix(NA, N, nSNPs)

for(i in 1:nSNPs){
MZdadA1[,i] <- rbinom(N, 1, .5)
MZdadA2[,i] <- rbinom(N, 1, .5)
}
MZdadG <- MZdadA1 + MZdadA2

# Parental Unique Environments 
MZmomE <- rnorm(N)
MZdadE <- rnorm(N)

#Parental Phenotypes
MZmom <- MZmomG%*%snpBETA + MZmomE * sqrt(.5)
MZdad <- MZdadG%*%snpBETA + MZdadE * sqrt(.5)

# Selecting a random allele from mom and dad for each twin
mom_mz <- matrix(NA, N, nSNPs)
dad_mz <- matrix(NA, N, nSNPs)

for(i in 1:nSNPs){
mom_mz[,i] <- rbinom(N, 1, .5)
dad_mz[,i] <- rbinom(N, 1, .5)
}

mz_A1 <- matrix(NA, N, nSNPs)
mz_A2 <- matrix(NA, N, nSNPs)

for(i in 1:nSNPs){
mz_A1[,i] <- ifelse(mom_mz[,i] == 0, MZmomA1[,i], MZmomA2[,i])
mz_A2[,i] <- ifelse(dad_mz[,i] == 0, MZdadA1[,i], MZdadA2[,i])

}

# MZ twin genotypes
mzG <- mz_A1 + mz_A2

# MZ twin moderators
MZmods <- as.data.frame(mvrnorm(N, c(0,0), matrix(c(1, MZr, MZr, 1),2,2), empirical = T))
colnames(MZmods) <- c("mz1M", "mz2M")

if (binMod == T) MZmods <- as.data.frame(ifelse(MZmods>0, 1, 0))

mz1M <- MZmods$mz1M
mz2M <- MZmods$mz2M

# MZ twin unique environments 
mz1E <- rnorm(N)
mz2E <- rnorm(N)

# define the MZ phenotypes
mz1 <- b0 + mzG%*%snpBETA + b2*mz1M + (mzG* mz1M)%*%modBETA + mz1E* sqrt(.5) 
mz2 <- b0 + mzG%*%snpBETA + b2*mz2M + (mzG* mz2M)%*%modBETA + mz2E* sqrt(.5) 

mzData <- as.data.frame(cbind(mz1, mz2, mz1M, mz2M))
dzData <- as.data.frame(cbind(dz1, dz2, dz1M, dz2M))

selVars <- paste("t", 1:2, sep = "")
modVars <- paste("mod", 1:2, sep = "")

colnames(mzData) <- colnames(dzData) <- c(selVars, modVars)

return(list(mzData = mzData, dzData = dzData, mz1G = mzG, mz2G = mzG, dz1G = dz1G,dz2G = dz2G, 
	        DZmomG = DZmomG, DZmom = DZmom, DZdadG = DZdadG, DZdad = DZdad, 
			MZmomG = MZmomG, MZmom = MZmom, MZdadG = MZdadG, MZdad = MZdad, 
			snpBETA = snpBETA, modBETA = modBETA))

}

##############################################################################################################
##############################################################################################################
##                                                                                                          ##
##  Function to simulate only genome-wide for PRS analyses when the twin data is too large                  ##
##                                                                                                          ##
##############################################################################################################
##############################################################################################################


gxeSimSmall <- function(N, nSNPs, nModSNPs, binMod = F, b0, b2) {

snpBETA <- as.matrix(rep(sqrt(1/nSNPs), nSNPs))
modBETA <- as.matrix(c(rep(sqrt(1/nSNPs), nModSNPs), rep(0, nSNPs-nModSNPs)))


###################################
###################################
## Simulate data for DZ twins 

#Maternal Genotype

DZmomA1 <- matrix(NA, N, nSNPs)
DZmomA2 <- matrix(NA, N, nSNPs)

for(i in 1:nSNPs){
DZmomA1[,i] <- rbinom(N, 1, .5)
DZmomA2[,i] <- rbinom(N, 1, .5)
}

#Paternal Genotype
DZdadA1 <- matrix(NA, N, nSNPs)
DZdadA2 <- matrix(NA, N, nSNPs)

for(i in 1:nSNPs){
DZdadA1[,i] <- rbinom(N, 1, .5)
DZdadA2[,i] <- rbinom(N, 1, .5)
}

# Selecting a random allele from mom and dad for each twin
mom_dz1 <- matrix(NA, N, nSNPs)
dad_dz1 <- matrix(NA, N, nSNPs)

for(i in 1:nSNPs){
mom_dz1[,i] <- rbinom(N, 1, .5)
dad_dz1[,i] <- rbinom(N, 1, .5)

}

dz1_A1 <- matrix(NA, N, nSNPs)
dz1_A2 <- matrix(NA, N, nSNPs)

for(i in 1:nSNPs){
dz1_A1[,i] <- ifelse(mom_dz1[,i] == 0, DZmomA1[,i], DZmomA2[,i])
dz1_A2[,i] <- ifelse(dad_dz1[,i] == 0, DZdadA1[,i], DZdadA2[,i])
}

# DZ twin genotypes
dz1G <- dz1_A1 + dz1_A2

# DZ twin moderators
DZmods <- as.data.frame(mvrnorm(N, c(0), matrix(c(1),1,1), empirical = T))
colnames(DZmods) <- c("dz1M")

if (binMod == T) DZmods <- as.data.frame(ifelse(DZmods>0, 1, 0))

dz1M <- DZmods$dz1M

# DZ twin unique environments 
dz1E <- rnorm(N)

# define the DZ phenotypes
dz1 <- b0 + dz1G%*%snpBETA + b2*dz1M + (dz1G* dz1M)%*%modBETA + dz1E* sqrt(.5) 

dzData <- as.data.frame(cbind(dz1, dz1M))
colnames(dzData) <- c("t1", "mod1")


return(list(dzData = dzData, dz1G = dz1G, snpBETA = snpBETA, modBETA = modBETA))

}

gxeACE <- function(mzData, dzData) {
	
	selVars <- paste("t", 1:2, sep = "")
	modVars <- paste("mod", 1:2, sep = "")

	colnames(mzData) <- colnames(dzData) <- c(selVars, modVars)

	MZdata  <-  mxData(mzData, type="raw" )
	DZdata  <-  mxData(dzData, type="raw" )
	nv <- 1

	#=======================================================================#
	#   PREPARE MODEL                                                       #
	#=======================================================================#

	# Matrices to store a, c, and e Path Coefficients
	pathA    <- mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=.6, label="a11", name="a" ) 
	pathC    <- mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=.6, label="c11", name="c" )
	pathE    <- mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=.6, label="e11", name="e" )

	# Matrices to store the moderating a, c, and e Path Coefficients
	modPathA <- mxMatrix( type='Full', nrow=nv, ncol=nv, free=TRUE, values=.1, label="aM11", name="aM" )
	modPathC <- mxMatrix( type='Full', nrow=nv, ncol=nv, free=TRUE, values=.1, label="cM11", name="cM" )
	modPathE <- mxMatrix( type='Full', nrow=nv, ncol=nv, free=TRUE, values=.1, label="eM11", name="eM" )

	# Matrix for the moderator variable
	mod1 <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.mod1"), name="D1")
	mod2 <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.mod2"), name="D2")

	# Matrices to compute the moderated A, C, and E variance components
	covAmod <- mxAlgebra( expression=(a+ D1%*%aM) %*% t(a+ D2%*%aM), name="A" )
	covCmod <- mxAlgebra( expression=(c+ D1%*%cM) %*% t(c+ D2%*%cM), name="C" )
	covEmod <- mxAlgebra( expression=(e+ D1%*%eM) %*% t(e+ D2%*%eM), name="E" )

	# Algebra for the expected mean vector and expected variance/covariance matrices 
	covMZ <- mxAlgebra( expression= rbind ( cbind(A+C+E , A+C),
	                                        cbind(A+C, A+C+E)), name="expCovMZ" )
	covDZ <- mxAlgebra( expression= rbind ( cbind(A+C+E , 0.5%x%A+C),
	                                        cbind(0.5%x%A+C, A+C+E)), name="expCovDZ" )


	# Matrices for the expected means vector
	intercept <- mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=.1, label="interc", name="int" )

	self     <- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=0, label=c("reg1"), name="self" )
	cotwin   <- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=0, label=c("reg2"), name="cotwin" )

	meanG <- mxAlgebra( expression= cbind((int + self*D1 + cotwin*D2),(int + self*D2 + cotwin*D1)), name="expMean" )


	obs <-  list(pathA, pathC, pathE, modPathA, modPathC, modPathE, mod1, mod2, 
		         covAmod, covCmod, covEmod, covMZ, covDZ, intercept, self, cotwin, meanG)
	
	fun <- mxFitFunctionML()
	mzExp <- mxExpectationNormal(covariance="expCovMZ", means="expMean", dimnames=selVars )
	dzExp <- mxExpectationNormal(covariance="expCovDZ", means="expMean", dimnames=selVars )

	MZ <- mxModel("MZ", obs, MZdata, fun, mzExp)
	DZ <- mxModel("DZ", obs, DZdata, fun, dzExp)

	modFun <- mxFitFunctionMultigroup(c("MZ","DZ"))
	aceMod <- mxModel("ACE", MZ, DZ, fun, modFun)

	aceModFit <- mxRun(aceMod, silent = T)
	if(aceModFit$output$status$code >1) aceModFit <- mxRun(aceModFit, silent = T)
	

	aceModFit
}

##############################################################################################################
##############################################################################################################
##                                                                                                          ##
##  Function to plot the twin results                                                                       ##
##                                                                                                          ##
##############################################################################################################
##############################################################################################################

gxePlot <- function(aceModFit = aceModFit, binMod = F, title = "Add Title", legend = T, legend.loc = "topleft"){

a <- aceModFit$output$matrices$MZ.a
c <- aceModFit$output$matrices$MZ.c
e <- aceModFit$output$matrices$MZ.e
aM <- aceModFit$output$matrices$MZ.aM
cM <- aceModFit$output$matrices$MZ.cM
eM <- aceModFit$output$matrices$MZ.eM

if (binMod == T) mod <- c(0,1)
	if (binMod == F) mod <- seq(-2.5,2.5, .25)

Va <- (c(a) + c(aM)*mod)^2
Vc <- (c(c) + c(cM)*mod)^2
Ve <- (c(e) + c(eM)*mod)^2
Vt <- Va + Vc + Ve
SVa <- Va/Vt
SVc <- Vc/Vt
SVe <- Ve/Vt
out <- as.matrix(cbind(Va,Vc,Ve,Vt,SVa,SVc,SVe))
rownames(out) <- paste0("mod=", mod)

matplot(mod, out[,1:4], type="l", lty=4:1, col=4:1, xlab="Moderator", 
        ylab="Variance Component", main= title, 
		lwd = 3, xaxt = "n")
if (legend == T) legend(legend.loc, c("Va","Vc","Ve","Vt"), lty=4:1, col=4:1, , bty = "n")

out
}

##############################################################################################################
##############################################################################################################
##                                                                                                          ##
##  Function to conduct the GxE GWAS (from the GxE sim data function)                                       ##
##                                                                                                          ##
##############################################################################################################
##############################################################################################################


GxEgwas <- function(Sims, person = "dz1"){

if (person == "dz1") phenoData <- as.data.frame(cbind(Sims$dzData$t1, Sims$dzData$mod1))
if (person == "dz2") phenoData <- as.data.frame(cbind(Sims$dzData$t2, Sims$dzData$mod2))
if (person == "mz1") phenoData <- as.data.frame(cbind(Sims$mzData$t1, Sims$mzData$mod1))
if (person == "mz2") phenoData <- as.data.frame(cbind(Sims$mzData$t2, Sims$mzData$mod2))

colnames(phenoData) <- c("DV", "mod")

if (person == "dz1") genoData <- Sims$dz1G
if (person == "dz2") genoData <- Sims$dz2G
if (person == "mz1") genoData <- Sims$mz1G
if (person == "mz2") genoData <- Sims$mz2G

nSnp <- dim(genoData)[2]
colnames(genoData) <- paste0("rs", 1:nSnp)

Res <- as.data.frame(matrix(NA, nSnp, 13))

for(i in 1:nSnp)	{

variant <- genoData[,i]
phenoData$snp <- variant
phenoData$inter <- phenoData$snp * phenoData$mod

fit <- lm(DV ~ snp + mod + inter, data = phenoData)
summ <- summary(fit)

genInfo <- 	names(genoData[1,i])
COEF <- summ$coef[,1]
vCov <- diag(vcov(fit))
offDiag <- vcov(fit)[c("inter"),c("snp", "mod")]
N <- fit$df.residual +  dim(summary(fit)$coef)[1]
mu <- mean(phenoData$snp, na.rm = T)/2
maf <- ifelse(mu< .5 , mu , 1- mu)


Res[i,] <- as.data.frame(t(c(genInfo, COEF, vCov, offDiag, N, maf)	))


	}
	colnames(Res) <- c("rsID", "b0", "b1_snp", "b2_mod", "b3_inter", "v0", "v_snp", "v_mod", "v_inter", "v_inter_snp", "v_inter_mod", "N", "maf")
	
	for(j in 2:13) Res[,j] <- as.numeric(Res[,j])
	Res

}

##############################################################################################################
##############################################################################################################
##                                                                                                          ##
##  Function to conduct a standard GWAS (from the GxE sim data function)                                    ##
##                                                                                                          ##
##############################################################################################################
##############################################################################################################

TypGWAS <- function(Sims, person = "dz1", incl_mod = F){

if (person == "dz1") phenoData <- as.data.frame(cbind(Sims$dzData$t1, Sims$dzData$mod1))
if (person == "dz2") phenoData <- as.data.frame(cbind(Sims$dzData$t2, Sims$dzData$mod2))
if (person == "mz1") phenoData <- as.data.frame(cbind(Sims$mzData$t1, Sims$mzData$mod1))
if (person == "mz2") phenoData <- as.data.frame(cbind(Sims$mzData$t2, Sims$mzData$mod2))

colnames(phenoData) <- c("DV", "mod")

if (person == "dz1") genoData <- Sims$dz1G
if (person == "dz2") genoData <- Sims$dz2G
if (person == "mz1") genoData <- Sims$mz1G
if (person == "mz2") genoData <- Sims$mz2G

nSnp <- dim(genoData)[2]
colnames(genoData) <- paste0("rs", 1:nSnp)

nPar <- ifelse(incl_mod == F, 7, 9)

Res <- as.data.frame(matrix(NA, nSnp, nPar))

for(i in 1:nSnp)	{

variant <- genoData[,i]
phenoData$snp <- variant

if(incl_mod == F) fit <- lm(DV ~ snp, data = phenoData)
if(incl_mod == T) fit <- lm(DV ~ snp + mod, data = phenoData)
	
summ <- summary(fit)

genInfo <- 	names(genoData[1,i])
COEF <- summ$coef[,1]
vCov <- diag(vcov(fit))
N <- fit$df.residual +  dim(summary(fit)$coef)[1]
mu <- mean(phenoData$snp, na.rm = T)/2
maf <- ifelse(mu< .5 , mu , 1- mu)


Res[i,] <- as.data.frame(t(c(genInfo, COEF, vCov, N, maf)	))


	}
	
if (incl_mod == F) colnames(Res) <- c("rsID", "b0", "b1_snp",           "v0", "v_snp",          "N", "maf")
if (incl_mod == T) colnames(Res) <- c("rsID", "b0", "b1_snp", "b2_mod", "v0", "v_snp", "v_mod", "N", "maf")
	
	if (incl_mod == F)	for(j in 2:7) Res[,j] <- as.numeric(Res[,j])
		if (incl_mod == T)	for(j in 2:9) Res[,j] <- as.numeric(Res[,j])
	Res

}

##############################################################################################################
##############################################################################################################
##                                                                                                          ##
##  Function to calculate marginal genetic effects (not used in the paper)                                  ##
##                                                                                                          ##
##############################################################################################################
##############################################################################################################

 
 marg <- function(result, level, mainEffect, inter, main_V, int_V, cov_main_int, lab){
 	
result$ME    <- result[[mainEffect]] + level * result[[inter]]
result$ME.SE <- sqrt(result[[main_V]] +
                level^2 * result[[int_V]] +
                2*level * result[[cov_main_int]])

colnames(result)[colnames(result) == "ME"] <- lab
colnames(result)[colnames(result) == "ME.SE"] <- paste0(lab, "_se")

	result
 }
 
 ##############################################################################################################
 ##############################################################################################################
 ##                                                                                                          ##
 ##  Function to calculate statistical significance of the GWAS parameters                                   ##
 ##                                                                                                          ##
 ##############################################################################################################
 ##############################################################################################################
 
 
 sig <- function(result, focus, focus_se, SE = F, lab){

	 if(SE == T) result$Z <-  result[[focus]] / result[[focus_se]]
		 	 if(SE == F) result$Z <-  result[[focus]] / sqrt(result[[focus_se]])

 result$P <- 2*pnorm(-abs(result$Z))

 colnames(result)[colnames(result) == "Z"] <- paste0(lab, "_Z")
 colnames(result)[colnames(result) == "P"] <- paste0(lab, "_P")

 result
 }

 ##############################################################################################################
 ##############################################################################################################
 ##                                                                                                          ##
 ##  Function to PRS                                                                                         ##
 ##                                                                                                          ##
 ##############################################################################################################
 ##############################################################################################################
 
prs <- function(disc, param, p_val, thresh, targ, person, lab){

	if (person == "dz1") phenoData <- as.data.frame(cbind(targ$dzData$t1, targ$dzData$mod1))
	if (person == "dz2") phenoData <- as.data.frame(cbind(targ$dzData$t2, targ$dzData$mod2))
	if (person == "mz1") phenoData <- as.data.frame(cbind(targ$mzData$t1, targ$mzData$mod1))
	if (person == "mz2") phenoData <- as.data.frame(cbind(targ$mzData$t2, targ$mzData$mod2))

	colnames(phenoData) <- c("DV", "mod")

	if (person == "dz1") genoData <- targ$dz1G
	if (person == "dz2") genoData <- targ$dz2G
	if (person == "mz1") genoData <- targ$mz1G
	if (person == "mz2") genoData <- targ$mz2G


	nSnp <- dim(genoData)[2]
	colnames(genoData) <- paste0("rs", 1:nSnp)
	
	disc[[param]][ disc[[p_val]] >  thresh ] <- 0
	
	BETA <- disc[[param]]
	
	
	phenoData$PRS <- genoData %*% BETA
	
    colnames(phenoData)[colnames(phenoData) == "PRS"] <- lab
	
	phenoData
}


