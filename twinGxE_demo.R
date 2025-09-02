require(OpenMx)
require(MASS)

source("twinGxEsim_Fun.R")



############################################################################################################
############################################################################################################
## Study 1
## SIMULATE DATA FOR A TWIN MODEL WITH a single variant and *NO* GxE
## 
############################################################################################################
############################################################################################################

set.seed(11111)
oneSNP_noGxE   <- gxeSim(N = 10000, DZr = 0, MZr = 0, nSNPs = 1, nModSNPs = 0, binMod = F, b0 = 0, b2 = 0, snpSIGN = F)
oneSNP_noGxEFit <- gxeACE(mzData = oneSNP_noGxE$mzData, dzData = oneSNP_noGxE$dzData)
summary(oneSNP_noGxEFit)
oneSNP_noGxEres <- gxePlot(aceModFit = oneSNP_noGxEFit, binMod = F, title = "Normal - 1 SNP and No GxE", legend = F, legend.loc = "left")


dzT1  <- summary(lm(oneSNP_noGxE$dzData$t1 ~ oneSNP_noGxE$dz1G[,1]))
dzT2  <- summary(lm(oneSNP_noGxE$dzData$t2 ~ oneSNP_noGxE$dz2G[,1]))
mzT1  <- summary(lm(oneSNP_noGxE$mzData$t1 ~ oneSNP_noGxE$mz1G[,1]))
mzT2  <- summary(lm(oneSNP_noGxE$mzData$t2 ~ oneSNP_noGxE$mz2G[,1]))
dzMom <- summary(lm(oneSNP_noGxE$DZmom ~ oneSNP_noGxE$DZmomG[,1]))
mzMom <- summary(lm(oneSNP_noGxE$MZmom ~ oneSNP_noGxE$MZmomG[,1]))
dzDad <- summary(lm(oneSNP_noGxE$DZdad ~ oneSNP_noGxE$DZdadG[,1]))
mzdad <- summary(lm(oneSNP_noGxE$MZdad ~ oneSNP_noGxE$MZdadG[,1]))



regRes <- cbind(
c(dzT1 $coef[,1], dzT1 $r.squared),
c(dzT2 $coef[,1], dzT2 $r.squared),
c(mzT1 $coef[,1], mzT1 $r.squared),
c(mzT2 $coef[,1], mzT2 $r.squared),
c(dzMom$coef[,1], dzMom$r.squared),
c(mzMom$coef[,1], mzMom$r.squared),
c(dzDad$coef[,1], dzDad$r.squared),
c(mzdad$coef[,1], mzdad$r.squared)
)
rownames(regRes) <- c("Intercept", "BetaSNP", "r2")
colnames(regRes) <- c("dzT1", "dzT2", "mzT1", "mzT2", "dzMom", "mzMom", "dzDad", "mzdad")

regRes

round(regRes, 2)


# PRS

# Run "discovery" gwas (on the DZ twins (MZ twins will be the target sample))
DZ1 <- GxEgwas(Sims = oneSNP_noGxE, person = "dz1")
DZ2 <- GxEgwas(Sims = oneSNP_noGxE, person = "dz2")

typDZ1 <- TypGWAS(Sims = oneSNP_noGxE, person = "dz1", incl_mod = T)
typDZ2 <- TypGWAS(Sims = oneSNP_noGxE, person = "dz2", incl_mod = T)

#################
# calculate statistical significance

DZ1 <-  sig(result = DZ1, focus = "b1_snp", focus_se = "v_snp", SE = F, lab = "b1")
DZ1 <-  sig(result = DZ1, focus = "b3_inter", focus_se = "v_inter", SE = F, lab = "b3")

DZ2 <-  sig(result = DZ2, focus = "b1_snp", focus_se = "v_snp", SE = F, lab = "b1")
DZ2 <-  sig(result = DZ2, focus = "b3_inter", focus_se = "v_inter", SE = F, lab = "b3")

typDZ1 <-  sig(result = DZ1, focus = "b1_snp", focus_se = "v_snp", SE = F, lab = "b1")
typDZ2 <-  sig(result = DZ2, focus = "b1_snp", focus_se = "v_snp", SE = F, lab = "b1")


# calculate PRS and run linear regression
mz1.prs.main  <- prs(disc = DZ1, param = "b1_snp", p_val = "b1_P", thresh = 1, targ = oneSNP_noGxE, person = "mz1", lab = "dz1_2_mz1_PRS")  
prs.main <- summary(lm(DV ~  dz1_2_mz1_PRS, data = mz1.prs.main))

mz1.prs.int  <- prs(disc = DZ1, param = "b3_inter", p_val = "b3_P", thresh = 1, targ = oneSNP_noGxE, person = "mz1", lab = "dz1_2_mz1_PRS")  
mz1.prs.int$dz1_2_mz1_PRS <- mz1.prs.int$dz1_2_mz1_PRS* mz1.prs.int$mod
prs.int <- summary(lm(DV ~  dz1_2_mz1_PRS, data = mz1.prs.int))

mz1.prs <- cbind(mz1.prs.main, mz1.prs.int[,3])
colnames(mz1.prs) <- c("DV", "mod", "main", "int")
prs.both <- summary(lm(DV ~  main + mod + int, data = mz1.prs))


mz1.prs.typ  <- prs(disc = typDZ1, param = "b1_snp", p_val = "b1_P", thresh = 1, targ = oneSNP_noGxE, person = "mz1", lab = "dz1_2_mz1_PRS")  
prs.typ <- summary(lm(DV ~  dz1_2_mz1_PRS, data = mz1.prs.typ))

prs.typMod <- summary(lm(DV ~  dz1_2_mz1_PRS + mod + dz1_2_mz1_PRS * mod, data = mz1.prs.typ))

# These results are directly in Table 2
prs.main
prs.int
prs.both

prs.typ
prs.typMod

############################################################################################################
############################################################################################################
## 
## Study 2: SIMULATE DATA FOR A TWIN MODEL WITH a single variant and GxE with a binary moderator
## 
############################################################################################################
############################################################################################################


set.seed(11111)
oneSNP_GxE   <- gxeSim(N = 10000, DZr = 0, MZr = 0, nSNPs = 1, nModSNPs = 1, binMod = T, b0 = 0, b2 = 0, snpSIGN = F, modSIGN = F)
oneSNP_GxEFit <- gxeACE(mzData = oneSNP_GxE$mzData, dzData = oneSNP_GxE$dzData)
summary(oneSNP_GxEFit)
# This is Figure 1
oneSNP_GxEres <- gxePlot(aceModFit = oneSNP_GxEFit, binMod = T, title = " ", legend = T, legend.loc = "topleft")
axis(1, 0:1, labels = c("No", "Yes"))


dzT1  <- summary(lm(oneSNP_GxE$dzData$t1 ~ oneSNP_GxE$dz1G[,1]   + oneSNP_GxE$dzData$mod1 + oneSNP_GxE$dz1G[,1] : oneSNP_GxE$dzData$mod1))
dzT2  <- summary(lm(oneSNP_GxE$dzData$t2 ~ oneSNP_GxE$dz2G[,1]   + oneSNP_GxE$dzData$mod2 + oneSNP_GxE$dz2G[,1] : oneSNP_GxE$dzData$mod2))
mzT1  <- summary(lm(oneSNP_GxE$mzData$t1 ~ oneSNP_GxE$mz1G[,1]   + oneSNP_GxE$mzData$mod1 + oneSNP_GxE$mz1G[,1] : oneSNP_GxE$mzData$mod1))
mzT2  <- summary(lm(oneSNP_GxE$mzData$t2 ~ oneSNP_GxE$mz2G[,1]   + oneSNP_GxE$mzData$mod2 + oneSNP_GxE$mz2G[,1] : oneSNP_GxE$mzData$mod2))

dz1_0 <- summary(lm(oneSNP_GxE$dzData$t1[oneSNP_GxE$dzData$mod1 == 0] ~ oneSNP_GxE$dz1G[,1][oneSNP_GxE$dzData$mod1 == 0] ))$r.squared
dz1_1 <- summary(lm(oneSNP_GxE$dzData$t1[oneSNP_GxE$dzData$mod1 == 1] ~ oneSNP_GxE$dz1G[,1][oneSNP_GxE$dzData$mod1 == 1] ))$r.squared

dz2_0 <- summary(lm(oneSNP_GxE$dzData$t2[oneSNP_GxE$dzData$mod2 == 0] ~ oneSNP_GxE$dz2G[,1][oneSNP_GxE$dzData$mod2 == 0] ))$r.squared
dz2_1 <- summary(lm(oneSNP_GxE$dzData$t2[oneSNP_GxE$dzData$mod2 == 1] ~ oneSNP_GxE$dz2G[,1][oneSNP_GxE$dzData$mod2 == 1] ))$r.squared

mz1_0 <- summary(lm(oneSNP_GxE$mzData$t1[oneSNP_GxE$mzData$mod1 == 0] ~ oneSNP_GxE$mz1G[,1][oneSNP_GxE$mzData$mod1 == 0] ))$r.squared
mz1_1 <- summary(lm(oneSNP_GxE$mzData$t1[oneSNP_GxE$mzData$mod1 == 1] ~ oneSNP_GxE$mz1G[,1][oneSNP_GxE$mzData$mod1 == 1] ))$r.squared

mz2_0 <- summary(lm(oneSNP_GxE$mzData$t2[oneSNP_GxE$mzData$mod2 == 0] ~ oneSNP_GxE$mz2G[,1][oneSNP_GxE$mzData$mod2 == 0] ))$r.squared
mz2_1 <- summary(lm(oneSNP_GxE$mzData$t2[oneSNP_GxE$mzData$mod2 == 1] ~ oneSNP_GxE$mz2G[,1][oneSNP_GxE$mzData$mod2 == 1] ))$r.squared


# These results are presented in Table 3
regRes <- cbind(
c(dzT1 $coef[,1], dzT1 $r.squared, dz1_0, dz1_1),
c(dzT2 $coef[,1], dzT2 $r.squared, dz2_0, dz2_1),
c(mzT1 $coef[,1], mzT1 $r.squared, mz1_0, mz1_1),
c(mzT2 $coef[,1], mzT2 $r.squared, mz2_0, mz2_1)
)
rownames(regRes) <- c("Intercept", "BetaSNP", "BetaMod", "BetaGxE", "r2", "R2_no", "r2_yes")
colnames(regRes) <- c("dzT1", "dzT2", "mzT1", "mzT2")

regRes 

# PRS

# Run "discovery" gwas on DZ twin 1 (MZ twin 1 will be the target sample))
DZ1 <- GxEgwas(Sims = oneSNP_GxE, person = "dz1")

typDZ1 <- TypGWAS(Sims = oneSNP_GxE, person = "dz1", incl_mod = T)

#################
# calculate statistical significance

DZ1 <-  sig(result = DZ1, focus = "b1_snp", focus_se = "v_snp", SE = F, lab = "b1")
DZ1 <-  sig(result = DZ1, focus = "b3_inter", focus_se = "v_inter", SE = F, lab = "b3")

typDZ1 <-  sig(result = DZ1, focus = "b1_snp", focus_se = "v_snp", SE = F, lab = "b1")


# calculate PRS and run linear regression
mz1.prs.main  <- prs(disc = DZ1, param = "b1_snp", p_val = "b1_P", thresh = 1, targ = oneSNP_GxE, person = "mz1", lab = "dz1_2_mz1_PRS")  
prs.main <- summary(lm(DV ~  dz1_2_mz1_PRS, data = mz1.prs.main))

mz1.prs.int  <- prs(disc = DZ1, param = "b3_inter", p_val = "b3_P", thresh = 1, targ = oneSNP_GxE, person = "mz1", lab = "dz1_2_mz1_PRS")  
mz1.prs.int$dz1_2_mz1_PRS <- mz1.prs.int$dz1_2_mz1_PRS* mz1.prs.int$mod
prs.mod <- summary(lm(DV ~  dz1_2_mz1_PRS, data = mz1.prs.int))

mz1.prs <- cbind(mz1.prs.main, mz1.prs.int[,3])
colnames(mz1.prs) <- c("DV", "mod", "main", "int")
prs.both <- summary(lm(DV ~  main + mod + int, data = mz1.prs))


mz1.prs.typ  <- prs(disc = typDZ1, param = "b1_snp", p_val = "b1_P", thresh = 1, targ = oneSNP_GxE, person = "mz1", lab = "dz1_2_mz1_PRS")  
prs.typ <- summary(lm(DV ~  dz1_2_mz1_PRS, data = mz1.prs.typ))

prs.typMod <- summary(lm(DV ~  dz1_2_mz1_PRS + mod + dz1_2_mz1_PRS*mod, data = mz1.prs.typ))

# These results are directly in Table 4
prs.main
prs.mod
prs.both

prs.typ
prs.typMod


############################################################################################################
############################################################################################################
## 
## Study 3: SIMULATE DATA FOR A TWIN MODEL WITH POLYGENIC VARIANTS and GxE
## 
############################################################################################################
############################################################################################################

set.seed(11111)
poly_GxE   <- gxeSim(N = 10000, DZr = .0, MZr = .0, nSNPs = 1000, nModSNPs = 1000, binMod = T, b0 = 0, b2 = 0, snpSIGN = F, modSIGN = F)
poly_GxEFit <- gxeACE(mzData = poly_GxE$mzData, dzData = poly_GxE$dzData)
#poly_GxEFit <- mxRun(poly_GxEFit)
summary(poly_GxEFit)
# This is Figure 2
poly_GxEres <- gxePlot(aceModFit = poly_GxEFit, binMod = T, title = " ", legend = T, legend.loc = "topleft")
axis(1, 0:1, labels = c("No", "Yes"))

# hist(poly_GxE$mzData$t1)

# Run "discovery" gwas on DZ twin 1 (MZ twin 1 will be the target sample))
DZ1 <- GxEgwas(Sims = poly_GxE, person = "dz1")
typDZ1 <- TypGWAS(Sims = poly_GxE, person = "dz1", incl_mod = T)

# calculate statistical significance

DZ1 <-  sig(result = DZ1, focus = "b1_snp", focus_se = "v_snp", SE = F, lab = "b1")
DZ1 <-  sig(result = DZ1, focus = "b3_inter", focus_se = "v_inter", SE = F, lab = "b3")
typDZ1 <-  sig(result = DZ1, focus = "b1_snp", focus_se = "v_snp", SE = F, lab = "b1")

table(DZ1 $b1_P < 5e-8)
table(DZ1 $b3_P < 5e-8)

table(DZ1 $b1_P < .05/1000)
table(DZ1 $b3_P < .05/1000)

mean(DZ1 $b1_snp)
mean(DZ1 $b3_inter)
#################

# calculate PRS and run linear regression
mz1.prs.main  <- prs(disc = DZ1, param = "b1_snp", p_val = "b1_P", thresh = 1, targ = poly_GxE, person = "mz1", lab = "dz1_2_mz1_PRS")  
prs.main <- summary(lm(DV ~  dz1_2_mz1_PRS, data = mz1.prs.main))

mz1.prs.int  <- prs(disc = DZ1, param = "b3_inter", p_val = "b3_P", thresh = 1, targ = poly_GxE, person = "mz1", lab = "dz1_2_mz1_PRS")  
mz1.prs.int$dz1_2_mz1_PRS <- mz1.prs.int$dz1_2_mz1_PRS* mz1.prs.int$mod
prs.mod <- summary(lm(DV ~  dz1_2_mz1_PRS, data = mz1.prs.int))

mz1.prs <- cbind(mz1.prs.main, mz1.prs.int[,3])
colnames(mz1.prs) <- c("DV", "mod", "main", "int")
prs.both <- summary(lm(DV ~  main + mod + int, data = mz1.prs))


mz1.prs.typ  <- prs(disc = typDZ1, param = "b1_snp", p_val = "b1_P", thresh = 1, targ = poly_GxE, person = "mz1", lab = "dz1_2_mz1_PRS")  
prs.typ <- summary(lm(DV ~  dz1_2_mz1_PRS, data = mz1.prs.typ))

prs.typMod <- summary(lm(DV ~  dz1_2_mz1_PRS + mod + dz1_2_mz1_PRS*mod, data = mz1.prs.typ))

# These results are directly in Table 7
prs.main
prs.mod
prs.both

prs.typ
prs.typMod



############################################################################################################
############################################################################################################
##
## Examine how the number of Main effect and GxE SNP relates to the twin model
## 
############################################################################################################
############################################################################################################

require(OpenMx)
require(MASS)
source("twinGxEsim_Fun.R")

#set.seed(11111)

Ns <- c(1, 10, 100, 1000)
Nm <- c(1, 10, 100, 1000)

PAR <- cbind(Ns, Nm)
RO <- 0

nReps <- 1000
nCond <- dim(PAR)[1]


res <- as.data.frame(matrix(NA, nCond*nReps, 39))
colnames(res) <- c("a11", "c11", "e11", "aM11", "cM11", "eM11", "interc", "reg1", "reg2",
                    "a11_se", "c11_se", "e11_se", "aM11_se", "cM11_se", "eM11_se", "interc_se", "reg1_se", "reg2_se",
                    "status", "Nsnps", "Nmods", 
					"new_Intercept", "new_main", "new_mod", "new_int", "new_Intercept_se", "new_main_se", "new_mod_se", "new_int_se", "new_r2",
					"typ_Intercept", "typ_main", "typ_mod", "typ_int", "typ_Intercept_se", "typ_main_se", "typ_mod_se", "typ_int_se", "new_r2")


res2 <- as.data.frame(matrix(NA, 1000, 39)) # changed
colnames(res2) <- c("a11", "c11", "e11", "aM11", "cM11", "eM11", "interc", "reg1", "reg2",
                    "a11_se", "c11_se", "e11_se", "aM11_se", "cM11_se", "eM11_se", "interc_se", "reg1_se", "reg2_se",
                    "status", "Nsnps", "Nmods", 
					"new_Intercept", "new_main", "new_mod", "new_int", "new_Intercept_se", "new_main_se", "new_mod_se", "new_int_se", "new_r2",
					"typ_Intercept", "typ_main", "typ_mod", "typ_int", "typ_Intercept_se", "typ_main_se", "typ_mod_se", "typ_int_se", "new_r2")
					


start <- Sys.time()
Sys.time()

#for (p in 1:nCond){
#	print(PAR[p,])
	for(i in 1:nReps){
		if(i == 250) print("250")
		if(i == 250) print(Sys.time())
			if(i == 500) print("500")
			if(i == 500) print(Sys.time())
				if(i == 750) print("750")
				if(i == 750) print(Sys.time())
					if(i == 999) print("999")
					if(i == 999) print(Sys.time())

RO <- RO + 1
#DAT   <- gxeSim(N = 10000, DZr = 0, MZr = 0, nSNPs = PAR[p,1], nModSNPs = PAR[p,2], binMod = T, b0 = 0, b2 = 0, snpSIGN = F, modSIGN = F)	
DAT   <- gxeSim(N = 20000, DZr = 0, MZr = 0, nSNPs = 100, nModSNPs = 100, binMod = T, b0 = 0, b2 = 0, snpSIGN = F, modSIGN = F)	

# run the twin model and save the estimates and standard errors
repFit <- gxeACE(mzData = DAT$mzData, dzData = DAT$dzData)
twin   <- c(repFit$output$estimate, repFit$output$standardErrors, repFit$output$status$code, PAR[3,])# changed

# PRS

# Run "discovery" gwas on DZ twin 1 (MZ twin 1 will be the target sample))
DZ1 <- GxEgwas(Sims = DAT, person = "dz1")
typDZ1 <- TypGWAS(Sims = DAT, person = "dz1", incl_mod = T)

# calculate statistical significance
DZ1 <-  sig(result = DZ1, focus = "b1_snp", focus_se = "v_snp", SE = F, lab = "b1")
DZ1 <-  sig(result = DZ1, focus = "b3_inter", focus_se = "v_inter", SE = F, lab = "b3")
typDZ1 <-  sig(result = DZ1, focus = "b1_snp", focus_se = "v_snp", SE = F, lab = "b1")

# calculate PRSs 
mz1.prs.main  <- prs(disc = DZ1, param = "b1_snp", p_val = "b1_P", thresh = 1, targ = DAT, person = "mz1", lab = "dz1_2_mz1_PRS")  
mz1.prs.int  <- prs(disc = DZ1, param = "b3_inter", p_val = "b3_P", thresh = 1, targ = DAT, person = "mz1", lab = "dz1_2_mz1_PRS")  
mz1.prs.int$dz1_2_mz1_PRS <- mz1.prs.int$dz1_2_mz1_PRS* mz1.prs.int$mod
mz1.prs <- cbind(mz1.prs.main, mz1.prs.int[,3])
colnames(mz1.prs) <- c("DV", "mod", "main", "int")


mz1.prs.typ  <- prs(disc = typDZ1, param = "b1_snp", p_val = "b1_P", thresh = 1, targ = DAT, person = "mz1", lab = "dz1_2_mz1_PRS")  


prs.both <- summary(lm(DV ~  main + mod + int, data = mz1.prs))
prs.typMod <- summary(lm(DV ~  dz1_2_mz1_PRS + mod + dz1_2_mz1_PRS*mod, data = mz1.prs.typ))


# These results are directly in the table
PRS <- c(coef(prs.both)[,1], coef(prs.both)[,2], prs.both$r.squared,
coef(prs.typMod)[,1], coef(prs.typMod)[,2], prs.typMod$r.squared)


res2[i,] <- c(twin, PRS)# changed

	}
	#}
end <- Sys.time()
end-start

#save.image("SimMainInt.Rdata")
#load("SimMainInt.Rdata")

############################################################################################
############################################################################################
############################################################################################

twina.1     <- abs(res$a11[res$Nsnps == 1 ])
twina.10    <- abs(res$a11[res$Nsnps == 10 ])
twina.100   <- abs(res$a11[res$Nsnps == 100 ])
twina.1000  <- abs(res$a11[res$Nsnps == 1000 ])

twinaM.1    <- abs(res$aM11[res$Nsnps == 1 ])
twinaM.10   <- abs(res$aM11[res$Nsnps == 10 ])
twinaM.100  <- abs(res$aM11[res$Nsnps == 100 ])
twinaM.1000 <- abs(res$aM11[res$Nsnps == 1000 ])

new_main.1    <- (res$new_main[res$Nsnps == 1 ])
new_main.10   <- (res$new_main[res$Nsnps == 10 ])
new_main.100  <- (res$new_main[res$Nsnps == 100 ])
new_main.1000 <- (res$new_main[res$Nsnps == 1000 ])

new_int.1     <- (res$new_int[res$Nsnps == 1 ])
new_int.10    <- (res$new_int[res$Nsnps == 10 ])
new_int.100   <- (res$new_int[res$Nsnps == 100 ])
new_int.1000  <- (res$new_int[res$Nsnps == 1000 ])

typ_main.1    <- (res$typ_main[res$Nsnps == 1 ])
typ_main.10   <- (res$typ_main[res$Nsnps == 10 ])
typ_main.100  <- (res$typ_main[res$Nsnps == 100 ])
typ_main.1000 <- (res$typ_main[res$Nsnps == 1000 ])

typ_int.1     <- (res$typ_int[res$Nsnps == 1 ])
typ_int.10    <- (res$typ_int[res$Nsnps == 10 ])
typ_int.100   <- (res$typ_int[res$Nsnps == 100 ])
typ_int.1000  <- (res$typ_int[res$Nsnps == 1000 ])



### Constructing Figure 3
require(vioplot)

par(xaxt = "n", mfrow = c(1, 3))

vioplot(twina.1, twina.10, twina.100, twina.1000, 
	    col = "blue", side = "left" , 
		plotCentre = "line", ylim = c(.4, .8), ylab = " ", xlab = " ")
				
		mtext(c(1, 10, 100, 1000), 1, line = .75, at = 1:4)
		mtext("Number of SNPs", 1, line =2.5, at = 2.5)
		mtext("Unstandardized Twin Parameters", 2, line =2.5, at = .6)
		
#ylab = "Unstandardized Twin Parameters", xlab = "Number of SNPs"
		
		
vioplot(twinaM.1, twinaM.10, twinaM.100, twinaM.1000, 
	    col = "red", side = "right" , 
		plotCentre = "line", ylim = c(.4, .8),  add = T)

abline(h = .71, col = "red", lwd = 3)

legend("bottomleft", fill = c("blue", "red"), 
       legend = c("Main Effect", "Interaction Effect"), title = " ",
				  cex = .7, box.col = "white", border = "white")


### Constructing Figure (GxE GWAS PRS)

vioplot(new_main.1, new_main.10, new_main.100, new_main.1000, 
	    col = "blue", side = "left" , 
		plotCentre = "line", ylim = c(.5, 1.1), ylab = " ", xlab = " ")
				
		mtext(c(1, 10, 100, 1000), 1, line = .75, at = 1:4)
		mtext("Number of SNPs", 1, line =2.5, at = 2.5)
		mtext("Unstandardized PRS (GxE GWAS)", 2, line =2.5, at = .85)
		
		
		
vioplot(new_int.1, new_int.10, new_int.100, new_int.1000, 
	    col = "red", side = "right" , 
		plotCentre = "line", ylim = c(.5, 1.1),  add = T)

abline(h = 1, col = "red", lwd = 3)

legend("bottomleft", fill = c("blue", "red"), 
       legend = c("Main Effect", "Interaction Effect"), title = " ",
				  cex = .7, box.col = "white", border = "white")


### Constructing Figure  (Standard GWAS PRS)

vioplot(typ_main.1, typ_main.10, typ_main.100, typ_main.1000, 
	    col = "blue", side = "left" , 
		plotCentre = "line", ylim = c(.5, 1.1), ylab = " ", xlab = " ")
				
		mtext(c(1, 10, 100, 1000), 1, line = .75, at = 1:4)
		mtext("Number of SNPs", 1, line =2.5, at = 2.5)
		mtext("Unstandardized PRS (Standard)", 2, line =2.5, at = .85)
		
		
		
vioplot(typ_int.1, typ_int.10, typ_int.100, typ_int.1000, 
	    col = "red", side = "right" , 
		plotCentre = "line", ylim = c(.5, 1.1),  add = T)

abline(h = 1, col = "red", lwd = 3)




legend("bottomleft", fill = c("blue", "red"), 
       legend = c("Main Effect", "Interaction Effect"), title = " ",
				  cex = .7, box.col = "white", border = "white")





############################################################################################################
############################################################################################################
##
## Study 4: SIMULATE GENOMIC DATA WITH POLYGENIC VARIANTS and GxE from a normal distribution
## 
############################################################################################################
############################################################################################################

set.seed(11111)
norm_GxE   <- gxeSim(N = 10000, DZr = 0, MZr = 0, nSNPs = 1000, nModSNPs = 1000, binMod = T, b0 = 0, b2 = 0, snpSIGN = T, modSIGN = T, effDist = "norm")
norm_GxEFit <- gxeACE(mzData = norm_GxE$mzData, dzData = norm_GxE$dzData)
#norm_GxEFit <- mxRun(norm_GxEFit)
summary(norm_GxEFit)
norm_GxEres <- gxePlot(aceModFit = norm_GxEFit, binMod = T, title = " ", legend = T, legend.loc = "topleft")


# Run "discovery" gwas on DZ twin 1 (MZ twin 1 will be the target sample))

DZ1 <- GxEgwas(Sims = norm_GxE, person = "dz1")
typDZ1 <- TypGWAS(Sims = norm_GxE, person = "dz1", incl_mod = T)

# calculate statistical significance

DZ1 <-  sig(result = DZ1, focus = "b1_snp", focus_se = "v_snp", SE = F, lab = "b1")
DZ1 <-  sig(result = DZ1, focus = "b3_inter", focus_se = "v_inter", SE = F, lab = "b3")
typDZ1 <-  sig(result = DZ1, focus = "b1_snp", focus_se = "v_snp", SE = F, lab = "b1")

mean(DZ1 $b1_snp)
mean(DZ1 $b3_inter)

table(DZ1 $b1_P < 5e-8)
table(DZ1 $b3_P < 5e-8)

DZ1[DZ1$b1_P < 5e-8,]

table(DZ1 $b1_P < .05)
table(DZ1 $b3_P < .05)

table(DZ1 $b1_P < .05/1000)
table(DZ1 $b3_P < .05/1000)

DZ1[DZ1$b1_P < .05/1000,]
DZ1[DZ1$b3_P < .05/1000,]

#################

# calculate PRS and run linear regression
mz1.prs.main  <- prs(disc = DZ1, param = "b1_snp", p_val = "b1_P", thresh = 1, targ = norm_GxE, person = "mz1", lab = "dz1_2_mz1_PRS")  
prs.main <- summary(lm(DV ~  dz1_2_mz1_PRS, data = mz1.prs.main))

mz1.prs.int  <- prs(disc = DZ1, param = "b3_inter", p_val = "b3_P", thresh = 1, targ = norm_GxE, person = "mz1", lab = "dz1_2_mz1_PRS")  
mz1.prs.int$dz1_2_mz1_PRS <- mz1.prs.int$dz1_2_mz1_PRS* mz1.prs.int$mod
prs.mod <- summary(lm(DV ~  dz1_2_mz1_PRS, data = mz1.prs.int))

mz1.prs <- cbind(mz1.prs.main, mz1.prs.int[,3])
colnames(mz1.prs) <- c("DV", "mod", "main", "int")
prs.both <- summary(lm(DV ~  main + mod + int, data = mz1.prs))


mz1.prs.typ  <- prs(disc = typDZ1, param = "b1_snp", p_val = "b1_P", thresh = 1, targ = norm_GxE, person = "mz1", lab = "dz1_2_mz1_PRS")  
prs.typ <- summary(lm(DV ~  dz1_2_mz1_PRS, data = mz1.prs.typ))

prs.typMod <- summary(lm(DV ~  dz1_2_mz1_PRS + mod + dz1_2_mz1_PRS*mod, data = mz1.prs.typ))

# These results are directly in the table
prs.main
prs.mod
prs.both

prs.typ
prs.typMod


