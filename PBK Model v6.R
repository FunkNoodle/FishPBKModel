### ><(((º> ><(((º> ><(((º> ><(((º> ><(((º> ><(((º> ><(((º> ><(((º> ><(((º> ### 

#Clearing the environment
rm(list = ls())

### Input compound name ###

# Example uses Carbamazepine as an input
target_name <- "Example"

### Input parameters ###

PBKinput <- read.csv("PBKinput.csv")

# Find the row with the target compound name
row <- PBKinput[PBKinput$Name == target_name, ]

# Extract values from the row and assign them to variables
if (nrow(row) == 1) {
  Name <- row$Name
  BM <- row$BM
  Sex <- row$Sex
  temp <- row$temp
  pHw <- row$pHw
  pHf <- row$pHf
  pKa <- row$pKa
  logKow <- (row$logKow)
  logKow_ion <- ( row$logKow_ion)
  ii <- row$ii
  HL <- row$HL
  Cl_int <- row$Cl_int
  UF <- row$UF
  C_water <- row$C_water
  Q_ingest <- row$Q_ingest
  
} else {
  stop("Name not found or multiple entries for the same name.")
}

#Amount of time between timesteps (days)
#If HL is used, use fraction of HL; if Cl_int is used, use 1 second

if (!is.null(Cl_int) && Cl_int != 0) {
  dt = (1/30000)/Cl_int
  
}else{ 
  dt = HL/10000 #Increase this number if you get oscillations
}

# Converting logKow to Kow and changing ionic Kow to Kow if it is missing

Kow <- 10^logKow

if(is.na(logKow_ion)) {
  Kow_ion = Kow
} else {
  Kow_ion <- 10^logKow_ion
  
}

#Setting internal pH of the Fish to 7.4 if no input is found

if(is.null(pHf)) {
  pHf = 7.4
}

### Tissue abbreviations ###

#Fat: fat
#Total blood: bl
#Arterial blood: art 
#Venous blood: ven 
#Brain: br 
#Kidney: k 
#Liver: l 
#Gastrointestinal tract: git 
#Gonads: go 
#Richly perfused tissues: rpt
#Poorly perfused tissues: ppt
#Skin: s
#Skeleton: sk
#Muscle: mu

### Parameters ###

#Weighting factors for blood flow per tissue
WF_fat = 0.626 
WF_br = 3.6 
WF_k = 7.005 
WF_l = 2.23 
WF_git = 4.846 
WF_go = 3.367 
WF_pp = 0.733 
WF_rp = 3.6 
WF_s = 0.57
#Fraction of ppt blood flow to kidney
a_fpp = 0.6  
#Fraction of skin blood flow to kidney
a_fs = 0.9 
#Temperature in Celsius 
Temp_C = temp 
#Temperature in Kelvin
Temp_K = (temp + 273.15) 

#Volume per tissue
V_fat = (0.020 * BM ^ 1.24) 
V_bl = (0.030 * BM ^ 1.02) 
V_art = ((0.030 * BM ^ 1.02)*(1/3)) 
V_ven = ((0.030 * BM ^ 1.02)*(2/3)) 
V_br = (0.017 * BM ^ 0.61) 
V_k = 0.014 * BM ^ 0.87
V_l = 0.018 * BM ^ 0.95
V_git = 0.068 * BM ^ 0.96
V_sk = 0.033 * BM ^ 1.03

#Gonads volume: Female, Male, or average
if (Sex == 2){ 
  V_go = 0.13 * BM ^ 1.04
}else if (Sex == 1){
  V_go = 0.012 * BM ^ 1.06
}else{
  V_go = ((0.13 * BM ^ 1.04) + (0.012 * BM ^ 1.06))/2
}

V_rp = 0.0025 * BM ^ 1.04
V_s = 0.10 * BM ^ 0.67
V_pp = BM - (V_fat + V_art + V_ven + V_br + V_k + V_l + V_git + V_go + V_rp +
               V_s)
V_mu = V_pp-V_sk 
#Fraction muscle of ppt
mu_ppt = V_mu/V_pp
#Cardiac output 
F_card = exp(25.5) * exp(-5790/Temp_K) * BM ^ 0.75 

### Blood flow per tissue ###

#Blood flow to fat
F_fat = ((WF_fat * (V_fat/BM))/((WF_br * (V_br/BM)) + (WF_k * (V_k/BM)) +
                                  (WF_l * (V_l/BM)) + (WF_git * (V_git/BM)) + (WF_go * (V_go/BM)) +
                                  (WF_rp * (V_rp/BM)) + (WF_pp * (V_pp/BM)) + (WF_s * (V_s/BM)) + 
                                  (WF_fat * (V_fat/BM)))) * F_card
#Blood flow to brain
F_br = ((WF_br * (V_br/BM))/((WF_br * (V_br/BM)) + (WF_k * (V_k/BM)) + 
                               (WF_l * (V_l/BM)) + (WF_git * (V_git/BM)) + (WF_go * (V_go/BM)) + 
                               (WF_rp * (V_rp/BM)) + (WF_pp * (V_pp/BM)) + (WF_s * (V_s/BM)) + 
                               (WF_fat * (V_fat/BM)))) * F_card
#Blood flow to kidney
F_k = ((WF_k * (V_k/BM))/((WF_br * (V_br/BM)) + (WF_k * (V_k/BM)) + 
                            (WF_l * (V_l/BM)) + (WF_git * (V_git/BM)) + (WF_go * (V_go/BM)) +
                            (WF_rp * (V_rp/BM)) + (WF_pp * (V_pp/BM)) + (WF_s * (V_s/BM)) + 
                            (WF_fat * (V_fat/BM)))) * F_card
#Blood flow to liver
F_l = ((WF_l * (V_l/BM))/((WF_br * (V_br/BM)) + (WF_k * (V_k/BM)) + 
                            (WF_l * (V_l/BM)) + (WF_git * (V_git/BM)) + (WF_go * (V_go/BM)) + 
                            (WF_rp * (V_rp/BM)) + (WF_pp * (V_pp/BM)) + (WF_s * (V_s/BM)) + 
                            (WF_fat * (V_fat/BM)))) * F_card
#Blood flow to gastrointestinal tract
F_git = ((WF_git * (V_git/BM))/((WF_br * (V_br/BM)) + (WF_k * (V_k/BM)) + 
                                  (WF_l * (V_l/BM)) + (WF_git * (V_git/BM)) + (WF_go * (V_go/BM)) +
                                  (WF_rp * (V_rp/BM)) + (WF_pp * (V_pp/BM)) + (WF_s * (V_s/BM)) + 
                                  (WF_fat * (V_fat/BM)))) * F_card
#Blood flow to gonads
F_go = ((WF_go * (V_go/BM))/((WF_br * (V_br/BM)) + (WF_k * (V_k/BM)) + 
                               (WF_l * (V_l/BM)) + (WF_git * (V_git/BM)) + (WF_go * (V_go/BM)) + 
                               (WF_rp * (V_rp/BM)) + (WF_pp * (V_pp/BM)) + (WF_s * (V_s/BM)) + 
                               (WF_fat * (V_fat/BM)))) * F_card
#Blood flow to richly perfused tissues
F_rp = ((WF_rp * (V_rp/BM))/((WF_br * (V_br/BM)) + (WF_k * (V_k/BM)) + 
                               (WF_l * (V_l/BM)) + (WF_git * (V_git/BM)) + (WF_go * (V_go/BM)) + 
                               (WF_rp * (V_rp/BM)) + (WF_pp * (V_pp/BM)) + (WF_s * (V_s/BM)) + 
                               (WF_fat * (V_fat/BM)))) * F_card
#Blood flow to poorly perfused tissues
F_pp = ((WF_pp * (V_pp/BM))/((WF_br * (V_br/BM)) + (WF_k * (V_k/BM)) + 
                               (WF_l * (V_l/BM)) + (WF_git * (V_git/BM)) + (WF_go * (V_go/BM)) + 
                               (WF_rp * (V_rp/BM)) + (WF_pp * (V_pp/BM)) + (WF_s * (V_s/BM)) + 
                               (WF_fat * (V_fat/BM)))) * F_card
#Blood flow to skin
F_s = ((WF_s * (V_s/BM))/((WF_br * (V_br/BM)) + (WF_k * (V_k/BM)) + 
                            (WF_l * (V_l/BM)) + (WF_git * (V_git/BM)) + (WF_go * (V_go/BM)) + 
                            (WF_rp * (V_rp/BM)) + (WF_pp * (V_pp/BM)) + (WF_s * (V_s/BM)) + 
                            (WF_fat * (V_fat/BM)))) * F_card

#Oxygen consumption rate 
VO2 = exp(19.1) * exp(-5040/Temp_K) * BM ^ 0.75 
#Dissolved oxygen concentration
C_ox = exp(-10.6) * exp(1743.3/Temp_K) 
#Fraction of the neutral:ionised compounds at gill surface and in fish body
if (ii == 0){ 
  Fn_fish = 1
  Fn_gill = 1
}else{ 
  Fn_fish = 1/(1 + (10 ^ (ii * (pHf-pKa)))) 
  Fn_gill = 1/(1 + (10 ^ (ii * (pHw-pKa)))) 
} 

#Plasma:water distribution ratio
D_plw = 10 ^ ((0.75 * Kow_ion) + 0.58)  
#Octanol-water distribution coefficient at fish body pH
Dow = (Fn_fish * Kow) + (1-Fn_fish) * Kow_ion  
#Octanol-water distribution coefficient at gill surface pH
Dow_gill = Fn_gill * Kow + (1-Fn_gill) * Kow_ion  
#Unbound fraction in blood
if (is.null(UF) || UF == 0) {
  UF <- 1 / (1 + D_plw)
}

#Volume divided by mass per tissue
ViM_fat = V_fat/BM
ViM_bl = V_bl/BM
ViM_br = V_br/BM
ViM_k = V_k/BM
ViM_l = V_l/BM
ViM_git = V_git/BM
ViM_go = V_go/BM
ViM_pp = V_pp/BM
ViM_rp = V_rp/BM
ViM_s = V_s/BM
#Fraction neutral lipids per tissue
ViM_nlfat = (88.4/100)* ViM_fat * BM
ViM_nlbl = (0.8/100)* ViM_bl * BM
ViM_nlbr = (4.7/100)* ViM_br * BM
ViM_nlk = (5.9/100)* ViM_k * BM
ViM_nll = (3.7/100)* ViM_l * BM
ViM_nlgit = (5.2/100)* ViM_git * BM
ViM_nlgo = (4.2/100)* ViM_go * BM
ViM_nlrp = (6.6/100)* ViM_rp * BM
ViM_nlpp = (2.4/100)* ViM_pp * BM
ViM_nls = (3.7/100)* ViM_s * BM
#Fraction polar lipids per tissue
ViM_plfat = (0.9/100)* ViM_fat * BM
ViM_plbl = (0.7/100)* ViM_bl * BM
ViM_plbr = (4.3/100)* ViM_br * BM
ViM_plk = (2.3/100)* ViM_k * BM
ViM_pll = (2.5/100)* ViM_l * BM
ViM_plgit = (2.3/100)* ViM_git * BM
ViM_plgo = (2.4/100)* ViM_go * BM
ViM_plrp = (0.8/100)* ViM_rp * BM
ViM_plpp = (0.9/100)* ViM_pp * BM
ViM_pls = (1.2/100)* ViM_s * BM
#Fraction non-lipid organic matter per tissue
ViM_nlomfat = (4.8/100)* ViM_fat * BM
ViM_nlombl = (13.4/100)* ViM_bl * BM
ViM_nlombr = (22.8/100)* ViM_br * BM
ViM_nlomk = (21.5/100)* ViM_k * BM
ViM_nloml = (23.4/100)* ViM_l * BM
ViM_nlomgit = (25.5/100)* ViM_git * BM
ViM_nlomgo = (26.2/100)* ViM_go * BM
ViM_nlomrp = (33.8/100)* ViM_rp * BM
ViM_nlompp = (19.9/100)* ViM_pp * BM
ViM_nloms = (25.9/100)* ViM_s * BM
#Fraction water per tissue
ViM_watfat = (5.9/100)* ViM_fat * BM
ViM_watbl = (85.1/100)* ViM_bl * BM
ViM_watbr = (68.2/100)* ViM_br * BM
ViM_watk = (70.3/100)* ViM_k * BM
ViM_watl = (70.4/100)* ViM_l * BM
ViM_watgit = (67/100)* ViM_git * BM
ViM_watgo = (67.2/100)* ViM_go * BM
ViM_watrp = (58.7/100)* ViM_rp * BM
ViM_watpp = (76.9/100)* ViM_pp * BM
ViM_wats = (69.2/100)* ViM_s * BM
#Total neutral lipids
nl_total = (ViM_nlfat + ViM_nlbl + ViM_nlbr + ViM_nlk + ViM_nll + ViM_nlgit + 
              ViM_nlgo + ViM_nlrp + ViM_nlpp + ViM_nls)/BM
#Total polar lipids
pl_total = (ViM_plfat + ViM_plbl + ViM_plbr + ViM_plk + ViM_pll + ViM_plgit + 
              ViM_plgo + ViM_plrp + ViM_plpp + ViM_pls)/BM
#Total non-lipid organic matter
nlom_total = (ViM_nlomfat + ViM_nlombl + ViM_nlombr + ViM_nlomk + ViM_nloml + 
                ViM_nlomgit + ViM_nlomgo + ViM_nlomrp + ViM_nlompp + ViM_nloms)/BM
#Total water                
wat_total = (ViM_watfat + ViM_watbl + ViM_watbr + ViM_watk + ViM_watl + 
               ViM_watgit + ViM_watgo + ViM_watrp + ViM_watpp + ViM_wats)/BM
#Partition coefficients per tissue and BCF: unionized compounds
if(ii == 0){
  P_bw = 0.008 * Kow + 0.007 * (Kow ^ 0.94) + 0.134 * (Kow ^ 0.63) + 0.851
  P_gob = (0.042 * Kow + 0.024 * (Kow ^ 0.94) + 0.262 * (Kow ^ 0.63) + 
             0.672) / P_bw
  P_kb = (0.059 * Kow + 0.023 * (Kow ^ 0.94) + 0.215 * (Kow ^ 0.63) + 
            0.703) / P_bw
  P_lb = (0.037 * Kow + 0.025 * (Kow ^ 0.94) + 0.234 * (Kow ^ 0.63) + 
            0.704) / P_bw
  P_ppb = (0.024 * Kow + 0.009 * (Kow ^ 0.94) + 0.199 * (Kow ^ 0.63) + 
             0.769) / P_bw
  P_rpb = (0.066 * Kow + 0.008 * (Kow ^ 0.94) + 0.338 * (Kow ^ 0.63) + 
             0.587) / P_bw
  P_sb = (0.037 * Kow + 0.012 * (Kow ^ 0.94) + 0.259 * (Kow ^ 0.63) + 
            0.692) / P_bw
  P_fatb = (0.884 * Kow + 0.009 * (Kow ^ 0.94) + 0.048 * (Kow ^ 0.63) + 
              0.059) / P_bw
  P_brb = (0.047 * Kow + 0.043 * (Kow ^ 0.94) + 0.228 * (Kow ^ 0.63) + 
             0.682) / P_bw
  P_gitb = (0.052 * Kow + 0.023 * (Kow ^ 0.94) + 0.255 * (Kow ^ 0.63) + 
              0.670) / P_bw
  BCF = ((nl_total)*Kow + (pl_total) *(Kow^0.94) + (nlom_total) * (Kow^0.63) + 
           (wat_total))
  
  #Partition coefficients per tissue and BCF: ionized compounds 
}else{
  P_bw = 0.008 * 0.3 * Dow + 0.007 * 2.0 * (Dow ^ 0.94) + 
    0.134 * 2.9 * (Dow^ 0.63) + 0.851
  P_gob = (0.042 * 0.3 * Dow + 0.024 * 2.0 * (Dow ^ 0.94) + 
             0.262 * 2.9 * (Dow^ 0.63) + 0.672)/P_bw
  P_kb = (0.059 * 0.3 * Dow + 0.023 * 2.0 * (Dow ^ 0.94) + 0.215 * 2.9 *
            (Dow^ 0.63) + 0.703)/P_bw
  P_lb = (0.037 * 0.3 * Dow + 0.025 * 2.0 * (Dow ^ 0.94) + 0.234 * 2.9 *
            (Dow^ 0.63) + 0.704)/P_bw
  P_ppb = (0.024 * 0.3 * Dow + 0.009 * 2.0 * (Dow ^ 0.94) + 0.199 * 2.9 *
             (Dow ^ 0.63) + 0.769)/P_bw
  P_rpb = (0.066 * 0.3 * Dow + 0.008 * 2.0 * (Dow ^ 0.94) + 0.338 * 2.9 *
             (Dow ^ 0.63) + 0.587)/P_bw
  P_sb = (0.037 * 0.3 * Dow + 0.012 * 2.0 * (Dow ^ 0.94) + 0.259 * 2.9 *
            (Dow^ 0.63) + 0.692)/P_bw
  P_fatb = (0.884 * 0.3 * Dow + 0.009 * 2.0 * (Dow ^ 0.94) + 0.048 * 2.9 *
              (Dow ^ 0.63) + 0.059)/P_bw
  P_brb = (0.047 * 0.3* Dow + 0.043 * 2.0 * (Dow ^ 0.94) + 0.228 * 2.9 *
             (Dow^ 0.63) + 0.682)/P_bw
  P_gitb = (0.052 * 0.3 * Dow + 0.023 * 2.0 * (Dow ^ 0.94) + 0.255 * 2.9 *
              (Dow ^ 0.63) + 0.670)/P_bw
  BCF = ((nl_total) * 0.3 * Dow + (pl_total) * 2.0 * (Dow^0.94) + (nlom_total) *
           2.9 * (Dow^0.63) + (wat_total))
}
#Gill ventilation coefficient 
Gamma_water = (VO2/(0.71 * C_ox)) * (1/((1000 ^ 0.25) * (BM ^ 0.75))) 
#Blood perfusion coefficient
Gamma_blood = F_card * P_bw * (1/((1000 ^ 0.25) * (BM ^ 0.75))) 
#Exchange coefficient between blood and water 
Kx = (((BM/1000) ^ 0.75)/((2.8 * 10 ^ -3) + (68/Dow_gill) + (1/Gamma_water) +
                            (1/Gamma_blood))) * 1000 
#Uptake rate constant 
Ku = (0.8/(1-0.8) * (1/(0.03 * ((BM/1000) ^ 0.04) * (Dow - 1) + 1)) *
        (((BM/1000) ^ -0.25)/((1.1 * 10 ^ -5) + (68/Dow) + (1/(0.03 *
                                                                 ((BM/1000) ^ 0.04) * Dow * (1-0.8) * 0.005)))))
#Excretion rateconstant 
Ke_feces = (1/((0.03 * (BM/1000) ^ 0.04) * (Dow - 1) + 1)) *
  ((BM/1000) ^ -0.25)/(1.1 * 10 ^ -5 + (68/Dow) + (1/(0.03 *
                                                        ((BM/1000) ^ 0.04) * Dow * (1-0.8) * 0.005)))
#Whole-body primary biotransformation rate constant 
#Checks if intrinsic hepatic clearence is provided to calculate hepatic clearance
#Otherwise normalised to a 10-g fish at 15 ℃ 
if (!is.null(Cl_int) && Cl_int != 0){ 
  Cl_hep = (F_l * UF * Cl_int)/(F_l + UF * Cl_int)
}else{ 
  #Normalised to a 10-g fish at 15 ℃ 
  Km_n = (log(2))/HL 
  #Corrected for body mass and temperature differences
  Km_x = Km_n * ((BM/10)^-0.25) * exp(-874*((1/Temp_K)-(1/288))) 
  #Apparent volume of distribution
  V_D = BCF/P_bw 
  #Hepatic clearance
  Cl_hep = Km_x * V_D 
} 

#Fraction assimilated from food
f_abs = Ku/(0.005*1*(BM/1000)^-0.25) 

### PBK Model over time ###

#Setting all tissue concentrations and quantities to 0
Q_fat = 0
Q_br = 0
Q_k = 0
Q_l = 0
Q_gitlumen = 0
Q_git = 0
Q_go = 0
Q_rp = 0
Q_pp = 0
Q_s = 0
Q_art = 0
Q_ven = 0
C_ven = 0
C_art = 0
C_fat = 0
C_br = 0
C_k = 0
C_l = 0
C_git = 0
C_go = 0
C_rp = 0
C_pp = 0
C_s = 0
C_total = 0
C_mu = 0

# Running the loop until steady state
i <- 1
max_iter <- 100000 # Prevent infinite loop
C_art_start = (Kx * (C_water - (C_ven * UF)/P_bw))/F_card #Used to find the steady state
# Define how close to zero is "steady state" based on starting arterial concentrations
# 0.01% of starting plasma concentration
tolerance <- (C_art_start/10000)

#Calculating tissue concentrations and quantities, and putting values into "PBKoutput"
PBKoutput <- data.frame(
  Timestep= seq(0, (dt*max_iter), by = dt), 
  C_ven= numeric(max_iter+1),
  C_art= numeric(max_iter+1),
  C_fat= numeric(max_iter+1),
  C_br= numeric(max_iter+1),
  C_k= numeric(max_iter+1),
  C_l= numeric(max_iter+1),
  C_git= numeric(max_iter+1),
  C_go= numeric(max_iter+1),
  C_rp= numeric(max_iter+1),
  C_pp= numeric(max_iter+1),
  C_mu= numeric(max_iter+1),
  C_s= numeric(max_iter+1),
  C_total= numeric(max_iter+1),
  
  dQ_fat= numeric(max_iter+1),
  dQ_br= numeric(max_iter+1),
  dQ_k= numeric(max_iter+1),
  dQ_l= numeric(max_iter+1),
  dQ_git= numeric(max_iter+1),
  dQ_go= numeric(max_iter+1),
  dQ_rp= numeric(max_iter+1),
  dQ_pp= numeric(max_iter+1),
  dQ_s= numeric(max_iter+1),
  
  Q_fat= numeric(max_iter+1),
  Q_br= numeric(max_iter+1),
  Q_k= numeric(max_iter+1),
  Q_l= numeric(max_iter+1),
  Q_git= numeric(max_iter+1),
  Q_go= numeric(max_iter+1),
  Q_rp= numeric(max_iter+1),
  Q_pp= numeric(max_iter+1),
  Q_s= numeric(max_iter+1),
  Q_art= numeric(max_iter+1),
  Q_ven= numeric(max_iter+1)
)


try({
  while (i <= max_iter) {
    
    dQ_fat = F_fat * (C_art - (Q_fat/(V_fat * P_fatb)))
    dQ_br = F_br * (C_art - (Q_br/(V_br * P_brb)))
    dQ_go = F_go * (C_art - (Q_go/(V_go * P_gob)))
    dQ_rp = F_rp * (C_art - (Q_rp/(V_rp * P_rpb))) 
    dQ_pp = F_pp * (C_art - (Q_pp/(V_pp * P_ppb))) 
    dQ_s = F_s * (C_art - (Q_s/(V_s * P_sb))) 
    
    dQ_gitlumen = (f_abs * Q_ingest) - ((Ku + Ke_feces) * Q_gitlumen)
    dQ_git = Ku * Q_gitlumen + F_git * (C_art - (Q_git/(V_git * P_gitb)))
    dQ_l = F_l * (C_art - (Q_l/(V_l * P_lb))) + F_rp * ((Q_rp/(V_rp * P_rpb))-
                                                          (Q_l/(V_l * P_lb))) + F_git * ((Q_git/(V_git * P_gitb))-
                                                                                           (Q_l/(V_l * P_lb))) + F_go * ((Q_go/(V_go * P_gob))-(Q_l/(V_l * P_lb))) -
      Cl_hep * BM * (Q_l/(V_l * P_lb))
    dQ_k = F_k * (C_art - (Q_k/(V_k * P_kb))) + (1- a_fpp) * F_pp * ((Q_pp/(V_pp *
                                                                              P_ppb))-(Q_k/(V_k * P_kb))) + (1 - a_fs) * F_s * ((Q_s/(V_s * P_sb))-
                                                                                                                                  (Q_k/(V_k * P_kb)))
    
    Q_fat = dQ_fat*dt + Q_fat
    Q_br = dQ_br*dt + Q_br
    Q_go = dQ_go*dt + Q_go
    Q_rp = dQ_rp*dt + Q_rp
    Q_pp = dQ_pp*dt + Q_pp
    Q_s = dQ_s*dt + Q_s
    Q_gitlumen = dQ_gitlumen*dt + Q_gitlumen
    Q_git = dQ_git*dt + Q_git
    Q_l = dQ_l*dt + Q_l
    Q_k = dQ_k*dt + Q_k
    
    C_ven = (((1-a_fs) * F_s * (Q_s/(V_s * P_sb))) + ((F_k + a_fpp*F_pp +
                                                         a_fs*F_s) * (Q_k/(V_k * P_kb))) + ((1-a_fpp) * F_pp * (Q_pp/(V_pp *
                                                                                                                        P_ppb))) + (F_fat * (Q_fat/(V_fat * P_fatb))) + (F_br * (Q_br/(V_br *
                                                                                                                                                                                         P_brb))) + ((F_l + F_rp + F_go + F_git) * (Q_l/(V_l * P_lb))))/ F_card
    C_art = ((Kx * (C_water - (C_ven * UF)/P_bw))/F_card) + C_ven
    C_fat = Q_fat / V_fat
    C_br = Q_br / V_br
    C_k = Q_k / V_k
    C_l = Q_l / V_l
    C_git = Q_git / V_git
    C_go = Q_go / V_go
    C_rp = Q_rp / V_rp
    C_pp = Q_pp / V_pp
    C_s = Q_s / V_s
    C_mu = C_pp * mu_ppt
    C_total = (Q_fat + Q_br + Q_go + Q_rp + Q_pp + Q_s + Q_gitlumen + Q_git +
                 Q_l + Q_k + Q_art + Q_ven)/BM
    
    AchievedSteadyState = tolerance/((C_art*V_bl/3) - Q_art) #Used to find the steady state
    
    Q_art = C_art*V_bl/3
    Q_ven = C_ven*V_bl*2/3
    
    PBKoutput$C_ven[i] <- C_ven
    PBKoutput$C_art[i] <- C_art
    PBKoutput$C_fat[i] <- C_fat
    PBKoutput$C_br[i] <- C_br
    PBKoutput$C_k[i] <- C_k
    PBKoutput$C_l[i] <- C_l
    PBKoutput$C_git[i] <- C_git
    PBKoutput$C_go[i] <- C_go
    PBKoutput$C_rp[i] <- C_rp
    PBKoutput$C_pp[i] <- C_pp
    PBKoutput$C_mu[i] <- C_mu
    PBKoutput$C_s[i] <- C_s
    PBKoutput$C_total[i] <- C_total
    
    PBKoutput$dQ_fat[i] <- dQ_fat
    PBKoutput$dQ_br[i] <- dQ_br
    PBKoutput$dQ_k[i] <- dQ_k
    PBKoutput$dQ_l[i] <- dQ_l
    PBKoutput$dQ_git[i] <- dQ_git
    PBKoutput$dQ_go[i] <- dQ_go
    PBKoutput$dQ_rp[i] <- dQ_rp
    PBKoutput$dQ_pp[i] <- dQ_pp
    PBKoutput$dQ_s[i] <- dQ_s
    
    PBKoutput$Q_fat[i] <- Q_fat
    PBKoutput$Q_br[i] <- Q_br
    PBKoutput$Q_k[i] <- Q_k
    PBKoutput$Q_l[i] <- Q_l
    PBKoutput$Q_git[i] <- Q_git
    PBKoutput$Q_go[i] <- Q_go
    PBKoutput$Q_rp[i] <- Q_rp
    PBKoutput$Q_pp[i] <- Q_pp
    PBKoutput$Q_s[i] <- Q_s
    PBKoutput$Q_art[i] <- Q_art
    PBKoutput$Q_ven[i] <- Q_ven
    
    # Stop loop once steady state is reached based on AchievedSteadyState, ignoring the first 5000 loops
    
    if (i > 5000 && !is.na(AchievedSteadyState) && abs(AchievedSteadyState) > 1) {
      message("Steady state reached at iteration ", i)
      break
    }
    
    i <- i + 1
  }
})

#Trimming the data frame based on the number of iterations
PBKoutput <- PBKoutput[1:i-1, ]

###Plot example###

png(paste0(Name,".png"),width = 800, height = 600)

# Find the max and minimum concentrations from the concentrations to adjust the graph axis
max_ylim <- 10^ceiling(log10(max(PBKoutput[, grep("^C_", names(PBKoutput))], na.rm = TRUE)))
min_ylim <- 10^floor(log10(min(PBKoutput[, grep("^C_", names(PBKoutput))][PBKoutput[, grep("^C_", names(PBKoutput))] > 0], na.rm = TRUE)))
max_xlim <- max(PBKoutput$Timestep, na.rm = TRUE)

plot(PBKoutput$C_art ~ PBKoutput$Timestep, xlim = c(0, max_xlim), ylim = c(min_ylim, max_ylim), 
     type = 'l', col = "red", lwd = 1.5, log = 'y', 
     xlab = "Time (d)", ylab = "Concentration (μg/g)", 
     main = paste("Concentration vs Time for", Name))

# Add additional lines
lines(PBKoutput$C_br ~ PBKoutput$Timestep, lwd = 1.5, col = "darkorange2")
lines(PBKoutput$C_mu ~ PBKoutput$Timestep, lwd = 1.5, col = "hotpink")
lines(PBKoutput$C_s ~ PBKoutput$Timestep, lwd = 1.5, col = "purple")
lines(PBKoutput$C_k ~ PBKoutput$Timestep, lwd = 1.5, col = "navy")
lines(PBKoutput$C_l ~ PBKoutput$Timestep, lwd = 1.5, col = "green")
lines(PBKoutput$C_art ~ PBKoutput$Timestep, lwd = 1.5, col = "blue")

# Add the legend
legend("bottomright", inset(0), xpd = TRUE, 
       legend = c("Plasma", "Brain", "Muscle", "Skin", "Kidney", "Liver"),
       col = c("red", "blue", "darkorange2", "hotpink", "purple", "navy"), 
       lty = 1, lwd = 1.5)
dev.off()

### Output Data ###

file_name <- paste0(Name, "_PBKoutput.csv")
write.csv(PBKoutput, file_name, row.names = FALSE)

### ><(((º> ><(((º> ><(((º> ><(((º> ><(((º> ><(((º> ><(((º> ><(((º> ><(((º> ### 
