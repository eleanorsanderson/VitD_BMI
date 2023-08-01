#Code for MR estimation of Vitamin D on BMI


library (TwoSampleMR)
library(writexl)

#Select the SNPs associated with 25 hydroxyvitamin D (id code: ieu-b-4812)
vitD_exp_dat <- extract_instruments(outcomes='ieu-b-4812')
#'clump' the data to remove correlated SNPs
vitD_exp_dat <- clump_data(vitD_exp_dat)
#Extract data on the SNP-BMI association for the SNPs associated with 25 hydroxyvitamin D
#BMI id ieu-b-40
bmi_out_dat <- extract_outcome_data (snps = vitD_exp_dat$SNP, outcomes = 'ieu-b-40')


#Harmonise the data for the exposure and outcome so the same allele is the reference allele
dat <- harmonise_data (exposure_dat = vitD_exp_dat, outcome_dat = bmi_out_dat)
dat$Fstat <- (dat$beta.exposure^2)/(dat$se.exposure^2)
Fstat<- mean(dat$Fstat)

#Run the main and robust MR analyses
res_bmi <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_weighted_mode", 
                                   "mr_egger_regression"))
res_bmi
#Run the per SNP MR analysis
res_single <- mr_singlesnp(dat)
res_single
write_xlsx(res_single,"res_single.xlsx")

# Visualization of data in scatter plot
dat$outcome <- "Body Mass Index"
dat$exposure <- "Vitamin D"
p1<- mr_scatter_plot(res_bmi, dat)
# Visualization of data in a forest plot
forest_plot <- mr_forest_plot(res_single)
forest_plot 



# To apply Steiger filtering
dat$samplesize.exposure <- 443734

dat$r_exp <- get_r_from_bsen(
  b = dat$beta.exposure, 
  se = dat$se.exposure, 
  n = dat$samplesize.exposure)

dat$samplesize.outcome <- 700000

dat$r_out <- get_r_from_bsen(
  b = dat$beta.outcome, 
  se = dat$se.outcome, 
  n = dat$samplesize.outcome)

# Steiger filtering
steiger_dat <- steiger_filtering(dat)
res_bmi <- mr(steiger_dat)


#BMI to Vitamin D
library(TwoSampleMR)

# Select the SNPs associated with BMI
BMI_exp_dat<-extract_instruments(outcomes= 'ieu-b-40')

#Clump the data to remove correlated SNPs
BMI_exp_dat<-clump_data(BMI_exp_dat)

#Extract data on SNP-Vitamin D associations for the SNPs associated with BMI
#Vitamin D ieu-b-4812
VitaminD_out_dat<-extract_outcome_data(snps = BMI_exp_dat$SNP,outcomes = 'ieu-b-4812')

#Harmonise the data for the exposure & outcome so that same data is reference allele
dat<-harmonise_data(exposure_dat = BMI_exp_dat,outcome_dat = VitaminD_out_dat)

#run the main & robust MR analysis
res_VitaminD<-mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_weighted_mode", 
                                      "mr_egger_regression"))
res_VitaminD

#run the per SNP MR analysis
res_single<-mr_singlesnp(dat)
res_single


