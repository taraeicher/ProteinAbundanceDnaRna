if (!require(corrplot)) {
  install.packages("corrplot")
  library(corrplot)
}
breast_dna = read.table("C:\\Users\\tarae\\OneDrive\\Documents\\PhD\\BMI7830_Integrative_Methods_of_Bioinformatics_for_Human_Diseases\\dream_challenge_2\\breast_dna_corr.rds", fill = TRUE)
corrplot(as.numeric(as.character(breast_dna))