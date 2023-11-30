

SNP_names_file <- paste0("/raid6/Tianyu/PRS/bert_sample/ReferencePopulation-Package/", file.title,"/CHR/",file.title,"-chr", chr,".bim")
paste0("/raid6/Tianyu/PRS/bert_sample/ReferencePopulation-Package/", file.title,"/CHR/",file.title,"-chr", chr,".bed")

####save the selected preliminary beta
beta_generate_data_file <- './temp_file_generate_data/beta_generate_data.txt'
A1 <- unlist(lapply(strsplit(SNP, ":"),`[[`,4))
effect_size_df <- data.frame(SNP = SNP, 
                             A1 = A1,
                             BETA = betaGenerateData)
write.table(effect_size_df, beta_generate_data_file, sep = "\t", quote = FALSE, row.names = FALSE)

##calculate the PRS for each individual in the reference panel
if(anc == 'CEU'){
  file.title <- 'CEU-20K'
}else{
  file.title <- 'YRI-4K'  
}
reference_genotype_file <- paste0("/raid6/Tianyu/PRS/bert_sample/ReferencePopulation-Package/", file.title, "/", file.title)
PRS_file <- "./temp_file_generate_data/PRS_generate_data"
plink2.command = paste(plink2,"--nonfounders","--allow-no-sex","--threads",24,"--memory",25000,
                       "--bfile", reference_genotype_file ,
                       "--score", beta_generate_data_file,"header-read",1,2,
                       "--score-col-nums",3,
                       "--out",PRS_file,
                       sep=" ")

system(plink2.command)

PRS_result <- read.table(paste0(PRS_file, '.sscore'))
PRS_result$PRS <- PRS_result$V5 * PRS_result$V6