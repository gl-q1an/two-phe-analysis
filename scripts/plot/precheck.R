# DESCRIPTION:This script is used to draw a Manhattan plot of the imported summary
#      file and output it as an html structure to check.
# Note:If the results of the Manhattan plot are not significant,  the results of 
#      subsequent analyses may also be less significant,and it may be necessary to 
#      change the phenotype or perform a meta-analysis on the Summary data
# Input: The file formated before.

library(optparse)
library(data.table)
library(dplyr)
library(CMplot)
library(ieugwasr)
library(htmltools)
library(parallel)

help_d <- "The dirctory where the input files are placed."
help_i <- "The input summary, more than one can be entered at a time."
help_o <- "The output directory."
help_p <- "If this parameter is included, the program is processed in parallel, 
followed by the number of parallel cpus to be attached."
help_c <- "If this parameter is specified, it should be followed by a reference 
template to check for chain unbalance relationships"
help_l <- "The clump is based on plink, so if -c is specified, this param is necessary"
help_description <- "This script is used to draw a Manhattan plot of the imported summary roughly, 
which to to perform checks before data processing."

option_list <- list(
    make_option(c("-d", "--directory"), type = "character", default = ".",  
              help = help_d),
    make_option(c("-i", "--input"), type = "character", default = NULL, 
              list=TRUE, required = TRUE, help = help_i),
    make_option(c("-o", "--output"), type = "character", default = ".",
              help = help_o)
    make_option(c("-p", "--parallel"), type = "integer", default = NULL,
              help = help_p)
    make_option(c("-c", "--clump"), type = "character", default = NULL,
              help = help_c)
    make_option(c("-l", "--plink"), type = "character", default = NULL,
              help = help_l)          
)

opt <- parse_args(OptionParser(usage = "Usage: %prog [options]", 
                               option_list = option_list,
                               description=help_description))

directory <- opt$dirctory
input_list <- opt$input
output <- opt$output
parallel <- opt$parallel
clump <- opt$clump
plink <- opt$plink
ldsc_sign <- NULL

if (!is.null(clump) && is.null(plink)) {
    stop("Error: When --clump is set, --plink must also be set.\n")
}

plot_manhattan <- function(input, output, clump, plink) {
    
    data <- fread(input) %>%
      select(SNP,CHR,BP,P) %>%
      filter(P < 1e-3)
    
    n_5e8 <- data %>%
      filter(P < 5e-8) %>%
      nrow(.)

    n_5e6 <- data %>%
      filter(P < 5e-6) %>%
      nrow(.)
    
    if (!is.null(clump)){
        data_clump <- data %>%
          select(SNP,P) %>%
          rename(rsid = SNP, pval = P) %>%
          ld_clump(., clump_kb = 10000, clump_r2 = 0.001, clump_p = 5e-8, 
                  plink_bin=plink, bfile=clump)
        ldsc_sign <- data_clump$rsid
    }

    plot <- CMplot(data, plot.type="m", LOG10=TRUE, ylim=NULL,col=c("grey30","grey60"),
               threshld=c(5e-8,5e-6), threshold.lty=c(2,2),threshold.lwd=c(1,1), 
               threshold.col=c("black","grey"), signal.cex=c(1.5,1.5),signal.pch=c(19,19),
               highlight=ldsc_sign, highlight.col="pink", highlight.cex=1,highlight.pch=19,
               file="jpg",file.name="",dpi=300,file.output=TRUE,verbose=TRUE,
               width=14,height=6)
    
    list_result <- list(n_5e8=n_5e8, n_5e6=n_5e6, ldsc_sign=ldsc_sign, plot=plot)

    ！！！有问题 需要看一下CMplot会不会直接输出

    return(list_result)
}

write_hmtl <- function(file ,list_result, html_report){
    n_5e8 <- list_result$n_5e8
    n_5e6 <- list_result$n_5e6
    ldsc_sign <- list_result$ldsc_sign
    plot <- list_result$plot

    result_5e8 <- paste("The number of SNPs with P < 5e-8 is",as.character(n_5e8),".\n")
    result_5e6 <- paste("The number of SNPs with P < 5e-6 is",as.character(n_5e6),".\n")
    result_clump <- ""
    if (!is.null(ldsc_sign)) {
        result_clump <- paste("Clumped P < 5e-8 is",as.character(len(ldsc_sign)),".\n")
    }

    sub_section <- tags$div(
        tags$h1(paste("The summary of",file)),
        tags$p(paste(result_5e8, result_5e6, result_clump))
        img(src=plot, width = 60, height = 40)
        tags$p("Please Check!!\n\n")
    )

    html_report <- html_report %>% 
      html_add_section(sub_section)
    
    return(html_report)
}

html_report <- html_document(
    title="The overview of the GWAS summary"
)

if (!is.null(parallel)) {
    ！！！并行处理的代码记得写
} else {
    for (file in input){
        plot_manhattan(file, output, clump, plink)
    }
}

for (file in input){
    html_report <- write_hmtl(file, list_result, html_report)
}

current_time <- Sys.time()
formatted_time <- format(current_time, format = "%Y%m%d%H%M")
save_html(html_report,paste0("result_",formatted_time,".html"))