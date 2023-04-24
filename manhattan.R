library(qqman)
library(ggplot2)


#reading file from plink
pulmo <- read.delim("C:/Users/brovo/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu18.04onWindows_79rhkp1fndgsc/LocalState/rootfs/home/olga/DKFZ/pulmo_new.assoc",sep = "" , header = T,  na.strings ="", stringsAsFactors= F)
hetero <- read.delim("C:/Users/brovo/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu18.04onWindows_79rhkp1fndgsc/LocalState/rootfs/home/olga/DKFZ/pulmo.het",sep = "" , header = T,  na.strings ="", stringsAsFactors= F)
pulmo_pca <- read.delim("C:/Users/brovo/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu18.04onWindows_79rhkp1fndgsc/LocalState/rootfs/home/olga/DKFZ/pulmo_pca.eigenvec",sep = "" , header = F,  na.strings ="", stringsAsFactors= F)


##managing data
pulmo$F_A <- as.numeric(pulmo$F_A)
pulmo <- na.omit(pulmo)
#creating vectors with snp of interest
snp_int <-c("rs111404182",
            "rs9282543",
            "rs111033363",
            "rs1800386",
            "rs112428886",
            "rs200631556",
            "rs35766612",
            "rs756241691",
            "rs199472880",
            "rs137854528",
            "rs879254678",
            "rs9282543",
            "rs144830740"
)

#performing manhattan plot
str(pulmo)
pulmo$CHR <- as.numeric(pulmo$CHR)
pulmo$P <- as.numeric(pulmo$P)
pulmo <- na.omit(pulmo)

snp_int <- snp_int[snp_int %in% pulmo$SNP]

#manhattan plot
png("C:/Bioinf/HSE/manh.png",width = 900, height = 500)
manhattan(pulmo, highlight = snp_int, annotatePval = 0.01)
dev.off()
##qqplot
png("C:/Bioinf/HSE/qq.png",width = 900, height = 500)
qq(pulmo$P)
dev.off()

# heterozygosity rate = (N(NM) - O(Hom))/N(NM)

hetero$rate <- (hetero$N.NM.-hetero$O.HOM.)/hetero$N.NM.
hetero$F_stat <- hetero$F
quantile_hetero <- quantile(hetero$rate)

png("C:/Bioinf/HSE/hetero_dist.png")
hist(hetero$rate, xlab = "Heterozygosity",
     main = "Distribution of heterozygosity rate per individual")
abline(v= c(0.08220352, 0.10076230), col="red")
dev.off()

## ploting PCA
rownames(pulmo_pca) <- pulmo_pca[,2]
pulmo_pca <- pulmo_pca[,3:ncol(pulmo_pca)]
colnames(pulmo_pca) <- paste('Principal Component ', c(1:20), sep = '')

population <- as.data.frame(rownames(pulmo_pca))
population$pop <- substr(population[,1],start=1,stop=2)
population$pop <- gsub("HG","European",population$pop)
population$pop <- gsub("NA","Utah",population$pop)
population$pop <- gsub("PH","Patients",population$pop)
population$pop <- as.factor(population$pop)

col <- colorRampPalette(c("yellow","forestgreen","royalblue"))(length(unique(population$pop)))[factor(population$pop)]


plot(pulmo_pca[,1], pulmo_pca[,2],
     type = 'n',
     adj = 0.5,
     xlab = 'PC1',
     ylab = 'PC2',
     font = 2,
     font.lab = 2)
points(pulmo_pca[,1], pulmo_pca[,2],pch = 20, cex = 2.25, col = col)
legend('topleft',
       title = 'Cohorts',
       c('European', 'Utah', 'Patients'),
       fill = c("yellow","forestgreen","royalblue"))

plot(pulmo_pca[,1], pulmo_pca[,3],
     type = 'n',
     adj = 0.5,
     xlab = 'PC1',
     ylab = 'PC2',
     font = 2,
     font.lab = 2)
points(pulmo_pca[,1], pulmo_pca[,3], pch = 20, cex = 2.25)

