library(VennDiagram)
#--------------------------------------------------------------------------------------------------
#                                       load land functions                                           |
#--------------------------------------------------------------------------------------------------
load(url("http://omicsoft.com/downloads/land/Rapi/Land_R_API.Rda"));
OshellDirectory = 'D:/Oshell'
BaseDir = 'C:/Users/binbin/Documents/Omicsoft'
TempDir = 'C:/Users/binbin/Documents/Omicsoft/Temp'
Land.InitiateOshell(OshellDirectory = OshellDirectory,BaseDirectory = BaseDir, TempDirectory = TempDir);

#--------------------------------------------------------------------------------------------------
#                                       Prepare setting                                           |
#--------------------------------------------------------------------------------------------------
main_folder <- "Z:/Users/Daniel/Test/20150702 Oncoland2015Q2 bench Mark"
setwd(main_folder)
Detail_folder <- "Details"
if (!file.exists(Detail_folder)){
  dir.create(file.path(main_folder, Detail_folder))
}
final_output_file = "summary_report.txt"
final_report = file(final_output_file)
#--------------------------------------------------------------------------------------------------
#                                         Basic setting                                            |
#--------------------------------------------------------------------------------------------------
old_server = "192.168.1.106:8065"
new_server = "192.168.1.106:7065"
land_name = "TCGA2015"
Primary_grouping = "Tissue"
lots_samples = F
Gene_list = "Z:/Users/Daniel/Test/20150702 Oncoland2015Q2 bench Mark/Gene_list/Oncogenes.txt"

#------------------------------------------------

# for the old land
old_land <- c(Server = old_server, UserID = "admin", Password = "xxxxxx", LandName = land_name)
# for the new land
new_land <- c(Server = new_server, UserID = "admin", Password = "xxxxxx", LandName = land_name)

#connect to the server
Land.InitiateLand(Server = old_server, UserID = "admin",Password = "xxxxxx", LandName = land_name)
Land.ConnectServer() 

#--------------------------------------------------------------------------------------------------
#                                   Check Data Availability                                         |
#--------------------------------------------------------------------------------------------------
Land.CurrentLand <- old_land
old_land_DataAvailability <- Land.ListDataAvailability()

Land.CurrentLand <- new_land
new_land_DataAvailability <- Land.ListDataAvailability()

#test if old_land_DataAvailability and new_land_DataAvailability are identical 
DataAvailability_check <- identical(old_land_DataAvailability, new_land_DataAvailability)
if (DataAvailability_check) {
  writeLines(c("The two lands have the same Data Availbility \n"), final_report)
}else{
  writeLines(c("The two lands have different Data Availbility, Please check the DataAvailability.txt under the Details subfolder for more information \n"),final_report)
}

write.table(old_land_DataAvailability,file="DataAvailability_old.txt", row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")
write.table(new_land_DataAvailability,file="DataAvailability_new.txt", row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")

#--------------------------------------------------------------------------------------------------
#                                   Check Meta Data                                               |
#--------------------------------------------------------------------------------------------------

Land.CurrentLand <- old_land
old_landMetaData <- Land.DownloadMetaData()

Land.CurrentLand <- new_land
new_landMetaData <- Land.DownloadMetaData()

Metadata_check <- identical(old_landMetaData, new_landMetaData)
#FALSE

output_line = paste('old land has ',nrow(old_landMetaData),' samples, new land has ',nrow(new_landMetaData),' samples\n',sep='')
writeLines(output_line, final_report)
cat(output_line)

output_line = paste('old land has ',ncol(old_landMetaData),' metadata, new land has ',ncol(new_landMetaData),' metadata\n',sep='')
writeLines(output_line, final_report)
cat(output_line)

#if not identical, you might use "sort"
sample_set_check <- identical(sort(old_landMetaData$ID),sort(new_landMetaData$ID))
if (sample_set_check) {
  writeLines("samples are the same between two lands", final_report)
}else{
  writeLines("samples are different between two lands", final_report)
}

#--------------------------------------------------------------------------------------------------
#                            Generate the Gene and samples                                        |
#--------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------
# build a function to randomly pick up a sample from each tumor type
random_select_sample<-function(metaData,Primary_grouping)
{ 
  total_primary_factors <- metaData[,Primary_grouping]
  index <- duplicated(total_primary_factors)
  unique_primary_factor <- total_primary_factors[!index]
  unique_primary_factor <- sort(unique_primary_factor) 
  primary.position <- which(colnames(metaData)==Primary_grouping)
  results=NULL
  for(a in unique_primary_factor)
  { 
    ids=metaData[primary.position,1]
    results[a]=ids[round(runif(1,min=1,max=length(ids)))]	
  }   
  return(results)  
}

if (lots_samples==T)
{
  test_samples = old_landMetaData[,1]
  genes = c('BRAF','TP53','KRAS','ERG','TMPRSS2')
}else{
  genes_temp = read.table(Gene_list,header=T)
  genes = as.vector(genes_temp$GeneSymbol)
  test_samples=random_select_sample(old_landMetaData,Primary_grouping=Primary_grouping)
}

#--------------------------------------------------------------------------------------------------
#                 retrieve data from land based on different data modes                           |
#--------------------------------------------------------------------------------------------------
#-------------------------------------
#retrieve Different data from old land
Land.CurrentLand = old_land
samples = test_samples
Expression_Ratio_old        <- Land.TextDumpArrayLandData(Genes = c("MET","egfr","braf","KRas"), Samples = "(all)", DataMode = "Expression_Ratio")
General_Expression_old      <- Land.TextDumpArrayLandData(Genes = genes, Samples = samples, DataMode = "General_Expression")
CNV_old                     <- Land.TextDumpArrayLandData(Genes = genes, Samples = samples, DataMode = "CNV")
RPPA_old                    <- Land.TextDumpArrayLandData(Genes = genes, Samples = samples, DataMode = "RPPA")
DnaSeq_Mutation_old         <- Land.TextDumpArrayLandData(Genes = genes, Samples = samples, DataMode = "DnaSeq_Mutation")
DnaSeq_SomaticMutation_old  <- Land.TextDumpArrayLandData(Genes = genes, Samples = samples, DataMode = "DnaSeq_SomaticMutation")
RnaSeq_Transcript_old       <- Land.TextDumpArrayLandData(Genes = genes, Samples = samples, DataMode = "RnaSeq_Transcript")
RnaSeq_GeneBas_old          <- Land.TextDumpArrayLandData(Genes = genes, Samples = samples, DataMode = "RnaSeq_GeneBas")
RnaSeq_Fusion_old           <- Land.TextDumpArrayLandData(Genes = genes, Samples = samples, DataMode = "RnaSeq_Fusion")
RnaSeq_Mutation_old         <- Land.TextDumpArrayLandData(Genes = genes, Samples = samples, DataMode = "RnaSeq_Mutation")
RnaSeq_SomaticMutation_old  <- Land.TextDumpArrayLandData(Genes = genes, Samples = samples, DataMode = "RnaSeq_SomaticMutation")
Methylation450_B37_old      <- Land.TextDumpArrayLandData(Genes = genes, Samples = samples, DataMode = "Methylation450_B37")

#---- OK, let's play with the new land data -------
Land.CurrentLand = new_land
samples = test_samples
Expression_Ratio_new        <- Land.TextDumpArrayLandData(Genes = c("MET","egfr","braf","KRas"), Samples = "(all)", DataMode = "Expression_Ratio")
General_Expression_new      <- Land.TextDumpArrayLandData(Genes = genes, Samples = samples, DataMode = "General_Expression")
CNV_new                     <- Land.TextDumpArrayLandData(Genes = genes, Samples = samples, DataMode = "CNV")
RPPA_new                    <- Land.TextDumpArrayLandData(Genes = genes, Samples = samples, DataMode = "RPPA")
DnaSeq_Mutation_new         <- Land.TextDumpArrayLandData(Genes = genes, Samples = samples, DataMode = "DnaSeq_Mutation")
DnaSeq_SomaticMutation_new  <- Land.TextDumpArrayLandData(Genes = genes, Samples = samples, DataMode = "DnaSeq_SomaticMutation")
RnaSeq_Transcript_new       <- Land.TextDumpArrayLandData(Genes = genes, Samples = samples, DataMode = "RnaSeq_Transcript")
RnaSeq_GeneBas_new          <- Land.TextDumpArrayLandData(Genes = genes, Samples = samples, DataMode = "RnaSeq_GeneBas")
RnaSeq_Fusion_new           <- Land.TextDumpArrayLandData(Genes = genes, Samples = samples, DataMode = "RnaSeq_Fusion")
RnaSeq_Mutation_new         <- Land.TextDumpArrayLandData(Genes = genes, Samples = samples, DataMode = "RnaSeq_Mutation")
RnaSeq_SomaticMutation_new  <- Land.TextDumpArrayLandData(Genes = genes, Samples = samples, DataMode = "RnaSeq_SomaticMutation")
Methylation450_B37_new      <- Land.TextDumpArrayLandData(Genes = genes, Samples = samples, DataMode = "Methylation450_B37")

#---- Now compare new with old -------
if (!is.null(Expression_Ratio_new)){write.table(Expression_Ratio_new,file="Expression_ratio_new.txt",row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}
if (!is.null(Expression_Ratio_old)){write.table(Expression_Ratio_old,file="Expression_ratio_old.txt",row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}

if (!is.null(Expression_Ratio_new)){write.table(Expression_Ratio_new,file="Expression_ratio_new.txt",row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}
if (!is.null(Expression_Ratio_old)){write.table(Expression_Ratio_old,file="Expression_ratio_old.txt",row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}

if (!is.null(General_Expression_new)){write.table(General_Expression_new,file="General_Expression_new.txt",row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}
if (!is.null(General_Expression_old)){write.table(General_Expression_old,file="General_Expression_old.txt",row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}

if (!is.null(CNV_new)){ write.table(CNV_new,file="CNV_new.txt",row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}
if (!is.null(CNV_old)){ write.table(CNV_old,file="CNV_old.txt",row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}


if (!is.null(RPPA_new)){ write.table(RPPA_new,file="RPPA_new.txt",row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}
if (!is.null(RPPA_old)){ write.table(RPPA_old,file="RPPA_old.txt",row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}


if (!is.null(DnaSeq_Mutation_new)){ write.table(DnaSeq_Mutation_new$DnaSeq_Mutation,file="DnaSeq_Mutation_new.txt",row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}
if (!is.null(DnaSeq_Mutation_old)){ write.table(DnaSeq_Mutation_old$DnaSeq_Mutation,file="DnaSeq_Mutation_old.txt",row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}


if (!is.null(DnaSeq_SomaticMutation_new)){ write.table(DnaSeq_SomaticMutation_new$DnaSeq_SomaticMutation,file="DnaSeq_SomaticMutation.txt_new",row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}
if (!is.null(DnaSeq_SomaticMutation_old)){ write.table(DnaSeq_SomaticMutation_old$DnaSeq_SomaticMutation,file="DnaSeq_SomaticMutation.txt_old",row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}


if (!is.null(RnaSeq_Transcript_new)){ write.table(RnaSeq_Transcript_new,file="RnaSeq_Transcript_new.txt",row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}
if (!is.null(RnaSeq_Transcript_old)){ write.table(RnaSeq_Transcript_old,file="RnaSeq_Transcript_old.txt",row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}

AllGene_RnaSeq_Transcript <- merge(RnaSeq_Transcript_old, RnaSeq_Transcript_new, by=c("SampleID","GeneID","Name"), all.x=TRUE, all.y=TRUE);
UnequalTest <- (AllGene_RnaSeq_Transcript$FPKM.x != AllGene_RnaSeq_Transcript$FPKM.y);
UnequalTest[is.na(UnequalTest)] <- TRUE;
AllGene_RnaSeq_Transcript_Dif <- AllGene_RnaSeq_Transcript[UnequalTest,]
AllGene_RnaSeq_Transcript_Dif$Difference = AllGene_RnaSeq_Transcript_Dif$FPKM.x - AllGene_RnaSeq_Transcript_Dif$FPKM.y
write.csv(AllGene_RnaSeq_Transcript_Dif, paste(version,"Transcript_Dif",new_land_name,".csv"),  row.names=FALSE, quote=FALSE);

if (!is.null(RnaSeq_GeneBas_new)){ write.table(RnaSeq_GeneBas_new,file=paste(new_land_name,version,"RnaSeq_GeneBas.txt",sep=''),row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}
if (!is.null(RnaSeq_GeneBas_old)){ write.table(RnaSeq_GeneBas_old,file=paste(old_land_name,version,"RnaSeq_GeneBas.txt",sep=''),row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}
AllGene_RnaSeq_GeneBas <- merge(RnaSeq_GeneBas_old, RnaSeq_GeneBas_new, by=c("SampleID","GeneID","Chromosome","Start","End"), all.x=TRUE, all.y=TRUE);
UnequalTest <- (AllGene_RnaSeq_GeneBas$Coverage.x != AllGene_RnaSeq_GeneBas$Coverage.y);
UnequalTest[is.na(UnequalTest)] <- TRUE;
AllGene_RnaSeq_GeneBas_Dif <- AllGene_RnaSeq_GeneBas[UnequalTest,]
AllGene_RnaSeq_GeneBas_Dif$Difference = AllGene_RnaSeq_GeneBas_Dif$Coverage.x - AllGene_RnaSeq_GeneBas_Dif$Coverage.y
write.csv(AllGene_RnaSeq_GeneBas_Dif, paste(version,"GeneBas_Dif",new_land_name,".csv"),  row.names=FALSE, quote=FALSE);


if (!is.null(RnaSeq_Fusion_new)){ write.table(RnaSeq_Fusion_new$RnaSeq_Fusion,file=paste(new_land_name,version,"RnaSeq_Fusion.txt",sep=''),row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}
if (!is.null(RnaSeq_Fusion_old)){ write.table(RnaSeq_Fusion_old$RnaSeq_Fusion,file=paste(old_land_name,version,"RnaSeq_Fusion.txt",sep=''),row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}
AllGene_RnaSeq_Fusion <- merge(RnaSeq_Fusion_old$RnaSeq_Fusion, RnaSeq_Fusion_new$RnaSeq_Fusion, by=c("SampleID","GeneID","FusionID","FusionStrand","Gene1","Gene2"), all.x=TRUE, all.y=TRUE);
UnequalTest <- (AllGene_RnaSeq_Fusion$FPKM.x != AllGene_RnaSeq_Fusion$FPKM.y & (AllGene_RnaSeq_Fusion$SeedCount.x>3 | AllGene_RnaSeq_Fusion$SeedCount.y>3));
UnequalTest[is.na(UnequalTest)] <- TRUE;
AllGene_RnaSeq_Fusion_Dif <- AllGene_RnaSeq_Fusion[UnequalTest,]
AllGene_RnaSeq_Fusion_Dif$Difference = AllGene_RnaSeq_Fusion_Dif$FPKM.x - AllGene_RnaSeq_Fusion_Dif$FPKM.y
write.csv(AllGene_RnaSeq_Fusion_Dif, paste(version,"Fusion_Dif",new_land_name,".csv"),  row.names=FALSE, quote=FALSE);


# 
# 
# identical(RnaSeq_Mutation_new,RnaSeq_Mutation_old)
# if (!is.null(RnaSeq_Mutation_new)){ write.table(RnaSeq_Mutation_new$RnaSeq_Mutation,file=paste(new_land_name,version,"RnaSeq_Mutation.txt",sep=''),row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}
# if (!is.null(RnaSeq_Mutation_old)){ write.table(RnaSeq_Mutation_old$RnaSeq_Mutation,file=paste(old_land_name,version,"RnaSeq_Mutation.txt",sep=''),row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}
# 
# AllGene_RnaSeq_Mutation <- merge(RnaSeq_Mutation_old$RnaSeq_Mutation, RnaSeq_Mutation_new$RnaSeq_Mutation, by=c("SampleID","GeneID","Chromosome","Start","End"), all.x=TRUE, all.y=TRUE);
# UnequalTest <- (AllGene_RnaSeq_Mutation$MutationID.x != AllGene_RnaSeq_Mutation$MutationID.y & (AllGene_RnaSeq_Mutation$Frequency.x >0.3 | AllGene_RnaSeq_Mutation$Frequency.y >0.3));
# UnequalTest[is.na(UnequalTest)] <- TRUE;
# AllGene_RnaSeq_Mutation_Dif <- AllGene_RnaSeq_Mutation[UnequalTest,]
# write.csv(AllGene_RnaSeq_Mutation_Dif, paste(version,"Mutation_Dif",new_land_name,".csv"),  row.names=FALSE, quote=FALSE);
# 
# #---- play with the venn diagam
# 
# old.temp = RnaSeq_Mutation_old$RnaSeq_Mutation
# cols = c('SampleID','MutationID')
# old.temp$Sample_Mutation <- apply( old.temp[ , cols ] , 1 , paste , collapse = "-" )
# old.temp = old.temp[(old.temp$Frequency>0.3),]
# old.all = dim(old.temp)[1]
# #write.csv(old.temp,"test.csv",  row.names=FALSE, quote=FALSE);
# 
# new.temp = RnaSeq_Mutation_new$RnaSeq_Mutation
# cols = c('SampleID','MutationID')
# new.temp$Sample_Mutation <- apply( new.temp[ , cols ] , 1 , paste , collapse = "-" )
# new.temp = new.temp[(old.temp$Frequency>0.3),]
# new.all = dim(new.temp)[1]
# overlap = length(which(old.temp$Sample_Mutation %in% new.temp$Sample_Mutation))
# grid.newpage()
# draw.pairwise.venn(area1 = old.all, area2 = new.all, cross.area = overlap, category = c(old_land_name, new_land_name),fill = c("light blue", "pink"), alpha = rep(0.5, 2))
# png(paste(version,"Mutation_Dif",new_land_name,".png"))
# grid.newpage()
# draw.pairwise.venn(area1 = old.all, area2 = new.all, cross.area = overlap, category = c(old_land_name, new_land_name),fill = c("light blue", "pink"), alpha = rep(0.5, 2))
# dev.off()
# 
# if (!is.null(RnaSeq_SomaticMutation_new)){ write.table(RnaSeq_SomaticMutation_new$RnaSeq_SomaticMutation,file=paste(new_land_name,version,"RnaSeq_SomaticMutation.txt",sep=''),row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}
# if (!is.null(RnaSeq_SomaticMutation_old)){ write.table(RnaSeq_SomaticMutation_old$RnaSeq_SomaticMutation,file=paste(old_land_name,version,"RnaSeq_SomaticMutation.txt",sep=''),row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}
# 
# 
# if (!is.null(Methylation450_B37_new)){ write.table(Methylation450_B37_new,file=paste(new_land_name,version,"Methylation450_B37.txt",sep=''),row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}
# if (!is.null(Methylation450_B37_old)){ write.table(Methylation450_B37_old,file=paste(old_land_name,version,"Methylation450_B37.txt",sep=''),row.names=FALSE,col.names=TRUE,quote = FALSE, sep = "\t")}
# 
# identical(Expression_Ratio_new, Expression_Ratio_old)
# identical(General_Expression_new, General_Expression_old)
# identical(CNV_new, CNV_old)
# identical(RPPA_new,RPPA_old)
# identical(DnaSeq_Mutation_new,DnaSeq_Mutation_old)
# identical(DnaSeq_SomaticMutation_new,DnaSeq_SomaticMutation_old)
# identical(RnaSeq_Transcript_new,RnaSeq_Transcript_old)
# identical(RnaSeq_GeneBas_new,RnaSeq_GeneBas_old)
# identical(RnaSeq_Fusion_new,RnaSeq_Fusion_old)
# identical(RnaSeq_SomaticMutation_new,RnaSeq_SomaticMutation_old)
# identical(Methylation450_B37_new,Methylation450_B37_old)
close(final_output_file)