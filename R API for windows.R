#Debug mode: IfDeleteResult=FALSE; then check log file under Temp folder.

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
main_folder <- "Z:/Users/Daniel/Test/Immunodown";
setwd(main_folder);
#Detail_folder <- "Details"
#if (!file.exists(Detail_folder)){
#  dir.create(file.path(main_folder, Detail_folder))
#}
#final_output_file = "summary_report.txt"
#final_report = file(final_output_file)

#--------------------------------------------------------------------------------------------------
#                                         Basic setting                                            |
#--------------------------------------------------------------------------------------------------
Land.CurrentLand <- c(Server = "xxxxxx", UserID = "admin", Password = "xxxxxx", LandName = "TCGA2015");

print(ListLands <- Land.ListLands())

ListDataAvailability <- Land.ListDataAvailability()

MetaData_Immuno2015 <- Land.DownloadMetaData();
write.table(MetaData_Immuno2015, "MetaData_Immuno2015.txt", sep="\t", row.names=FALSE, quote=FALSE);

Land.InitiateLand("xxxxxx", "admin", "xxxxxx", "ImmunoLand2015");

Land.CheckInstall();

Land.CheckVersion();

Land.ConnectServer();

ListLands = Land.ListLands();

ListDataAvailability = Land.ListDataAvailability();

Land.CurrentLand <- c(Server = "xxxxxx", UserID = "admin", Password = "xxxxxx", LandName = "TCGA2015");

textDump <- Land.TextDumpArrayLandGeneData(Genes="(all)",SampleSet="GSE54456 Normal Non Lesional Samples", DataMode="Rna.Seq.FPKM", IfDeleteResult=TRUE,OutputFolder=main_folder)

genes <- c("BRAF", "KRAS", "TP53");
Expression_Ratio_TCGA2015A <- Land.TextDumpArrayLandData(Genes=genes, Samples="(all)", DataMode="Expression_Ratio");
Expression_Ratio_TCGA2015A <- Land.TextDumpArrayLandData(Genes=genes, Samples="(all)", DataMode="RnaSeq_Transcript");
Expression_Ratio_TCGA2015A <- Land.TextDumpArrayLandData(Genes=genes,  SampleSet="liver cancer", DataMode="RnaSeq_Transcript");

textDump = Land.TextDumpArrayLandGeneData(Genes=genes,GeneSet=NULL,Samples="(all)", DataMode="RnaSeq_Transcript", IfDeleteResult=TRUE,OutputFolder="/ngs-data1/datafolder1/")



textDump = Land.TextDumpArrayLandGeneData(Genes=genes,GeneSet=NULL,Samples=NULL,SampleSet="GSE54456 Normal Non Lesional Samples", DataMode="RnaSeq_Transcript", IfDeleteResult=TRUE,OutputFolder="/ngs-data1/datafolder1/")


textDump = Land.TextDumpArrayLandGeneData(Genes=genes,  SampleSet="Lesional in RNASeq", DataMode="RnaSeq_Transcript", IfDeleteResult=T,OutputFolder=main_folder)

> genes = c("MET","egfr","braf","KRas")
> textDump = Land.TextDumpArrayLandGeneData(Genes=genes,GeneSet=NULL,Samples=NULL,SampleSet="GSE54456 Normal Non Lesional Samples", DataMode="RnaSeq_Transcript", IfDeleteResult=TRUE,OutputFolder="/ngs-data1/datafolder1/")



FinalSampleSet=Land.TermSets_Terms(TermSets="GSE54456 Normal Non Lesional Samples", TermSetType="SampleSet", Terms=NULL,TermType="Samples")
