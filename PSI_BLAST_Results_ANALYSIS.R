#12.5.20
#PSI-BLAST Results
setwd('/Volumes/Protein_Passport/MSA_PAH')

#The following command was run to generate the input file for this script
# psiblast -db nr_db/nr -query MSA_PAH/Human_PAH_J3KND8_114-444.fasta \
# -out MSA_PAH/5iter_psiBLAST_HumanPAH_J3kND8_114_444.full_nr.txt \
# -outfmt 6 -max_target_seqs 100023 -evalue 1e-4 -num_iterations 5 -num_threads 2

#Load Libraries
library(dplyr)
library(tidyr)
library(tidyverse)
library(RCurl)
library(magrittr)

PSI_Results <- as.data.frame(read_delim('5iter_psiBLAST_HumanPAH_J3kND8_114_444.full_nr.txt','\t',col_names = F))
#Add column names to the PSI-BLAST Results
colnames(PSI_Results)[1] <- 'Query_ID'
colnames(PSI_Results)[2] <- 'Subject_ID'
colnames(PSI_Results)[3] <- 'Percent_Match'
colnames(PSI_Results)[4] <- 'Alignment_Length'
colnames(PSI_Results)[5] <- 'Mismatch'
colnames(PSI_Results)[6] <- 'Gap_Open'
colnames(PSI_Results)[7] <- 'Q_Start'
colnames(PSI_Results)[8] <- 'Q_End'
colnames(PSI_Results)[9] <- 'S_Start'
colnames(PSI_Results)[10] <- 'S_End'
colnames(PSI_Results)[11] <- 'evalue'
colnames(PSI_Results)[12] <- 'bit_score'
PSI_Results <- filter(PSI_Results,PSI_Results$Query_ID != 'Search has CONVERGED!')

#Find out how many unique IDs there are
length(unique(PSI_Results$Subject_ID))

#For each unique ID, trim the PSI-BLAST Results to find:
# 1.) highest % Match
# 2.) longest match length
# 3.) lowest S_Start
# 4.) highest S_end
# 5.) lowest evalue

#List of unique IDs
Unique_IDs <- unique(PSI_Results$Subject_ID)

#1.)
#######################################
# Data frame of highest percent match #
#######################################
Top_Percent_Match <- data.frame(matrix(ncol=1))
holder <- data.frame()
for(i in Unique_IDs){
  holder <- filter(PSI_Results,Subject_ID == i)
  Max_match <- max(holder$Percent_Match)
  Top_Percent_Match <- rbind2(Top_Percent_Match,Max_match)
  holder <- data.frame()
}
Top_Percent_Match <- Top_Percent_Match[-1,]
Top_Percent_Match <- as.data.frame(Top_Percent_Match)

#2.)
#######################################
# Data frame of longest match lengths #
#######################################
Longest_Unique_Read <- data.frame(matrix(ncol=1))
holder <- data.frame()
for(i in Unique_IDs){
  holder <- filter(PSI_Results,Subject_ID == i)
  Max_len <- max(holder$Alignment_Length)
  Longest_Unique_Read <- rbind2(Longest_Unique_Read,Max_len)
  holder <- data.frame()
}
Longest_Unique_Read <- Longest_Unique_Read[-1,]
Longest_Unique_Read <- as.data.frame(Longest_Unique_Read)

#3.)
################################
# Data frame of lowest S_Start #
################################
Unique_Low_Sstart <- data.frame(matrix(ncol=1))
holder <- data.frame()
for(i in Unique_IDs){
  holder <- filter(PSI_Results,Subject_ID == i)
  LSs <- min(holder$S_Start)
  Unique_Low_Sstart <- rbind2(Unique_Low_Sstart,LSs)
  holder <- data.frame()
}
Unique_Low_Sstart <- Unique_Low_Sstart[-1,]
Unique_Low_Sstart <- as.data.frame(Unique_Low_Sstart)

#4.)
###############################
# Data frame of highest S_End #
###############################
Unique_Max_Send <- data.frame(matrix(ncol=1))
holder <- data.frame()
for(i in Unique_IDs){
  holder <- filter(PSI_Results,Subject_ID == i)
  HSe <- max(holder$S_End)
  Unique_Max_Send <- rbind2(Unique_Max_Send,HSe)
  holder <- data.frame()
}
Unique_Max_Send <- Unique_Max_Send[-1,]
Unique_Max_Send <- as.data.frame(Unique_Max_Send)

#5.)
###############################
# Data frame of lowest evalue #
###############################
Unique_Low_evalue <- data.frame(matrix(ncol=1))
holder <- data.frame()
for(i in Unique_IDs){
  holder <- filter(PSI_Results,Subject_ID == i)
  ev <- min(holder$evalue)
  Unique_Low_evalue <- rbind2(Unique_Low_evalue,ev)
  holder <- data.frame()
}
Unique_Low_evalue <- Unique_Low_evalue[-1,]
Unique_Low_evalue <- as.data.frame(Unique_Low_evalue)

#Now I need to combine the 5 data frames I just made
Unique_IDs <- as.data.frame(Unique_IDs)
Unique_IDs <- cbind(Unique_IDs,Top_Percent_Match)
Unique_IDs <- cbind(Unique_IDs,Longest_Unique_Read)
Unique_IDs <- cbind(Unique_IDs,Unique_Low_Sstart)
Unique_IDs <- cbind(Unique_IDs,Unique_Max_Send)
Unique_IDs <- cbind(Unique_IDs,Unique_Low_evalue)
Unique_Seq_Stats <- Unique_IDs
Unique_Seq_Stats <- mutate(Unique_Seq_Stats, Sstart_Send_Domain = ((Unique_Max_Send- Unique_Low_Sstart) + 1))

# Make a histogram of Max Read Lengths to help determine where to cut
mean(Longest_Unique_Read$Longest_Unique_Read) #276
ggplot(Longest_Unique_Read, aes(x=Longest_Unique_Read)) + 
  geom_histogram(color="black", fill="darkgreen") +
  geom_vline(aes(xintercept=mean(Longest_Unique_Read),color="mean: 276"), #Enter mean value here
             linetype="solid", size=1)+
  ggtitle("Distribution of Max Read Length for PSI-BLAST Output")+
  xlab("Length of Unique Read Match")+
  ylab("Number of Sequences")

# Make a histogram of Maximum Extent of Alignment to help determine where to cut
mean(Unique_Seq_Stats$Sstart_Send_Domain) #281
ggplot(Unique_Seq_Stats, aes(x=Sstart_Send_Domain)) + 
  geom_histogram(color="black", fill="darkgreen") +
  geom_vline(aes(xintercept=mean(Sstart_Send_Domain),color="mean: 281"), #Enter mean value here
             linetype="solid", size=1)+
  ggtitle("Distribution of lengths for Maximum Extent of Alignment for PSI-BLAST Output")+
  xlab("Length of Sequence in Aligned")+
  ylab("Number of Sequences")

#Select a domain cut off based on the histogram 
CUT_OFF = 200 #Enter cut off value
Domain_Trimmed <- filter(Unique_Seq_Stats,Longest_Unique_Read >= CUT_OFF)
nrow(Domain_Trimmed) #Prints how many sequences made it past the cut
nrow(Unique_IDs) - nrow(Domain_Trimmed) #Prints how many sequences didn't make the cut

#Write gi numbers to a new file and fetch the sequences 
write.table(Domain_Trimmed$Unique_IDs,'IDs_to_fetch_seqs_for.txt',
            row.names = FALSE,col.names=FALSE,quote=FALSE)

##############
# NEXT STEPS #
##############

#Run a blastdbcmd using the file that was just generated 
#
# blastdbcmd -entry_batch IDs_to_fetch_seqs_for.txt -out PSI_FASTAS.fasta \
# -db ../nr_db/nr

#After running that command, it skipped 1306389B and 1306389C in the output
#Need to remove them from Domain_Trimmed
Updated_Domain_Trimmed <- subset(Domain_Trimmed,
                                 Unique_IDs!= '1306389B')
Updated_Domain_Trimmed <- subset(Updated_Domain_Trimmed,
                                 Unique_IDs!= '1306389C')

#Make blastdbcmd output fasta a single line fasta
# python /Users/seanquigley/Documents/PAH/toSingleLineFasta.py PSI_FASTAS.fasta --out singleLine_PSI_FASTAS.fa

#Make dataframes for sequences and headers only
# jupyter notebook Separate_Sequences_and_Headers.ipynb

#########################################################
# Read in and trim seqs from the above command's output #
#########################################################

#Read in headers and sequences from Separate_Sequences_and_Headers.ipynb output
PSI_Headers <- as.data.frame(read_lines('PSI_Headers_Only.txt'))
PSI_Sequences <- as.data.frame(read_lines('PSI_Sequences_Only.txt'))
colnames(PSI_Headers)[1] <- 'Header'
colnames(PSI_Sequences)[1] <- 'Sequence'

#Bind the two data frames 
DBCMD_Out <- cbind(PSI_Headers,PSI_Sequences)
#Bind the updated Domain trim to the DBCMD_Out
DBCMD_Annotated_PSI_Results <- cbind(Updated_Domain_Trimmed,DBCMD_Out)

#Find Length of Full Length Sequences
DBCMD_Annotated_PSI_Results <- mutate(DBCMD_Annotated_PSI_Results,
                                      Full_Seq_Len=str_length(Sequence))
#Trim Seqs to the Sstart and Send
DBCMD_Annotated_PSI_Results <- mutate(DBCMD_Annotated_PSI_Results,Sstart_Send_Trimmed=substr(DBCMD_Annotated_PSI_Results$Sequence,
                                                Unique_Low_Sstart,Unique_Max_Send))
#Find Sstar_Send Domain Lengths
DBCMD_Annotated_PSI_Results <- mutate(DBCMD_Annotated_PSI_Results,
                                      Sstart_Send_Length=str_length(Sstart_Send_Trimmed))
#Make a histogram of the SstartSend Lengths and cut where appropriate 
mean(DBCMD_Annotated_PSI_Results$Sstart_Send_Length) #290
ggplot(DBCMD_Annotated_PSI_Results, aes(x=Sstart_Send_Length)) + 
  geom_histogram(color="black", fill="darkgreen") +
  geom_vline(aes(xintercept=mean(Sstart_Send_Length),color="mean: 290"), #Enter mean value here
             linetype="solid", size=1)+
  ggtitle("Distribution of lengths for Maximum Extent of Alignment for PSI-BLAST Output")+
  xlab("Length of Sequence in Aligned")+
  ylab("Number of Sequences")

#Trim SstartSend at 450 and view new histogram
DBCMD_Annotated_PSI_Results <- subset(DBCMD_Annotated_PSI_Results,
                                      Sstart_Send_Length <= 450)
#Histogram
mean(DBCMD_Annotated_PSI_Results$Sstart_Send_Length) #290
ggplot(DBCMD_Annotated_PSI_Results, aes(x=Sstart_Send_Length)) + 
  geom_histogram(color="black", fill="darkgreen") +
  geom_vline(aes(xintercept=mean(Sstart_Send_Length),color="mean: 290"), #Enter mean value here
             linetype="solid", size=1)+
  ggtitle("Distribution of lengths for Maximum Extent of Alignment for PSI-BLAST Output\n450 Trim")+
  xlab("Length of Sequence in Aligned")+
  ylab("Number of Sequences")

###################################################################
# Write a file for the trimmed seqs and run either famsa or pasta #
###################################################################
#Isolate Headers and Sequences
#Headers
famsa_Headers <- select(DBCMD_Annotated_PSI_Results,Header)
colnames(famsa_Headers)[1] <- 'famsa'
#Sequences
famsa_Sequences <- select(DBCMD_Annotated_PSI_Results,Sstart_Send_Trimmed)
colnames(famsa_Sequences)[1] <- 'famsa'
#Interleave isolated data frames in a fasta formated style
fasta_for_famsa <- data.frame()
for(i in 1:nrow(famsa_Headers)){
  fasta_for_famsa <- rbind2(fasta_for_famsa,famsa_Headers[i,])
  fasta_for_famsa <- rbind2(fasta_for_famsa,famsa_Sequences[i,])
}
#Specify Column name of output file
colnames(fasta_for_famsa)[1] <- 'Seq'

#Trimmed Seqs
write.table(fasta_for_famsa$Seq,'Sequences_for_famsa_input.fa',
            row.names = FALSE,col.names=FALSE,quote=FALSE)

#############################
# Run famsa in command line #
#############################

# famsa -v Sequences_for_famsa_input.fa FAMSA_Alignment.aln

