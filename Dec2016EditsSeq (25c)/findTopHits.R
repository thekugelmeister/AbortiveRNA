#Examples for histograms
require("grid")

# Import files
library(readr)

seq4_RBSMatches <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq4.RBSMatches.dG", 
                              " ", escape_double = FALSE, trim_ws = TRUE)
seq5_RBSMatches <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq5.RBSMatches.dG", 
                              " ", escape_double = FALSE, trim_ws = TRUE)
seq6_RBSMatches <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq6.RBSMatches.dG", 
                             " ", escape_double = FALSE, trim_ws = TRUE)
seq7_RBSMatches <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq7.RBSMatches.dG", 
                             " ", escape_double = FALSE, trim_ws = TRUE)
seq8_RBSMatches <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq8.RBSMatches.dG", 
                             " ", escape_double = FALSE, trim_ws = TRUE)
seq9_RBSMatches <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq9.RBSMatches.dG", 
                             " ", escape_double = FALSE, trim_ws = TRUE)
seq10_RBSMatches <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq10.RBSMatches.dG", 
                             " ", escape_double = FALSE, trim_ws = TRUE)
seq11_RBSMatches <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq11.RBSMatches.dG", 
                             " ", escape_double = FALSE, trim_ws = TRUE)
seq12_RBSMatches <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq12.RBSMatches.dG", 
                             " ", escape_double = FALSE, trim_ws = TRUE)
seq13_RBSMatches <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq13.RBSMatches.dG", 
                             " ", escape_double = FALSE, trim_ws = TRUE)
seq14_RBSMatches <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq14.RBSMatches.dG", 
                             " ", escape_double = FALSE, trim_ws = TRUE)
seq15_RBSMatches <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq15.RBSMatches.dG", 
                             " ", escape_double = FALSE, trim_ws = TRUE)

seq4_trmout <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq4.trmout.dG", 
                          " ", escape_double = FALSE, trim_ws = TRUE)
seq5_trmout <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq5.trmout.dG", 
                          " ", escape_double = FALSE, trim_ws = TRUE)
seq6_trmout <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq6.trmout.dG", 
                         " ", escape_double = FALSE, trim_ws = TRUE)
seq7_trmout <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq7.trmout.dG", 
                         " ", escape_double = FALSE, trim_ws = TRUE)
seq8_trmout <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq8.trmout.dG", 
                         " ", escape_double = FALSE, trim_ws = TRUE)
seq9_trmout <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq9.trmout.dG", 
                         " ", escape_double = FALSE, trim_ws = TRUE)
seq10_trmout <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq10.trmout.dG", 
                         " ", escape_double = FALSE, trim_ws = TRUE)
seq11_trmout <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq11.trmout.dG", 
                         " ", escape_double = FALSE, trim_ws = TRUE)
seq_12trmout <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq12.trmout.dG", 
                         " ", escape_double = FALSE, trim_ws = TRUE)
seq_13trmout <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq13.trmout.dG", 
                         " ", escape_double = FALSE, trim_ws = TRUE)
seq_14trmout <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq14.trmout.dG", 
                         " ", escape_double = FALSE, trim_ws = TRUE)
seq_15trmout <- read_delim("~/NairLab/Dec2016EditsSeq (25c)/seq15.trmout.dG", 
                         " ", escape_double = FALSE, trim_ws = TRUE)


currentTheme = theme(axis.text = element_text(size = rel(6)), title = element_text(size = rel(6)), panel.grid.major = element_line(size = rel(4)), panel.grid.minor = element_line(size = rel(4)), panel.border = element_rect(size = rel(5), colour = "black"), legend.key.size = unit(units = "cm", 3), legend.key = element_rect(size = rel(4), colour = "black"), legend.text = element_text(size = rel(5)))

ggplot(allSeq, aes(x=DeltaGBinding, fill=factor(faRNALength))) + geom_histogram(binwidth=1, color="black", size = 2) + scale_y_sqrt() + xlab(expression(paste(Delta, "G"["Binding"], sep = ""))) + ylab("Count") + theme_bw() + scale_fill_discrete(name = "faRNA\nLength\n") + currentTheme

ggplot(allWob, aes(x=DeltaGBinding, fill=factor(faRNALength))) + geom_histogram(binwidth=1, color="black", size = 2) + scale_y_sqrt() + xlab(expression(paste(Delta, "G"["Binding"], sep = ""))) + ylab("Count") + theme_bw() + scale_fill_discrete(name = "faRNA\nLength\n") + currentTheme

ggplot(uRBSWob, aes(x=DeltaGBinding, fill=factor(faRNALength))) + geom_histogram(binwidth=1, color="black", size = 2) + scale_y_sqrt() + xlab(expression(paste(Delta, "G"["Binding"], sep = ""))) + ylab("Count") + theme_bw() + scale_fill_discrete(name = "faRNA\nLength\n") + currentTheme

ggplot(uTermWob, aes(x=DeltaGBinding, fill=factor(faRNALength))) + geom_histogram(binwidth=1, color="black", size = 2) + scale_y_sqrt()+ xlab(expression(paste(Delta, "G"["Binding"], sep = ""))) + ylab("Count") + theme_bw() + scale_fill_discrete(name = "faRNA\nLength\n") + currentTheme


allSeq = rbind(seq4_trmout, seq5_trmout, seq6_trmout, seq7_trmout, seq8_trmout, seq9_trmout, seq10_trmout, seq11_trmout, seq_12trmout, seq_13trmout, seq_14trmout, seq_15trmout)
allWob = rbind(wob4_trmout, wob5_trmout, wob6_trmout, wob7_trmout, wob8_trmout, wob9_trmout, wob10_trmout, wob11_trmout, wob12_trmout, wob13_trmout, wob14_trmout, wob15_trmout)
# allSeq$DeltaGBinding = -allSeq$DeltaGBinding
# allWob$DeltaGBinding = -allWob$DeltaGBinding


#SEQUENCE TERMINATOR ANALYSIS
#separate forward and reverse strand matches
forwardSeq = allSeq[allSeq$Strand == "forward",]
reverseSeq = allSeq[allSeq$Strand == "reverse",]
#order forward hits on TSS, matchRB, and BindingDG in that order
forwardSeq = forwardSeq[order(forwardSeq$TSS, forwardSeq$matchRB, forwardSeq$DeltaGBinding),]
#order reverse hits on TSS, matchLB, and BindingDG in that order
reverseSeq = reverseSeq[order(reverseSeq$TSS, reverseSeq$matchLB, reverseSeq$DeltaGBinding),]
#select uniuque hits from both lists based on the TSS and the "static" boundary of the match site across all abortive lengths
uTermSeq = rbind(forwardSeq[!duplicated(forwardSeq[, c("TSS", "matchRB"),]),], reverseSeq[!duplicated(reverseSeq[, c("TSS", "matchLB"),]),])
uTermSeq = uTermSeq[order(uTermSeq$DeltaGBinding),]

#WOBBLE TERMINATOR ANALYSIS
#separate forward and reverse strand matches
forwardWob = allWob[allWob$Strand == "forward",]
reverseWob = allWob[allWob$Strand == "reverse",]
#order forward hits on TSS, matchRB, and BindingDG in that order
forwardWob = forwardWob[order(forwardWob$TSS, forwardWob$matchRB, forwardWob$DeltaGBinding),]
#order reverse hits on TSS, matchLB, and BindingDG in that order
reverseWob = reverseWob[order(reverseWob$TSS, reverseWob$matchLB, reverseWob$DeltaGBinding),]
#select uniuque hits from both lists based on the TSS and the "static" boundary of the match site across all abortive lengths
uTermWob = rbind(forwardWob[!duplicated(forwardWob[, c("TSS", "matchRB"),]),], reverseWob[!duplicated(reverseWob[, c("TSS", "matchLB"),]),])
uTermWob = uTermWob[order(uTermWob$DeltaGBinding),]





allSeq_RBSMatches = rbind(seq4_RBSMatches, seq5_RBSMatches, seq6_RBSMatches, seq7_RBSMatches, seq8_RBSMatches, seq9_RBSMatches, seq10_RBSMatches, seq11_RBSMatches, seq12_RBSMatches, seq13_RBSMatches, seq14_RBSMatches, seq15_RBSMatches)
allWob_RBSMatches = rbind(wob4_RBSMatches, wob5_RBSMatches, wob6_RBSMatches, wob7_RBSMatches, wob8_RBSMatches, wob9_RBSMatches, wob10_RBSMatches, wob11_RBSMatches, wob12_RBSMatches, wob13_RBSMatches, wob14_RBSMatches, wob15_RBSMatches)
# allSeq_RBSMatches$DeltaGBinding = -allSeq_RBSMatches$DeltaGBinding
# allWob_RBSMatches$DeltaGBinding = -allWob_RBSMatches$DeltaGBinding

#SEQUENCE RBS ANALYSIS
#separate forward and reverse strand matches
forwardSeq = allSeq_RBSMatches[allSeq_RBSMatches$Strand == "forward",]
reverseSeq = allSeq_RBSMatches[allSeq_RBSMatches$Strand == "reverse",]
#order forward hits on TSS, matchRB, and BindingDG in that order
forwardSeq = forwardSeq[order(forwardSeq$TSS, forwardSeq$matchRB, forwardSeq$DeltaGBinding),]
#order reverse hits on TSS, matchLB, and BindingDG in that order
reverseSeq = reverseSeq[order(reverseSeq$TSS, reverseSeq$matchLB, reverseSeq$DeltaGBinding),]
#select uniuque hits from both lists based on the TSS and the "static" boundary of the match site across all abortive lengths
uRBSSeq = rbind(forwardSeq[!duplicated(forwardSeq[, c("TSS", "matchRB"),]),], reverseSeq[!duplicated(reverseSeq[, c("TSS", "matchLB"),]),])
uRBSSeq = uRBSSeq[order(uRBSSeq$DeltaGBinding),]

#WOBBLE RBS ANALYSIS
#separate forward and reverse strand matches
forwardWob = allWob_RBSMatches[allWob_RBSMatches$Strand == "forward",]
reverseWob = allWob_RBSMatches[allWob_RBSMatches$Strand == "reverse",]
#order forward hits on TSS, matchRB, and BindingDG in that order
forwardWob = forwardWob[order(forwardWob$TSS, forwardWob$matchRB, forwardWob$DeltaGBinding),]
#order reverse hits on TSS, matchLB, and BindingDG in that order
reverseWob = reverseWob[order(reverseWob$TSS, reverseWob$matchLB, reverseWob$DeltaGBinding),]
#select unique hits from both lists based on the TSS and the "static" boundary of the match site across all abortive lengths
uRBSWob = rbind(forwardWob[!duplicated(forwardWob[, c("TSS", "matchRB"),]),], reverseWob[!duplicated(reverseWob[, c("TSS", "matchLB"),]),])
uRBSWob = uRBSWob[order(uRBSWob$DeltaGBinding),]
