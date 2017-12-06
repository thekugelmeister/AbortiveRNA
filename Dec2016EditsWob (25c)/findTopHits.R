#Examples for histograms
require("grid")

# Import files
library(readr)

wob4_RBSMatches <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob4.RBSMatches.dG", 
                              " ", escape_double = FALSE, trim_ws = TRUE)
wob5_RBSMatches <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob5.RBSMatches.dG", 
                              " ", escape_double = FALSE, trim_ws = TRUE)
wob6_RBSMatches <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob6.RBSMatches.dG", 
                             " ", escape_double = FALSE, trim_ws = TRUE)
wob7_RBSMatches <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob7.RBSMatches.dG", 
                             " ", escape_double = FALSE, trim_ws = TRUE)
wob8_RBSMatches <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob8.RBSMatches.dG", 
                             " ", escape_double = FALSE, trim_ws = TRUE)
wob9_RBSMatches <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob9.RBSMatches.dG", 
                             " ", escape_double = FALSE, trim_ws = TRUE)
wob10_RBSMatches <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob10.RBSMatches.dG", 
                             " ", escape_double = FALSE, trim_ws = TRUE)
wob11_RBSMatches <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob11.RBSMatches.dG", 
                             " ", escape_double = FALSE, trim_ws = TRUE)
wob12_RBSMatches <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob12.RBSMatches.dG", 
                             " ", escape_double = FALSE, trim_ws = TRUE)
wob13_RBSMatches <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob13.RBSMatches.dG", 
                             " ", escape_double = FALSE, trim_ws = TRUE)
wob14_RBSMatches <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob14.RBSMatches.dG", 
                             " ", escape_double = FALSE, trim_ws = TRUE)
wob15_RBSMatches <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob15.RBSMatches.dG", 
                             " ", escape_double = FALSE, trim_ws = TRUE)

wob4_trmout <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob4.trmout.dG", 
                          " ", escape_double = FALSE, trim_ws = TRUE)
wob5_trmout <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob5.trmout.dG", 
                          " ", escape_double = FALSE, trim_ws = TRUE)
wob6_trmout <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob6.trmout.dG", 
                         " ", escape_double = FALSE, trim_ws = TRUE)
wob7_trmout <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob7.trmout.dG", 
                         " ", escape_double = FALSE, trim_ws = TRUE)
wob8_trmout <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob8.trmout.dG", 
                         " ", escape_double = FALSE, trim_ws = TRUE)
wob9_trmout <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob9.trmout.dG", 
                         " ", escape_double = FALSE, trim_ws = TRUE)
wob10_trmout <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob10.trmout.dG", 
                         " ", escape_double = FALSE, trim_ws = TRUE)
wob11_trmout <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob11.trmout.dG", 
                         " ", escape_double = FALSE, trim_ws = TRUE)
wob12_trmout <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob12.trmout.dG", 
                         " ", escape_double = FALSE, trim_ws = TRUE)
wob13_trmout <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob13.trmout.dG", 
                         " ", escape_double = FALSE, trim_ws = TRUE)
wob14_trmout <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob14.trmout.dG", 
                         " ", escape_double = FALSE, trim_ws = TRUE)
wob15_trmout <- read_delim("~/NairLab/Dec2016EditsWob (25c)/wob15.trmout.dG", 
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
