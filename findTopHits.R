#Examples for histograms
require("grid")

currentTheme = theme(axis.text = element_text(size = rel(6)), title = element_text(size = rel(6)), panel.grid.major = element_line(size = rel(4)), panel.grid.minor = element_line(size = rel(4)), panel.border = element_rect(size = rel(5), colour = "black"), legend.key.size = unit(units = "cm", 3), legend.key = element_rect(size = rel(4), colour = "black"), legend.text = element_text(size = rel(5)))

ggplot(allSeq, aes(x=DeltaGBinding, fill=factor(faRNALength))) + geom_histogram(binwidth=1, color="black", size = 2) + scale_y_sqrt() + xlab(expression(paste(Delta, "G"["Binding"], sep = ""))) + ylab("Count") + theme_bw() + scale_fill_discrete(name = "faRNA\nLength\n") + currentTheme

ggplot(allWob, aes(x=DeltaGBinding, fill=factor(faRNALength))) + geom_histogram(binwidth=1, color="black", size = 2) + scale_y_sqrt() + xlab(expression(paste(Delta, "G"["Binding"], sep = ""))) + ylab("Count") + theme_bw() + scale_fill_discrete(name = "faRNA\nLength\n") + currentTheme

ggplot(uRBSWob, aes(x=DeltaGBinding, fill=factor(faRNALength))) + geom_histogram(binwidth=1, color="black", size = 2) + scale_y_sqrt() + xlab(expression(paste(Delta, "G"["Binding"], sep = ""))) + ylab("Count") + theme_bw() + scale_fill_discrete(name = "faRNA\nLength\n") + currentTheme

ggplot(uTermWob, aes(x=DeltaGBinding, fill=factor(faRNALength))) + geom_histogram(binwidth=1, color="black", size = 2) + scale_y_sqrt()+ xlab(expression(paste(Delta, "G"["Binding"], sep = ""))) + ylab("Count") + theme_bw() + scale_fill_discrete(name = "faRNA\nLength\n") + currentTheme


allSeq = rbind(seq4_trmout, seq5_trmout, seq6_trmout, seq7_trmout, seq8_trmout, seq9_trmout, seq10_trmout)
allWob = rbind(wob4_trmout, wob5_trmout, wob6_trmout, wob7_trmout, wob8_trmout, wob9_trmout, wob10_trmout, wob11_trmout)
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





allSeq_RBSMatches = rbind(seq4_RBSMatches, seq5_RBSMatches)
allWob_RBSMatches = rbind(wob4_RBSMatches, wob5_RBSMatches, wob6_RBSMatches, wob7_RBSMatches, wob8_RBSMatches, wob9_RBSMatches, wob10_RBSMatches)
allSeq_RBSMatches$DeltaGBinding = -allSeq_RBSMatches$DeltaGBinding
allWob_RBSMatches$DeltaGBinding = -allWob_RBSMatches$DeltaGBinding

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
