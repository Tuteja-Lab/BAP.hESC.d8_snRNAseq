v.clus1 <- ggplot(data=clus1.five.twenty, aes(x=log2fc, y=-log10(p_val), col=diffexpressed, label=delabel)) +
    geom_point(alpha = 0.5) +
    theme_classic() +
    scale_color_manual(name = "Expression", values=c("#4d4d4d", "#0571b0", "#ca0020"), labels = c(expression(paste("Higher in 20% ", "O"[2])), "All other genes", expression(paste("Higher in 5% ", "O"[2])) )) +
    geom_text_repel(show.legend = F) +
    ggtitle(expression(paste("Cluster 1: 20% ", "O"[2], " vs. 5% ", "O"[2]))) +
    xlab(expression(paste("log"[2], " fold change"))) +
    ylab(expression(paste("-log"[10], " pvalue"))) +
    theme(legend.text.align = 0)
ggsave("Figure_6_A.svg", plot = v.clus1, dpi=900, width = 8, height = 6)
v.clus2 <- ggplot(data=clus2.five.twenty, aes(x=log2fc, y=-log10(p_val), col=diffexpressed, label=delabel)) +
    geom_point(alpha = 0.5) +
    theme_classic() +
    scale_color_manual(name = "Expression", values=c("#4d4d4d", "#0571b0", "#ca0020"), labels = c(expression(paste("Higher in 20% ", "O"[2])), "All other genes", expression(paste("Higher in 5% ", "O"[2])) )) +
    geom_text_repel(show.legend = F) +
    ggtitle(expression(paste("Cluster 2: 20% ", "O"[2], " vs. 5% ", "O"[2]))) +
    xlab(expression(paste("log"[2], " fold change"))) +
    ylab(expression(paste("-log"[10], " pvalue"))) +
    theme(legend.text.align = 0)
ggsave("Figure_6_B.svg", plot = v.clus2, dpi=900, width = 8, height = 6)
v.clus3 <- ggplot(data=clus3.five.twenty, aes(x=log2fc, y=-log10(p_val), col=diffexpressed, label=delabel)) +
    geom_point(alpha = 0.5) +
    theme_classic() +
    scale_color_manual(name = "Expression", values=c("#4d4d4d", "#0571b0", "#ca0020"), labels = c(expression(paste("Higher in 20% ", "O"[2])), "All other genes", expression(paste("Higher in 5% ", "O"[2])) )) +
    geom_text_repel(show.legend = F) +
    ggtitle(expression(paste("Cluster 3: 20% ", "O"[2], " vs. 5% ", "O"[2]))) +
    xlab(expression(paste("log"[2], " fold change"))) +
    ylab(expression(paste("-log"[10], " pvalue"))) +
    theme(legend.text.align = 0)
ggsave("Figure_6_C.svg", plot = v.clus3, dpi=900, width = 8, height = 6)
v.clus4 <- ggplot(data=clus4.five.twenty, aes(x=log2fc, y=-log10(p_val), col=diffexpressed, label=delabel)) +
    geom_point(alpha = 0.5) +
    theme_classic() +
    scale_color_manual(name = "Expression", values=c("#4d4d4d", "#0571b0", "#ca0020"), labels = c(expression(paste("Higher in 20% ", "O"[2])), "All other genes", expression(paste("Higher in 5% ", "O"[2])) )) +
    geom_text_repel(show.legend = F) +
    ggtitle(expression(paste("Cluster 4: 20% ", "O"[2], " vs. 5% ", "O"[2]))) +
    xlab(expression(paste("log"[2], " fold change"))) +
    ylab(expression(paste("-log"[10], " pvalue"))) +
    theme(legend.text.align = 0)
ggsave("Figure_6_D.svg", plot = v.clus4, dpi=900, width = 8, height = 6)
v.clus5 <- ggplot(data=clus5.five.twenty, aes(x=log2fc, y=-log10(p_val), col=diffexpressed, label=delabel)) +
    geom_point(alpha = 0.5) +
    theme_classic() +
    scale_color_manual(name = "Expression", values=c("#4d4d4d", "#0571b0", "#ca0020"), labels = c(expression(paste("Higher in 20% ", "O"[2])), "All other genes", expression(paste("Higher in 5% ", "O"[2])) )) +
    geom_text_repel(show.legend = F) +
    ggtitle(expression(paste("Cluster 5: 20% ", "O"[2], " vs. 5% ", "O"[2]))) +
    xlab(expression(paste("log"[2], " fold change"))) +
    ylab(expression(paste("-log"[10], " pvalue"))) +
    theme(legend.text.align = 0)
ggsave("Figure_6_E.svg", plot = v.clus5, dpi=900, width = 8, height = 6)
v.clus6 <- ggplot(data=clus6.five.twenty, aes(x=log2fc, y=-log10(p_val), col=diffexpressed, label=delabel)) +
    geom_point(alpha = 0.5) +
    theme_classic() +
    scale_color_manual(name = "Expression", values=c("#4d4d4d", "#0571b0", "#ca0020"), labels = c(expression(paste("Higher in 20% ", "O"[2])), "All other genes", expression(paste("Higher in 5% ", "O"[2])) )) +
    geom_text_repel(show.legend = F) +
    ggtitle(expression(paste("Cluster 6: 20% ", "O"[2], " vs. 5% ", "O"[2]))) +
    xlab(expression(paste("log"[2], " fold change"))) +
    ylab(expression(paste("-log"[10], " pvalue"))) +
    theme(legend.text.align = 0)
ggsave("Figure_6_F.svg", plot = v.clus6, dpi=900, width = 8, height = 6)
v.clus7 <- ggplot(data=clus7.five.twenty, aes(x=log2fc, y=-log10(p_val), col=diffexpressed, label=delabel)) +
    geom_point(alpha = 0.5) +
    theme_classic() +
    scale_color_manual(name = "Expression", values=c("#4d4d4d", "#0571b0", "#ca0020"), labels = c(expression(paste("Higher in 20% ", "O"[2])), "All other genes", expression(paste("Higher in 5% ", "O"[2])) )) +
    geom_text_repel(show.legend = F) +
    ggtitle(expression(paste("Cluster 7: 20% ", "O"[2], " vs. 5% ", "O"[2]))) +
    xlab(expression(paste("log"[2], " fold change"))) +
    ylab(expression(paste("-log"[10], " pvalue"))) +
    theme(legend.text.align = 0)
ggsave("Figure_6_G.svg", plot = v.clus7, dpi=900, width = 8, height = 6)
v.clus8 <- ggplot(data=clus8.five.twenty, aes(x=log2fc, y=-log10(p_val), col=diffexpressed, label=delabel)) +
    geom_point(alpha = 0.5) +
    theme_classic() +
    scale_color_manual(name = "Expression", values=c("#4d4d4d", "#0571b0", "#ca0020"), labels = c(expression(paste("Higher in 20% ", "O"[2])), "All other genes", expression(paste("Higher in 5% ", "O"[2])) )) +
    geom_text_repel(show.legend = F) +
    ggtitle(expression(paste("Cluster 8: 20% ", "O"[2], " vs. 5% ", "O"[2]))) +
    xlab(expression(paste("log"[2], " fold change"))) +
    ylab(expression(paste("-log"[10], " pvalue"))) +
    theme(legend.text.align = 0)
ggsave("Figure_6_H.svg", plot = v.clus8, dpi=900, width = 8, height = 6)
v.clus9 <- ggplot(data=clus9.five.twenty, aes(x=log2fc, y=-log10(p_val), col=diffexpressed, label=delabel)) +
    geom_point(alpha = 0.5) +
    theme_classic() +
    scale_color_manual(name = "Expression", values=c("#4d4d4d", "#0571b0", "#ca0020"), labels = c(expression(paste("Higher in 20% ", "O"[2])), "All other genes", expression(paste("Higher in 5% ", "O"[2])) )) +
    geom_text_repel(show.legend = F) +
    ggtitle(expression(paste("Cluster 9: 20% ", "O"[2], " vs. 5% ", "O"[2]))) +
    xlab(expression(paste("log"[2], " fold change"))) +
    ylab(expression(paste("-log"[10], " pvalue"))) +
    theme(legend.text.align = 0)
ggsave("Figure_6_I.svg", plot = v.clus9, dpi=900, width = 8, height = 6)

