library(ggplot2)
library(ggpubr)
library(ggrepel)

cut_off_logFC = 1
cut_off_pvalue = 0.05

all_compare$change = ifelse(all_compare$T2vsT0_padj < cut_off_pvalue & abs(all_compare$T2vsT0_log2FoldChange) >= cut_off_logFC, 
                            ifelse(all_compare$T2vsT0_log2FoldChange> cut_off_logFC ,'Up','Down'),'Stable')
all_compare_na <- na.omit(all_compare)
all_compare_na$label <- ifelse(all_compare_na$gene_name == 'Klf4', as.character(all_compare_na$gene_name),"")
all_compare_na$label <- ifelse(all_compare_na$gene_name == 'Klf2', as.character(all_compare_na$gene_name),all_compare_na$label)
p1 <- ggplot(all_compare_na, aes(x = T2vsT0_log2FoldChange, 
                                 y = -log10(T2vsT0_padj), 
                                 colour=change)) +
  geom_point(alpha=0.4, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  geom_vline(xintercept=c(-1,1),
             lty=4,
             col="black",
             lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),
             lty=4,col="black",
             lwd=0.8) +
  theme_bw()+
  coord_cartesian(ylim = c(0,5))+
  labs(x="log2(fold change)",)+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
p2 <- ggplot(all_compare_na, aes(x = T2vsT0_log2FoldChange, 
                                 y = -log10(T2vsT0_padj), 
                                 colour=change)) +
  geom_point(alpha=0.4, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  geom_vline(xintercept=c(-1,1),
             lty=4,
             col="black",
             lwd=0.8) +
  #geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  geom_label_repel( aes(x = T2vsT0_log2FoldChange, 
                        y = -log10(T2vsT0_padj), 
                        label = label),
                    size = 3)+
  labs(,y="-log10 (p.adj)")+
  theme_bw()+coord_cartesian(ylim = c(5,35))


ggarrange(p2,p1,heights=c(4/6, 2/6),ncol = 1, nrow = 2,common.legend = TRUE,legend="right",align = "v") 
