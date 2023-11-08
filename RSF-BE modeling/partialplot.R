extract_partial<-function(i,result,time,npts){
  gene_i.year<-subset(result,gene.name== gene_list[i]&Time_point==time)
  output <- as.data.frame(cbind(gene.name=rep(gene_list[i],npts), 
                                Gene.exp=sapply(1:npts, function(j) mean(gene_i.year$Gene.exp[seq(j, npts*B, by = npts)])),
                                survprob=sapply(1:npts, function(j)  mean(gene_i.year$survprob[seq(j,  npts*B, by = npts)])),
                                survprob_l=as.numeric(sapply(1:npts, function(j)  mean(gene_i.year$survprob[seq(j,  npts*B, by = npts)])))-2*as.numeric(sapply(1:npts, function(j)  sd(gene_i.year$survprob[seq(j,  npts*B, by = npts)])/5)),
                                survprob_u=as.numeric(sapply(1:npts, function(j)  mean(gene_i.year$survprob[seq(j,  npts*B, by = npts)])))+2*as.numeric(sapply(1:npts, function(j)  sd(gene_i.year$survprob[seq(j,  npts*B, by = npts)])/5)),
                                Time_point=rep(time,npts)))
  return(output)
  
}
partialplot.func<-function(rs,subtitle,surv_label,ncol){ggplot(rs, aes(x = Gene.exp, y = survprob)) + theme_bw()+
    theme_minimal() + 
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
    theme(panel.grid.major = element_blank(), 
          # panel.grid.minor = element_blank(),
          legend.position = "top",
          strip.background = element_blank(), strip.placement = "outside",
          text = element_text(size = 23),panel.spacing = unit(1, "lines"))+
    geom_point(aes(colour = factor(Time_point),fill=factor(Time_point),shape = factor(Time_point)), size=1.5)+
    geom_smooth(aes(x = Gene.exp, y = survprob_u,colour = factor(Time_point)),data=rs, linetype = 2, se = F, linewidth=0.6) +
    geom_smooth(aes(x = Gene.exp, y = survprob_l,colour = factor(Time_point)),data=rs, linetype = 2, se = F, linewidth=0.6)+
    geom_smooth(aes(x = Gene.exp, y = survprob,colour = factor(Time_point)),data=rs, linetype = 1, se = F, linewidth=0.7)+
    labs(title = "Partial dependence plots",
         subtitle = subtitle,
         y = surv_label,
         x = "") + 
    facet_wrap(~ gene.name,ncol =ncol, scales = "free_x", strip.position = "bottom")+
    scale_shape_manual(values=c(21,24))+  
    scale_fill_manual(values=c("#FF0033","#0066FF"))
}

