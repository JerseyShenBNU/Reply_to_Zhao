analyze_correlation_climate_ato3_to_tws_across_eurasia_include_tibet <- function(
){
  
  font_size = 7
  text_theme = theme(
    panel.grid = element_blank(),
    plot.background = element_blank(),
    axis.text = element_text(size = font_size,color = 'black'),
    legend.text = element_text(size = font_size,color = 'black'),
    axis.title = element_text(size = font_size+2,color = 'black',face = 'bold'),
    legend.title = element_text(size = font_size+2,color = 'black',face = 'bold')
  )
  
  
  library(ggplot2)
  library(data.table)
  library(ggsci)

  prato = as.data.frame(fread('data/prato.csv'))
  evaato = as.data.frame(fread('data/evaato.csv'))
  pmeato = as.data.frame(fread('data/pmeato.csv'))
  tws_tibet = as.data.frame(fread('data/tws_tibet.csv'))
  tws_sr = as.data.frame(fread('data/tws_srs.csv'))
  
  stand_trend_fun<-function(x){
    x = ts(x,start = c(2003,1),frequency = 12)
    x = decompose(x)$trend
    x = x[-which(is.na(x))]
    x = (x-mean(x))/sd(x)
    return(x)
  }
  # Whole gap date is from 2017-07 to 2018-09
  date = seq(as.Date('2003-01-01'),as.Date('2023-10-01'),'1 month')
  # analyze 2006,2007,2008 to 201706
  loc_201706 = 174
  loc_200601 = which(date == '2006-01-01')
  loc_200701 = which(date == '2007-01-01')
  loc_200801 = which(date == '2008-01-01')
  loc_200901 = which(date == '2009-01-01')
  loc_201901 = which(date == '2019-01-01')
  loc_202310 = which(date == '2023-10-01')
  
  sdate = c(loc_200601,loc_200701,loc_200801,loc_200901,loc_201901)
  edate = c(loc_201706,loc_201706,loc_201706,loc_201706,loc_202310)
  
  i = 1:5
  analyze_fun_distribution_cor <- function(i){
    tmpsdate = sdate[i]
    tmpedate = edate[i]
    
    tmp_prato3 = prato[tmpsdate:tmpedate,3]
    tmp_evaato3 = evaato[tmpsdate:tmpedate,3]
    tmp_pmeato3 = pmeato[tmpsdate:tmpedate,3]
    
    tmp_twsdf = cbind(tws_sr,tws_tibet[,1:2])
    tmp_twsdf = tmp_twsdf[tmpsdate:tmpedate,]
    
    colnames(tmp_twsdf) = c(paste0('SR',1:14),'TPM1','TPM2')
    
    tmp_sprato3 = stand_trend_fun(tmp_prato3)
    tmp_sevaato3 = stand_trend_fun(tmp_evaato3)
    tmp_spmeato3 = stand_trend_fun(tmp_pmeato3)
    tmp_tws = apply(tmp_twsdf,2,stand_trend_fun)
    
    tmp_tws = tmp_tws[,-c(2,3)]
    
    tmp_twssr_mean = apply(tws_sr[tmpsdate:tmpedate,-c(2,3)],1,mean,na.rm = T)
    tmp_twssr_min = apply(tws_sr[tmpsdate:tmpedate,-c(2,3)],1,min,na.rm = T)
    tws_twssr_max = apply(tws_sr[tmpsdate:tmpedate,-c(2,3)],1,max,na.rm = T)
    
    tmp_twssr_mean = stand_trend_fun(tmp_twssr_mean)
    
    cor_box_pme = 1
    p_box_pme = 1
    for(k in 1:ncol(tmp_tws)){
      #tmpcor_pr = max(ccf(tmp_tws[,i],tmp_sprato3,plot = F)$acf)
      #tmpcor_eva = max(ccf(tmp_tws[,i],tmp_sevaato3,plot = F)$acf)
      tmpcor_pme = max(ccf(tmp_tws[,k],tmp_spmeato3,plot = F)$acf)
      tmpp_pme = cor.test(tmp_tws[,k],tmp_spmeato3)$p.value
      cor_box_pme = c(cor_box_pme,tmpcor_pme)
      p_box_pme = c(p_box_pme,tmpp_pme)
    }
    cor_box_pme = cor_box_pme[-1]
    p_box_pme = p_box_pme[-1]
    
    cor_box_spr = 1
    p_box_spr = 1
    for(k in 1:ncol(tmp_tws)){
      #tmpcor_pr = max(ccf(tmp_tws[,i],tmp_sprato3,plot = F)$acf)
      #tmpcor_eva = max(ccf(tmp_tws[,i],tmp_sevaato3,plot = F)$acf)
      tmpcor_spr = max(ccf(tmp_tws[,k],tmp_sprato3,plot = F)$acf)
      tmpp_spr = cor.test(tmp_tws[,k],tmp_sprato3)$p.value
      cor_box_spr = c(cor_box_spr,tmpcor_spr)
      p_box_spr = c(p_box_spr,tmpp_spr)
    }
    cor_box_spr = cor_box_spr[-1]
    p_box_spr = p_box_spr[-1]
    
    cor_box_seva = 1
    p_box_seva = 1
    for(k in 1:ncol(tmp_tws)){
      #tmpcor_eva = max(ccf(tmp_tws[,i],tmp_sevaato3,plot = F)$acf)
      #tmpcor_eva = max(ccf(tmp_tws[,i],tmp_sevaato3,plot = F)$acf)
      tmpcor_seva = max(ccf(tmp_tws[,k],tmp_sevaato3,plot = F)$acf)
      tmpp_seva = cor.test(tmp_tws[,k],tmp_sevaato3)$p.value
      cor_box_seva = c(cor_box_seva,tmpcor_seva)
      p_box_seva = c(p_box_seva,tmpp_seva)
    }
    cor_box_seva = cor_box_seva[-1]
    p_box_seva = p_box_seva[-1]
    
    cor_matrix = data.frame(region = c(paste0('SR',1:12),'TPM1','TPM2'),
                            cor_pr = cor_box_spr,
                            cor_eva = cor_box_seva,
                            cor_pme = cor_box_pme)
    
    p_matrix = data.frame(region = c(paste0('SR',1:12),'TPM1','TPM2'),
                          p_pr = p_box_spr,
                          p_eva = p_box_seva,
                          p_pme = p_box_pme)
    
    return(cor_matrix)
  }
  
  ret_cordis = lapply(i,analyze_fun_distribution_cor)
  # distribution result 
  plist = list()
  for(i in 1:length(ret_cordis)){
    tmpcormat = ret_cordis[[i]]
    tmpcormat$region[1:12] = paste0('H',tmpcormat$region[1:12] )
    
    point_df = data.frame(
      y = rep(c(1,2,3),each = nrow(tmpcormat)),
      value = c(tmpcormat[,2],tmpcormat[,3],tmpcormat[,4]),
      variable = rep(tmpcormat$region,3)
    )
    
    #cols = colorRampPalette(brewer.pal(9,'RdYlBu'))(length(unique(tmpcormat$region)))
    cols = pal_npg()(9)
    cols = colorRampPalette(cols)(length(unique(tmpcormat$region)))
    
    
    min_pr1 = min(tmpcormat[,2])
    max_pr1 = max(tmpcormat[,2])
    
    min_eva = min(tmpcormat[,3])
    max_eva = max(tmpcormat[,3])
    
    min_pme = min(tmpcormat[,4])
    max_pme = max(tmpcormat[,4])
    
    hlinedf = data.frame(
      y = rep(c(1,2,3),each = 1),
      x = c(min_pr1,min_eva,min_pme),
      xend = c(max_pr1,max_eva,max_pme)
    )
    
    vlinedf = data.frame(
      y = rep(c(1,2,3)-0.1,2),
      x = c(min_pr1,min_eva,min_pme,max_pr1,max_eva,max_pme),
      yend = rep(c(1,2,3)+0.1, 2)
    )
  
    
    mean_cor = apply(tmpcormat[,2:4],2,mean)
    
    mean_cordf = data.frame(
      x = mean_cor,
      y = c(1,2,3)+0.15,
      label = round(mean_cor,1)
    )
    
    plot_label = data.frame(
      x = 0.4,
      y = 3.4,
      label = paste0('(',letters[i+2],' ',substr(date[sdate[i]],1,7),'~',
                     substr(date[edate[i]],1,7),
                     ')')
    )
    
    ppoint = ggplot()+
      geom_segment(data = hlinedf,aes(x = x,y = y,xend = xend,yend = y),
                   color = 'black',linewidth = 0.5)+
      geom_segment(data = vlinedf,aes(x = x,y = y,xend = x,yend = yend),
                   color = 'black',linewidth = 0.5)+
      geom_point(data = point_df,
                 aes(x = value,y = y,color= variable),
                 size = 2,shape = 1)+
      geom_point(data = point_df,
                 aes(x = value,y = y,color= variable),
                 size = 2.2,shape = 1)+
      geom_point(data = point_df,
                 aes(x = value,y = y,color= variable),
                 size = 1.9,shape = 1)+
      geom_text(data=mean_cordf,
                aes(x = x,y = y,label = label),
                size= 2.5,color = 'black')+
      #geom_text(data=plot_label,
      #          aes(x = x,y = y,label = label),
      #          size= 2.5,color = 'black',hjust = 1)+
      scale_color_manual(values = cols)+
      theme_bw()+
      text_theme+
      theme(axis.text.y = element_text(angle = 90,hjust = 0.5))+
      scale_y_continuous(limits = c(0.5,3.5),
                         breaks = c(1,2,3),
                         labels = c('Pr_NATO3','Eva_NATO3','PME_NATO3'))+
      theme(legend.position = 'none')+
      ylab('Hydroclimate condition')+
      xlab("Cross correlation")
    
    
    if(i == 2){
      ppoint = ppoint+
        scale_x_continuous(breaks = c(0.2,0.5,0.8),labels = c(0.2,0.5,0.8))
    }else if(i == 3){
      ppoint = ppoint+
        scale_x_continuous(breaks = c(0.2,0.5,0.8),labels = c(0.2,0.5,0.8))
    }else if(i == 4){
      ppoint = ppoint+
        scale_x_continuous(breaks = c(0.2,0.5,0.8),labels = c(0.2,0.5,0.8))
    }else if(i == 5){
      ppoint = ppoint+
        scale_x_continuous(breaks = c(0.2,0.5,0.8),labels = c(0.2,0.5,0.8))
    }
    
    
    if(i>1){
      ppoint = ppoint+
        theme(axis.title.y = element_blank())
    }
    
    plist[[i]] = ppoint
    
  }
  
  
  
  
  date_id = 1:5
  analyze_fun_tws1_pmeato3_twssrs <- function(i){
    
    tmpsdate = sdate[i]
    tmpedate = edate[i]
    
    tmp_prato3 = prato[tmpsdate:tmpedate,3]
    tmp_evaato3 = evaato[tmpsdate:tmpedate,3]
    tmp_pmeato3 = pmeato[tmpsdate:tmpedate,3]
    
    tmp_twsdf = cbind(tws_sr,tws_tibet[,1:2])
    tmp_twsdf = tmp_twsdf[tmpsdate:tmpedate,]
    
    colnames(tmp_twsdf) = c(paste0('SR',1:14),'TPM1','TPM2')
    
    tmp_sprato3 = stand_trend_fun(tmp_prato3)
    tmp_sevaato3 = stand_trend_fun(tmp_evaato3)
    tmp_spmeato3 = stand_trend_fun(tmp_pmeato3)
    tmp_tws = apply(tmp_twsdf,2,stand_trend_fun)
    
    tmp_tws = tmp_tws[,-c(2,3)]
    
    tmp_twssr_mean = apply(tws_sr[tmpsdate:tmpedate,-c(2,3)],1,mean,na.rm = T)
    tmp_twssr_min = apply(tws_sr[tmpsdate:tmpedate,-c(2,3)],1,min,na.rm = T)
    tws_twssr_max = apply(tws_sr[tmpsdate:tmpedate,-c(2,3)],1,max,na.rm = T)
    
    tmp_twssr_mean = stand_trend_fun(tmp_twssr_mean)
    
    
    cor_sr_tpm1 = max(ccf(tmp_twssr_mean,tmp_tws[,13],plot =F)$acf)
    cor_pr_tpm1 = max(ccf(tmp_sprato3,tmp_tws[,13],plot = F)$acf)
    cor_eva_tpm1 = max(ccf(tmp_sevaato3,tmp_tws[,13],plot = F)$acf)
    cor_pme_tpm1 = max(ccf(tmp_spmeato3,tmp_tws[,13],plot = F)$acf)
    
    p_sr_tpm1 = cor.test(tmp_twssr_mean,tmp_tws[,13])$p.value
    p_pr_tpm1 = cor.test(tmp_sprato3,tmp_tws[,13])$p.value
    p_eva_tpm1 = cor.test(tmp_sevaato3,tmp_tws[,13])$p.value
    p_pme_tpm1 = cor.test(tmp_spmeato3,tmp_tws[,13])$p.value
    
    
    
    
    cor_sr_tpm2 = max(ccf(tmp_twssr_mean,tmp_tws[,14],plot =F)$acf)
    cor_pr_tpm2 = max(ccf(tmp_sprato3,tmp_tws[,14],plot = F)$acf)
    cor_eva_tpm2 = max(ccf(tmp_sevaato3,tmp_tws[,14],plot = F)$acf)
    cor_pme_tpm2 = max(ccf(tmp_spmeato3,tmp_tws[,14],plot = F)$acf)
    
    p_sr_tpm2 = cor.test(tmp_twssr_mean,tmp_tws[,14])$p.value
    p_pr_tpm2 = cor.test(tmp_sprato3,tmp_tws[,14])$p.value
    p_eva_tpm2 = cor.test(tmp_sevaato3,tmp_tws[,14])$p.value
    p_pme_tpm2 = cor.test(tmp_spmeato3,tmp_tws[,14])$p.value
    
    
    dfcor_tpm1 = data.frame(
      date = date_id[i],
      CCor_TWS_SRs = cor_sr_tpm1,
      CCor_Pr_NATO3 = cor_pr_tpm1,
      CCor_Eva_NATO3 = cor_eva_tpm1,
      CCor_PME_NATO3 = cor_pme_tpm1
    )
    
    
    dfcor_tpm2 = data.frame(
      date = date_id[i],
      CCor_TWS_SRs = cor_sr_tpm2,
      CCor_Pr_NATO3 = cor_pr_tpm2,
      CCor_Eva_NATO3 = cor_eva_tpm2,
      CCor_PME_NATO3 = cor_pme_tpm2
    )
    
    dfp_tpm1 = data.frame(
      date = date_id[i],
      P_TWS_SRs = p_sr_tpm1,
      P_Pr_NATO3 = p_pr_tpm1,
      P_Eva_NATO3 = p_eva_tpm1,
      P_PME_NATO3 = p_pme_tpm1
    )
    dfp_tpm2 = data.frame(
      date = date_id[i],
      P_TWS_SRs = p_sr_tpm2,
      P_Pr_NATO3 = p_pr_tpm2,
      P_Eva_NATO3 = p_eva_tpm2,
      P_PME_NATO3 = p_pme_tpm2
    )
    
    
    
    if(i !=5){
      return(list(dfcor_tpm1,dfcor_tpm2,dfp_tpm1,dfp_tpm2))
    }else{
      dfcor_tpm1_negeva = max(ccf(tmp_sevaato3*-1,tmp_tws[,13],plot = F)$acf)
      dfcor_tpm2_negeva = max(ccf(tmp_sevaato3*-1,tmp_tws[,14],plot = F)$acf)
      p_tpm1_negeva = cor.test(tmp_sevaato3*-1,tmp_tws[,13])$p.value
      p_tpm2_negeva = cor.test(tmp_sevaato3*-1,tmp_tws[,14])$p.value
      
      dfcor_tpm1$CCor_NegEva_NATO3 = dfcor_tpm1_negeva
      dfcor_tpm2$CCor_NegEva_NATO3 = dfcor_tpm2_negeva
      dfp_tpm1$P_NegEva_NATO3 = p_tpm1_negeva
      dfp_tpm2$P_NegEva_NATO3 = p_tpm2_negeva
      
      return(list(dfcor_tpm1,dfcor_tpm2,dfp_tpm1,dfp_tpm2))
    }
    
  }
  
  
  ret_list2 = lapply(date_id,analyze_fun_tws1_pmeato3_twssrs)
  
  cor_tpm1box = 1
  cor_tpm2box = 1
  p_tpm1box = 1
  p_tpm2box = 1
  for(i in 1:5){
    print(i)
    tmp = ret_list2[[i]]
    
    tmp1 = tmp[[1]]
    tmp2 = tmp[[2]]
    tmp3 = tmp[[3]]
    tmp4 = tmp[[4]]
    
    tmp1 = reshape2::melt(tmp1,'date')
    tmp2 = reshape2::melt(tmp2,'date')
    tmp3 = reshape2::melt(tmp3,'date')
    tmp4 = reshape2::melt(tmp4,'date')
    
    cor_tpm1box = rbind(cor_tpm1box,tmp1)
    cor_tpm2box = rbind(cor_tpm2box,tmp2)
    p_tpm1box = rbind(p_tpm1box,tmp3)
    p_tpm2box = rbind(p_tpm2box,tmp4)
  }
  cor_tpm1box = cor_tpm1box[-1,]
  cor_tpm2box = cor_tpm2box[-1,]
  p_tpm1box = p_tpm1box[-1,]
  p_tpm2box = p_tpm2box[-1,]
  
  
  cor_tpm1box$value[c(1:16,18:21)] = round(cor_tpm1box$value[c(1:16,18:21)],1)
  cor_tpm1box$value[c(17)] = round(cor_tpm1box$value[c(17)],2)
  cor_tpm2box$value[c(1:16,18:21)] = round(cor_tpm2box$value[c(1:16,18:21)],1)
  cor_tpm2box$value[c(17)] = round(cor_tpm2box$value[c(17)],2)
  
  df1 = cor_tpm1box 
  df2 = cor_tpm2box
  
  tell_significance <- function(x){
    if(x <0.01){
      label = "**"
    }else if(x <0.05){
      label = "*"
    }else if(x >=0.1){
      label = ' '
    }
    return(label)
  }
  
  plabel1 = sapply(p_tpm1box$value,tell_significance)
  df1$plabel = paste0(df1$value,'\n',plabel1)
  
  plabel2 = sapply(p_tpm2box$value,tell_significance)
  df2$plabel = paste0(df2$value,'\n',plabel2)
  
  
 
  #fils = pal_lancet(alpha = 0.7)(5)
  fils = pal_npg()(5)
  #fils = colorRampPalette(brewer.pal(9,'RdYlBu'))(5)
  
  date_label = paste0(substr(date[sdate],1,7),'~',substr(date[edate],1,7))
  
  p1 = ggplot()+
    geom_bar(data = df1,aes(x = date,y = value,fill = variable),
             width = 0.8,
             position = position_dodge2(0.8),
             stat = 'identity')+
    geom_text(data = df1,aes(x = date,y = value *1,label = plabel),
              size =2.5,color = 'black',position = position_dodge2(0.8),hjust = 0.5)+
    theme_bw()+
    text_theme+
    #geom_hline(yintercept = 0.4)+
    theme(legend.position = 'none')+
    scale_x_continuous(breaks = c(1,2,3,4,5),
                       labels = date_label)+
    scale_y_continuous(breaks = c(0,0.5,1),
                       labels = c(0,0.5,1),
                       limits = c(0,1.1))+
    scale_fill_manual(values = fils)+
    theme(axis.text.y = element_text(angle= 90,hjust= 0.5))+
    xlab('Time period')+
    ylab('Cross correlation')
  
  p2 = ggplot()+
    geom_bar(data = df2,aes(x = date,y = value,fill = variable),
             width = 0.8,
             position = position_dodge2(0.8),
             stat = 'identity')+
    geom_text(data = df2,aes(x = date,y = value *1,label = plabel),
              size =2.5,color = 'black',position = position_dodge2(0.8),hjust = 0.5)+
    theme_bw()+
    text_theme+
    #geom_hline(yintercept = 0.4)+
    theme(legend.position = 'none')+
    scale_x_continuous(breaks = c(1,2,3,4,5),
                       labels = date_label)+
    scale_y_continuous(breaks = c(0,0.5,1),
                       labels = c(0,0.5,1),
                       limits = c(0,1.1))+
    scale_fill_manual(values = fils)+
    theme(axis.text.y = element_text(angle= 90,hjust= 0.5))+
    xlab('Time period')+
    ylab('Cross correlation')
  
  
  # distribution plot 
  
  library(cowplot)
  
  pcom = plot_grid(p1,p2,ncol = 1,
                   rel_widths = c(1,1),
                   rel_heights = c(1,1))
  
  p3 = plot_grid(plist[[1]],plist[[2]],plist[[3]],
                 plist[[4]],plist[[5]],nrow = 1,
                 rel_widths = c(1,1,1,1,1),
                 rel_heights = c(1,1,1,1,1))
  
  pcom2 = plot_grid(pcom,p3,ncol = 1,
                    rel_heights = c(1,1),
                    rel_widths = c(1,1))
  
  library(ggpubr)
  plabel1 = p1+theme(legend.position = 'bottom')+
    guides(fill = guide_legend(title = '(a-b) Cross correlation between TWS in TPM1-2 and external factors',
                              nrow = 1,title.position = 'top'))
  
  plabel1 = as_ggplot(get_legend(plabel1))
  
  plabel2 = plist[[1]]+theme(legend.position = 'bottom')+
    guides(color = guide_legend(title = '(c-g) Cross correlation between TWS across HSRs and TPMs and hydroclimate conditions over NATO3',
                               nrow = 2,title.position = 'top'))
  plabel2 = as_ggplot(get_legend(plabel2))
  
  
  plabels = plot_grid(plabel1,plabel2,
                      ncol = 1,
                      rel_heights = c(1,1),
                      rel_widths =  c(1,1))
  
  output = 'fig'
  dir.create(output)
  output1 = paste0(output,'/fig_correaltion.pdf')
  output2 = paste0(output,'/fig_correaltion_label.pdf')
  ggsave(output1,plot = pcom2,
         width = 18 ,height = 18,units = 'cm',dpi = 800)
  
  
  ggsave(output2,plot = plabels,
         width = 18 ,height = 18,units = 'cm',dpi = 800)
  
  
}