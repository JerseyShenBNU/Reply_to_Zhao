analyza_correaltion_vw_snowcover <- function(
){
  
  input_w1 = 'output_rebuttal/uploads/data/vertical_velocity_tpm1.csv'
  input_w2 = 'output_rebuttal/uploads/data/vertical_velocity_tpm2.csv'

  stpm1_w_box = as.data.frame(fread(input_w1))
  stpm2_w_box = as.data.frame(fread(input_w2))
  
  tpmcor1 = 1
  tpmcor2 = 1
  pbox1 = 1
  pbox2 = 1
  for(i in 1:3){
    print(i)
    stpm1vw = stpm1_w_box[,i]
    stpm2vw = stpm2_w_box[,i]
    
  
    snowcover = as.data.frame(fread('output_rebuttal/uploads/data/snowcover.csv'))
    snowcover = snowcover[,1:2]
    
    
    tmpcor1 = min(ccf(stpm1vw,snowcover[,1],plot = F)$acf)
    tmpcor2 = min(ccf(stpm2vw,snowcover[,2],plot = F)$acf)
    
    
    print(hpa[i])
    print(tmpcor1)
    print(tmpcor2)
    matplot(cbind(stpm1vw,snowcover[,1]),type = 'l')
    matplot(cbind(stpm2vw,snowcover[,2]),type = 'l')
    plot(stpm1vw,snowcover[,1])
    plot(stpm2vw,snowcover[,2])
    
    tpmcor1 = c(tpmcor1,tmpcor1)
    tmpp1 = cor.test(stpm1vw,snowcover[,1])$p.value
    tmpp1 = scientific(tmpp1,2)
    pbox1 = c(pbox1,tmpp1)
    
    tpmcor2 = c(tpmcor2,tmpcor2)
    tmpp2 = cor.test(stpm2vw,snowcover[,2])$p.value
    tmpp2 = scientific(tmpp2,2)
    pbox2 = c(pbox2,tmpp2)
    
  }
  tpmcor1 = tpmcor1[-1]
  tpmcor2 = tpmcor2[-1]
  pbox1 = pbox1[-1]
  pbox2 = pbox2[-1]
  
  colnames(stpm1_w_box) = c('700hpa','650hpa','600hpa')
  colnames(stpm2_w_box) = c('400hpa','350hpa','300hpa')
  
 
  
  dflabel1 = data.frame(
    x = 10,
    y = -3,
    label = paste0('ccr: ',round(tpmcor1,1),'(',pbox1,')'),
    type = paste0('(',letters[1:3],') ',hpa[1:3],'hpa at TPM1')
  )
  
  
  dflabel2 = data.frame(
    x = 10,
    y = -3,
    label = paste0('ccr: ',round(tpmcor2,1),'(',pbox2,')'),
    type = paste0('(',letters[4:6],') ',hpa[7:9],'hpa at TPM2')
  )
  
  dflabel2$x[3] = 80
  
  
  datelabel = seq(as.Date('2003-07-01'),as.Date('2016-12-01'),'1 month')
  datelabel = substr(datelabel,1,4)
  
  colnames(stpm1_w_box) = c('TPM1_W1','TPM1_W2','TPM1_W3')
  colnames(stpm2_w_box) = c('TPM2_W1','TPM2_W2','TPM2_W3')
  
  plinedf1 = data.frame(
    Time = 1:162,
    #TPM1_SC = snowcover[,1],
    stpm1_w_box
    #TPM1_W = c(stpm1_w_box[,1],stpm1_w_box[,2],stpm1_w_box[,3])
  )
  
  plinedf1m = reshape2::melt(plinedf1,'Time')
  plinedf1m$type = rep(paste0('(',letters[1:3],') ',hpa[1:3],'hpa at TPM1'),each = 162)
  plinedf1m$cols = 'Vertical velocity of wind'
  
  snow_df1 = data.frame(
    Time = rep(1:162,3),
    TPM1_SC = rep(snowcover[,1],3),
    type  = rep(paste0('(',letters[1:3],') ',hpa[1:3],'hpa at TPM1'),each = 162),
    cols = 'Snow cover area'
  )
  
  plinedf2 = data.frame(
    Time = 1:162,
    stpm2_w_box
  )
  
  plinedf2m = reshape2::melt(plinedf2,'Time')
  plinedf2m$type = rep(paste0('(',letters[4:6],') ',hpa[7:9],'hpa at TPM2'),each = 162)
  plinedf2m$cols = 'Vertical velocity of wind'
  
  
  snow_df2 = data.frame(
    Time = rep(1:162,3),
    TPM2_SC = rep(snowcover[,2],3),
    type  = rep(paste0('(',letters[4:6],') ',hpa[7:9],'hpa at TPM2'),each = 162),
    cols = 'Snow cover area'
  )
  
  library(ggsci)
  cols = pal_lancet()(9)[c(1,7)]
  cols = c('Snow cover area'=cols[1],
           'Vertical velocity of wind' = cols[2])
  
  library(ggplot2)
  library(ggpubr)
  fontsize = 8
  textsize = 3
  theme_text = theme(
    axis.text = element_text(size = fontsize,color = 'black'),
    strip.text = element_text(size = fontsize,color = 'black'),
    legend.text = element_text(size = fontsize,color = 'black'),
    axis.title = element_text(size = fontsize,color = 'black',face = 'bold'),
    legend.title = element_text(size = fontsize,color = 'black',face = 'bold'),
    panel.grid = element_blank(),
    plot.background = element_blank()
  )

  
  
  
  pline = ggplot()+
    geom_line(data = plinedf1m,
              aes(x = Time,y = value,color = cols),
              linewidth = 1)+
    geom_line(data = plinedf2m,
              aes(x = Time,y = value,color = cols),
              linewidth = 1)+
    geom_line(data = snow_df1,
              aes(x = Time,y = TPM1_SC,color = cols),
              linewidth = 1)+
    geom_line(data = snow_df2,
              aes(x = Time,y = TPM2_SC,color = cols),
              linewidth = 1)+
    geom_text(data = dflabel1,
              aes(x = x,y = y,label = label),
              size = textsize,color = 'black',hjust = 0)+
    geom_text(data = dflabel2,
              aes(x = x,y = y,label = label),
              size = textsize,color = 'black',hjust = 0)+
    facet_wrap(~type,nrow = 2)+
    theme_bw()+
    theme_text+
    scale_color_manual(values = cols)+
    scale_x_continuous(breaks = c(1,81,162),
                       labels = datelabel[c(1,81,162)])+
    theme(legend.position = 'bottom')+
    guides(color = guide_legend(title = 'Standardized indices',title.position = 'top',nrow = 1))+
    xlab("Time")+
    ylab('Standardized indices')
    
  output_fig = 'output_rebuttal/fig/relation_w_snowcover2.pdf'
  ggsave(output_fig,plot = pline,dpi = 800,
         width = 18,height = 12,units = 'cm')
    
    
}