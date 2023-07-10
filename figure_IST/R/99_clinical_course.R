library(ggplot2)

IST1.clinical.data = readxl::read_excel('./data/IST/clinical_data/20220119_CBC_values.xlsx', sheet = 1)
IST1.clinical.data$Patient = 'IST1'
IST3.clinical.data = readxl::read_excel('./data/IST/clinical_data/20220119_CBC_values.xlsx', sheet = 2)
IST3.clinical.data$Patient = 'IST2'
IST4.clinical.data = readxl::read_excel('./data/IST/clinical_data/20220119_CBC_values.xlsx', sheet = 3)
IST4.clinical.data$Patient = 'IST3'
IST5.clinical.data = readxl::read_excel('./data/IST/clinical_data/20220119_CBC_values.xlsx', sheet = 4)
IST5.clinical.data$Patient = 'IST4'
IST.clinical.data = rbind(IST1.clinical.data, IST3.clinical.data)
IST.clinical.data = rbind(IST.clinical.data, IST4.clinical.data)
IST.clinical.data = rbind(IST.clinical.data, IST5.clinical.data)

Taper.dates = data.frame(Start = c(min(IST1.clinical.data$Day[which(IST1.clinical.data$IST != 'NA')]),
                                   min(IST3.clinical.data$Day[which(IST3.clinical.data$IST != 'NA')]),
                                   min(IST4.clinical.data$Day[which(IST4.clinical.data$IST != 'NA')]),
                                   min(IST5.clinical.data$Day[which(IST5.clinical.data$IST != 'NA')])),
                         End = c(max(IST1.clinical.data$Day[which(IST1.clinical.data$IST != 'NA')]),
                                 max(IST3.clinical.data$Day[which(IST3.clinical.data$IST != 'NA')]),
                                 max(IST4.clinical.data$Day[which(IST4.clinical.data$IST != 'NA')]),
                                 max(IST5.clinical.data$Day[which(IST5.clinical.data$IST != 'NA')])),
                         Y = c(1.1*max(IST1.clinical.data$WBC),
                               1.1*max(IST3.clinical.data$WBC),
                               1.1*max(IST4.clinical.data$WBC),
                               1.1*max(IST5.clinical.data$WBC[which(IST5.clinical.data$WBC != 'NA')])),
                         Patient = c('IST1', 'IST3', 'IST4', 'IST5'))

for (patient in unique(IST.clinical.data$Patient)) {
  p1=ggplot() + 
    geom_line(data=IST.clinical.data[which(IST.clinical.data$Patient == patient & IST.clinical.data$WBC != 'NA'),], aes(x=Day, y=WBC)) +
    geom_point(data=IST.clinical.data[which(IST.clinical.data$Patient == patient & IST.clinical.data$WBC != 'NA'),], aes(x=Day, y=WBC), size=0.5) +
    geom_segment(data=Taper.dates[which(Taper.dates$Patient == patient),], aes(x=Start, xend=End, y=Y, yend=Y), size=3, color='magenta') + 
    scale_x_continuous(limits = c(0, max(IST.clinical.data$Day[which(IST.clinical.data$Patient == patient)]))) + 
    scale_y_continuous('WBC / nl') + 
    theme_classic() +
    theme(axis.title = element_text('Arial', size=10, color='black'),
          axis.text = element_text('Arial', size=10, color='black'),
          axis.title.x = element_blank())
  ggsave(paste0('./figure_IST/figures/clinical_plots/20220120_',patient,'_WBC.svg'), width = 2.5, height = 1.5, plot=p1)
  
  p2=ggplot() + 
    geom_line(data=IST.clinical.data[which(IST.clinical.data$Patient == patient),], aes(x=Day, y=as.numeric(Platelets)), color='orange') +
    geom_point(data=IST.clinical.data[which(IST.clinical.data$Patient == patient),], aes(x=Day, y=as.numeric(Platelets)), size=0.5) +
    scale_x_continuous(limits = c(0, max(IST.clinical.data$Day[which(IST.clinical.data$Patient == patient)]))) + 
    scale_y_continuous('Platelets / nl') + 
    theme_classic() +
    theme(axis.title = element_text('Arial', size=10, color='black'),
          axis.text = element_text('Arial', size=10, color='black'),
          axis.title.x = element_blank())
  ggsave(paste0('./figure_IST/figures/clinical_plots/20220120_',patient,'_Platelets.svg'), width = 2.5, height = 1.5, plot=p2)
  
  p3=ggplot() + 
    geom_line(data=IST.clinical.data[which(IST.clinical.data$Patient == patient & IST.clinical.data$Chimerism != 'NA'),], 
              aes(x=Day, y=as.numeric(Chimerism)), color='red') +
    geom_point(data=IST.clinical.data[which(IST.clinical.data$Patient == patient & IST.clinical.data$Chimerism != 'NA'),], 
               aes(x=Day, y=as.numeric(Chimerism)), size=0.5) +
    scale_x_continuous(limits = c(0, max(IST.clinical.data$Day[which(IST.clinical.data$Patient == patient)]))) + 
    scale_y_continuous('% Chimerim') + 
    theme_classic() +
    theme(axis.title = element_text('Arial', size=10, color='black'),
          axis.text = element_text('Arial', size=10, color='black'),
          axis.title.x = element_blank())
  ggsave(paste0('./figure_IST/figures/clinical_plots/20220120_',patient,'_Chimerism.svg'), width = 2.5, height = 1.5, plot=p3)
  
  p4=ggplot() + 
    geom_line(data=IST.clinical.data[which(IST.clinical.data$Patient == patient & IST.clinical.data$ALC != 'NA'),], 
              aes(x=Day, y=as.numeric(ALC)), color='blue') +
    geom_point(data=IST.clinical.data[which(IST.clinical.data$Patient == patient & IST.clinical.data$ALC != 'NA'),], 
               aes(x=Day, y=as.numeric(ALC)), size=0.5) +
    scale_x_continuous(limits = c(0, max(IST.clinical.data$Day[which(IST.clinical.data$Patient == patient)]))) + 
    scale_y_continuous('Lymphocytes / nl') + 
    theme_classic() +
    theme(axis.title = element_text('Arial', size=10, color='black'),
          axis.text = element_text('Arial', size=10, color='black'),
          axis.title.x = element_blank())
  ggsave(paste0('./figure_IST/figures/clinical_plots/20220120_',patient,'_Lymphocytes.svg'), width = 2.5, height = 1.5, plot=p4)

  p5=ggplot() + 
    geom_line(data=IST.clinical.data[which(IST.clinical.data$Patient == patient & IST.clinical.data$Mono != 'NA'),], 
              aes(x=Day, y=as.numeric(Mono)), color='darkgreen') +
    geom_point(data=IST.clinical.data[which(IST.clinical.data$Patient == patient & IST.clinical.data$Mono != 'NA'),], 
               aes(x=Day, y=as.numeric(Mono)), size=0.5) +
    scale_x_continuous(limits = c(0, max(IST.clinical.data$Day[which(IST.clinical.data$Patient == patient)]))) + 
    scale_y_continuous('Monocytes / nl') + 
    theme_classic() +
    theme(axis.title = element_text('Arial', size=10, color='black'),
          axis.text = element_text('Arial', size=10, color='black'),
          axis.title.x = element_blank())
  ggsave(paste0('./figure_IST/figures/clinical_plots/20220120_',patient,'_Monocytes.svg'), width = 2.5, height = 1.5, plot=p5)
  
  p=cowplot::plot_grid(plotlist = list(p1, p3, p2, p4, p5), ncol = 1, align='v')
  ggsave(paste0('./figure_IST/figures/clinical_plots/20220120_',patient,'_all.svg'), width = 2.5, height = 7.5, plot=p)
  
}