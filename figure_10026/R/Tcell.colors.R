
Tcell.colors = c('CD4.naive' = RColorBrewer::brewer.pal(6, name = 'Blues')[2],
                 'CD4.CM' = RColorBrewer::brewer.pal(6, name = 'Blues')[3],
                 'CD4.EM' = RColorBrewer::brewer.pal(6, name = 'Blues')[4],
                 'CD4.TEMRA' = RColorBrewer::brewer.pal(6, name = 'Blues')[5],
                 'CD4.senescent' = RColorBrewer::brewer.pal(6, name = 'Blues')[6],
                 'CD8' = RColorBrewer::brewer.pal(4, name = 'Reds')[2],
                 'CD8.CM' = RColorBrewer::brewer.pal(4, name = 'Reds')[3],
                 'CD8.senescent' = RColorBrewer::brewer.pal(4, name = 'Reds')[4],
                 'NK' = 'black')
                 
AML.1010.mito.colors = c('HSCT' = RColorBrewer::brewer.pal(n=5, name='Reds')[5], 
                         'GMP' = RColorBrewer::brewer.pal(n=5, name='Reds')[4], 
                         'Mono' = RColorBrewer::brewer.pal(n=5, name='Reds')[3],
                         'erythroid' = RColorBrewer::brewer.pal(n=5, name='Reds')[2], 
                         'CD4' = RColorBrewer::brewer.pal(n=3, name='Blues')[2],
                         'CD8' = RColorBrewer::brewer.pal(n=3, name='Blues')[3])

AML.1012.mito.colors = c('HSCT' = RColorBrewer::brewer.pal(n=6, name='Reds')[6], 
                         'GMP' = RColorBrewer::brewer.pal(n=6, name='Reds')[5], 
                         'Mono1' = RColorBrewer::brewer.pal(n=6, name='Reds')[4],
                         'Mono2' = RColorBrewer::brewer.pal(n=6, name='Reds')[3],
                         'erythroid' = RColorBrewer::brewer.pal(n=6, name='Reds')[2], 
                         'DC' = 'yellow',
                         'CD4' = RColorBrewer::brewer.pal(n=3, name='Blues')[2],
                         'CD8' = RColorBrewer::brewer.pal(n=3, name='Blues')[3],
                         'none' = 'grey')

AML.1026.mito.colors = c('HSCT' = RColorBrewer::brewer.pal(n=5, name='Reds')[5], 
                         'GMP' = RColorBrewer::brewer.pal(n=5, name='Reds')[4], 
                         'Mono' = RColorBrewer::brewer.pal(n=5, name='Reds')[3],
                         'erythroid' = RColorBrewer::brewer.pal(n=5, name='Reds')[2], 
                         'T cell' = RColorBrewer::brewer.pal(n=3, name='Blues')[3])