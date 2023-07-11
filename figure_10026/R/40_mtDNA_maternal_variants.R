# maternal mtDNA variants used for donor-recipient deconvolution

donor.variants.1010 = gtools::mixedsort(c('195T>C', '150C>T', '6146A>G', '1811A>G', '10907T>C', '6047A>G', '5999T>C', '14866C>T', '11009T>C', '11467A>G', '9070T>G', '15693T>C', '12308A>G',
                                     '12372G>A', '4646T>C', '14620C>T', '15530T>C', '4811A>G', '499G>A', '16356T>C', '11332C>T', '16179C>T'))
recipient.variants.1010 = gtools::mixedsort(c('16126T>C', '4216T>C', '10084T>C', '489T>C', '462C>T', '11251A>G', '12612A>G', '15452C>A', '3010G>A', '14798T>C', '13708G>A', 
                                         '16069C>T', '16319G>A', '295C>T', '55T>C', '56A>G', '185G>A', '228G>A'))

donor.variants.1026 = gtools::mixedsort(c('10463T>C', '16519T>C', '15928G>A', '16126T>C', '16153G>A', '73A>G', '15607A>G', '15452C>A', '7028C>T', '4917A>G',
                                     '13368G>A', '4216T>C', '8697G>A', '11812A>G', '14766C>T', '11719G>A', '709G>A', '8269G>A', '14905G>A', '2706A>G',
                                     '11251A>G', '14233A>G', '1888G>A', '9947G>A', '150C>T', '16294C>T', '16296C>T'))
recipient.variants.1026 = gtools::mixedsort(c('16304T>C', '456C>T', '8433T>C', '15833C>T', '4336T>C', '9722T>C', '4011C>T', '93A>G'))

germline.variants = data.frame()
germline.variants = rbind(germline.variants, data.frame(
  variant = c(donor.variants.1010, recipient.variants.1010),
  individual = c(rep('donor', length(donor.variants.1010)), rep('recipient', length(recipient.variants.1010))),
  sample = rep('AML1010', length(c(donor.variants.1010, recipient.variants.1010)))
))
germline.variants = rbind(germline.variants, data.frame(
  variant = c(donor.variants.1026, recipient.variants.1026),
  individual = c(rep('donor', length(donor.variants.1026)), rep('recipient', length(recipient.variants.1026))),
  sample = rep('AML1026', length(c(donor.variants.1026, recipient.variants.1026)))
))
write.csv2(file='./data/10026/mtDNA/20220117_AML_germline_variants.csv', germline.variants, quote = F, row.names = F)

length(donor.variants.1010)
length(donor.variants.1026)

length(recipient.variants.1010)
length(recipient.variants.1026)


variant.positions = c(as.numeric(gsub("([0-9]+).*$", "\\1", donor.variants.1010)), 
                      as.numeric(gsub("([0-9]+).*$", "\\1", donor.variants.1026)), 
                      as.numeric(gsub("([0-9]+).*$", "\\1", recipient.variants.1010)), 
                      as.numeric(gsub("([0-9]+).*$", "\\1", recipient.variants.1026)))
                      