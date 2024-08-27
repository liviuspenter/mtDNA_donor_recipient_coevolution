clinical.data <- readxl::read_excel("./data/CRC/20231127_CRC_mtDNA_meta.xlsx", sheet = 3) %>% as.data.frame()
clinical.data <- clinical.data[-which(clinical.data$Day < 0), ]

stable.patients <- c(15, 36, 14, 17, 18, 2, 40, 41)
natural.patients <- c(12, 3, 7, 1, 9, 19, 4)
dynamic.patients <- c(11, 20, 21, 39, 6, 10, 37, 13, 5, 16, 38)
patients.order <- c(dynamic.patients, natural.patients, stable.patients)

clinical.data <- clinical.data[which(clinical.data$Patient %in% c(stable.patients, natural.patients, dynamic.patients)), ]
clinical.data$Patient <- factor(clinical.data$Patient,
  levels = as.character(patients.order)
)

sample.data <- readxl::read_excel("./data/CRC/20231127_CRC_mtDNA_meta.xlsx", sheet = 1) %>%
  filter(Patient %in% c(stable.patients, natural.patients, dynamic.patients)) %>%
  as.data.frame()
sample.data$Patient <- factor(sample.data$Patient,
  levels = as.character(patients.order)
)

last.day <- data.frame()
for (p in unique(clinical.data$Patient)) {
  date.last.sample <- max(sample.data$Day[which(sample.data$Patient == p)])
  clinical.data <- clinical.data[-which(clinical.data$Day > date.last.sample & clinical.data$Patient == p), ]
  last.day <- rbind(last.day, data.frame(Patient = p, Day = max(clinical.data$Day[which(clinical.data$Patient == p)])))
}

ggplot() +
  geom_bar(
    data = clinical.data[order(-clinical.data$Day), ],
    aes(x = Patient, y = Day, fill = log10(WBC)),
    stat = "identity", position = "identity"
  ) +
  geom_bar(data = last.day, aes(x = Patient, y = Day), stat = "identity", position = "identity", color = "black", fill = NA) +
  geom_point(data = sample.data, aes(x = Patient, y = Day, color = Description)) +
  # geom_col(data=clinical.data[order(-clinical.data$Day),], aes(x=as.character(Patient), y=max(as.numeric(Day))/100), fill=NA, color='black') +
  # scale_fill_gradientn(colors = BuenColors::jdb_palette(name = 'brewer_spectra')) +
  scale_fill_gradient(low = "white", high = "orange") +
  scale_color_manual(values = c("pretreatment" = "blue", "posttreatment" = "red", "posttreatment2" = "firebrick")) +
  scale_x_discrete("CLL") +
  scale_y_continuous("Days from baseline sample") +
  coord_flip() +
  theme_classic() +
  theme( # legend.position = 'none',
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./figure_CLL_WES_RNA/figures/CRC/20240118_swimmers_plot.svg", width = 4, height = 4)
