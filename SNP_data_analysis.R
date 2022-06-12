interesting %>%
  group_by(SNP, Magic) %>%
  summarise(Trait = str_c(unique(Trait), collapse = " | "),
            n_trials = length(unique(Trial)),
            Method = str_c(unique(str_c(Phenotype_transformation, GWAS_method, sep = "_")), collapse = " | "), n = n(),
            Project = str_c(unique(Project), collapse = " | "),
            n_projects = length(unique(Project))) %>%
  write_csv("interesting_SNPs_summarized.csv")

SNP %>%
  filter(Trait %like% "ilum") %>%
  ggplot(mapping = aes(x = Position, y = -log10(P.value), color = as_factor(Chromosome), size = -log10(P.value))) +
  geom_point() 
  scale_x_continuous(label = Chromosome#, breaks = axis_set$center
                     ) +
  #  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  #  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  #  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chr)))) +
   # scale_size_continuous(range = c(0.5,3)) +
  #  labs(x = NULL, 
  #       y = "-log<sub>10</sub>(p)") + 
  #  theme_minimal() +
  #  theme( 
   #   legend.position = "none",
  #    panel.border = element_blank(),
  #    panel.grid.major.x = element_blank(),
  #    panel.grid.minor.x = element_blank(),
  #    axis.title.y = element_markdown(),
  #    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
   # )
NULL
