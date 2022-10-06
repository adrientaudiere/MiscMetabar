data("GlobalPatterns")
GP <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table, 1 == "Archaea")
expect_s3_class()

accu_plot(GP, "SampleType", nb_seq = TRUE, by.fact = TRUE)
