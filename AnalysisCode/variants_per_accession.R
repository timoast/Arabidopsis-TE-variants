library(readr)
library(dplyr)
library(ggplot2)

# deletions
dat <- read_tsv("../RawData/TEPID_TEPAV.tsv.gz") %>%
        filter(Absence_classification == "True deletion") %>%
        mutate(acc_list = ifelse(Minor_allele == "absence", Accessions_TE_absent, Accessions_TE_present)) %>%
        select(acc_list) %>%
        rowwise() %>%
        mutate(n_accessions = length(unlist(strsplit(acc_list, ","))))

accessions <- lapply(dat$acc_list, strsplit, split = ",", fixed = TRUE)
accessions <- lapply(accessions, unlist)

# make a dataframe
df <- data_frame(accessions = character())
for(i in seq_along(accessions)) {
        t_df <- data_frame(accession = accessions[[i]])
        df <- rbind(df, t_df)
}

df %>%
        group_by(accession) %>%
        mutate(true_deletions = n()) %>%
        unique() %>%
        ggplot(., aes(true_deletions)) + geom_histogram(color = "black", bins = 50) +
        theme_bw() +
        theme(text = element_text(size=20)) +
        geom_vline(xintercept = 8) +
        annotate("text", x = 30, y = 12, label = "Col-0") +
        xlab("True deleltions") + ylab("Number of accessions") +
        ggtitle("True deletions per accession") +
        ggsave("../Plots/deletions_per_accession.pdf", height = 10, width = 15, units = "cm")

# insertions
dat <- read_tsv("../RawData/TEPID_TEPAV.tsv.gz") %>%
        filter(Absence_classification == "No insertion") %>%
        mutate(acc_list = ifelse(Minor_allele == "presence", Accessions_TE_present, Accessions_TE_absent)) %>%
        select(acc_list) %>%
        rowwise() %>%
        mutate(n_accessions = length(unlist(strsplit(acc_list, ","))))

insertions <- lapply(dat$acc_list, strsplit, split = ",", fixed = TRUE)
insertions <- lapply(insertions, unlist)
# make a dataframe
ins_df <- data_frame(accessions = character())
for(i in seq_along(insertions)) {
        t_df <- data_frame(accession = insertions[[i]])
        ins_df <- rbind(ins_df, t_df)
}

ins_df %>%
        group_by(accession) %>%
        mutate(insertions = n()) %>%
        unique() %>%
        ggplot(., aes(insertions)) + geom_histogram(color = "black", bins = 50) +
        theme_bw() +
        theme(text = element_text(size=20)) +
        geom_vline(xintercept = 388) +
        annotate("text", x = 410, y = 15, label = "Col-0") +
        xlab("TE insertions") + ylab("Number of accessions") +
        ggtitle("True insertions per accession") + 
        ggsave("../Plots/insertions_per_accession.pdf", height = 10, width = 15, units = "cm")
