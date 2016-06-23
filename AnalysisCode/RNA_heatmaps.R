library(dplyr)
library(RColorBrewer)
library(readr)

color <- colorRampPalette(brewer.pal(9,"Reds"))(100)

expression <- read.table("../RawData/gene_expression.tsv.gz", header = T, row.names = 1)
exon <- read_tsv('../ProcessedData/GeneFeatures/exon_intersections.bed.gz', col_names = F)

lookupExp <- function(l, gene, exp) {
  x <- vector()
  accessions <- unlist(strsplit(gsub("-", "_", l), ","))
  for(acc in accessions) {
    x <- c(x, as.numeric(exp[gene,acc]))
  }
  return(x)
}

pos <- c("Db-1,Uk-1,Pu2-23,Bl-1,La-0,Bor-4,Bd-0,Ta-0,Krot-0,Bch-1,Amel-1,Fi-0,In-0,Knox-18,Kyoto,En-D,Zdr-1,Dra-0,RRS-7,Ting-1,Di-G,Jl-3,Uod-1,Ak-1,Gu-0,Mz-0,Et-0,Ove-0,Jm-0,Sq-8,Hs-0,Wl-0,Appt-1,Durh-1,Pu2-7,Ba-1")
neg <- c("Ty-0,Sp-0,Su-0,Gre-0,Ei-2,Nok-3,Wa-1,Gel-1,Wt-5,Ra-0,Tu-0,Gifu-2,Np-0,Bs-1,Com-1,Ha-0,Old-1,Bu-0,Qar-8a,Rennes-1,Kro-0,Kelsterbach-4,Cnt-1,Anholt-1,Cerv-1,Kondara,Rome-1,Etna-2,Nw-0,Lan-0,Co-1,Seattle-0,Gie-0,Ann-1,Gr-1,An-1,Sus-1,Sg-1,Utrecht,Aa-0,Got-7,Bsch-0,Stw-0,Pla-0,Kil-0,Kas-1,Per-1,Anz-0,Chi-0,Rhen-1,Ts-1,El-0,Er-0,Chat-1,Fr-2,Tamm-2,Je-0,Bla-1,Vind-1,Tol-0,Dja-1,Ca-0,Gy-0,Col-0,Wc-1,Ang-0,Si-0,Kin-0,Abd-0,Hey-1,Alst-1,HR-5,Bik-1,Tul-0,Boot-1,Neo-6,Sei-0,Ven-1,Est,Cal-0,Baa-1,Rld-1,Hau-0,Pna-17,Ga-0,Kl-5,Ob-0,RRS-10,Yo-0,Do-0,Pog-0,Hh-0,Kar-1,Kz-9,Is-0,Nc-1,Ws-2,Van-0,Se-0,Lp2-2,Altai-5,Benk-1,Tscha-1,Es-0")


AT2G15040.pos.exp <- as.data.frame(lookupExp(pos, "AT2G15040", expression))
AT2G15040.neg.exp <- as.data.frame(lookupExp(neg, "AT2G15040", expression))

AT2G15042.pos.exp <- as.data.frame(lookupExp(pos, "AT2G15042", expression))
AT2G15042.neg.exp <- as.data.frame(lookupExp(neg, "AT2G15042", expression))

AT2G15040.all <- rbind( AT2G15040.neg.exp , setNames( AT2G15040.pos.exp , names( AT2G15040.neg.exp ) ) )
AT2G15042.all <- rbind( AT2G15042.neg.exp , setNames( AT2G15042.pos.exp , names( AT2G15042.neg.exp ) ) )

all <- cbind(AT2G15040.all, AT2G15042.all)

all.floor <- all > 2500000
all[all.floor] <- 2500000

pdf("../Plots/atg215040_atg215042.pdf", height=5, width=3)
image(as.matrix(t(all)), col = color)
abline(h = nrow(AT2G15040.neg.exp)/nrow(all), lwd=2, col="blue")
dev.off()

# upstream insertion causing activation

up.pos <- c("Gre-0,Tol-0,Tul-0,RRS-10,Gifu-2,Pog-0,Pna-17,Gel-1")
up.neg <- c("Uk-1,La-0,Sp-0,Krot-0,Amel-1,Fi-0,In-0,Knox-18,Ei-2,Nok-3,Wa-1,Etna-2,Su-0,Jl-3,Wt-5,Ra-0,Gu-0,Tu-0,Mz-0,Np-0,Bs-1,Com-1,Ha-0,Rennes-1,Ty-0,Jm-0,Qar-8a,Kro-0,Kelsterbach-4,Ba-1,Cnt-1,Anholt-1,Cerv-1,Kondara,Rome-1,Nw-0,Lan-0,Bu-0,Co-1,Seattle-0,Gie-0,Ann-1,Gr-1,An-1,Di-G,Sus-1,Sg-1,Utrecht,Aa-0,Got-7,Bsch-0,Stw-0,Pla-0,Kil-0,Kas-1,Per-1,Chi-0,Rhen-1,Ts-1,El-0,Er-0,Pu2-7,Ven-1,Chat-1,Bl-1,Fr-2,Bd-0,Tamm-2,Ta-0,Je-0,Sq-8,Bla-1,Vind-1,Dra-0,Uod-1,Dja-1,Ca-0,Gy-0,Col-0,Wc-1,Ang-0,Si-0,Kin-0,Abd-0,Hey-1,Alst-1,Ove-0,HR-5,Bik-1,Boot-1,Neo-6,Wl-0,Anz-0,Db-1,Est,Cal-0,Baa-1,Rld-1,Hau-0,Bor-4,Sei-0,Bch-1,Zdr-1,Ga-0,Kyoto,En-D,Kl-5,Ob-0,Yo-0,Do-0,Old-1,Pu2-23,RRS-7,Hh-0,Ting-1,Kar-1,Kz-9,Ak-1,Is-0,Nc-1,Ws-2,Van-0,Se-0,Et-0,Lp2-2,Altai-5,Hs-0,Benk-1,Appt-1,Durh-1,Tscha-1,Es-0")

AT2G01360.pos.exp <- as.data.frame(lookupExp(up.pos, "AT2G01360", expression))
AT2G01360.neg.exp <- as.data.frame(lookupExp(up.neg, "AT2G01360", expression))

AT2G01360.all <- rbind( AT2G01360.neg.exp , setNames( AT2G01360.pos.exp , names( AT2G01360.neg.exp ) ) )

# set max value for heatmap
AT2G01360.floor <- AT2G01360.all > 590000
AT2G01360.all[AT2G01360.floor] <- 590000

pdf("../Plots/AT2G01360_heatmap.pdf", height=5, width=3)
image(as.matrix(t(AT2G01360.all)), col = color)
abline(h = nrow(AT2G01360.neg.exp)/nrow(AT2G01360.all), lwd=2, col="blue")
dev.off()

# AT2G01350

AT2G01350.pos.exp <- as.data.frame(lookupExp(up.pos, "AT2G01350", expression))
AT2G01350.neg.exp <- as.data.frame(lookupExp(up.neg, "AT2G01350", expression))

AT2G01350.all <- rbind( AT2G01350.neg.exp , setNames( AT2G01350.pos.exp , names( AT2G01350.neg.exp ) ) )

# max(AT2G01350.all)

# set max value for heatmap
AT2G01350.floor <- AT2G01350.all > 3100000
AT2G01350.all[AT2G01350.floor] <- 3100000

pdf("../Plots/AT2G01350_heatmap.pdf", height=5, width=3)
image(as.matrix(t(AT2G01350.all)), col = color)
abline(h = nrow(AT2G01350.neg.exp)/nrow(AT2G01350.all), lwd=2, col="blue")
dev.off()