load("C:/Temp/20251217/taxid.hierarchy.Rdata")
load("C:/Temp/20251217/cid.taxid.Rdata")

levs.all.taxa <- sapply(1:ncol(taxid.hierarchy), FUN = function(x) length(unique(taxid.hierarchy[,x])))
names(levs.all.taxa) <- names(taxid.hierarchy)

cid.taxid <- cid.taxid[cid.taxid$data.source %in% "LOTUS - the natural products occurrence database", ]


glc <- 5793
glc.taxids <- cid.taxid$taxid[cid.taxid$cid %in% glc]

glc.tax.h <- lapply(1:ncol(taxid.hierarchy), FUN = function(x) {
  which(taxid.hierarchy[,x] %in% glc.taxids)
}
)
glc.tax.h <- unlist(glc.tax.h)
glc.tax.h <- sort(unique(glc.tax.h))
glc.tax.h <- taxid.hierarchy[glc.tax.h,]

levs.glc.taxa <- sapply(1:ncol(glc.tax.h), FUN = function(x) length(unique(glc.tax.h[,x])))
names(levs.glc.taxa) <- names(glc.tax.h)

par(mfrow = c(2,1))
use <- !grepl('sub', names(levs.all.taxa))
plot(log(levs.all.taxa[use]))
plot(log(levs.glc.taxa[use]))

fc <- 5280371
cid.taxid$taxid[cid.taxid$cid %in% fc]
fc.taxids <- cid.taxid$taxid[cid.taxid$cid %in% fc]

fc.tax.h <- lapply(1:ncol(taxid.hierarchy), FUN = function(x) {
  which(taxid.hierarchy[,x] %in% fc.taxids)
}
)
fc.tax.h <- unlist(fc.tax.h)
fc.tax.h <- sort(unique(fc.tax.h))
fc.tax.h <- taxid.hierarchy[fc.tax.h,]

levs.fc.taxa <- sapply(1:ncol(fc.tax.h), FUN = function(x) length(unique(fc.tax.h[,x])))
names(levs.fc.taxa) <- names(fc.tax.h)

par(mfrow = c(3,1))
use <- !grepl('sub', names(levs.all.taxa))
# use <- 1:length(levs.all.taxa)
plot(log(levs.all.taxa[use]), type = 'l', main = "all taxa")
plot(log(levs.glc.taxa[use]), type = 'l', main = "glucose occurance")
plot(log(levs.fc.taxa[use]), type = 'l', main = "Bergaptol occurance")

occ <- data.frame(
  "all.taxa" = log(levs.all.taxa[use]),
  "glucose" = log(levs.glc.taxa[use]),
  "bergaptol" = log(levs.fc.taxa[use])
)
for(i in 1:ncol(occ)) {
  occ[,i] <- occ[,i]/max(occ[,i], na.rm = TRUE)
}
occ$level <- factor(row.names(occ), levels = row.names(occ))
occ <- data.frame(
  taxa.level = rep(occ$level, 3),
  compound = c(rep("all.taxa", nrow(occ)), rep("glucose", nrow(occ)), rep("bergaptol", nrow(occ))),
  taxa.count = c(occ$all.taxa, occ$glucose, occ$bergaptol)
  )

ggplot(occ, aes(x = taxa.level, y = taxa.count, group = compound, color = compound)) +
  geom_line() +
  labs(title = "")

levels
library(ggplot2)
ggplot(occ, aes(x = level)) +
  geom_line(aes(y = all.taxa, group = 1), color = "black") + # group=1 ensures lines are connected
  geom_line(aes(y = glucose, group = 1), color = "blue") +
  geom_line(aes(y = bergaptol, group = 1), color = "red") +
  theme(axis.title.y = element_blank()) +
  scale_color_manual(
    name = "",        # Custom legend title
    values = c("black", "blue", "red"),     # Custom colors
    labels = c("all.taxa", "glucose", "bergaptol")      # Custom labels
  ) +
  theme_minimal() +
  guides(color = guide_legend(nrow = 3))

