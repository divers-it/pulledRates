rm(list = ls())

library(ape)
library(phytools)
library(RColorBrewer)

palette(brewer.pal(3, "Set2"))
#LTT only

pdf("figures/fig1b.pdf", height = 5, width = 5)

par(mfrow = c(1, 1))

dips_phy <- read.tree("data/dipsidae.tre")
ltt.plot(
  dips_phy,
  cex = 0.5,
  ylab = "No. of lineages",
  log = "y",
  xlim = c(-300, 0),
  col = 1,
  lwd = 2,
  lty = 1,
  xlab = bquote("Age (" * tau * ")"),
  xaxt = "n"
)
axis(1, at=pretty(range(c(-300, 0))), labels=rev(pretty(range(-c(-300, 0)))))

cycad_phy <- read.tree("data/Cycadales_MCC.tre")
ltt.lines(
  cycad_phy,
  cex = 0.5,
  xlab = "Time (Ma)",
  ylab = "No. of lineages",
  col = 2,
  lwd = 2
)

whales_phy <- read.tree("data/Cetaceae-timetree.tre")
ltt.lines(
  force.ultrametric(whales_phy),
  cex = 0.5,
  xlab = "Time (Ma)",
  ylab = "No. of lineages",
  col = 3,
  lwd = 2
)

legend(
  "topleft",
  legend = c("Campanulids", "Cycadales", "Cetaceae"),
  bty = "n",
  col = c(1, 2, 3),
  lty = 1,
  lwd = 2
)

dev.off()
