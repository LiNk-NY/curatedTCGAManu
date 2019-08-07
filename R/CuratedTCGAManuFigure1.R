library(ggplot2)
library(ggrepel)
xvec <- c(2, 4, 3)
yvec <- c(2, 2, 4)

triangle <- data.frame(x = xvec, y = yvec,
    r = c("ease-of-use", "completeness" ,"integration"))

midpoint <- function(x, y) {
    xa <- Reduce(`+`, x) / 2
    yb <- Reduce(`+`, y) / 2
    c(xa, yb)
}

cmat <- cbind(xvec, yvec)
mids <- apply(combn(1:3, 2), 2L, function(x) {
    midpoint(cmat[x, 1L], cmat[x, 2L])
})

mids2 <- t(mids)
cmat2 <- cmat[rev(seq_len(nrow(cmat))), ]

m <- mids2[1, ]
n <- cmat2[1, ]

midxy <- c(n[1]+0.666666*(m[1]-n[1]),  n[2]+0.666666*(m[2]-n[2]))

# x11()

# ggplot() + xlim(1.5, 4.5) + ylim(1.5, 4.5) +
#     geom_polygon(data = triangle, mapping = aes(x = x, y = y), fill = "white",
#         size = 2, color = "black") +
#     geom_text(data = triangle, aes(x = x, y = y, label = r),
#         hjust = c(0.5, 0.5, 0.5), vjust = c(1.2, 1.2, -0.5), size = 8) +
#     theme_bw() +
#     geom_segment(aes(x = cmat[, 1], y = cmat[, 2], xend = 3, yend = 2.666))

# library(ggmap)
# gglocator(n = 7)

ratings <- data.frame(
        x = c(2.87957519262172, 2.62101863265239, 3.23234225831191,
            3.21648660231498, 3.15130223877203, 3.36447272495301,
            3.00331611613399),
        y = c(2.71593776544351,  2.23753058101305, 3.23145436778364,
            2.90589275552184, 2.85162307159464, 2.95337659029189,
            2.19002009185595),
        resource = c("curatedTCGAData", "cBioPortal", "cgdsr", "firebrowseR",
    "TCGAbiolinks", "GenomicDataCommons", "GDAC Firehose")
)

png("ggtriangle_TCGAresources.png", width = 9, height = 8, units = "in", res = 300)

ggplot() + xlim(1.75, 4.25) + ylim(1.75, 4.25) +
    geom_polygon(data = triangle, mapping = aes(x = x, y = y), fill = "white",
        size = 1, color = "black") +
    geom_segment(aes(x = cmat[, 1], y = cmat[, 2], xend = 3, yend = 2.666),
        lty = 3, color = "grey80", lwd = 2) +
    geom_text(data = triangle, aes(x = x, y = y, label = r),
        hjust = c(0.5, 0.5, 0.5), vjust = c(1.2, 1.2, -0.5), size = 8,
        fontface = "bold") +
    theme(
        axis.line=element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.border=element_blank(), plot.background=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()
    ) +
    geom_point(data = ratings, mapping = aes(x = x, y = y), color = "red",
        size = 3) +
    geom_text_repel(data = ratings, mapping = aes(x, y, label = resource),
        size = 6, fontface = "bold")

dev.off()
