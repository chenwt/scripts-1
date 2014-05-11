install.packages("plot3D")
require(plot3D)
data(Oxsat)
head(Oxsat)
par(mfrow = c(1, 1))

datafile = ""
z 
hist3D( x = seq(0, 1, length.out = nrow(z)),
        y = seq(0, 1, length.out = ncol(z)), z, 
        colvar = z, phi = 45, theta = 20,
        xlab = "ceRNA driver", ylab ="selected tumor samples", 
        col = NULL, NAcol = "white", border = NA, facets = TRUE,
        colkey = NULL, image = FALSE, contour = FALSE,
        panel.first = NULL, clim = NULL, clab = NULL, bty = "b",
        lighting = T, shade = F, ltheta = -135, lphi = 0,
        space = 0.2, add = F, plot = TRUE)
