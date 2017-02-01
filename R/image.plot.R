# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2
"image.plot" <- function(..., add = FALSE,
    breaks= NULL, nlevel = 64, col = NULL,
    horizontal = FALSE, legend.shrink = 0.9, legend.width = 1.2,
    legend.mar = ifelse(horizontal, 3.1, 5.1), legend.lab = NULL,
    legend.line= 2,
    graphics.reset = FALSE, bigplot = NULL, smallplot = NULL,
    legend.only = FALSE,  lab.breaks = NULL,
    axis.args = NULL, legend.args = NULL, legend.cex=1.0, midpoint = FALSE, border = NA,
    lwd = 1, verbose=FALSE) {
    # Thanks to S. Koehler and  S. Woodhead
    # for comments on making this a better function
    #
    # save current graphics settings
    old.par <- par(no.readonly = TRUE)
    # set defaults for color scale
    # note this works differently than the image function.
    if( is.null(col))  {
    	col<-  tim.colors(nlevel)}
    	else{
    		nlevel<- length( col)
    		}
    #  figure out zlim from passed arguments
    #  also set the breaks for colors if they have not been passed,
    info <- imagePlotInfo(..., breaks=breaks, nlevel=nlevel)
    # breaks have been computed if not passed in the call
    breaks<- info$breaks
    if( verbose){
    	print(info)
    }
    if (add) {
        big.plot <- old.par$plt
    }
    if (legend.only) {
        graphics.reset <- TRUE
    }
    if (is.null(legend.mar)) {
        legend.mar <- ifelse(horizontal, 3.1, 5.1)
    }
    # figure out how to divide up the plotting real estate
    temp <- imageplot.setup(add = add, legend.shrink = legend.shrink,
        legend.width = legend.width, legend.mar = legend.mar,
        horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
    # bigplot has plotting region coordinates for image
    # smallplot has plotting coordinates for legend strip
    smallplot <- temp$smallplot
    bigplot <- temp$bigplot
    # draw the image in bigplot, just call the R base function
    # or poly.image for polygonal cells
    # note the logical switch
    # for poly.grid is parsed out of call from image.plot.info
    if (!legend.only) {
        if (!add) {
            par(plt = bigplot)
        }
        if (!info$poly.grid) {
            image(..., breaks=breaks, add = add, col = col)
        }
        else {
            poly.image(..., add = add, col = col, midpoint = midpoint,
                border = border, lwd.poly = lwd)
        }
        big.par <- par(no.readonly = TRUE)
    }
    ##
    ## check dimensions of smallplot
    if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
        par(old.par)
        stop("plot region too small to add legend\n")
    }
    # Following code draws the legend using the image function
    # and a one column image.
    # What might be confusing is the values of the "image" are the same
    # as the locations on the legend axis.
    # Moreover the image values are in the middle of each breakpoint category
    # thanks to Tobias Nanu Frechen and Matthew Flickinger
    # for sorting out some problems with the breaks position in the legend.
    ix <- 1:2
    iy<- breaks
    nBreaks<- length( breaks)
    midpoints<- (breaks[1:(nBreaks-1)] +  breaks[2:nBreaks] )/2
    iz <- matrix(midpoints, nrow = 1, ncol = length(midpoints))
    if( verbose){print(breaks)
    	print( midpoints)
    	print( ix)
    	print( iy)
    	print( iz)
    	print( col)}
        #
    # next par call sets up a new plotting region just for the legend strip
    # at the smallplot coordinates
    par(new = TRUE, pty = "m", plt = smallplot, err = -1)
    # draw color scales the two  cases are horizontal/vertical
    # add a label if this is passed.
    if (!horizontal) {
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "",
                ylab = "", col = col, breaks=breaks)
    }
    else {
            image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "",
                ylab = "", col = col, breaks=breaks)
    }
    # create the argument list to draw the axis
    #  this avoids 4 separate calls to axis and allows passing extra
    # arguments.
    if (!is.null(lab.breaks)) {
        # axis with labels at break points
        axis.args <- c(list(side = ifelse(horizontal, 1, 4),
            mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2),
            at = breaks, labels = lab.breaks), axis.args)
    }
    else {
        # If lab.breaks is not specified ( with or without breaks), pretty
        # tick mark locations and labels are computed internally,
        # or as specified in axis.args at the function call
        axis.args <- c(list(side = ifelse(horizontal, 1, 4),
            mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)),
            axis.args)
    }
    #
    # now add the axis to the legend strip.
    # notice how all the information is in the list axis.args
    do.call("axis", axis.args)
    # add a box around legend strip
    box()
    #
    # add a label to the axis if information has been  supplied
    # using the mtext function. The arguments to mtext are
    # passed as a list like the drill for axis (see above)
    #
    if (!is.null(legend.lab)) {
        legend.args <- list(text = legend.lab, side = ifelse(horizontal,
            1, 4), line = legend.line, cex=legend.cex)
        #                    just guessing at a good default for line argument!
    }
    # add the label using mtext function
    if (!is.null(legend.args)) {
        do.call(mtext, legend.args)
    }
    #
    # clean up graphics device settings
    # reset to larger plot region with right user coordinates.
    mfg.save <- par()$mfg
    if (graphics.reset | add) {
        par(old.par)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
    else {
        par(big.par)
        par(plt = big.par$plt, xpd = FALSE)
        par(mfg = mfg.save, new = FALSE)
        # Suggestion from Karline Soetaert <Karline.Soetaert@nioz.nl>
        # this is to reset margins to be based on the mar arguments
  #      par(mar = par("mar"))  or
  #      par(mar = big.par$mar)
  # unfortunately this causes problems by allowing plotting outside of the
  # original plot region.
        invisible()
    }
}

#' @export
"tim.colors" <- function(n = 64, alpha = 1) {
  # tims original 64 color definition definition:
  orig <- c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF",
            "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF",
            "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF",
            "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF",
            "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF",
            "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F",
            "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50",
            "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00",
            "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00",
            "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000",
            "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000",
            "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000",
            "#AF0000", "#9F0000", "#8F0000", "#800000")
  if (n == 64 & alpha == 1)
    return(orig)
  rgb.tim <- t(col2rgb(orig))
  temp <- matrix(NA, ncol = 3, nrow = n)
  x <- seq(0, 1, , 64)
  xg <- seq(0, 1, , n)
  for (k in 1:3) {
    hold <- splint(x, rgb.tim[, k], xg)
    hold[hold < 0] <- 0
    hold[hold > 255] <- 255
    temp[, k] <- round(hold)
  }
  if (alpha == 1) {
    rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
  }
  else {
    rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255,
        alpha = alpha)
  }
}

#' @export
imagePlotInfo=function (..., breaks = NULL, nlevel)
{
  temp <- list(...)
  xlim <- NA
  ylim <- NA
  zlim <- NA
  poly.grid <- FALSE
  if (is.list(temp[[1]])) {
    xlim <- range(temp[[1]]$x, na.rm = TRUE)
    ylim <- range(temp[[1]]$y, na.rm = TRUE)
    zlim <- range(temp[[1]]$z, na.rm = TRUE)
    if (is.matrix(temp[[1]]$x) & is.matrix(temp[[1]]$y) &
        is.matrix(temp[[1]]$z)) {
      poly.grid <- TRUE
    }
  }
  if (length(temp) >= 3) {
    if (is.matrix(temp[[1]]) & is.matrix(temp[[2]]) & is.matrix(temp[[3]])) {
      poly.grid <- TRUE
    }
  }
  if (is.matrix(temp[[1]]) & !poly.grid) {
    xlim <- c(0, 1)
    ylim <- c(0, 1)
    zlim <- range(temp[[1]], na.rm = TRUE)
  }
  if (length(temp) >= 3) {
    if (is.matrix(temp[[3]])) {
      xlim <- range(temp[[1]], na.rm = TRUE)
      ylim <- range(temp[[2]], na.rm = TRUE)
      zlim <- range(temp[[3]], na.rm = TRUE)
    }
  }
  if (!is.na(zlim[1])) {
    if (zlim[1] == zlim[2]) {
      if (zlim[1] == 0) {
        zlim[1] <- -1e-08
        zlim[2] <- 1e-08
      }
      else {
        delta <- 0.01 * abs(zlim[1])
        zlim[1] <- zlim[1] - delta
        zlim[2] <- zlim[2] + delta
      }
    }
  }
  if (is.matrix(temp$x) & is.matrix(temp$y) & is.matrix(temp$z)) {
    poly.grid <- TRUE
  }
  xthere <- match("x", names(temp))
  ythere <- match("y", names(temp))
  zthere <- match("z", names(temp))
  if (!is.na(zthere))
    zlim <- range(temp$z, na.rm = TRUE)
  if (!is.na(xthere))
    xlim <- range(temp$x, na.rm = TRUE)
  if (!is.na(ythere))
    ylim <- range(temp$y, na.rm = TRUE)
  if (!is.null(temp$zlim))
    zlim <- temp$zlim
  if (!is.null(temp$xlim))
    xlim <- temp$xlim
  if (!is.null(temp$ylim))
    ylim <- temp$ylim
  if (is.null(breaks)) {
    midpoints <- seq(zlim[1], zlim[2], , nlevel)
    delta <- (midpoints[2] - midpoints[1])/2
    breaks <- c(midpoints[1] - delta, midpoints + delta)
  }
  list(xlim = xlim, ylim = ylim, zlim = zlim, poly.grid = poly.grid,
       breaks = breaks)
}

