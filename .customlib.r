########################################################################
# .customlib.r
#
# Copyright (C) 2019-2020  Mahmoud M Ibrahim
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see http://www.gnu.org/licenses.
#
# Contact: mmibrahim@pm.me (Mahmoud M Ibrahim)
########################################################################


#custom functions
makeTransparent = function(..., alpha=0.5) {
	alpha = floor(255*alpha)  
	newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
	.makeTransparent = function(col, alpha) {
		rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
	}
	newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
	return(newColor)
}
color.gradient = function(x, colors=c("red","yellow","green"), colsteps=100) {
	colgrad = colorRampPalette(colors, bias = 10) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ]
	return(colgrad)
}
legend.col = function(col, lev){
	opar = par
	n = length(col)
	bx = par("usr")
	box.cx = c(bx[2] + (bx[2] - bx[1]) / 1000,
	bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
	box.cy = c(bx[3], bx[3])
	box.sy = (bx[4] - bx[3]) / n
	xx = rep(box.cx, each = 2)
	par(xpd = TRUE)
	for(i in 1:n){
		yy = c(box.cy[1] + (box.sy * (i - 1)),
		box.cy[1] + (box.sy * (i)),
		box.cy[1] + (box.sy * (i)),
		box.cy[1] + (box.sy * (i - 1)))
		polygon(xx, yy, col = col[i], border = col[i])
	}
	par(new = TRUE)
	plot(0, 0, type = "n",
	ylim = c(min(lev), max(lev)),
	yaxt = "n", ylab = "",
	xaxt = "n", xlab = "",
	frame.plot = FALSE)
	axis(side = 4, las = 2, tick = FALSE, line = .25)
	par	= opar
}
minor.ticks.axis = function(ax,n,t.ratio=0.5,mn,mx,...){
	lims = par("usr")
	if(ax %in%c(1,3)) lims = lims[1:2] else lims[3:4]

	major.ticks = pretty(lims,n=5)
	if(missing(mn)) mn = min(major.ticks)
	if(missing(mx)) mx = max(major.ticks)
	
	major.ticks = major.ticks[major.ticks >= mn & major.ticks <= mx]
	labels = sapply(major.ticks,function(i) as.expression(bquote(10^ .(i))))
	axis(ax,at=major.ticks,labels=labels,...)

	n = n+2
	minors = log10(pretty(10^major.ticks[1:2],n))-major.ticks[1]
	minors = minors[-c(1,n)]
	
	minor.ticks = c(outer(minors,major.ticks,`+`))
	minor.ticks = minor.ticks[minor.ticks > mn & minor.ticks < mx]
	
	axis(ax,at=minor.ticks,tcl=par("tcl")*t.ratio,labels=FALSE)
}
score = function(x){
	x  = ((x-min(x))/(max(x)-min(x)))
	return(x)
}
dist2d = function(a,b,c) {
	v1 = b - c
	v2 = a - b
	m = cbind(v1,v2)
	d = abs(det(m))/sqrt(sum(v1*v1))
	return(d)
}
ReadAlevin = function( base.path = NULL ){
    if (! dir.exists(base.path )){
      stop("Directory provided does not exist")
    }
	barcode.loc = paste0( base.path, "/quants_mat_rows.txt" )
	gene.loc = paste0( base.path, "/quants_mat_cols.txt" )
	matrix.loc = paste0( base.path, "/quants_mat.csv" )
	if (!file.exists( barcode.loc )){
		stop("Barcode file missing")
	}
	if (! file.exists(gene.loc) ){
		stop("Gene name file missing")
	}
	if (! file.exists(matrix.loc )){
		stop("Expression matrix file missing")
	}
	matrix <- as.matrix(read.csv( matrix.loc, header=FALSE))
	matrix <- t(matrix[,1:ncol(matrix)-1])

	cell.names <- readLines( barcode.loc )
	gene.names <- readLines( gene.loc )

	colnames(matrix) <- cell.names
	rownames(matrix) <- gene.names
	matrix[is.na(matrix)] <- 0
	return(matrix)
}


#general color palette
col_palette = c("#d64e3c","#7dd554","#8641c6","#cfca45","#7778cb","#59803d","#d04d9c","#73d6a8","#492f60","#ccc497","#7f343b","#72acc0","#b97d40","#c796b5","#45483a","purple","green","yellow")
col_palette_trans = makeTransparent(col_palette, alpha = 0.5)

col_palette_short = c("#999999","#56B4E9","#009E73","#0072B2")
col_palette_short_trans = makeTransparent(col_palette_short, alpha = 0.4)
