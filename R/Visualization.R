PoisUMI_plot_fit <- function(output, genes=NA) {
	if (!is.logical(genes)) {
		logi = names(output$s) %in% genes;
		if (sum(logi) < length(genes)) {
			warning(paste(length(genes)-sum(logi),"/",length(genes)," genes could not be matched to data.", sep=""))
		}
		genes = logi;
	}
	# Check input
	if (sum(c("s","p_obs","p_exp") %in% names(output)) != 3) {stop("First argument does not contain required items: s, p_obs and p_exp!")}

	arrangement = order(output$s)
	xes = log(output$s)/log(10);
	plot(xes, output$p_obs, xlab="log10(expression)", ylab="Number of Dropouts");
	lines(xes[arrangement], output$p_exp[arrangement], col="green", lwd=3)
	if (genes[1] != NA & sum(genes) > 0) {
		points(xes[genes], output$p_obs[genes], col="purple", pch=16)
	}
}

PoisUMI_dropout_plot <- function(fit) {
	plot(fit$s, fit$p_obs/length(diradj[1,]), log="x", xlab="Expression",ylab="Dropout Rate")
	stddev = sqrt(fit$p_exp_var)
	arrows(fit$s,(fit$p_exp-stddev*2)/length(diradj[1,]), fit$s, (fit$p_exp+stddev*2)/length(diradj[1,]), len=0, col="forestgreen")
	points(fit$s, fit$p_exp/length(diradj[1,]), col="green", pch=16, cex=0.9)
	legend("topright", paste("alpha =", round(fit$alpha, digits=2)), bty="n")
}

PoisUMI_plot_distance<-function(distance_matrix, labels) {
	require("RColorBrewer")
	require("gplots")
	keep = rowSums(is.na(cell_cell_dist)) != length(cell_cell_dist[1,])-1
	heat_data = cell_cell_dist[keep,keep]
	Side_Colours=brewer.pal("Set1",n=length(unique(labels)))[labels[keep]]
	heatmap.2(heat_data, trace="n", col=rainbow(30),key.title="", key.xlab="Distance", RowSideColors=Side_Colours, ColSideColors=Side_Colours, hclustfun = function(x){hclust(x,method="ward.D2")})
}
