PoisUMI_plot_fit <- function(output, genes=NA) {
	if (!is.logical(genes)) {
		logi = names(output$s) %in% genes;
		if (sum(logi) < length(genes)) {
			warn(paste(length(genes)-sum(logi),"/",length(genes)," genes could not be matched to data.", sep=""))
		}
		genes = logi;
	}
	# Check input

	arrangement = order(output$s)
	xes = log(output$s)/log(10);
	plot(xes, output$p_obs, xlab="log10(expression)", ylab="Number of Dropouts");
	lines(xes[arrangement], output$p_exp[arrangement], col="green", lwd=3)
	if (sum(genes) > 0 | genes[1] != NA) {
		points(xes[genes], output$p_obs[genes], col="purple", pch=16)
	}
}
