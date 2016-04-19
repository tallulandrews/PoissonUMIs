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
