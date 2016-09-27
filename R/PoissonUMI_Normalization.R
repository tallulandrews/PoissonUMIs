#Copyright (c) 2016 Genome Research Ltd .
#Author : Tallulah Andrews <tallulandrews@gmail.com>
#This file is part of PoissonUMIs.

#PoissonUMIs is free software : you can redistribute it and/or modify it under
#the terms of the GNU General Public License as published by the Free Software
#Foundation; either version 2 of the License, or (at your option) any later
#version.

#This program is distributed in the hope that it will be useful, but WITHOUT
#ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

#You should have received a copy of the GNU General Public License along with
#this program . If not , see <http://www.gnu.org/licenses/>.


PoisUMI_Normalize_Data <- function(expr_mat, fit=NA) {
	if (is.na(fit)) {
		fit = PoisUMI_Fit_Full_Poisson(expr_mat)
	}
	return((expr_mat*fit$alpha-fit$lambdas)/sqrt(fit$lambdas))
#	return((expr_mat*fit$alpha-fit$lambdas))
}

# New Correlation idea: Fisher's method to combine p-values + individual p-values based on skellam

hidden_PoisUMI_getAUC <- function(gene, labels) {
        ranked=gene;
	ranked[gene>=0] = rank(gene[gene>=0])
	ranked[gene < 0] = 0;
        ms = aggregate(ranked~unlist(labels),FUN=mean); #Get average score for each cluster
        posgroup = as.character(unlist(ms[which(ms[,2]==max(ms[,2])),1])); #Get cluster with highest average score
        if (length(posgroup) > 1 | max(ms) < 0) {return (c(-1,-1,-1))} # Return negatives if there is a tie for cluster with highest average score (by definition this is not cluster specific), or if gene is detected in few cells in the cluster.

        # Create 1/0 vector of truths for predictions, cluster with highest average score vs everything else
        truth = labels == posgroup

        #Make predictions & get auc using RCOR package.
        pred=ROCR::prediction(ranked,as.numeric(truth))
        val = unlist(ROCR::performance(pred,"auc")@y.values)
        pval = wilcox.test(gene[truth],gene[!truth])$p.value
        if (!exists("pval")) {pval=NA}

        return(c(val,posgroup,pval))
}

PoisUMI_getmarkers <- function(expr_mat, labels) {
	if (length(labels) != length(expr_mat[1,])) {
		stop("Length of labels does not match number of cells.")
	}
        aucs = apply(expr_mat,1,hidden_PoisUMI_getAUC,labels=labels)
        auc_df <- data.frame(matrix(unlist(aucs), ncol=3, byrow=T))
        rownames(auc_df) = rownames(expr_mat)
        colnames(auc_df) = c("AUC","Group", "pval")
        auc_df[,1] = as.numeric(as.character(auc_df[,1]))
        auc_df[,3] = as.numeric(as.character(auc_df[,3]))
        auc_df = auc_df[auc_df[,1] > 0,]
	auc_df = auc_df[order(-auc_df$AUC),]
        return(auc_df);
}
