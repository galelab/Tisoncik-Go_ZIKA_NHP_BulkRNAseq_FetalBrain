##############-- Libraries--##############
options(mc.cores = 8)
library(SetRank)
library(GeneSets.Homo.sapiens)
library(stringr)
library(factoextra)
library(umap)
library(amap)
library(ggplot2)
library(mclust)
library(dendextend)
library(limma)
library(edgeR)
library(Rtsne)
library(corrplot)
library(fossil)
library(svglite)
library(data.table)
source("./heatmap3LW_function.r")


# TITLE: NHP Fetal-Maternal study with Jenny Go.  Previous analyses has been done by Dan Newhouse and Rich Green
# AUTHOR: LEANNE WHITMORE
###################### -- FILES/OBJECTS--######################

pigtail2human <- read.csv(
    file = "./pigtail2human.csv",
    header = TRUE,
    stringsAsFactors = FALSE
)

# When set to true code will run enrichment analysis (takes a bit of time so have added this option to easily turn off)
SetRankRun <- FALSE
if (isTRUE(SetRankRun)) {
    # converters
    symbol2EntrezID <- createIDConverter(
        "Homo.sapiens", "SYMBOL",
        "ENTREZID"
    )
    IDConverter <- createIDConverter(
        "Homo.sapiens", "ENTREZID",
        "SYMBOL"
    )
}

############# -- Functions --#############

theme_Publication <- function(base_size = 14, base_family = "arial") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size = base_size)
    + theme(
            plot.title = element_text(
                face = "bold",
                size = rel(1.2), hjust = 0.5
            ),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(1)),
            axis.title.y = element_text(angle = 90, vjust = 2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour = "#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size = unit(0.6, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face = "italic"),
            plot.margin = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(face = "bold")
        ))
}

generate_folder <- function(foldername) {
    workDir <- getwd()
    subDir <- foldername
    results_path <- file.path(workDir, subDir)
    if (file.exists(subDir)) {
    } else {
        dir.create(results_path)
    }
    return(results_path)
}

pca_fun <- function(exprs, labels, results_path,
                    base_file_name, target_columns,
                    figres = 100, size = 1, pca = FALSE) {
    # Run PCA/SVD reduction
    if (isFALSE(pca)) {
        pca <- prcomp(t(exprs))
    }
    E <- get_eig(pca)
    cx <- sweep(t(exprs), 2, colMeans(t(exprs)), "-")
    sv <- svd(cx)


    vizualize_pca(
        file.path(results_path, paste0("svd_", base_file_name)),
        sv$u, labels[, target_columns[1]],
        labels[, target_columns[2]], figres, E, size
    )
    vizualize_pca(
        file.path(results_path, paste0("pca_", base_file_name)),
        pca$x, labels[, target_columns[1]],
        labels[, target_columns[2]],
        figres, E, size
    )
    vizualize_scree_plot(
        file.path(
            results_path,
            paste0("scree_", base_file_name)
        ), pca, figres
    )

    loadingscores <- as.data.frame(pca$rotation)
    is_pc1_0 <- loadingscores$PC1 > 0
    is_pc2_0 <- loadingscores$PC2 > 0

    loadingscores <- loadingscores[is_pc1_0, ]
    loadingscores <- loadingscores[with(loadingscores, order(-PC1)), ]
    save_loading_scores(
        file.path(results_path, paste0("loadingscores_pc1", base_file_name, ".txt")),
        loadingscores["PC1"], figres
    )

    loadingscores <- as.data.frame(pca$rotation)
    loadingscores <- loadingscores[is_pc2_0, ]
    loadingscores <- loadingscores[with(loadingscores, order(-PC2)), ]
    save_loading_scores(
        file.path(results_path, paste0("loadingscores_pc2", base_file_name, ".txt")),
        loadingscores["PC2"], figres
    )
    return(pca)
}

umap_fun <- function(exprs, labels, results_path,
                     base_file_name, target_columns,
                     figres = 100, size = 1, UMAP = FALSE) {
    # Runs default paramaters of umap
    if (isFALSE(UMAP)) {
        UMAP <- umap(t(exprs))
    }
    vizualize_umap(
        file.path(results_path, paste0("umap_", base_file_name)),
        UMAP$layout, labels[, target_columns[1]],
        labels[, target_columns[2]], figres, size
    )

    return(UMAP)
}

vizualize_umap <- function(plot_file, U, class1, class2, figres, size) {
    # Vizualize umap reduction
    library(Polychrome)
    minx <- min(U[, 1])
    maxx <- max(U[, 1])
    miny <- min(U[, 2])
    maxy <- max(U[, 2])
    if (length(levels(factor(class2))) <= 3) {
        if (length(levels(factor(class1))) <= 6) {
            qplot(U[, 1], U[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                xlab("UMAP 1") +
                ylab("UMAP 2") +
                xlim(minx, maxx) + ylim(miny, maxy) +
                scale_color_manual(values = c("Prot" = "red", "NonProt" = "black", "Uninfected" = "gray")) +
                scale_fill_manual(values = c("Prot" = "red", "NonProt" = "black", "Uninfected" = "gray")) +
                theme(legend.position = "right")
        } else {
            qplot(U[, 1], U[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                xlab("UMAP 1") +
                ylab("UMAP 2") +
                xlim(minx, maxx) + ylim(miny, maxy) +
                scale_color_manual(values = c("Prot" = "red", "NonProt" = "black", "Uninfected" = "gray")) +
                scale_fill_manual(values = c("Prot" = "red", "NonProt" = "black", "Uninfected" = "gray")) +
                theme(legend.position = "right") +
                scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        }
    } else {
        P36 <- createPalette(length(levels(factor(class2))), c("#ff0000", "#00ff00", "#0000ff"))
        if (length(levels(factor(class1))) <= 6) {
            qplot(U[, 1], U[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                scale_color_manual(values = as.character(P36)) +
                scale_fill_manual(values = as.character(P36)) +
                xlab("UMAP 1") +
                ylab("UMAP 2") +
                xlim(minx, maxx) + ylim(miny, maxy) +
                theme(legend.position = "right")
        } else if (length(levels(factor(class1))) > 6) {
            qplot(U[, 1], U[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                xlab("UMAP 1") +
                ylab("UMAP 2") +
                scale_color_manual(values = as.character(P36)) +
                scale_fill_manual(values = as.character(P36)) +
                xlim(minx, maxx) + ylim(miny, maxy) +
                theme(legend.position = "right") +
                scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        }
    }
    ggsave(plot_file, width = 6, height = 4, units = "in", dpi = figres)
}

vizualize_pca <- function(plot_file, PCA, class1, class2, figres, E, size) {
    # Vizualize PCA  results
    library(Polychrome)
    minx <- min(PCA[, 1])
    maxx <- max(PCA[, 1])
    miny <- min(PCA[, 2])
    maxy <- max(PCA[, 2])
    if (length(levels(factor(class2))) <= 3) {
        if (length(levels(factor(class1))) <= 6) {
            qplot(PCA[, 1], PCA[, 2], color = factor(class2), shape = factor(class1), size = I(size)) +
                theme_Publication() +
                theme(legend.title = element_blank()) +
                scale_color_manual(values =c("ZIKV" = "black", "mock" = "gray")) +
                scale_fill_manual(values =c("ZIKV" = "black", "mock" = "gray")) +
                xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
                ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
                theme(legend.position = "right")
        } else {
            qplot(PCA[, 1], PCA[, 2], color = factor(class2), shape = factor(class1), size = I(size)) +
                theme_Publication() +
                theme(legend.title = element_blank()) +
                scale_color_manual(values = c("ZIKV" = "black", "mock" = "gray")) +
                scale_fill_manual(values = c("ZIKV" = "black", "mock" = "gray")) +
                xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
                ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
                theme(legend.position = "right") + scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        }
    } else {
        P36 <- createPalette(length(levels(factor(class2))), c("#ff0000", "#00ff00", "#0000ff"))
        if (length(levels(factor(class1))) <= 6) {
            qplot(PCA[, 1], PCA[, 2], color = factor(class2), shape = factor(class1), size = I(size)) +
                theme_Publication() +
                theme(legend.title = element_blank()) +
                xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
                ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
                theme(legend.position = "right") +
                scale_color_manual(values = as.character(P36)) +
                scale_fill_manual(values = as.character(P36))
        } else if (length(levels(factor(class1))) > 6) {
            qplot(PCA[, 1], PCA[, 2], color = factor(class2), shape = factor(class1), size = I(size)) +
                theme_Publication() +
                theme(legend.title = element_blank()) +
                xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
                ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
                theme(legend.position = "right") +
                scale_color_manual(values = as.character(P36)) +
                scale_fill_manual(values = as.character(P36)) +
                scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        }
    }
    ggsave(plot_file, width = 6, height = 4, units = "in", dpi = 500)
    plotfilesvg <- str_replace_all(plot_file, ".png", ".svg")
    ggsave(plotfilesvg, width = 5.5, height = 4, units = "in", dpi = 500)

}

tsne_fun <- function(exprs, labels, results_path,
                     base_file_name, target_columns, figres = 300, size = 1, T = FALSE) {
    # Runs default paramaters of umap
    if (isFALSE(T)) {
        T <- Rtsne(t(exprs), perplexity = 1)
    }
    vizualize_tSNE(
        file.path(results_path, paste0("tsne_", base_file_name)),
        T$Y, labels[, target_columns[1]],
        labels[, target_columns[2]], figres, size
    )
    return(T)
}

vizualize_tSNE <- function(plot_file, T, class1, class2, figres, size) {
    # Vizualize tsne reduction
    library(Polychrome)
    minx <- min(T[, 1])
    maxx <- max(T[, 1])
    miny <- min(T[, 2])
    maxy <- max(T[, 2])
    if (length(levels(factor(class2))) <= 3) {
        if (length(levels(factor(class1))) <= 6) {
            qplot(T[, 1], T[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                xlab("tSNE 1") +
                ylab("tSNE 2") +
                xlim(minx, maxx) + ylim(miny, maxy) +
                scale_color_manual(values = c("Prot" = "red", "NonProt" = "black", "Uninfected" = "gray")) +
                scale_fill_manual(values = c("Prot" = "red", "NonProt" = "black", "Uninfected" = "gray")) +
                theme(legend.position = "right")
        } else {
            qplot(T[, 1], T[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                xlab("tSNE 1") +
                ylab("tSNE 2") +
                xlim(minx, maxx) + ylim(miny, maxy) +
                scale_color_manual(values = c("Prot" = "red", "NonProt" = "black", "Uninfected" = "gray")) +
                scale_fill_manual(values = c("Prot" = "red", "NonProt" = "black", "Uninfected" = "gray")) +
                theme(legend.position = "right") + scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        }
    } else {
        P36 <- createPalette(length(levels(factor(class2))), c("#ff0000", "#00ff00", "#0000ff"))
        if (length(levels(factor(class1))) <= 6) {
            qplot(T[, 1], T[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                scale_color_manual(values = as.character(P36)) +
                scale_fill_manual(values = as.character(P36)) +
                xlab("tSNE 1") +
                ylab("tSNE 2") +
                xlim(minx, maxx) + ylim(miny, maxy) +
                theme(legend.position = "right")
        } else if (length(levels(factor(class1))) > 6) {
            qplot(T[, 1], T[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                xlab("tSNE 1") +
                ylab("tSNE 2") +
                scale_color_manual(values = as.character(P36)) +
                scale_fill_manual(values = as.character(P36)) +
                xlim(minx, maxx) + ylim(miny, maxy) +
                theme(legend.position = "right") +
                scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        }
    }
    ggsave(plot_file, width = 6, height = 4, units = "in", dpi = figres)
}

generate_boxplots_voom <- function(data, labels, filename, figres, maintitle, ylabtitle) {
    png(filename, width = 10, height = 8, units = "in", res = figres)
    # par(mar=c(1,1,1,1))
    minvalue <- min(data)
    maxvalue <- max(data)
    boxplot(data,
        labels = labels, ylim = c(minvalue - 1, maxvalue + 1),
        ylab = ylabtitle, main = maintitle, cex.axis = .6, las = 2,
        frame = FALSE
    )
    dev.off()
}

generate_density_plot <- function(data, labels, filename, figres) {
    png(filename, res = figres)
    par(xpd = TRUE)
    if (length(labels) > 10) {
        plotDensities(data, legend = FALSE)
    } else {
        plotDensities(data,
            legend = "topright",
            inset = c(-0.2, 0), levels(labels)
        )
    }
    dev.off()
}

normalize_data <- function(CM2, targetfile, geneTrans, species="rhesus") {

    # order target and count matrix so they are the same (THIS IS IMPORTANT)
    CM2 <- CM2[, rownames(targetfile)]

    # CHECK IF ORDER IS THE SAME
    if (all.equal(colnames(CM2), rownames(targetfile)) != TRUE) {
        print("MASSIVE WARNING: RESULTS WILL BE WRONG IF THIS IS NOT EQUAL!!!!!!!!")
        print(rownames(targetfile))
        print(colnames(CM2))
    }

    # design matrix without outcome
    Xid <- factor(targetfile[, "Animal.ID"])

    mm <- model.matrix(~ 0 + Xid)

    rownames(mm) <- colnames(CM2)
    colnames(mm) <- make.names(colnames(mm))
    mm <- mm[, colnames(mm)[order(tolower(colnames(mm[, ])))]]
    mm <- mm[, colSums(mm) > 0]
    excludeO <- nonEstimable(mm)
    if ("ti" %in% excludeO) {
        return("interactions term non estimable")
    }
    mm <- mm[, !colnames(mm) %in% excludeO]
    if (!is.fullrank(mm)) {
        return("not full rank")
    }
    # normalize
    CM2 <- DGEList(counts = CM2)
    CM2 <- calcNormFactors(CM2, method = "TMM") # TMM normalization
    png(file.path(norm_results, paste0("mean_variance_norm_", species, ".png")))
    Pi.CPM <- voom(counts = CM2, design = mm, normalize.method = "none", plot = T, span = 0.1)
    dev.off()
    write.csv(Pi.CPM$E, file.path(norm_results, paste0("1.norm_matrix_", species, ".csv")))
    # message("STATUS: getting the corfit")
    # corfit <- duplicateCorrelation(CM2$counts,
    #     block = factor(targetfile$Animal.ID)
    # ) # account for repeated sampling of individuals

    # message("STATUS: renormalizing with corfit")
    # png(file.path(norm_results, "mean_variance_norm_corfit.png"))
    # Pi.CPM <- voom(
    #     counts = CM2, normalize.method = "none",
    #     correlation = corfit$consensus.correlation,
    #     plot = T, span = 0.1, save.plot = T
    # )
    # dev.off()

    # message("STATUS: recalculating corfit")
    # corfit <- duplicateCorrelation(Pi.CPM,
    #     block = factor(targetfile$Animal.ID)
    # )
    # write.csv(Pi.CPM$E, file.path(norm_results, paste0("1.norm_matrix_", species, "_corfit.csv")))


    sig_HGNC <- merge(geneTrans, Pi.CPM$E,
        by.x = "Gene.stable.ID",
        by.y = "row.names",
        all.X = T, all.Y = T
    )

    sig_HGNC <- sig_HGNC[, !(names(sig_HGNC) %in% c("Gene.stable.ID"))]
    sig_HGNC <- avereps(sig_HGNC,
        ID = sig_HGNC$HGNC.symbol
    )
    rownames(sig_HGNC) <- sig_HGNC[, "HGNC.symbol"]
    sig_HGNC <- sig_HGNC[, !(colnames(sig_HGNC) %in% c("HGNC.symbol"))]
    sig_HGNC <- as.matrix(data.frame(sig_HGNC))
    class(sig_HGNC) <- "numeric"
    write.csv(sig_HGNC, file.path(norm_results, paste0("1.norm_matrix_HGNC_", species, ".csv")), quote = FALSE)
    write.table(sig_HGNC, file.path(norm_results, paste0("1.norm_matrix_HGNC_", species, ".txt")), sep="\t", quote = FALSE)
    return(list("norm"=Pi.CPM))#, "corfit"=corfit))
}

generate_clusterdendograms <- function(hc, plotfilename1, adjvalue, labelsize = 0.7) {
    library(Polychrome)

    counter <- 0
    labelsf <- c()
    colors <- c()
    dend <- as.dendrogram(hc)
    dend_labels <- labels(dend)
    P36 <- createPalette(length(levels(factor(dend_labels))), c("#ff0000", "#00ff00", "#0000ff"))
    names(P36) <- unique(dend_labels)
    for (i in dend_labels) {
        labelsf <- c(labelsf, i)
        colors <- c(colors, P36[[i]])
    }
    labels_colors(dend) <- colors
    labels_cex(dend) <- labelsize
    png(plotfilename1,
        units = "in", # bg = "transparent",
        width = 14.5, height = 5, res = 300
    )
    par(mar = c(6, 3, 2, 0.5), xpd = TRUE)
    plot(dend, xlab = "", main = "")
    mtext(paste0("Adj Rand index ", round(adjvalue, 3)))
    dev.off()
}

filter_read_counts_mean <- function(cm, filter_cutoff) {
    # Filter value was calculated by:
    #### Filters by row means usually set at 10 reads per gene across all samples

    A <- rowMeans(cm)
    isexpr <- A >= filter_cutoff
    cmfl <- cm[isexpr, ]
    return(cmfl)
}

vizualize_scree_plot <- function(plot_file, PCA, figres) {
    # Vizualize principle component variation results
    scree.plot <- fviz_eig(PCA, addlabels = TRUE, hjust = -0.3)
    png(plot_file, width = 7, height = 6, units = "in", res = figres)
    print(scree.plot)
    dev.off()
}

save_loading_scores <- function(write_file, df, figres) {
    # Save list of genes that have a positive effect on variation of principle
    # component 1 and 2 sorted from most influential
    write.table(df, file = write_file)
}

generate_design_matrix <- function(normmatrix, target) {
    TR <- factor(target$Strain)
    Xid <- factor(target$Animal.ID)
    BR <- factor(target$Brain.Region)

    mm <- model.matrix(~ 0 + TR:BR)
    # mm <- model.matrix(~ 0 + TR:BR + Xid)
    rownames(mm) <- colnames(normmatrix)
    colnames(mm) <- make.names(colnames(mm))
    mm <- mm[, colnames(mm)[order(tolower(colnames(mm[, ])))]]
    mm <- mm[, colSums(mm) > 0]

    excludeAll <- nonEstimable(mm)
    if (length(excludeAll) > 0) {
        message("WARNING: These samples are nonEstimatable, design matrix ", excludeAll)
    }

    if ("ti" %in% excludeAll) {
        return("interactions term non estimable")
    }
    mm <- mm[, !colnames(mm) %in% excludeAll]
    if (!is.fullrank(mm)) {
        return("not full rank")
    }
    return(mm)
}

generate_design_matrix2 <- function(normmatrix, target) {
    ZV <- factor(target$zika)
    Xid <- factor(target$Animal.ID)
    BR <- factor(target$Brain.Region)

    mm <- model.matrix(~ 0 + ZV:BR)
    rownames(mm) <- colnames(normmatrix)
    colnames(mm) <- make.names(colnames(mm))
    mm <- mm[, colnames(mm)[order(tolower(colnames(mm[, ])))]]
    mm <- mm[, colSums(mm) > 0]

    excludeAll <- nonEstimable(mm)
    if (length(excludeAll) > 0) {
        message("WARNING: These samples are nonEstimatable, design matrix ", excludeAll)
    }

    if ("ti" %in% excludeAll) {
        return("interactions term non estimable")
    }
    mm <- mm[, !colnames(mm) %in% excludeAll]
    if (!is.fullrank(mm)) {
        return("not full rank")
    }
    return(mm)
}

generate_design_matrix3 <- function(normmatrix, target) {
    TR <- factor(target$Strain)
    Xid <- factor(target$Animal.ID)
    BR <- factor(target$Brain.Region)

    # mm <- model.matrix(~ 0 + TR:BR)
    mm <- model.matrix(~ 0 + TR:BR + Xid)
    rownames(mm) <- colnames(normmatrix)
    colnames(mm) <- make.names(colnames(mm))
    mm <- mm[, colnames(mm)[order(tolower(colnames(mm[, ])))]]
    mm <- mm[, colSums(mm) > 0]

    excludeAll <- nonEstimable(mm)
    if (length(excludeAll) > 0) {
        message("WARNING: These samples are nonEstimatable, design matrix ", excludeAll)
    }

    if ("ti" %in% excludeAll) {
        return("interactions term non estimable")
    }
    mm <- mm[, !colnames(mm) %in% excludeAll]
    if (!is.fullrank(mm)) {
        return("not full rank")
    }
    return(mm)
}

vizualize_DE_genes_bp <- function(results, plot_file) {
    print("STATUS: Generating bar plot of number of DE genes...")
    results_t <- t(summary(results))
    results_t <- results_t[, -2]

    for (i in 1:(length(row.names(results_t)))) {
        results_t[i, 1] <- results_t[i, 1] * -1
    }

    DE <- as.data.frame(results_t)
    DE <- setnames(DE,
        old = c("Var1", "Var2", "Freq"),
        new = c("Time_Point", "group", "DE_genes")
    )

    # Create plot
    ggplot(DE, aes(
        x = Time_Point, y = DE_genes, fill = group,
        label = DE$DE_genes
    )) +
        geom_bar(stat = "identity", position = "identity") +
        # geom_text(size = 5, position = position_stack(vjust = 0) )+
        # theme_light() +
        theme_minimal() +
        scale_fill_manual(values = c("#0808c4", "#da9618")) +
        # xlab("Time point")
        ylab("Number of Differentially Expressed Genes") +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
            axis.text.y = element_text(size = 15)
        )
    ggsave(plot_file, width = 6, height = 4, units = "in", dpi = 300)
}

gene_enrichment <- function(genes, results_folder, cluster) {
    inputGenesTrans <- pigtail2human[pigtail2human$Gene.stable.ID %in% genes, ]
    inputGenesHGNC <- unique(unlist(inputGenesTrans$HGNC.symbol))
    inputGenes <- symbol2EntrezID(inputGenesHGNC)
    network <- setRankAnalysis(inputGenes, collection,
        use.ranks = FALSE,
        setPCutoff = 0.01,
        fdrCutoff = 0.05
    )

    generate_folder(results_folder)
    #### IMPORTANT OUTPUT INFORMATION###
    # SetRank value -  value reflects the prominence of a gene set in a gene set network (based on PageRank algorithm developed by google)
    # pSetRank value - significance of the SetRank value
    exportSingleResult(network, inputGenes,
        collection, paste0(results_folder, "/de_unranked_", cluster),
        IDConverter = IDConverter
    )
    # png(file.path(results_folder, paste0("de_unranked_", cluster,".png")), res=100)
    # plot(network, layout = layout.spring)
    # dev.off()
    return(network)
}

convertensembl2hgnc <- function(df) {
    df_HGNC <- merge(pigtail2human, df,
        by.x = "Gene.stable.ID",
        by.y = "row.names",
        all.X = T, all.Y = T
    )

    df_HGNC <- df_HGNC[, !(names(df_HGNC) %in% c("Gene.stable.ID"))]
    df_HGNC <- avereps(df_HGNC,
        ID = df_HGNC$HGNC.symbol
    )
    rownames(df_HGNC) <- df_HGNC[, "HGNC.symbol"]
    df_HGNC <- df_HGNC[, !(colnames(df_HGNC) %in% c("HGNC.symbol"))]
    df_HGNC <- as.matrix(data.frame(df_HGNC, check.names = FALSE))
    class(df_HGNC) <- "numeric"
    return(df_HGNC)
}

convertensembl2hgncsupp <- function(df) {
    df_HGNC <- merge(pigtail2human, df,
        by.x = "Gene.stable.ID",
        by.y = "row.names",
        all.X = T, all.Y = F
    )

    # df_HGNC <- df_HGNC[, !(names(df_HGNC) %in% c("Gene.stable.ID"))]
    # df_HGNC <- avereps(df_HGNC,
    #     ID = df_HGNC$Gene.stable.ID
    # )
    # rownames(df_HGNC) <- df_HGNC[, "Gene.stable.ID"]
    # df_HGNC <- df_HGNC[, !(colnames(df_HGNC) %in% c("Gene.stable.ID"))]
    # df_HGNC <- as.matrix(data.frame(df_HGNC, check.names = FALSE))
    # class(df_HGNC) <- "numeric"
    return(df_HGNC)
}
############# --Load Data --#############
message("STATUS: Loading data ...")
## -- Load target --##
target <- read.csv("metadata.csv", header = T, row.names = 1)
target <- target[target$Sequencing.Type != "Virus", ]
targetFetBrain <- target[target$Maternal.or.Fetal == "Fet" & target$Tissue == "Brain", ]

## --Pigtail alignment-- ##
countpigtail <- read.table("count_matrix.txt", header = T, check.names = FALSE, row.names = 1)
colnames(countpigtail) <- str_remove_all(colnames(countpigtail), "_D\\d+_\\S+")
colnames(countpigtail) <- str_remove_all(colnames(countpigtail), "_RNA_\\S+")
colnames(countpigtail) <- str_remove_all(colnames(countpigtail), "_run01_\\S+")

countpigtail <- countpigtail[, rownames(targetFetBrain)]
colnames(countpigtail) <- paste(targetFetBrain$Animal.ID, targetFetBrain$Brain.Region, sep = "_")

############# --Processing count data aligned --#############
message("STATUS: Processing count data ...")

count_results <- "1.count_data/"
generate_folder(count_results)

## -- Pigtail

targetFetBrain$SampleName <- rownames(targetFetBrain)
rownames(targetFetBrain) <- paste(targetFetBrain$Animal.ID, targetFetBrain$Brain.Region, sep = "_")
write.csv(targetFetBrain, file.path(count_results, "target_renamed.csv"))

write.csv(countpigtail, file.path(count_results, "count_matrix_pigtail_renamed.csv"))
generate_density_plot(
    countpigtail, rownames(targetFetBrain), file.path(count_results, "de_intensities_pigtail_raw_counts.png"), 100
)
generate_boxplots_voom(countpigtail, colnames(countpigtail),
    file.path(count_results, "boxplot_pigtail_counts.png"),
    100,
    maintitle = "count matrix",
    ylabtitle = "counts"
)

## -- Rhesus v10

write.csv(countrhesus, file.path(count_results, "count_matrix_rhesus_renamed.csv"))
generate_density_plot(
    countrhesus, rownames(targetFetBrain), file.path(count_results, "de_intensities_rhesus_raw_counts.png"), 100
)
generate_boxplots_voom(countrhesus, colnames(countrhesus),
    file.path(count_results, "boxplot_rhesus_counts.png"),
    100,
    maintitle = "count matrix",
    ylabtitle = "counts"
)

############# --Normalize data--#############
message("STATUS: Normalizing data ...")

## -- Pigtail --## 
norm_results <- "1.norm_data/"
generate_folder(norm_results)
pigtail <- normalize_data(countpigtail, targetFetBrain, pigtail2human, species="pigtail")
Pi.CPMP <- pigtail$norm
# corfitP <- pigtail$corfit
saveRDS(Pi.CPMP, file.path(norm_results,"normobject_pigtail.rds"))
# saveRDS(corfitP, file.path(norm_results,"corfitobject_pigtail.rds"))
generate_density_plot(
    Pi.CPMP$E, rownames(targetFetBrain), file.path(norm_results, "de_intensities_pigtail_raw_counts.png"), 100
)

generate_boxplots_voom(Pi.CPMP$E, colnames(Pi.CPMP$E),
    file.path(norm_results, "boxplot_vnorm_pigtail.png"),
    100,
    maintitle = "Normalized count matrix",
    ylabtitle = "voom normalized expression"
)

## Filter low expressing gene (hook of mean variance trend )
countpigtailfl <- filter_read_counts_mean(countpigtail, 3)
write.csv(countpigtailfl, file.path(count_results, "count_matrix_renamed_fl_pigtail.csv"))
norm_results <- "1.norm_data_fl/"
generate_folder(norm_results)
pigtail <- normalize_data(countpigtailfl, targetFetBrain, pigtail2human, species = "pigtail")
Pi.CPMP <- pigtail$norm
# corfitP <- pigtail$corfit
saveRDS(Pi.CPMP, file.path(norm_results, "normobject_pigtail.rds"))
# saveRDS(corfitP, file.path(norm_results, "corfitobject_pigtail.rds"))

generate_density_plot(
    Pi.CPMP$E, rownames(targetFetBrain), file.path(norm_results, "de_intensities_pigtail_raw_counts.png"), 100
)

generate_boxplots_voom(Pi.CPMP$E, colnames(Pi.CPMP$E),
    file.path(norm_results, "boxplot_vnorm_pigtail.png"),
    100,
    maintitle = "Normalized count matrix",
    ylabtitle = "voom normalized expression"
)


############# --Feature Reduction --#############

targetFetBrain$Brain.Region2 <- targetFetBrain$Brain.Region
targetFetBrain$Brain.Region2 <- str_replace_all(targetFetBrain$Brain.Region2, "Br2", "BR1")
targetFetBrain$Brain.Region2 <- str_replace_all(targetFetBrain$Brain.Region2, "Br4", "BR3")
targetFetBrain$Brain.Region2 <- str_replace_all(targetFetBrain$Brain.Region2, "Br6", "BR5")
targetFetBrain$Brain.Region2 <- str_replace_all(targetFetBrain$Brain.Region2, "Br8", "BR7")
targetFetBrain$Brain.Region2 <- str_replace_all(targetFetBrain$Brain.Region2, "Br10", "BR9")

message("STATUS: Run Feature reduction")
feature_results <- "1.feature_red_without_5267"
generate_folder(feature_results)

targetFetBraintmp <- targetFetBrain[targetFetBrain$Animal.ID!="5627", ]
normtmp <- Pi.CPMP$E[, rownames(targetFetBraintmp)]

pca <- pca_fun(
    normtmp, targetFetBraintmp,
    feature_results, "_pigtail_norm.png",
    c("Strain", "Animal.ID"), 300, 4
)

pca <- pca_fun(
    normtmp, targetFetBraintmp,
    feature_results, "_pigtail_normBR.png",
    c("Strain", "Brain.Region"), 300, 4, pca=pca
)

pca <- pca_fun(
    normtmp, targetFetBraintmp,
    feature_results, "_pigtail_normBRZIKV.png",
    c("Treatment", "Brain.Region"), 300, 4,
    pca = pca
)
pca <- pca_fun(
    normtmp, targetFetBraintmp,
    feature_results, "_pigtail_normBRZIKVflipped.png",
    c("Brain.Region", "Treatment"), 300, 4,
    pca = pca
)
message("STATUS: Run Feature reduction")
feature_results <- "1.feature_red"
generate_folder(feature_results)

## --Pigtail
pca <- pca_fun(
    Pi.CPMP$E, targetFetBrain,
    feature_results, "_pigtail_norm.png",
    c("Strain", "Animal.ID"), 300, 3
)
pca <- pca_fun(
    Pi.CPMP$E, targetFetBrain,
    feature_results, "_pigtail_normBrainRegion.png",
    c("Strain", "Brain.Region"), 300, 3
)
pca <- pca_fun(
    Pi.CPMP$E, targetFetBrain,
    feature_results, "_pigtail_AnimalIDrainRegion.png",
    c("Brain.Region2", "Animal.ID"), 300, 3
)

umap <- umap_fun(
    Pi.CPMP$E, targetFetBrain,
    feature_results, "_pigtail_norm.png",
   c("Strain", "Animal.ID"), 300, 2
)

umap <- umap_fun(
    Pi.CPMP$E, targetFetBrain,
    feature_results, "_pigtail_normBrainRegion.png",
    c("Strain", "Brain.Region"), 300, 2
)

tsne <- tsne_fun(
    Pi.CPMP$E, targetFetBrain,
    feature_results, "_pigtail_norm.png",
  c("Strain", "Animal.ID"), 300, 3,
)

tsne <- tsne_fun(
    Pi.CPMP$E, targetFetBrain,
    feature_results, "_pigtail_normBrainRegion.png",
     c("Strain", "Brain.Region"), 300, 3,
    T = tsne
)

############# --Hierarchal clustering--#############

message("STATUS: Running hierarchal clustering...")
hi_cluster_results <- "1.hi_results"
generate_folder(hi_cluster_results)

## -- Pigtail
d <- Pi.CPMP$E
## Cluster plot with brain region labels 
colnames(d) <- targetFetBrain$Brain.Region
hc <- hcluster(t(d), method = "pearson", link = "average")
x <- cutree(hc, k = length(unique(targetFetBrain$Brain.Region)))
m <- adj.rand.index(as.integer(x), as.integer(factor(as.character(targetFetBrain$Brain.Region))))
generate_clusterdendograms(hc,
    file.path(hi_cluster_results, "dendogram_pigtail_brainRegion.png"), m, labelsize=1)

## animal ID
colnames(d) <- targetFetBrain$Animal.ID
hc <- hcluster(t(d), method = "pearson", link = "average")
x <- cutree(hc, k = length(unique(targetFetBrain$Animal.ID)))
m <- adj.rand.index(as.integer(x), as.integer(factor(as.character(targetFetBrain$Animal.ID))))
generate_clusterdendograms(
    hc,
    file.path(hi_cluster_results, "dendogram_pigtail_AnimalID.png"), m, labelsize=1
)

## animal ID + brain region
colnames(d) <- rownames(targetFetBrain)
hc <- hcluster(t(d), method = "pearson", link = "average")
x <- cutree(hc, k = length(unique(targetFetBrain$Animal.ID)))
m <- adj.rand.index(as.integer(x), as.integer(factor(as.character(rownames(targetFetBrain)))))
generate_clusterdendograms(
    hc,
    file.path(hi_cluster_results, "dendogram_pigtail_AnimalID_BR.png"), m,
    labelsize = 1
)

############# -- DE Analysis without animal 5267 Edelfsen --#############
message("STATUS: DE without animal 5267")
deresults_path <- "1.de"
generate_folder(deresults_path)

mm_all <- generate_design_matrix3(Pi.CPMP$E, targetFetBrain)

Pi.lm <- lapply(rownames(Pi.CPMP$E), function(.gene) { ### Here Pi.CPMO is normalized, with weights assigned by voom()
    # print(.gene)
    X <- Pi.CPMP$E[.gene, ]
    W <- Pi.CPMP$weights[rownames(Pi.CPMP$E) == .gene, ]
    lm(X ~ 0 + mm_all, weights = W) ### Okay, this results in the same coefficients as the corresponding row of Pi.lmfit:
    # unname( Pi.lm$coefficients == Pi.lmfit$coefficients[.gene, ] )
    #
    # ### Calculate coefficients and vcov for contrasts
    # .coefs <- Pi.lm$coefficients %*% contr
    # .vcov <- t(contr) %*% vcov(Pi.lm) %*% contr
    # .stdev <- sqrt(diag(.vcov))
    # .tstats <- .coefs/.stdev
    # .pt <- pt(abs(.tstats), lower.tail = F, df = Pi.lm$df.residual) * 2 #2-sided t test
    # array( c(.coefs,.stdev,.tstats,.pt), dim = c(length(.coefs),4), dimnames = list( colnames(contr), c("coef","stdev","t","p") ))
})
names(Pi.lm) <- rownames(Pi.CPMP$E)

contrastsmatrix <- c(
    "TRFSS.BRBR1 -TRmock.BRBR1",
    "TRFSS.BRBR3 -TRmock.BRBR3",
    "TRFSS.BRBR5 -TRmock.BRBR5",
    "TRFSS.BRBR7 -TRmock.BRBR7",
    "TRFSS.BRBR9 -TRmock.BRBR9",
    "TRBrazil.BRBR1 -TRmock.BRBR1",
    "TRBrazil.BRBR3 -TRmock.BRBR3",
    "TRBrazil.BRBR5 -TRmock.BRBR5",
    "TRBrazil.BRBR7 -TRmock.BRBR7",
    "TRBrazil.BRBR9 -TRmock.BRBR9"
)

contrasts.fit.lm <- function(fit, contrasts) {
    out <- sapply(fit, function(.lm) {
        ### Calculate coefficients and vcov for contrasts
        .coefs <- .lm$coefficients %*% contrasts
        .vcov <- t(contrasts) %*% vcov(.lm) %*% contrasts
        .stdev <- sqrt(diag(.vcov))
        .tstats <- .coefs / .stdev
        .pt <- pt(abs(.tstats), lower.tail = F, df = .lm$df.residual) * 2 # 2-sided t test
        array(c(.coefs, .stdev, .tstats, .pt), dim = c(length(.coefs), 4), dimnames = list(colnames(contr), c("coef", "stdev", "t", "p")))
    }, simplify = "array")
}
contr <- makeContrasts(contrasts = contrastsmatrix, levels = mm_all)

Pi.contrasts.lm <- contrasts.fit.lm(Pi.lm, contr) # This results in a 3-dimensional array, which is not entirely analogous to Pi.contrasts from above. The next two lines extract the coefficients and p-value matrices from the array, which is all you need to use the decideTests() function.
.coefs <- t(Pi.contrasts.lm[, "coef", ])
.tstats <- t(Pi.contrasts.lm[, "t", ])
.pt <- t(Pi.contrasts.lm[, "p", ])
.pt.adjusted <- apply(.pt, 2, p.adjust, method = "BH")

results.lmde <- decideTests(object = .pt, coefficients = .coefs, lfc = (.58),
            method = "separate", adjust.method = "BH", p.value = 0.05)

write.csv(.coefs, file = file.path(deresults_path, "coefficients.csv"), quote = F)
write.csv(.tstats, file = file.path(deresults_path, "t_stats.csv"), quote = F)
write.csv(.pt, file = file.path(deresults_path, "p_value.csv"), quote = F)
write.csv(.pt.adjusted, file = file.path(deresults_path, "p_value_adj.csv"), quote = F)

dataMatrixde <- .coefs
sigMask <- dataMatrixde * (results.lmde**2) # 1 if significant, 0 otherwise
ExpressMatrixde <- subset(dataMatrixde, rowSums(sigMask) != 0)

write.csv(ExpressMatrixde, file = file.path(deresults_path, "expression_matrix_de_lmtmp.csv"), quote = F)
write.csv(results.lmde, file = file.path(deresults_path, "results_de_lmtmp.csv"), quote = F)
write.csv(dataMatrixde, file = file.path(deresults_path, "full_expression_matrix_de_lmtmp.csv"), quote = F)

ExpressMatrixde_HGNC <- convertensembl2hgnc(ExpressMatrixde)
write.csv(ExpressMatrixde_HGNC, file.path(deresults_path, "expression_matrix_de_lm_hgnc.csv"), quote = F)

dataMatrixde_HGNC <- convertensembl2hgnc(dataMatrixde)
write.csv(dataMatrixde_HGNC, file.path(deresults_path, "full_expression_matrix_de_lm_hgnc.csv"), quote = F)

message(paste0("Dimensionality of DE genes ", dim(ExpressMatrixde)[1]))

colcolorlistgroup <- c(rep("#0b4e1a", 5), rep("#530d7e", 5))
colcolormatrix <- as.matrix(colcolorlistgroup)
colnames(colcolormatrix) <- c("Group")

# --barplot
vizualize_DE_genes_bp(results.lmde, file.path(deresults_path, "barplot.png"))

new_colnames <- rep(c("BR1", "BR3", "BR5", "BR7", "BR9"), 2)
# --heatmap
# png(file.path(deresults_path, "heatmap.png"), width = 8, height = 10, units = "in", res = 300)
svglite(file.path(deresults_path, "heatmap.svg"))
# par(mar = c(4, 4, -1, 2))
global_modulesde <- heatmap.L.4(ExpressMatrixde,
    figmargins = c(7, 5),
    cutoff = 1, distmethod = "euclidean", cexcol = 2, labCol = new_colnames,
    clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9,
    colsep = c(5), colcolorlist = colcolormatrix
)
dev.off()
if (isTRUE(all.equal(names(global_modulesde$modulesrows), rownames(global_modulesde$clustermatrix)))) {
    supplementaloutputavg <- global_modulesde$clustermatrix
    supplementaloutputavg <- convertensembl2hgncsupp(supplementaloutputavg)
    cluster <- as.data.frame(global_modulesde$modulesrows)
    supplementaloutputavg <- merge(supplementaloutputavg, cluster,
        by.x = "Gene.stable.ID", by.y = "row.names", all.x = TRUE, all.y = FALSE
    )
    colnames(supplementaloutputavg) <- c("ensembl.gene.id", "HGNC.symbol", "FSSlfcBR1", "FSSlfcBR3", "FSSlfcBR5", "FSSlfcBR7", "FSSlfcBR9", 
    "BrazillfcBR1", "BrazillfcBR3", "BrazillfcBR5", "BrazillfcBR7", "BrazillfcBR9", "cluster")
    pval <- .pt[supplementaloutputavg$ensembl.gene.id, ]
    colnames(pval) <- c(
        "Pval_FSSBR1-MockBR1", "Pval_FSSBR3-MockBR3", "Pval_FSSBR5-MockBR5", "Pval_FSSBR7-MockBR7", "Pval_FSSBR9-MockBR9",
        "Pval_BrazilBR1-MockBR1", "Pval_BrazilBR3-MockBR3", "Pval_BrazilBR5-MockBR5", "Pval_BrazilBR7-MockBR7", "Pval_BrazilBR9-MockBR9"
    )
    padj <- .pt.adjusted[supplementaloutputavg$ensembl.gene.id, ]
    colnames(padj) <- c(
        "Pvaladj_FSSBR1-MockBR1", "Pvaladj_FSSBR3-MockBR3", "Pvaladj_FSSBR5-MockBR5", "Pvaladj_FSSBR7-MockBR7", "Pvaladj_FSSBR9-MockBR9",
        "Pvaladj_BrazilBR1-MockBR1", "Pvaladj_BrazilBR3-MockBR3", "Pvaladj_BrazilBR5-MockBR5", "Pvaladj_BrazilBR7-MockBR7", "Pvaladj_BrazilBR9-MockBR9"
    )
    supplementaloutputavg <- merge(supplementaloutputavg, pval,
        by.x = "ensembl.gene.id", by.y = "row.names", all.x = TRUE, all.y = FALSE
    )
    supplementaloutputavg <- merge(supplementaloutputavg, padj,
        by.x = "ensembl.gene.id", by.y = "row.names", all.x = TRUE, all.y = FALSE
    )
    write.csv(supplementaloutputavg, file.path(deresults_path, "supplde.csv"))
} else {
    stop("WARNING: gene lists not in the same order")
}

colcolorlistgroup <- c(rep("white", 5))
colcolormatrix <- as.matrix(colcolorlistgroup)
colnames(colcolormatrix) <- c("Group")
avgBR1 <- rowMeans(ExpressMatrixde[, c("TRFSS.BRBR1 -TRmock.BRBR1", "TRBrazil.BRBR1 -TRmock.BRBR1")])
avgBR3 <- rowMeans(ExpressMatrixde[, c("TRFSS.BRBR3 -TRmock.BRBR3", "TRBrazil.BRBR3 -TRmock.BRBR3")])
avgBR5 <- rowMeans(ExpressMatrixde[, c("TRFSS.BRBR5 -TRmock.BRBR5", "TRBrazil.BRBR5 -TRmock.BRBR5")])
avgBR7 <- rowMeans(ExpressMatrixde[, c("TRFSS.BRBR7 -TRmock.BRBR7", "TRBrazil.BRBR7 -TRmock.BRBR7")])
avgBR9 <- rowMeans(ExpressMatrixde[, c("TRFSS.BRBR9 -TRmock.BRBR9", "TRBrazil.BRBR9 -TRmock.BRBR9")])
new_colnames <- rep(c("BR1", "BR3", "BR5", "BR7", "BR9"), 1)
df <- cbind(avgBR1, avgBR3, avgBR5, avgBR7, avgBR9)
png(file.path(deresults_path, "heatmap_avg.png"), width = 8, height = 10, units = "in", res = 300)
svglite(file.path(deresults_path, "heatmap_avg.svg"))
# par(mar = c(4, 4, -1, 2))
global_modulesdeavg <- heatmap.L.4(df,
    figmargins = c(7, 5),
    cutoff = 2, distmethod = "euclidean", cexcol = 2, labCol = new_colnames,
    clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9,
   colcolorlist = colcolormatrix
)
dev.off()
if (isTRUE(all.equal(names(global_modulesdeavg$modulesrows), rownames(global_modulesdeavg$clustermatrix)))) {
    supplementaloutputavg <- global_modulesdeavg$clustermatrix
    supplementaloutputavg <- convertensembl2hgncsupp(supplementaloutputavg)
    cluster <- as.data.frame(global_modulesdeavg$modulesrows)
    supplementaloutputavg <- merge(supplementaloutputavg, cluster,
             by.x="Gene.stable.ID", by.y="row.names", all.x=TRUE, all.y=FALSE)
    colnames(supplementaloutputavg) <- c("ensembl.gene.id", "HGNC.symbol", "avglfcBR1", "avglfcBR3", "avglfcBR5", "avglfcBR7", "avglfcBR9", "cluster")
    pval <- .pt[supplementaloutputavg$ensembl.gene.id, ]
    colnames(pval) <- c("Pval_FSSBR1-MockBR1", "Pval_FSSBR3-MockBR3", "Pval_FSSBR5-MockBR5", "Pval_FSSBR7-MockBR7", "Pval_FSSBR9-MockBR9",
                         "Pval_BrazilBR1-MockBR1", "Pval_BrazilBR3-MockBR3", "Pval_BrazilBR5-MockBR5", "Pval_BrazilBR7-MockBR7", "Pval_BrazilBR9-MockBR9")
    padj <- .pt.adjusted[supplementaloutputavg$ensembl.gene.id, ]
    colnames(padj) <- c("Pvaladj_FSSBR1-MockBR1", "Pvaladj_FSSBR3-MockBR3", "Pvaladj_FSSBR5-MockBR5", "Pvaladj_FSSBR7-MockBR7", "Pvaladj_FSSBR9-MockBR9",
                         "Pvaladj_BrazilBR1-MockBR1", "Pvaladj_BrazilBR3-MockBR3", "Pvaladj_BrazilBR5-MockBR5", "Pvaladj_BrazilBR7-MockBR7", "Pvaladj_BrazilBR9-MockBR9")
    supplementaloutputavg <- merge(supplementaloutputavg, pval,
        by.x = "ensembl.gene.id", by.y = "row.names", all.x = TRUE, all.y = FALSE
    )
    supplementaloutputavg <- merge(supplementaloutputavg, padj,
        by.x = "ensembl.gene.id", by.y = "row.names", all.x = TRUE, all.y = FALSE
    )
    write.csv(supplementaloutputavg, file.path(deresults_path, "suppldeavg.csv"))
} else {
    stop("WARNING: gene lists not in the same order")
}
# ### barplot combining the two strains for publication 

# #Brazil DE genes
dataMatrixdebr <- .coefs[,6:10]
sigMaskbr <- dataMatrixdebr * (results.lmde[, 6:10]**2) # 1 if significant, 0 otherwise
ExpressMatrixdebr <- subset(dataMatrixdebr, rowSums(sigMaskbr) != 0)

dataMatrixdefss <- .coefs[, 1:5]
sigMaskfss <- dataMatrixdefss * (results.lmde[, 1:5]**2) # 1 if significant, 0 otherwise
ExpressMatrixdefss <- subset(dataMatrixdefss, rowSums(sigMaskfss) != 0)

braziluniqueDEgenes <- setdiff(rownames(ExpressMatrixdebr), rownames(ExpressMatrixdefss))
fssuniqueDEgenes <- setdiff(rownames(ExpressMatrixdefss), rownames(ExpressMatrixdebr))
# number of unique genes in brazil comparison is 0 therefore combining all into one bar plot 
results.lmdecombined <- results.lmde
colnames(results.lmdecombined) <- rep(c("BR1", "BR3", "BR5", "BR7", "BR9"), 2)
vizualize_DE_genes_bp(results.lmdecombined, file.path(deresults_path, "barplotcombined.svg"))

df <- as.data.frame(matrix(ncol=3, nrow=1))
colnames(df) <- c("type", "count", "region")
generate_data_frame<-function(df, comp1, comp2, region) {
    fssbr1up <- results.lmde[results.lmde[, comp1]==1, ]
    fssbr1dwn <- results.lmde[results.lmde[, comp1] == -1, ]
    brazilbr1up <- results.lmde[results.lmde[, comp2] == 1, ]
    brazilbr1dwn <- results.lmde[results.lmde[, comp2] == -1, ]
    upintersect <- length(intersect(rownames(brazilbr1up), rownames(fssbr1up)))
    upbruniq <- length(setdiff(rownames(brazilbr1up), rownames(fssbr1up)))
    fssbruniq <- length(setdiff(rownames(fssbr1up), rownames(brazilbr1up)))
    dwnintersect <- length(intersect(rownames(brazilbr1dwn), rownames(fssbr1dwn)))
    dwnbruniq <- length(setdiff(rownames(brazilbr1dwn), rownames(fssbr1dwn)))
    dwnfssuniq <- length(setdiff(rownames(fssbr1dwn), rownames(brazilbr1dwn)))
    df <- rbind(df, c("FSS & brazil up", upintersect, region))
    df <- rbind(df, c("brazil uniq. up", upbruniq, region))
    df <- rbind(df, c("FSS uniq. up", fssbruniq, region))
    df <- rbind(df, c("FSS & brazil dwn", -dwnintersect, region))
    df <- rbind(df, c("brazil uniq. dwn", -dwnbruniq, region))
    df <- rbind(df, c("FSS uniq. dwn", -dwnfssuniq, region))

    return(df)
}
df <- generate_data_frame(df, "TRFSS.BRBR1 -TRmock.BRBR1", "TRBrazil.BRBR1 -TRmock.BRBR1", "BR1")
df <- generate_data_frame(df, "TRFSS.BRBR3 -TRmock.BRBR3", "TRBrazil.BRBR3 -TRmock.BRBR3", "BR3")
df <- generate_data_frame(df, "TRFSS.BRBR5 -TRmock.BRBR5", "TRBrazil.BRBR5 -TRmock.BRBR5", "BR5")
df <- generate_data_frame(df, "TRFSS.BRBR7 -TRmock.BRBR7", "TRBrazil.BRBR7 -TRmock.BRBR7", "BR7")
df <- generate_data_frame(df, "TRFSS.BRBR9 -TRmock.BRBR9", "TRBrazil.BRBR9 -TRmock.BRBR9", "BR9")
df <- df[-1, ]
df$count <- as.numeric(df$count)

df$type <- factor(df$type, levels = c("FSS uniq. up", "brazil uniq. up", "FSS & brazil up",
          "FSS uniq. dwn", "brazil uniq. dwn", "FSS & brazil dwn"))
write.csv(t(df), file.path(deresults_path, "barplottable.csv"), quote=FALSE, col.names=FALSE)
ggplot(df, aes(
    x = region, y = count, fill = type,
    label = df$count
)) +
    geom_bar(stat = "identity", position = "stack") +
    # geom_text(size = 5, position = position_stack(vjust = 0) )+
    # theme_light() +
    theme_minimal() +
    scale_fill_manual(values = c("FSS & brazil dwn"="#04046c",  "brazil uniq. dwn"="#3030de", 
                "FSS uniq. dwn"="#7b7bcd",  "FSS uniq. up"="#e27777", "brazil uniq. up"="#e61a1a",
                "FSS & brazil up"="darkred" )) +
    # xlab("Time point")
    ylab("# of Differentially Expressed Genes") +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10)
    )
ggsave(file.path(deresults_path, "barplot_combinedv2.png"), width = 4, height =6, units = "in", dpi = 300, bg="white")
ggsave(file.path(deresults_path, "barplot_combinedv2.svg"), width = 4, height = 6, units = "in", dpi = 300, bg = "white")
