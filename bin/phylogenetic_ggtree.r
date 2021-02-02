#!/usr/bin/env Rscript

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(ggtree)


################################################
################################################
## VALIDATE COMMAND-LINE PARAMETERS           ##
################################################
################################################


# Run some tests with  ggtree package
plotTree <- function(inputDir, outputFile){

    setwd(inputDir)

    info <- read.csv("info.csv")
    tree <- read.tree("tree.nwk")
    cols <- c(HCMC='black', Hue='purple2', KH='skyblue2')

    p <- ggtree(tree, layout='circular') %<+% info + 
        geom_tippoint(aes(color=location)) + 
        scale_color_manual(values=cols) + 
        geom_tiplab2(aes(label=name), align=T, linetype=NA, size=2, offset=4, hjust=0.5) +
        geom_tiplab2(aes(label=year), align=T, linetype=NA, size=2, offset=8, hjust=0.5)

    heatmapData=read.csv("res_genes.csv", row.names=1)
    rn <- rownames(heatmapData)
    heatmapData <- as.data.frame(sapply(heatmapData, as.character))
    rownames(heatmapData) <- rn

    heatmap.colours <- c("white","grey","seagreen3","darkgreen",
                        "green","brown","tan","red","orange",
                        "pink","magenta","purple","blue","skyblue3",
                        "blue","skyblue2")

    names(heatmap.colours) <- 0:15

    gheatmap(p, heatmapData, offset = 10, color=NULL, colnames_position="top", colnames_angle=90, colnames_offset_y = 1, hjust=0, font.size=2) +
            scale_fill_manual(values=heatmap.colours, breaks=0:15)

    ggsave(outputFile, width = 50, height = 50, units = "cm", limitsize = FALSE)
}

plotTree("/home/hnguyen/Documents/phylogenetic/data/plotTree/tree_example_april2015/","/home/hnguyen/Documents/RunNextflow/RScript/PhylogeneticTree.pdf")

