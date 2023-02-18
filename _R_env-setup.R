options(menu.graphics=FALSE) # disable graphics in menu
# install.packages(c('devtools', 'gam', 'RColorBrewer', 'BiocManager'))
update.packages(ask=F)
install.packages(c('Seurat','sctransform'))
BiocManager::install(c("scran",'BiocParallel', "MAST","monocle","ComplexHeatmap","slingshot"))