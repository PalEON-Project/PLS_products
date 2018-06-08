## code for aggregation of raw data, without statistical modeling
## to be used for exploration

library(raster)
library(dplyr)

if(!'cell' %in% names(mw)) 
    mw <- mw %>% add_cells_to_dataset()

xy <- xyFromCell(base_raster, getValues(base_raster))
grid <- tibble(x = xy[,1], y = xy[,2], cell = getValues(base_raster))

## stem density

cell_dens <- mw %>% group_by(cell) %>% summarize(dens = mean(density, na.rm = TRUE))


cell_dens <- cell_dens %>% right_join(grid, by = c("cell" = "cell"))


if(FALSE) {
    cell_dens <- cell_dens %>% mutate(dens_tr = ifelse(dens > 500, 500, dens))
    image.plot(sort(unique(cell_dens$x)), sort(unique(cell_dens$y)), matrix(cell_dens$dens_tr, 296, 180)[,180:1],xlim=c(-70000,1100000))

    tmp <- cell_dens$dens
    tmp[tmp>500] <- 500
    pmap(tmp, as.matrix(dplyr::select(cell_dens, x, y)))
}

## biomass

biomass_avg <- mw %>% dplyr::select(biomass1, biomass2) %>% as.matrix(.)  %>% apply(1, mean, na.rm = TRUE)

mw2 <- mw %>% mutate(biomass = ifelse(num_trees == 0, 0,
                                        ifelse(num_trees == 1, biomass1, biomass_avg)))


cell_biom <- mw2 %>% filter(!is.na(biomass) & !is.na(density_for_biomass))  %>% group_by(cell) %>% summarize(biom = mean(biomass*density_for_biomass / kg_per_Mg, na.rm = TRUE))
cell_biom <- cell_biom %>% right_join(grid, by = c("cell" = "cell"))

if(FALSE) {
    cell_biom <- cell_biom %>% mutate(biom_tr = ifelse(biom > 300, 300, biom))
    image.plot(sort(unique(cell_biom$x)), sort(unique(cell_biom$y)), matrix(cell_biom$biom_tr, 296, 180)[,180:1],xlim=c(-70000,1100000))

    image.plot(sort(unique(cell_biom$x)), sort(unique(cell_biom$y)), matrix(cell_biom$biom_tr, 296, 180)[,180:1],xlim=c(500000,700000), ylim=c(1000000, 1200000))

    

}

if(FALSE) {
    load('cleaned_point.Rda')
    tmp = mw[mw$diam1 == 0 & mw$diam2 == 0 & !is.na(mw$L3_tree1) & !is.na(mw$L3_tree2),]
    
    gr <- tmp %>% group_by(cell) %>% summarize(cnt = n())
    gr <- gr %>% right_join(grid, by = c("cell" = "cell"))
    image.plot(sort(unique(gr$x)), sort(unique(gr$y)), matrix(gr$cnt, 296, 180)[,180:1],xlim=c(-70000,1100000))

    tmp = adf(mw[mw$cell%in%cl,])
    tmp3=cbind(tmp$cell, round(tmp$biomass*tmp$density/1000))
    tmp3[order(tmp3[,1]),]
}
