

RedRibbonColoc <- function(data, .ensg, half = 6300, algorithm="ea", niter=96)
{
    dt <-  as.data.table(data)

    ## id, chr, position, nea, ea, gwas.eaf, gwas.pvalue, gwas.or, gwas.n, eqtl.eaf, eqtl.pvalue, eqtl.n, eqtl.direction
    deps <- rrho_ldfit_prediction(half, dt$pval.eQTL, dt$pos)
    rr <- RedRibbon(dt$pval.GWAS, dt$pval.eQTL, correlation=newLDFIT(dt$pos, deps, half=half) )

    whole.fraction=0.6
    if ( algorithm == "ea" )
        quad <- quadrants(rr, whole=TRUE, whole.fraction=whole.fraction, permutation=TRUE, algo
rithm="ea", niter=niter)
    else
        quad <- quadrants(rr, whole=TRUE, whole.fraction=whole.fraction, permutation=TRUE, m=75
0, n=750, niter=niter)


    structure(list(data = dt,
                   rr = rr,
                   quadrants = quad),
              class = "RedRibbonColoc")
        
}

coloc.RedRibbonColoc  <- function(self)
{
    ## Coloc
    dt.rr <- if ( quadrants$whole$log_padj >= -log(0.05) ) self$dt[self$quadrants$whole$positions] else dt
        
    mylist.eQTL <- list(pvalues=dt.rr$pval.eQTL,
                        N=dt.rr$n.eQTL,
                        MAF=dt.rr$maf,
                        snp=dt.rr$rsid,
                        type="quant"
                        )
    mylist.GWAS <- list(pvalues=dt.rr$pval.GWAS,
                        N=dt.rr$n.GWAS,
                        MAF=dt.rr$maf,
                        snp=dt.rr$rsid,
                        type="quant"
                        )

    coloc.abf.res <- coloc.abf(mylist.eQTL, mylist.GWAS)
    results <- as.data.table(coloc.abf.res$results)
    
    setorder(results, -SNP.PP.H4)
    
    results[,SNP.PP.H4.cumsum:=cumsum(SNP.PP.H4)]
    
    bestSnp <- results[1,snp]
    SNP.PP.H4 <- results[1, SNP.PP.H4]
    
    ## 99% cdredible set
    ncredibleSet99 <- sum(results$SNP.PP.H4.cumsum < 0.99) + 1
    credibleSet99 <- results[1:ncredibleSet99, snp]
    PP.H4.abf <- coloc.abf.res$summary["PP.H4.abf"]
    
    coloc.res <- list(bestSnp = bestSnp, PP.H4.abf = PP.H4.abf, SNP.PP.H4 = SNP.PP.H4, ncredibleSet99 = ncredibleSet99, credibleSet99 = credibleSet99)
    
}

ggplot.RedRibbonColoc <- function(self, plot.order=1:4, show.title=TRUE)
{
    gg_quad <- ggplot(self$rr, labels=c("GWAS", "eQTL"), quadrants=self$quadrants, base_size = 14, show.quadrants=FALSE) +
        coord_fixed(ratio=1) + theme(legend.position = "none")
    
    gg_manh <-  ggplot(self$dt, aes(x=-log(pval.eQTL), y=-log(pval.GWAS))) +
        geom_point() +
        geom_point(data=self$dt[self$quadrants$whole$positions],
                   mapping=aes(x=-log(pval.eQTL), y=-log(pval.GWAS)), col="steelblue", size=2)
    
    if (! is.null(self$coloc) )
    {
        gg_manh <- gg_manh +
            geom_point(data=self$dt[rsid %in%  self$coloc$credibleSet99,],
                       mapping=aes(x=-log(pval.eQTL), y=-log(pval.GWAS)), col="green", size=2) +
            geom_point(data=self$dt[rsid ==  self$coloc$bestSnp,],
                       mapping=aes(x=-log(pval.eQTL), y=-log(pval.GWAS)), col="red", size=2) +
            geom_text_repel(data=self$dt[rsid ==  self$coloc$bestSnp,],
                            mapping=aes(x=-log(pval.eQTL), y=-log(pval.GWAS), label=rsid),
                            force = 1, nudge_x = 10, nudge_y = 10, color="red") +
            theme(legend.position = "none")
    }

    
    gg_manh_gwas <-  ggplot(self$dt, aes(x=pos/1000000, y=-log(pval.GWAS))) + theme_bw()+
        geom_point() +
        geom_point(data=self$dt[self$quadrants$whole$positions],
                   mapping=aes(x=pos/1000000, y=-log(pval.GWAS)), col="steelblue", size=2)  +
        geom_vline(xintercept = self$tss / 1000000, linetype="dotted", color = "red") + 
        xlab("Position (MBp)") +
        scale_x_continuous(breaks= pretty_breaks())
    
    if (! is.null(self$coloc) )
    {
        gg_manh_gwas <- gg_manh_gwas +
            geom_point(data=self$dt[rsid %in%  self$coloc$credibleSet99,],
                       mapping=aes(x=pos/1000000, y=-log(pval.GWAS)), col="green", size=2) +
            geom_point(data=self$dt[rsid ==  self$coloc$bestSnp,],
                       mapping=aes(x=pos/1000000, y=-log(pval.GWAS)), col="red", size=2) 
    }

    gg_manh_eqtl <-  ggplot(self$dt, aes(x=pos/1000000, y=-log(pval.eQTL))) + theme_bw()+
        geom_point() +
        geom_point(data=self$dt[self$quadrants$whole$positions],
                   mapping=aes(x=pos/1000000, y=-log(pval.eQTL)), col="steelblue", size=2)  +
        geom_vline(xintercept = self$tss / 1000000, linetype="dotted", color = "red")  +
        xlab("Position (MBp)") +
        scale_x_continuous(breaks= pretty_breaks())

    if (! is.null(self$coloc) )
    {
        gg_manh_eqtl <- gg_manh_eqtl +
            geom_point(data=self$dt[rsid %in%  self$coloc$credibleSet99,],
                       mapping=aes(x=pos/1000000, y=-log(pval.eQTL)), col="green", size=2) +
            geom_point(data=self$dt[rsid ==  self$coloc$bestSnp,],
                       mapping=aes(x=pos/1000000, y=-log(pval.eQTL)), col="red", size=2) 
    }

    title  <- self$shortid
    if (! is.null(self$coloc) )
    {
        title  <- paste0(title,
                         " - ", self$coloc$bestSnp,
                         " (PP.H4.abf = ", formatC(self$coloc$PP.H4.abf, digit=2),
                          " ; SNP.PP.H4 = ", formatC(self$coloc$SNP.PP.H4, digit=2), ")")
    }

    list.of.plots <- list(gg_quad, gg_manh_eqtl, gg_manh, gg_manh_gwas)[plot.order]
    square.gg <- do.call(ggpubr::ggarrange, c(list.of.plots,
                                              ncol = 2, nrow = 2))
    gg_merge <- if (show.title) annotate_figure(square.gg, top = title) else square.gg 
    
    return(gg_merge)
}
