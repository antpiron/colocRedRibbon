
## TODO: Add function to run over all genes


#' Compute an enrichment based colocalization
#'
#'
#' @param data a data.frame with columns id, a, b, position. For respectively, the name
#' of the SNP, the p-value of the first analysis, the p-value of the second analysis, and
#' the position on the chromosome
#' @param algorithm the algorithm to use for minimal hypergeometric p-value searching
#' @param half the linkage desiquilibrium fitting function parameter for permutation
#' @param niter the number of iteration for adjusted p-value computation
#' @return RedRibbonColoc object
#' @export
RedRibbonColoc <- function(data, algorithm=c("ea", "classic"), half = 6300, niter=96,
                           columns=c(id="id", position="position", a="a", b="b",
                                     a.n="a.n", a.eaf="a.eaf",
                                     b.n="b.n", b.eaf="b.eaf"))
{
    ## TODO: add a parameter to force GWAS risk increase and compute eQTL in only one direction like c(a="RiskIncrease", b="EffectIncrease")
    if ( is.null(columns) )
        columns=c(id="id", position="position", a="a", b="b",
                  a.n="a.n", a.eaf="a.eaf",
                  b.n="b.n", b.eaf="b.eaf")
    
    dt <-  as.data.table(data)

    a.pval <- dt[[ columns[["a"]] ]]
    b.pval <- dt[[ columns[["b"]] ]]
    pos <- dt[[ columns[["position"]] ]]
    ## id, chr, position, nea, ea, gwas.eaf, gwas.pvalue, gwas.or, gwas.n, eqtl.eaf, eqtl.pvalue, eqtl.n, eqtl.direction
    deps <- rrho_ldfit_prediction(half, b.pval, pos)
    rr <- RedRibbon(a.pval, b.pval, correlation=newLDFIT(pos, deps, half=half) )

    whole.fraction=0.6
    if ( algorithm[[1]] == "ea" )
        quad <- quadrants(rr, whole=TRUE, whole.fraction=whole.fraction, permutation=TRUE, algorithm="ea", niter=niter)
    else
        quad <- quadrants(rr, whole=TRUE, whole.fraction=whole.fraction, permutation=TRUE, m=750, n=750, niter=niter)


    structure(list(data = dt,
                   rr = rr,
                   quadrants = quad,
                   columns = columns),
              class = "RedRibbonColoc")
        
}

#' @export
coloc <- function (self, ...)
{
    UseMethod("coloc", self)
}

#' @export
coloc.RedRibbonColoc  <- function(self, ...)
{
    ## keep the RRHO enrichment SNP if significant. Run on subset if enriched, otherwise classic coloc.
    dt.rr <- if ( self$quadrants$whole$log_padj >= -log(0.05) ) self$data[self$quadrants$whole$positions] else data

    a.eaf <- dt.rr[[ self$columns[[ "a.eaf" ]]  ]]
    mylist.a <- list(pvalues=dt.rr[[ self$columns[[ "a" ]]  ]],
                     N=dt.rr[[ self$columns[[ "a.n" ]]  ]],
                     MAF=ifelse(a.eaf > 0.5, 1-a.eaf, a.eaf),
                     snp=dt.rr[[ self$columns[[ "id" ]]  ]],
                     type="quant"
                     )
    
    b.eaf <- dt.rr[[ self$columns[[ "b.eaf" ]]  ]]
    mylist.b <- list(pvalues=dt.rr[[ self$columns[[ "b" ]]  ]],
                     N=dt.rr[[ self$columns[[ "b.n" ]]  ]],
                     MAF=ifelse(b.eaf > 0.5, 1-b.eaf, b.eaf),
                     snp=dt.rr[[ self$columns[[ "id" ]]  ]],
                     type="quant"
                     )

    coloc.abf.res <- coloc.abf(mylist.a, mylist.b)
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

    self$coloc <- coloc.res

    return(self)
}

#' @export
ggplot.RedRibbonColoc <- function(self, plot.order=1:4, show.title=TRUE, labels=c("a", "b"), tss = NULL, shortid = NULL)
{
    gg_quad <- ggplot(self$rr, labels=labels, quadrants=self$quadrants, base_size = 14, show.quadrants=FALSE) +
        coord_fixed(ratio=1) + theme(legend.position = "none")

    log.a <- paste0("-log(", self$columns[[ "a" ]],")")
    log.b <- paste0("-log(", self$columns[[ "b" ]],")")
    position.mb <- paste0(self$columns[[ "position" ]], "/1000000")
    id.col <- self$columns[[ 'id' ]]
    
    gg_manh <-  ggplot(self$data, aes_string(x=log.a, y=log.b)) +
        geom_point() +
        geom_point(data=self$data[self$quadrants$whole$positions],
                   mapping=aes_string(x=log.a, y=log.b),
                   col="steelblue", size=2) +
        xlab(paste0("-log ", labels[[1]])) +
        ylab(paste0("-log ", labels[[2]]))
    

    if (! is.null(self$coloc) )
    {
        gg_manh <- gg_manh +
            geom_point(data=self$data[get(self$columns[[ "id" ]])  %in%  self$coloc$credibleSet99,],
                       mapping=aes_string(x=log.a, y=log.b),
                                          col="green", size=2) +
            geom_point(data=self$data[get(self$columns[[ "id" ]]) ==  self$coloc$bestSnp,],
                       mapping=aes_string(x=log.a, y=log.b), col="red", size=2) +
            geom_text_repel(data=self$data[get(self$columns[[ "id" ]]) ==  self$coloc$bestSnp,],
                            mapping=aes_string(x=log.a,
                                               y=log.b,
                                               label=id.col),
                            force = 1, nudge_x = 10, nudge_y = 10, color="red") +
            theme(legend.position = "none")
    }


    
    gg_manh_b <-  ggplot(self$data,
                            aes_string(x=paste0(self$columns[[ "position" ]], "/1000000"),
                                       y=paste0("-log(", self$columns[[ "b" ]],")"))) + theme_bw()+
        geom_point() +
        geom_point(data=self$data[self$quadrants$whole$positions],
                   mapping=aes_string(x=position.mb,
                                      y=log.b), col="steelblue", size=2)  +
        geom_vline(xintercept = tss / 1000000, linetype="dotted", color = "red") + 
        xlab("Position (MBp)") +
        ylab(paste0("-log ", labels[[2]])) +
        scale_x_continuous(breaks= pretty_breaks())
    
    
    if (! is.null(self$coloc) )
    {
        gg_manh_b <- gg_manh_b +
            geom_point(data=self$data[get(self$columns[[ "id" ]]) %in%  self$coloc$credibleSet99,],
                       mapping=aes_string(x=position.mb,
                                          y=log.b), col="green", size=2) +
            geom_point(data=self$data[get(self$columns[[ "id" ]]) ==  self$coloc$bestSnp,],
                       mapping=aes_string(x=position.mb,
                                          y=log.b), col="red", size=2)
    }

    
    gg_manh_a <-  ggplot(self$data,
                         aes_string(x=paste0(self$columns[[ "position" ]], "/1000000"),
                                    y=paste0("-log(", self$columns[[ "a" ]],")"))) + theme_bw()+
        geom_point() +
        geom_point(data=self$data[self$quadrants$whole$positions],
                   mapping=aes_string(x=position.mb,
                                      y=log.a), col="steelblue", size=2)  +
        geom_vline(xintercept = tss / 1000000, linetype="dotted", color = "red") + 
        xlab("Position (MBp)") +
        ylab(paste0("-log ", labels[[1]])) +
        scale_x_continuous(breaks= pretty_breaks())

    
    if (! is.null(self$coloc) )
    {
        gg_manh_a <- gg_manh_a +
            geom_point(data=self$data[get(self$columns[[ "id" ]]) %in%  self$coloc$credibleSet99,],
                       mapping=aes_string(x=position.mb,
                                          y=log.a), col="green", size=2) +
            geom_point(data=self$data[get(self$columns[[ "id" ]]) ==  self$coloc$bestSnp,],
                       mapping=aes_string(x=position.mb,
                                          y=log.a), col="red", size=2)
    }

 
    title  <- shortid
    if (! is.null(self$coloc) )
    {
        title  <- paste0(title,
                         " - ", self$coloc$bestSnp,
                         " (PP.H4.abf = ", formatC(self$coloc$PP.H4.abf, digit=2),
                          " ; SNP.PP.H4 = ", formatC(self$coloc$SNP.PP.H4, digit=2), ")")
    }

    list.of.plots <- list(gg_quad, gg_manh_a, gg_manh, gg_manh_b)[plot.order]
    square.gg <- do.call(ggpubr::ggarrange, c(list.of.plots,
                                              ncol = 2, nrow = 2))
    gg_merge <- if (show.title) annotate_figure(square.gg, top = title) else square.gg 
    
    return(gg_merge)
}


