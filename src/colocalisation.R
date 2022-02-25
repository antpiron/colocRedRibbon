
## TODO: Add function to run over all genes


are.cols <- function (dt, cols)
{
    for (col in cols)
        if ( is.null(dt[[col]]) )
            stop(paste0("'", col, "' does not exist."))
}

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
                           columns=NULL, risk=NULL, effect=`>=`)
{
    ## TODO: add a parameter to force GWAS risk increase and compute eQTL in only one direction like c(a="or.increase", b="beta.increase")
    .columns <- c(id="id", position="position", a="a", b="b",
                  a.n="a.n", a.eaf="a.eaf", a.or="a.or", a.beta="a.beta", 
                  b.n="b.n", b.eaf="b.eaf", b.or="b.or", b.beta="b.beta")
    .columns[names(columns)]  <- columns
    columns <- .columns
    
    
    dt <-  as.data.table(data)
    .intersect <- intersect(columns, names(dt))
    dt <- dt[, .intersect, with=FALSE]
    .swap_columns <- setNames(names(columns), columns)
    colnames(dt) <- sapply(colnames(dt), function(x) .swap_columns[x])
    
    if (! is.null(risk) )
    {
        if (! is.function(effect) )
            stop("'effect' should be a function.")
 
        
        if ( "a" == risk)
        {
            print(dt)
            are.cols(dt, c("a.or", "a.eaf", "b.beta", "b.eaf"))
            
            dt[a.or < 1.0, c("a.or", "a.eaf", "b.beta", "b.eaf") := list(1.0 / a.or, 1 - a.eaf,
                                                                         -b.beta, 1 - b.eaf) ]
            dt <- dt[effect(b.beta, 0),]
        } else if ( "b" == risk )
        {
            are.cols(dt, c("b.or", "b.eaf", "a.beta", "a.eaf"))

            dt[b.or < 1.0, c("b.or", "b.eaf", "a.beta", "a.eaf") := list(1.0 / b.or, 1 - b.eaf,
                                                                 -a.beta, 1 - a.eaf) ]
            dt <- dt[effect(a.beta, 0),]
        } else
            stop("risk should be either 'a' or 'b'.")
    }

    a.pval <- dt$a
    b.pval <- dt$b
    pos <- dt$position
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

    a.eaf <- dt.rr$a.eaf
    mylist.a <- list(pvalues=dt.rr$a,
                     N=dt.rr$a.n,
                     MAF=ifelse(a.eaf > 0.5, 1-a.eaf, a.eaf),
                     snp=dt.rr$id,
                     type="quant"
                     )
    
    b.eaf <- dt.rr$b.eaf
    mylist.b <- list(pvalues=dt.rr$b,
                     N=dt.rr$b.n,
                     MAF=ifelse(b.eaf > 0.5, 1-b.eaf, b.eaf),
                     snp=dt.rr$id,
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
ggplot.RedRibbonColoc <- function(self, plot.order=1:4, show.title=TRUE, labels=NULL, tss = NULL, shortid = NULL)
{
    if (is.null(labels) )
        labels  <- c(self$columns[["a"]], self$columns[["b"]])
    
    gg_quad <- ggplot(self$rr, labels=labels, quadrants=self$quadrants, base_size = 14, show.quadrants=FALSE) +
        coord_fixed(ratio=1) +
        theme(legend.position = c(0.8, 0.7),
              legend.background = element_rect(fill=alpha('white', 0.2)))

    ## theme(legend.position = "none")

    
    gg_manh <-  ggplot(self$data, aes(x=-log(a), y=-log(b)) ) +
        geom_point() +
        geom_point(data=self$data[self$quadrants$whole$positions],
                   mapping=aes(x=-log(a), y=-log(b)),
                   col="steelblue", size=2) +
        xlab(paste0("-log ", labels[[1]])) +
        ylab(paste0("-log ", labels[[2]]))
    

    if (! is.null(self$coloc) )
    {
        gg_manh <- gg_manh +
            geom_point(data=self$data[id  %in%  self$coloc$credibleSet99,],
                       mapping=aes(x=-log(a), y=-log(b)),
                                          col="green", size=2) +
            geom_point(data=self$data[id ==  self$coloc$bestSnp,],
                       mapping=aes(x=-log(a), y=-log(b)), col="red", size=2) +
            geom_text_repel(data=self$data[id ==  self$coloc$bestSnp,],
                            mapping=aes(x=-log(a), y=-log(b), label=id),
                            force = 1, nudge_x = 10, nudge_y = 10, color="red") +
            theme(legend.position = "none")
    }


    ggmanhatan <- function (axis="a", label="")
    {        
        gg <-  ggplot(self$data, aes(x=position / 1000000, y=-log(get(axis))) ) +
            theme_bw() +
            geom_point() +
            geom_point(data=self$data[self$quadrants$whole$positions],
                       mapping=aes(x=position / 1000000, y=-log(get(axis))), col="steelblue", size=2)  +
            geom_vline(xintercept = tss / 1000000, linetype="dotted", color = "red") + 
            xlab("Position (MBp)") +
            ylab(paste0("-log ", label)) +
            scale_x_continuous(breaks= pretty_breaks())

        if (! is.null(self$coloc) )
        {
            gg <- gg +
                geom_point(data=self$data[id %in%  self$coloc$credibleSet99,],
                           mapping=aes(x=position / 1000000, y=-log(b)), col="green", size=2) +
                geom_point(data=self$data[id ==  self$coloc$bestSnp,],
                           mapping=aes(x=position / 1000000, y=-log(b)), col="red", size=2)
        }

        return(gg)
    }

    gg_manh_b <- ggmanhatan("b", labels[[2]])
    gg_manh_a <- ggmanhatan("a", labels[[1]])

 
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


