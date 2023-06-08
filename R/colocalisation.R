
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
#' the position on the chromosome. Additionnal columns are for colocalisation:
#' [ab].n (number of samples of the analysis), [ab].eaf (effect allele frequency),
#' [ab].or (odd-ratio) and/or [ab].beta (effect slope).
#' @param algorithm the algorithm to use for minimal hypergeometric p-value searching
#' @param half the linkage desiquilibrium fitting function parameter for permutation
#' @param niter the number of iteration for adjusted p-value computation
#' @param risk the GWAS dataset with an odd-ratio (e.g. a.or or b.or) either NULL, 'a' or 'b'
#' @param effect an operator like `>=` or `<=` indicating the effect direction for the non-risk
#' dataset
#' @param columns a named vector with column names in the `data' data frame
#' @param shortlist run RedRibbon before coloc (default = TRUE)
#' 
#' @return RedRibbonColoc object
#' @export
RedRibbonColoc <- function(data, algorithm=c("ea", "classic"), half = 6300, niter=96,
                           risk=NULL, effect=`>=`, columns=NULL, shortlist=TRUE)
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

    are.cols(dt, c("id", "a", "b"))
    
    
    if (! is.null(risk) )
    {
        if (! is.function(effect) )
            stop("'effect' should be a function.")
 
        
        if ( "a" == risk)
        {
            ## print(dt)
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

    if (shortlist)
    {
        pos <- dt$position
        deps <- rrho_ldfit_prediction(half, b.pval, pos)
        rr <- RedRibbon(a.pval, b.pval, correlation=newLDFIT(pos, deps, half=half) )
        
        whole.fraction=0.6
        if ( algorithm[[1]] == "ea" )
            quad <- quadrants(rr, whole=TRUE, whole.fraction=whole.fraction, permutation=TRUE, algorithm="ea", niter=niter)
        else
            quad <- quadrants(rr, whole=TRUE, whole.fraction=whole.fraction, permutation=TRUE, m=750, n=750, niter=niter)
    } else {
        rr <- RedRibbon(a.pval, b.pval, NULL)
        quad <- NULL
    }

    structure(list(data = dt,
                   rr = rr,
                   quadrants = quad,
                   columns = columns),
              class = "RedRibbonColoc")
        
}

#' Compute a colocalisation
#'
#' @param self a colocRedRibbon object
#' 
#' @return RedRibbonColoc object
#' @export
coloc <- function (self, ...)
{
    UseMethod("coloc", self)
}


#' Compute a colocalisation
#'
#' @param self a colocRedRibbon object
#' @param n.reduce function to reduce the number of sample columns for a and b into a number (Default: min)
#' This is needed as some eQTL/GWAS toolchains output an effective number of samples by SNP.
#' 
#' @return RedRibbonColoc object
#' @method coloc RedRibbonColoc
#' @export
coloc.RedRibbonColoc  <- function(self, n.reduce = min)
{
    ## keep the RRHO enrichment SNP if significant. Run on subset if enriched, otherwise classic coloc.
    dt.rr <- if ( ! is.null(self$quadrants) &&  self$quadrants$whole$log_padj >= -log(0.05) ) self$data[self$quadrants$whole$positions] else self$data

    a.n <- ceiling(n.reduce(dt.rr$a.n, na.rm=TRUE))
    if (a.n < 2)
        stop(paste0("coloc.RedRibbonColoc(): number of samples for `a` is abnormaly low (", a.n, " < 2)"))
    a.eaf <- dt.rr$a.eaf
    mylist.a <- list(pvalues=dt.rr$a,
                     N=a.n,
                     MAF=ifelse(a.eaf > 0.5, 1-a.eaf, a.eaf),
                     snp=dt.rr$id,
                     type="quant"
                     )
    
    b.n <- ceiling(n.reduce(dt.rr$b.n, na.rm=TRUE))
    if (b.n < 2)
        stop(paste0("coloc.RedRibbonColoc(): number of samples for `b` is abnormaly low (", b.n, " < 2)"))
    b.eaf <- dt.rr$b.eaf
    mylist.b <- list(pvalues=dt.rr$b,
                     N=b.n,
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

#' Plot a colocalisation with ggplot
#'
#' @param self a colocRedRibbon object
#' @param plot.order a vector specifying the plot order (default =1:4,  1 = RedRibbon plot, 2 =  manhantan plot for `a`, 3 = plot for `a`and 'b`, 4 = manhantan plot for `b`).
#' @param show.title shows the title (default = TRUE)
#' @param labels axis labels
#' @param tss transcription start site
#' @param shortid name of the gene
#' @param title the title of the plot
#' @param highlight a list of SNPs to highlight
#' @param .log10 output log10 pval (default = FALSE)
#' 
#' @return ggplot object
#' 
#' @export
ggRedRibbonColoc <- function (self, plot.order=1:4, show.title=TRUE,
                              labels=NULL, tss = NULL, shortid = NULL, title = NULL,
                              highlight=c(), .log10 = FALSE)
{
     UseMethod("ggRedRibbonColoc")
}


#' Plot a colocalisation with ggplot
#'
#' @param self a colocRedRibbon object
#' @param plot.order a vector specifying the plot order (default =1:4,  1 = RedRibbon plot, 2 =  manhantan plot for `a`, 3 = plot for `a`and 'b`, 4 = manhantan plot for `b`).
#' @param show.title shows the title (default = TRUE)
#' @param labels axis labels
#' @param tss transcription start site
#' @param shortid name of the gene
#' @param title the title of the plot
#' @param highlight a list of SNPs to highlight
#' @param .log10 output log10 pval (default = FALSE)
#' 
#' @return ggRedRibbonColoc RedRibbonColoc
#' @export
ggRedRibbonColoc.RedRibbonColoc <- function(self, plot.order=1:4, show.title=TRUE,
                                            labels=NULL, tss = NULL, shortid = NULL, title = NULL,
                                            highlight=c(), .log10 = FALSE)
{
    if (is.null(labels) )
        labels  <- c(self$columns[["a"]], self$columns[["b"]])
    
    gg_quad <- ggRedRibbon(self$rr, labels=labels, quadrants=self$quadrants, base_size = 14, show.quadrants=FALSE, .log10 = .log10) +
        coord_fixed(ratio=1) +
        theme(legend.position = c(0.8, 0.7),
              legend.background = element_rect(fill=alpha('white', 0.2)))

    ## theme(legend.position = "none")
    mylog <- if (.log10) log10 else log
    log.label <- if (.log10) "log10" else "log"
    
    gg_manh <-  ggplot(self$data, aes(x=-mylog(a), y=-mylog(b)) ) +
        geom_point() +
        xlab(paste0("-", log.label, " ", labels[[1]])) +
        ylab(paste0("-", log.label, " ", labels[[2]])) +
        theme_bw()

    if (! is.null(self$quadrants) )
        gg_manh <- gg_manh + geom_point(data=self$data[self$quadrants$whole$positions],
                                        mapping=aes(x=-mylog(a), y=-mylog(b)),
                                        col="steelblue", size=2) 

    if (! is.null(self$coloc) )
    {
        gg_manh <- gg_manh +
            geom_point(data=self$data[id  %in%  self$coloc$credibleSet99,],
                       mapping=aes(x=-mylog(a), y=-mylog(b)),
                                          col="green", size=2) +
            geom_point(data=self$data[id ==  self$coloc$bestSnp,],
                       mapping=aes(x=-mylog(a), y=-mylog(b)), col="red", size=2) +
            geom_text_repel(data=self$data[id ==  self$coloc$bestSnp,],
                            mapping=aes(x=-mylog(a), y=-mylog(b), label=id),
                            force = 1, nudge_x = 10, nudge_y = 10, color="red") +
            ## geom_point(data=self$data[id %in%  highlight,],
            ##            mapping=aes(x=-mylog(a), y=-mylog(b)), col="yellow", size=2) +
            geom_text_repel(data=self$data[id %in%  highlight,],
                            mapping=aes(x=-mylog(a), y=-mylog(b), label=id),
                            force = 1, nudge_x = 3, nudge_y = 3, color="darkgray") +
            theme(legend.position = "none")
    }


    ggmanhatan <- function (axis="a", label="")
    {        
        gg <-  ggplot(self$data, aes(x=position / 1000000, y=-mylog(get(axis))) ) +
            theme_bw() +
            geom_point()  +
            geom_vline(xintercept = tss / 1000000, linetype="dotted", color = "red") + 
            xlab("Position (MBp)") +
            ylab(paste0("-", log.label, " ", label)) +
            scale_x_continuous(breaks= pretty_breaks())

        if (! is.null(self$quadrants) )
            gg <- gg + geom_point(data=self$data[self$quadrants$whole$positions],
                                  mapping=aes(x=position / 1000000, y=-mylog(get(axis))), col="steelblue", size=2)
        
        if (! is.null(self$coloc) )
        {
            gg <- gg +
                geom_point(data=self$data[id %in%  self$coloc$credibleSet99,],
                           mapping=aes(x=position / 1000000, y=-mylog(get(axis))), col="green", size=2) +
                geom_point(data=self$data[id ==  self$coloc$bestSnp,],
                           mapping=aes(x=position / 1000000, y=-mylog(get(axis))), col="red", size=2) +
                geom_text_repel(data=self$data[id %in%  highlight,],
                                mapping=aes(x=position / 1000000, y=-mylog(get(axis)), label=id),
                                force = 1, nudge_x = 3, nudge_y = 3, color="darkgray")
        }

        return(gg)
    }

    gg_manh_b <- ggmanhatan("b", labels[[2]])
    gg_manh_a <- ggmanhatan("a", labels[[1]])


    if (is.null(title))
    {
        title  <- shortid
        if (! is.null(self$coloc) )
        {
            title  <- paste0(title,
                             " - ", self$coloc$bestSnp,
                             "\n(PP.H4.abf = ", formatC(self$coloc$PP.H4.abf, digit=2),
                             " ; SNP.PP.H4 = ", formatC(self$coloc$SNP.PP.H4, digit=2), ")")
        }
    }

    list.of.plots <- list(gg_quad, gg_manh_a, gg_manh, gg_manh_b)[plot.order]
    .ncol <- .nrow <- 1
    if (1 < length(list.of.plots))
    {
        .ncol <- if (length(list.of.plots) <= 2) 1 else 2
        .nrow <- if (length(list.of.plots) <= 2) length(list.of.plots) else 2
    }
    square.gg <- do.call(ggpubr::ggarrange, c(list.of.plots,
                                              ncol = .ncol, nrow = .nrow))
    gg_merge <- if (show.title) annotate_figure(square.gg, top = title) else square.gg 

    gg_merge <- gg_merge + bgcolor("#FFFFFF") + border(size = 0)
    
    return(gg_merge)
}

#' Plot a colocalisation with ggplot
#'
#' @param self a colocRedRibbon object
#' @param plot.order a vector specifying the plot order (default =1:4,  1 = RedRibbon plot, 2 =  manhantan plot for `a`, 3 = plot for `a`and 'b`, 4 = manhantan plot for `b`).
#' @param show.title shows the title (default = TRUE)
#' @param labels axis labels
#' @param tss transcription start site
#' @param shortid name of the gene
#' @param title the title of the plot
#' @param highlight a list of SNPs to highlight
#' 
#' @return ggplot object
#'
#' @method ggplot RedRibbonColoc
#' @export
ggplot.RedRibbonColoc <- function(self, plot.order=1:4, show.title=TRUE,
                                  labels=NULL, tss = NULL, shortid = NULL, title = NULL,
                                  highlight=c())
{
    return(ggcolocRedRibbon(self, plot.order, show.title, labels, tss, shortid, title, highlight))
}
