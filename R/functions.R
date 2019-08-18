
createModelTable = function(m,
	dir=NA,
	filename= NA,
	filenameInfix="",
	predictor.names=NULL,
	caption="",
	label=NA,
	digits.coef=2,
	digits.coefse=digits.coef,
	digits.statistic=1,
	...
) {
	require(Hmisc)
	if (class(m) %in% c('lrm','ols')) { require(rms) }
	coefs <- as.data.frame(summary(m)@coefs)

	coefs$Stars <- ifelse(coefs[,4] > .1,
				"",
				ifelse(coefs[,4] < .0001,
					"\\textbf{****}",
				ifelse(coefs[,4] < 0.001,
					"\\textbf{***}",
				ifelse(coefs[,4] < 0.01,
					"\\textbf{**}",
				ifelse(coefs[,4] < 0.05,
					"\\textbf{*}",
					"\\textbf{$.$}"
				)))))
	coefs[,1] <- round(coefs[,1],digits=digits.coef)
	coefs[,2] <- round(coefs[,2],digits=digits.coefse)
	coefs[,3] <- round(coefs[,3],digits=digits.statistic)
	coefs[,4] <- ifelse(round(coefs[,4],digits=2) == .1,
				"=0.1",
				ifelse(round(coefs[,4],digits=1) == 1,
					">0.9",
				ifelse(coefs[,4] > 0.1,
					paste(">",round(coefs[,4],digits=1),sep=""),
				ifelse(coefs[,4] < .0001,
					"<0.0001",
				ifelse(coefs[,4] < 0.001,
					"<0.001",
				ifelse(coefs[,4] < 0.01,
					"<0.01",
				ifelse(coefs[,4] < 0.05,
					"<0.05",
					"<0.1"
				))))))
	)
	colnames(coefs) <- c("$\\hat{\\beta}$","$\\hat{SE}_{\\hat{\\beta}}$", "\\textbf{z}","\\textbf{p}","\ ")

	if (is.na(label)) label= paste("table:", filename, sep="")
	if (is.na(filename)) filename= paste(row.names(coefs), collapse="_")
	filename= paste(dir, "table-", filename, filenameInfix, ".tex", sep="")
	if (!is.null(predictor.names)) row.names(coefs) = predictor.names


	sink(file=filename)
		latex(coefs,file=filename,
			title="",
			table.env=FALSE,
			booktabs=TRUE,
			caption=caption,
			label=label, ...
		)
	sink()
}

my.lmer.nagelkerke <- function(f, d) {
	# expects as first argument a character vector with
	#     1) the fixed effect part of the formula and 2) the random effect parts
	# as well as the data set

	require(rms)
	require(lme4)

	lmer.full= lmer(formula= as.formula(paste(f, collapse="+")), d, family="binomial")
	logLik.lmer.full= as.numeric(logLik(lmer.full))
	N.lmer.full= nrow(lmer.full@X)
	cat(paste("Full mixed model: L=", logLik.lmer.full, ", N=", N.lmer.full, "\n", sep=""))

	lmer.intercept= lmer(formula= as.formula(paste(unlist(strsplit(f[1], "~"))[1], paste("1", f[2], sep=" + "), sep="~ ")), data= d, family="binomial")
	logLik.lmer.intercept= as.numeric(logLik(lmer.intercept))
	N.lmer.intercept= nrow(lmer.intercept@X)
	cat(paste("Intercept mixed model: L=", logLik.lmer.intercept, ", N=", N.lmer.intercept, "\n", sep=""))

	lrm.full= lrm(formula= as.formula(f[1]), data= d)
	logLik.lrm.full= as.numeric(deviance(lrm.full)[2] / - 2)
	N.lrm.full= as.numeric(lrm.full$stats[1])
	cat(paste("Full ordinary model: L=", logLik.lrm.full, ", N=", N.lrm.full, "\n", sep=""))

	# lrm does not allow intercept model. So, I extract deviance and N from full model
	# (which also contains the intercept model stats)
	logLik.lrm.intercept= as.numeric(deviance(lrm.full)[1] / - 2)
	N.lrm.intercept= as.numeric(lrm.full$stats[1])
	cat(paste("Intercept ordinary model: L=", logLik.lrm.intercept, ", N=", N.lrm.intercept, "\n", sep=""))

	coxsnell.lmer= 1 - exp((logLik.lmer.intercept - logLik.lmer.full) * (2/N.lmer.full))
	nagelkerke.lmer= coxsnell.lmer / (1 - exp(logLik.lmer.intercept * (2/N.lmer.full)))
	cat(paste("Full mixed model evaluated against mixed intercept model:\n\tCoxSnell R2: ", coxsnell.lmer, "\n\tNagelkerke R2: ", nagelkerke.lmer,"\n", sep=""))
	coxsnell.lrm= 1 - exp((logLik.lrm.intercept - logLik.lmer.full) * (2/N.lmer.full))
	nagelkerke.lrm= coxsnell.lrm / (1 - exp(logLik.lrm.intercept * (2/N.lmer.full)))
	cat(paste("Full mixed model evaluated against ordinary intercept model:\n\tCoxSnell R2: ", coxsnell.lrm, "\n\tNagelkerke R2: ", nagelkerke.lrm,"\n", sep=""))
	coxsnell.lrm.full= 1 - exp((logLik.lrm.intercept - logLik.lrm.full) * (2/N.lrm.full))
	nagelkerke.lrm.full= coxsnell.lrm.full / (1 - exp(logLik.lrm.intercept * (2/N.lrm.full)))
	cat(paste("Full ordinary model evaluated against ordinary intercept model:\n\tCoxSnell R2: ", coxsnell.lrm.full, "\n\tNagelkerke R2: ", nagelkerke.lrm.full,"\n", sep=""))
	coxsnell.lrm.2= 1 - exp((logLik.lrm.intercept - logLik.lmer.intercept) * (2/N.lmer.full))
	nagelkerke.lrm.2= coxsnell.lrm.2 / (1 - exp(logLik.lrm.intercept * (2/N.lmer.full)))
	cat(paste("Mixed intercept model evaluated against ordinary intercept model:\n\tCoxSnell R2: ", coxsnell.lrm.2, "\n\tNagelkerke R2: ", nagelkerke.lrm.2,"\n", sep=""))
}
myUpdate <- function(x) { if (is.factor(x)) { return(as.factor(as.character(x))) } }

multiplot <- function(row,col){
     op <- par(mfrow=c(row,col),pty="s")
}

NAtoX <- function(v, x) # replace all NA in v by x
  {
    if (is.factor(v)) {
              as.factor(ifelse(as.character(v) == "", as.character(x), as.character(v)))
            }
    else if (is.numeric(v)) {
              as.numeric(ifelse(is.na(v), x, v))
            }
    else if (is.logical(v)) {
              as.logical(ifelse(is.na(v), x, v))
            }
    else {
        ifelse(v == "", x, v)
      }
  }

XtoNA <- function(v, x) # replace all x in v by NA
{
      if (is.factor(v)) {
        as.factor(ifelse(v == x, NA, as.character(v)))
      }
        else if (is.numeric(v)) {
        as.numeric(ifelse(v == x, NA, v))
      }
      else {
          ifelse(v == x, NA, v)
      }
}


predict.lmerBin <- function(object, X){
        # BY SPENCER GRAVES
        # object has class "lmer"
        # X = model matrix with columns matching object@X

        if(missing(X))
        X <- object@X
        b <- fixef(object)
        X %*% b
}

odds <- function(p) {
        p / (1-p)
}



myRug <- function (x, y, ticksize = 0.03, lwd = 0.5, col = par("fg"), stacked = F, ...)
{
    x <- as.vector(x)
    ok <- is.finite(x)

    if (stacked == T) {
		d <- density(x)
		maxy <- max(d$y)
	  	lines((d$x), y - .005 - ticksize - (d$y / (maxy * 20)), lwd = lwd, lty=2, col= "black")
    }

    x <- x[ok]
    segments(x0 = x, y0 = y, x1 = x, , y1 = y - ticksize,  lwd = lwd,
        col = col)
}



## --------------------------------------------------
# some general graphics functions
## --------------------------------------------------
myGplot.defaults = function(
	type = c("paper","poster","slides")[1],
	base_size= if (type == "paper") 11 else 18,
	margin=c(0.4,0.2,0.5,0.1)
)
{
	require(ggplot2)

	theme_set(theme_bw(base_size=base_size))
	theme_update(
		axis.text.x = element_text(size=base_size-2, vjust=1),
		axis.text.y = element_text(size=base_size-2, hjust=1, vjust=.5),
		axis.title.x = element_text(size=base_size, vjust=0, hjust=0.6, face = "bold"),
		axis.title.y = element_text(angle=90, size=base_size, hjust=0.6, vjust=0.5, face = "bold"),
		legend.title = element_text(size=base_size, face = "bold", hjust= 0),
		legend.text = element_text(size=base_size-2)
		#  1) top, 2) right, 3) bottom, 4) left
		# plot.margin = unit(margin, "lines")
	)
}


panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    if(!is.numeric(x) | !is.numeric(y)) {
	r2 = spearman2(x, y)["Adjusted rho2"]
    } else {
      r2 <- cor(x, y, use="complete.obs")^2
    }
    txt <- format(c(r2, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex * sqrt(r2) , col="#000099")
}

panel.smooth <- function (x, y, col = par("col"), bg = NA, pch = par("pch"),
    cex = 1, col.smooth = "#000099", span = 2/3, iter = 3, ...)
{
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok))
        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),
            col = col.smooth, lwd=2, ...)
}

my.pairs = function(...) {
	pairs(lower.panel=panel.smooth, upper.panel=panel.cor, ...)
}
# how to use:
#par(mar=c(1,.5,.5,1))
#pairs(d[,c('Ppreceding','CndP_preceding','Pfollowing','RC')],
#	labels=c('log P(w-1)','log P(RC)','log P(w+1)', 'REL present\nx Original'),
#	lower.panel=panel.smooth, upper.panel=panel.cor)



## --------------------------------------------------
# some functions to plot lmer logit outputs
## --------------------------------------------------
my.plotCoefficients = function(
	model,
	mcmc,
	nsim= if(!is.null(mcmc)) { 2000 } else { NULL },
	indices= 1:length(fixef(model)),
	labels= fixef(model)[indices], ...
) {
	if (is.null(mcmc)) {
		require(languageR)
		mcmc = pvals.fnc(model, nsim=nsim)
	}

	c <- as.numeric(fixef(model))[indices]
	lc <- as.vector(as.numeric(mcmc[["fixed"]]$HPD95lower))[indices]
	hc <- as.vector(as.numeric(mcmc[["fixed"]]$HPD95upper))[indices]
	s <-  seq(0.7, 0.7 + (length(indices) - 1) * 1.2, 1.2)

	par(mar=c(10,5.5,3.5,.5), cex= 1.3, cex.main=1.5, cex.lab=.9, cex.axis=.7)
	barplot(c,
	        ylim = c(min(c), max(c) + max(hc)),
        	axisnames = F,
	        las = 2, border=NA
	)
	text(0.35, max(c) + max(hc) * .6,
		paste("N=",nrow(model@frame),"\nk=",length(fixef(model))),
		adj=0,
		cex=1
	)
	axis(	side=1,
		at= s,
        	labels= labels,
	        tick=F,
	        cex.axis=1,
	        las=2
	)
	my.errbar(s, c, hc, lc, lwd=2, add=TRUE)
}

my.glmergplot <- function(
# version 0.43
# written by tiflo@csli.stanford.edu
# code contributions from Austin Frank, Ting Qian, and Harald Baayen
# remaining errors are mine (tiflo@csli.stanford.edu)
#
# last modified 12/15/10
#
# now also supports linear models
# backtransforms centering and standardization
#
# known bugs:
#   too simple treatment of random effects
#
	model,
	name.predictor,
	name.outcome= "outcome",
	predictor= NULL,

	# is the predictor centered IN THE MODEL?
	# is the predictor transformed before
	# (centered and) ENTERED INTO THE MODEL?
	predictor.centered= if(!is.null(predictor)) { T } else { F },
	predictor.standardized= F,
	predictor.transform= NULL,
	fun= NULL,

	type= "hex",
	main= NA,
	xlab= NA,
	ylab= NA,
	xlim= NA,
	ylim= NA,
	legend.position="right",
	fontsize=16,
	col.line= "#333333",
	col.ci= col.line,
	lwd.line= 1.2,
	lty.line= "solid",
	alpha.ci= 3/10,
	hex.mincnt= 1,
	hex.maxcnt= nrow(model@frame) / 2,
	hex.limits = c(round(hex.mincnt), round(hex.maxcnt)),
	hex.limits2 = c(round(match.fun(hex.trans)(hex.mincnt)), round(match.fun(hex.trans)(hex.maxcnt))),
	hex.midpoint = (max(hex.limits) - (min(hex.limits) - 1)) / 2,
	hex.nbreaks = min(5, round(match.fun(hex.trans)(max(hex.limits)) - match.fun(hex.trans)(min(hex.limits))) + 1),
	hex.breaks = round(seq(min(hex.limits), max(hex.limits), length.out=hex.nbreaks)),
	hex.trans = "log10",
	...
)
{
    	if (!is(model, "mer")) {
     		stop("argument should be a mer model object")
    	}
    	if ((length(grep("^glmer", as.character(model@call))) == 1) &
          (length(grep("binomial", as.character(model@call))) == 1)) {
		model.type = "binomial"
	} else {
		if (length(grep("^lmer", as.character(model@call))) == 1) {
			model.type = "gaussian"
		}
	}
	if (!(model.type %in% c("binomial","gaussian"))) {
		stop("argument should be a glmer binomial or gaussian model object")
	}
    	if (!is.na(name.outcome)) {
     	   	if (!is.character(name.outcome))
            	stop("name.outcome should be a string\n")
    	}
	if (!is.na(xlab[1])) {
        	if (!is.character(xlab))
            	stop("xlab should be a string\n")
    	}
    	if (!is.na(ylab)) {
     	   	if (!is.character(ylab))
            	stop("ylab should be a string\n")
    	}
	# load libaries
	require(lme4)
	require(rms)
	require(ggplot2)

	if (predictor.standardized) { predictor.centered = T }
    	if (is.null(fun)) {
		if (is.na(ylab)) {
		 	if (model.type == "binomial") { ylab= paste("Predicted log-odds of", name.outcome) }
			if (model.type == "gaussian") { ylab= paste("Predicted ", name.outcome) }
		}
		fun= I
	} else {
		if (is.na(ylab)) {
		 	if (model.type == "binomial") { ylab= paste("Predicted probability of", name.outcome) }
			if (model.type == "gaussian") { ylab= paste("Predicted ", name.outcome) }
		}
		fun= match.fun(fun)
	}
    	if (!is.null(predictor.transform)) {
		predictor.transform= match.fun(predictor.transform)
    	} else { predictor.transform= I }

	indexOfPredictor= which(names(model@fixef) == name.predictor)

	# get predictor
	if (is.null(predictor)) {
		# simply use values from model matrix X
 		predictor= model@X[,indexOfPredictor]

		# function for predictor transform
		fun.predictor= I

		if (is.na(xlab)) { xlab= name.predictor }
    	} else {
		# make sure that only defined cases are included
 	predictor = predictor[-na.action(model@frame)]
	# function for predictor transform
		trans.pred = predictor.transform(predictor)
		m= mean(trans.pred, na.rm=T)
		rms = sqrt(var(trans.pred, na.rm=T) / (sum(ifelse(is.na(trans.pred),0,1)) - 1))
		fun.predictor <- function(x) {
			x= predictor.transform(x)
			if (predictor.centered == T) { x= x - m }
			if (predictor.standardized == T) { x= x / rms }
			return(x)
		}

		if ((is.na(xlab)) & (label(predictor) != "")) { xlab= label(predictor) }
   	}

	# get outcome for binomial or gaussian model
	if (model.type == "binomial") {
		outcome= fun(qlogis(fitted(model)))
	} else {
		outcome= fun(fitted(model))
	}


     ## calculate grand average but exclude effect to be modeled
     ## (otherwise it will be added in twice!)
     ## random effects are all included, even those for predictor (if any).
     ## should random slope terms for the predictor be excluded?

     ## prediction from fixed effects
	if (ncol(model@X) > 2) {
	     Xbeta.hat = model@X[, -indexOfPredictor] %*% model@fixef[-indexOfPredictor]
	} else {
		Xbeta.hat = model@X[, -indexOfPredictor] %*% t(model@fixef[-indexOfPredictor])
	}

     ## adjustment from random effects
     Zb = crossprod(model@Zt, model@ranef)@x

     ## predicted value using fixed and random effects
     Y.hat = Xbeta.hat + Zb

     ## intercept is grand mean of predicted values
     ## (excluding fixed effect of predictor)
     ## (including random effects of predictor, if any)
     int = mean(Y.hat)

	# slope
	slope <- fixef(model)[name.predictor]

	## error and confidence intervals
	stderr <- sqrt(diag(vcov(model)))
	names(stderr) <- names(fixef(model))
	slope.se <- stderr[name.predictor]
	lower <- -1.96 * slope.se
	upper <- 1.96 * slope.se

	# setting graphical parameters
	if (is.na(ylim)) { ylim= c(min(outcome) - 0.05 * (max(outcome) - min(outcome)), max(outcome) + 0.05 * (max(outcome) - min(outcome)) ) }
   	if (is.na(xlim)) { xlim= c(min(predictor) - 0.05 * (max(predictor) - min(predictor)), max(predictor) + - 0.05 * (max(predictor) - min(predictor))) }

	print("Printing with ...")
	print(paste("   int=", int))
	print(paste("   slope=", slope))
	print(paste("   centered=", predictor.centered))
	print("   fun:")
	print(fun.predictor)

	pdata= data.frame( 	predictor=predictor, outcome=outcome	)
	x= seq(xlim[1], xlim[2], length=1000)
	fit= int + slope * fun.predictor(x)
	ldata= data.frame(
				predictor= x,
				outcome= fun(fit),
				transformed.lower= fun(fit + lower),
				transformed.upper= fun(fit + upper)
	)
	theme_set(theme_grey(base_size=fontsize))
	theme_update(axis.title.y=theme_text(angle=90, face="bold", size=fontsize, hjust=.5, vjust=.5))
	theme_update(axis.title.x=theme_text(angle=0, face="bold", size=fontsize, hjust=.5, vjust=.5))
	p <- ggplot(data=pdata, aes(x=predictor, y=outcome)) +
		xlab(xlab) +
		ylab(ylab) +
		xlim(xlim) +
		ylim(ylim) +
		opts(legend.position=legend.position, aspect.ratio=1)

	     	# for degbugging:
		# panel.lines(rep(mean(x),2), c(min(y),max(y)))
		# panel.lines(c(min(x),max(x)), c(mean(y),mean(y)))

	if (type == "points") {
		p <- p + geom_point(alpha=3/10)
	} else if (type == "hex") {
		p <- p + geom_hex(bins = 30) +
			scale_fill_gradient2(low= "lightyellow",
							mid="orange",
							high=muted("red"),
							midpoint= hex.midpoint,
							space="rgb",
							name= "Count",
							limits= hex.limits,
							breaks= hex.breaks,
							trans = hex.trans
							)
	}
	p + 	geom_ribbon(data=ldata,
			aes(	x= predictor,
				ymin=transformed.lower,
				ymax=transformed.upper
			),
			fill= col.ci,
			alpha= alpha.ci
		) +
		geom_line(data=ldata,
			aes(x= predictor,
				y=outcome
			),
			colour= col.line,
			size= lwd.line,
			linetype= lty.line,
			alpha=1
		)
}

my.rcs.glmergplot <- function(
# version 0.3
# written by tiflo@csli.stanford.edu
# based on code from Austin and Ting
# last modified 12/30/09
#
# known bugs:
#   too simple treatment of random effects
#   crashes if model has only one fixed effect (needs an if-statement below)
	model,
	name.predictor,
	name.outcome= "outcome",
	predictor= NULL,

	# is the predictor centered IN THE MODEL?
	# is the predictor transformed before
	# (centered and) ENTERED INTO THE MODEL?
	predictor.centered= if(!is.null(predictor)) { T } else { F },
	fun= NULL,

	type= "hex",
	main= NA,
	xlab= NA,
	ylab= NA,
	xlim= NULL,
	ylim= NULL,
	legend.position="right",
	fontsize=16,
	col.line= "#333333",
	lwd.line= 1.2,
	lty.line= "solid",
	col.ci= col.line,
	alpha.ci = 3/10,
	mincnt= 1,

	...
)
{
    	if (!is(model, "mer")) {
     		stop("argument should be a mer model object")
    	}
    	if ((length(grep("^glmer", as.character(model@call))) != 1) |
         (length(grep("binomial", as.character(model@call))) != 1)) {
		stop("argument should be a glmer binomial model object")
	}
    	if (!is.na(name.outcome)) {
     	   	if (!is.character(name.outcome))
            	stop("name.outcome should be a string\n")
    	}
	if (!is.na(xlab[1])) {
        	if (!is.character(xlab))
            	stop("xlab should be a string\n")
    	}
    	if (!is.na(ylab)) {
     	   	if (!is.character(ylab))
            	stop("ylab should be a string\n")
    	}
    	if (is.null(fun)) {
		if (is.na(ylab)) { ylab= paste("Predicted log-odds of", name.outcome) }
		fun= I
	} else {
		fun= match.fun(fun)
		if (is.na(ylab)) { ylab= paste("Predicted probability of", name.outcome) }
	}

	# load libaries
	require(lme4)
	require(ggplot2)

	# determine which predictor is to be modeled
	lab.predictor= sub("^rcs\\(([A-Za-z][^,]*)\\)$", "\\1", name.predictor)

	# the following will be problematic for interactions with rcs
	name.predictor= names(model@fixef)[grep(name.predictor, names(model@fixef))]
	indexOfPredictor= which(names(model@fixef) %in% name.predictor)

	# get predictor
	if (is.null(predictor)) {
		# simply use values from model matrix X
 		predictor= model@X[,indexOfPredictor]

		# function for predictor transform
		fun.predictor= I

		if (is.na(xlab)) { xlab= lab.predictor }
    	} else {
		# function for predictor transform
		m= mean(predictor)
		fun.predictor <- function(x) {
			if (predictor.centered == T) { x= x - m }
			x= rcs(x)
			return(x)
		}

		if (is.na(xlab)) { xlab= label(predictor) }
   	}

	# get outcome
	outcome= fun(qlogis(fitted(model)))


     ## calculate grand average but exclude effect to be modeled
     ## (otherwise it will be added in twice!)
     ## random effects are all included, even those for predictor (if any).
     ## should random slope terms for the predictor be excluded?

     ## prediction from fixed effects
     Xbeta.hat = model@X[, -indexOfPredictor] %*% model@fixef[-indexOfPredictor]

     ## adjustment from random effects
     Zb = crossprod(model@Zt, model@ranef)@x

     ## predicted value using fixed and random effects
     Y.hat = Xbeta.hat + Zb

     ## intercept is grand mean of predicted values
     ## (excluding fixed effect of predictor)
     ## (including random effects of predictor, if any)
     int = mean(Y.hat)

	# slope
	slope <- fixef(model)[name.predictor]

	## error and confidence intervals
	stderr <- sqrt(diag(vcov(model)))
	names(stderr) <- names(fixef(model))
	slope.se <- stderr[name.predictor]
	lower <- -1.96 * slope.se
	upper <- 1.96 * slope.se

	# setting graphical parameters
	if (is.null(ylim)) { ylim= c(min(outcome) - 0.05 * (max(outcome) - min(outcome)), max(outcome) + 0.05 * (max(outcome) - min(outcome)) ) }
   	if (is.null(xlim)) { xlim= c(min(predictor) - 0.05 * (max(predictor) - min(predictor)), max(predictor) + 0.05 * (max(predictor) - min(predictor))) }

	print("Printing with ...")
	print(paste("   int=", int))
	print(paste("   slope=", slope))
	print(paste("   centered=", predictor.centered))
	print("   fun:")
	print(fun.predictor)

	pdata= data.frame( 	predictor=predictor, outcome=outcome	)
	x= sort(predictor)
	ldata= data.frame(
				predictor= x,
				outcome= as.numeric(fun(int + fun.predictor(x) %*% slope)),
				transformed.lower= as.numeric(fun(int + fun.predictor(x) %*% (slope + lower))),
				transformed.upper= as.numeric(fun(int + fun.predictor(x) %*% (slope + upper)))
	)
	ldata$transformed.lower <- ifelse(ldata$transformed.lower < ylim[1], ylim[1], ldata$transformed.lower)
	ldata$transformed.upper <- ifelse(ldata$transformed.upper > ylim[2], ylim[2], ldata$transformed.upper)
	theme_set(theme_grey(base_size=fontsize))
	theme_update(axis.title.y=theme_text(angle=90, face="bold", size=fontsize, hjust=.5, vjust=.5))
	theme_update(axis.title.x=theme_text(angle=0, face="bold", size=fontsize, hjust=.5, vjust=.5))
	p <- ggplot(data=pdata, aes(x=predictor, y=outcome)) +
		xlab(xlab) +
		ylab(ylab) +
		xlim(xlim) +
		ylim(ylim) +
		opts(legend.position=legend.position, aspect.ratio=1)

	     	# for degbugging:
		# panel.lines(rep(mean(x),2), c(min(y),max(y)))
		# panel.lines(c(min(x),max(x)), c(mean(y),mean(y)))

	if (type == "points") {
		p <- p + geom_point(alpha=3/10)
	} else if (type == "hex") {
		p <- p + geom_hex(bins = 30) +
			scale_fill_gradient2(low= "lightyellow", mid="orange", high=muted("red"), midpoint= 1.5, space="rgb", name= "Count", limits=c(log(mincnt, base=10),3), breaks=c(0,1,2,3), trans="log10")
	}
	if (is.null(col.ci)) {
		p + 	geom_line(data=ldata,
			aes(x= predictor,
				y=outcome
			),
			colour= col.line,
			size= lwd.line,
			linetype= lty.line,
			alpha=1
		)
	} else {
		p + 	geom_ribbon(data=ldata,
			aes(	x= predictor,
				ymin=transformed.lower,
				ymax=transformed.upper
			),
			fill= col.ci,
			alpha=3/10
		) +
		geom_line(data=ldata,
			aes(x= predictor,
				y=outcome
			),
			colour= col.line,
			size= lwd.line,
			linetype= lty.line,
			alpha=1
		)
	}
}

my.glmerplot <- function(
# version 0.22
# written by tiflo@csli.stanford.edu
# based on code from Austin
# last modified 02/25/09
#
# known bugs:
#   too simple treatment of random effects
#   crashes if model has only one fixed effect (needs an if-statement below)
	model,
	name.predictor,
	name.outcome= "outcome",
	predictor= NULL,

	# is the predictor centered IN THE MODEL?
	# is the predictor transformed before
	# (centered and) ENTERED INTO THE MODEL?
	predictor.centered= if(!is.null(predictor)) { T } else { F },
	predictor.transform= NULL,
	fun= NULL,
	main= NA,
	xlab= NA,
	ylab= NA,
	xlim= NA,
	ylim= NA,
	lty.line= "solid",
	lty.int= "dotted",
	col.line= "orange",
	col.int= col.line,
	lwd.line= 6,
	lwd.int=  3,
	mincnt= 1,

	# BTC means "blue to cyan"... defines the color
     # gradient that is used for the hexes
	colramp = function (n) { BTC(n, beg=225, end=1) },
	...
)
{
	# based on code by Austin Frank
	# incorporates code from plotLMER.fnc (languageR)
    	if (!is(model, "mer")) {
     		stop("argument should be a mer model object")
    	}
    	if ((length(grep("^glmer", as.character(model@call))) != 1) |
         (length(grep("binomial", as.character(model@call))) != 1)) {
		stop("argument should be a glmer binomial model object")
	}
    	if (!is.na(name.outcome)) {
     	   	if (!is.character(name.outcome))
            	stop("name.outcome should be a string\n")
    	}
	if (!is.na(xlab[1])) {
        	if (!is.character(xlab))
            	stop("xlab should be a string\n")
    	}
    	if (!is.na(ylab)) {
     	   	if (!is.character(ylab))
            	stop("ylab should be a string\n")
    	}
    	if (is.null(fun)) {
		if (is.na(ylab)) { ylab= paste("Predicted log-odds of", name.outcome) }
		fun= I
	} else {
		fun= match.fun(fun)
		if (is.na(ylab)) { ylab= paste("Predicted probability of", name.outcome) }
	}
    	if (!is.null(predictor.transform)) {
		predictor.transform= match.fun(predictor.transform)
    	} else { predictor.transform= I }

	# load libaries
	library(lme4)
	library(latticeExtra)
	library(hexbin)

	# determine which predictor is to be modeled
	indexOfPredictor= which(names(model@fixef) == name.predictor)

	# get predictor
	if (is.null(predictor)) {
		# simply use values from model matrix X
 		predictor= model@X[,indexOfPredictor]

		# function for predictor transform
		fun.predictor= I

		if (is.na(xlab)) { xlab= name.predictor }
    	} else {
		# function for predictor transform
		m= mean(predictor.transform(predictor))
		fun.predictor <- function(x) {
			x= predictor.transform(x)
			if (predictor.centered == T) { x= x - m }
			return(x)
		}

		if (is.na(xlab)) { xlab= label(predictor) }
   	}

	# get outcome
	outcome= fun(qlogis(fitted(model)))


     ## calculate grand average but exclude effect to be modeled
     ## (otherwise it will be added in twice!)
     ## random effects are all included, even those for predictor (if any).
     ## should random slope terms for the predictor be excluded?

     ## prediction from fixed effects
     Xbeta.hat = model@X[, -indexOfPredictor] %*% model@fixef[-indexOfPredictor]

     ## adjustment from random effects
     Zb = crossprod(model@Zt, model@ranef)@x

     ## predicted value using fixed and random effects
     Y.hat = Xbeta.hat + Zb

     ## intercept is grand mean of predicted values
     ## (excluding fixed effect of predictor)
     ## (including random effects of predictor, if any)
     int = mean(Y.hat)

	# slope
	slope <- fixef(model)[name.predictor]

	## error and confidence intervals
	stderr <- sqrt(diag(vcov(model)))
	names(stderr) <- names(fixef(model))
	slope.se <- stderr[name.predictor]
	lower <- -1.96 * slope.se
	upper <- 1.96 * slope.se

	# setting graphical parameters
	if (is.na(main)) { par(mar=c(3.1,4.1,0.1,2.1)) }
	lattice.options(default.args = list(as.table=TRUE))
	trellis.par.set(fontsize = list(text = 18), font = list(lab = 2) )
	if (is.na(ylim)) { ylim= c(min(outcome) - 0.05 * (max(outcome) - min(outcome)), max(outcome) + 0.05 * (max(outcome) - min(outcome)) ) }
   	if (is.na(xlim)) { xlim= c(min(predictor), max(predictor)) }

	print("Printing with ...")
	print(paste("   int=", int))
	print(paste("   slope=", slope))
	print(paste("   centered=", predictor.centered))
	print("   fun:")
	print(fun.predictor)

	hexbinplot(outcome ~ predictor,
		# define some extra variables to pass into the panel
          	int = int, slope = slope, lower = lower, upper = upper,
		fun= fun,

		# setting up general aspects of the plot
          	main = main,
		xlab = xlab,
		ylab = ylab,
		xlim = xlim,
         	ylim = ylim,
		aspect = 1,
         	colramp = colramp,
		mincnt = mincnt,

	   	panel = function (x, y, ...) {
          		panel.hexbinplot(x, y, ...)

                # actually plot the curves
                #   fun: determines whether graph is plotted in probability or logit space
                #   fun.predictor: redoes any transform function (e.g. log) and centering, if necessary.

			# for degbugging:
			# panel.lines(rep(mean(x),2), c(min(y),max(y)))
			# panel.lines(c(min(x),max(x)), c(mean(y),mean(y)))

                panel.curve(fun(int + slope * fun.predictor(x) + lower), lty = lty.int, col.line = col.int, lwd = lwd.int)
                panel.curve(fun(int + slope * fun.predictor(x) + upper), lty = lty.int, col.line = col.int, lwd = lwd.int)
                panel.curve(fun(int + slope * fun.predictor(x)), lty= lty.line, col.line = col.line, lwd = lwd.line)
          	}
	)
}


my.errbar <- function (x, y, yplus, yminus, cap = 0.025, xlab = as.character(substitute(x)),
    ylab = if (is.factor(x) || is.character(x)) "" else as.character(substitute(y)),
    add = FALSE, lty = 1, ylim, lwd = 1, Type = rep(1, length(y)), pch= NA_integer_,
    ...)
{
    if (missing(ylim))
        ylim <- range(y[Type == 1], yplus[Type == 1], yminus[Type ==
            1], na.rm = TRUE)
    if (is.factor(x) || is.character(x)) {
        x <- as.character(x)
        n <- length(x)
        t1 <- Type == 1
        t2 <- Type == 2
        n1 <- sum(t1)
        n2 <- sum(t2)
        omai <- par("mai")
        mai <- omai
        mai[2] <- max(strwidth(x, "inches")) + 0.25 * .R.
        par(mai = mai)
        on.exit(par(mai = omai))
        plot(0, 0, xlab = ylab, ylab = "", xlim = ylim, ylim = c(1,
            n + 1), axes = FALSE, ...)
        axis(1)
        w <- if (any(t2))
            n1 + (1:n2) + 1
        else numeric(0)
        axis(2, at = c(1:n1, w), labels = c(x[t1], x[t2]), las = 1,
            adj = 1)
        points(y[t1], 1:n1, pch = pch, ...)
        segments(yplus[t1], 1:n1, yminus[t1], 1:n1, ...)
        if (any(Type == 2)) {
            abline(h = n1 + 1, lty = 2, ...)
            offset <- mean(y[t1]) - mean(y[t2])
            if (min(yminus[t2]) < 0 & max(yplus[t2]) > 0)
                lines(c(0, 0) + offset, c(n1 + 1, par("usr")[4]),
                  lty = 2, ...)
            points(y[t2] + offset, w, pch = pch, ...)
            segments(yminus[t2] + offset, w, yplus[t2] + offset,
                w, ...)
            at <- pretty(range(y[t2], yplus[t2], yminus[t2]))
            axis(3, at = at + offset, label = format(round(at,
                6)))
        }
        return(invisible())
    }
    if (add)
        points(x, y, pch= pch, ...)
    else plot(x, y, ylim = ylim, xlab = xlab, ylab = ylab, ...)
    xcoord <- par()$usr[1:2]
    segments(x, yminus, x, yplus, lty = lty, lwd = lwd, ...)
    smidge <- cap * (xcoord[2] - xcoord[1])/2
    segments(x - smidge, yminus, x + smidge, yminus, lwd = lwd,
        ...)
    segments(x - smidge, yplus, x + smidge, yplus, lwd = lwd,
        ...)
    invisible()
}




my.plot.logistic.fit <- function (x, data,
	method = "cut",
	type = "points",
	lty = 2,
	where = seq(0, 1, by = 0.1),
	line.col= "grey", line.width=1.5, line.lty = 1,
	rug= T, rug.col= "#555555", rug.stacked = T, rug.y = 1.04,
	legend.cex= 1, legend.font= 2, legend.col= "black",
	print.minbinsize= T, MinBin = NA,
   	...)
{
    	require(lme4, quietly = TRUE)
   	require(rms, quietly = TRUE)
    	if (class(x)[1] == "glmer" | class(x)[1] == "mer") {
		data=data[-na.action(x@frame),]
        	depvar = as.character(formula(attr(x@frame, "terms")))[2]
        	probs = fitted(x)
    	}
	else {
        	if (class(x)[1] == "lrm") {
            	depvar = as.character(formula(x$call))[2]
            	probs = predict(x, type = "fitted")
        	}
        	else {
            	stop("first argument is not an lmer or lrm model")
      	}
    	}
    	if (method == "cut") {
        	classes = cut2(probs, where, levels.mean = TRUE)
		normdepvar = (as.numeric(data[, depvar]) - min(as.numeric(data[, depvar]))) / (max(as.numeric(data[, depvar])) - min(as.numeric(data[, depvar])))
		if (length(classes) != length(normdepvar)) {
			cat("Predicted values (fitted) and dependent variable must be of equal length.\n")
			cat("Probably the data frame has more cases than used in the model (e.g. because\n")
   			cat("of missing values for predictors or the dependent variable.\n")
		}
        	means = tapply(normdepvar, classes,
            		mean)
        	lengths = tapply(normdepvar, classes,
            		length)
    	}
    	else {
		if (method == "shingle") {
            	sh = equal.count(probs)
            	means = rep(0, length(levels(sh)))
            	midpoints = rep(0, length(means))
			lengths = rep(NA, length(levels(sh)))
            	for (i in 1:length(levels(sh))) {
                	means[i] = mean(probs[probs > levels(sh)[[i]][1] &
                  	probs < levels(sh)[[i]][2]])
                	midpoints[i] = as.character(mean(levels(sh)[[i]]))
				lengths[i] = length(probs[probs > levels(sh)[[i]][1] &
                  	probs < levels(sh)[[i]][2]])
            	}
            	names(means) = as.character(midpoints)
        	}
   	}
    	plot(as.numeric(names(means)), means, xlab = "Mean predicted probabilities",
        	ylab = "Observed proportions", type = "n", ...)
    	if (rug == T) {
		myRug(x= probs, y= rug.y, ticksize= 0.02, col = rug.col, stacked = rug.stacked)
	}
    	abline(0, 1, col = line.col, lwd = line.width, lty = line.lty)


	if (is.na(MinBin)) { minbin = min(lengths)
	} else { minbin = MinBin }

    	if (type == "sunflower") {
	    	sunflowerplot(as.numeric(names(means)), means, number= (lengths / minbin), add=T, ...)
	}
	else {
		if (type == "squares") {
			# to print squares sized by amount of data
			points(as.numeric(names(means)), means, pch = 46, cex = lengths / minbin * 2, ...)
		}
		else {
    			# to print standard points
			points(as.numeric(names(means)), means, pch = 19, cex = 1, ...)
		}
	}

#	lines(smooth(means ~ as.numeric(names(means))), lty = 3)
#	lines(loess(means ~ as.numeric(names(means)), weights=lengths, span=0.05), lty = 3)

    	if(print.minbinsize == T) {
		text(locator(1), pos= 4, paste(
			paste("R-squared: ", round(cor(as.numeric(names(means)), means)^2, 2), sep = ""),
			paste("Min. bin size: ", minbin, sep= ""), sep="\n"),
			cex= legend.cex, font= legend.font, col= legend.col)
	}
	else {
		text(locator(1), pos= 4,
			paste("R-squared: ", round(cor(as.numeric(names(means)), means)^2, 2), sep = ""),
			cex= legend.cex, font= legend.font, col= legend.col)
	}
}









## --------------------------------------------------
my.plot.lmer <- function (model, xlabel = NA, ylabel = NA, ylimit = NA, fun = NA,
    pred = NA, n = 100, intr = NA, mcmcMat = NA, lockYlim = TRUE,
    addlines = FALSE, withList = FALSE, cexsize = 0.5, ...)
{
    library(languageR)
    if (!(is(model, "lmer") | is(model, "glmer") | is(model,
        "mer"))) {
        stop("argument should be an lmer, glmer or mer model object")
    }
    if (!is.na(xlabel)) {
        if (!is.character(xlabel))
            stop("xlabel should be a string\n")
    }
    if (!is.na(ylabel)) {
        if (!is.character(ylabel))
            stop("ylabel should be a string\n")
    }
    if (!is.na(ylimit[1])) {
        if ((!is.numeric(ylimit)) | (length(ylimit) != 2))
            stop("ylimit should be a two-element numeric vector\n")
    }
    if (!is.na(mcmcMat[1])) {
        if (!is.matrix(mcmcMat))
            stop("mcmcMat should be a matrix\n")
    }
    if (!is.na(intr[1])) {
        if (!is.list(intr))
            stop("intr should be a list\n")
    }
    if (!is.numeric(n)) {
        stop("n should be an integer\n")
    }
    if (!is.na(pred)) {
        if (!is.character(pred))
            stop("pred should be a string\n")
    }
    if (!is.function(fun)) {
        if (!is.na(fun)) {
            stop("fun should be a function (not the name of a function)\n")
        }
    }
    else {
        if (is(model, "glmer"))
            fun = plogis
    }
    conditioningPred = ""
    conditioningVals = NULL
    conditioningPos = NA
    conditioningColors = 1
    conditioningLines = 1
    if (!is.na(intr[[1]])) {
        conditioningPred = intr[[1]]
        conditioningVals = intr[[2]]
        conditioningPos = intr[[3]]
        if (length(intr) == 4) {
            conditioningColors = intr[[4]][[1]]
            if (length(conditioningColors) != length(conditioningVals)) {
                stop("number of colors and number of conditioning values mismatch")
            }
            conditioningLines = intr[[4]][[2]]
            if (length(conditioningLines) != length(conditioningLines)) {
                stop("number of line types and number of conditioning values mismatch")
            }
        }
    }
    if (length(ylimit) > 1) {
        lockYlim = FALSE
    }
    if (is.na(ylabel))
        ylabel = as.character(model@call[2]$formula)[2]
    if (is.na(pred)) {
        predictors = colnames(model@frame)
        ranefnames = unique(names(ranef(model)))
        depvar = as.character(model@call[2]$formula)[2]
        predictors = predictors[1:(which(predictors == ranefnames[1]) -
            1)]
        predictors = predictors[!predictors %in% c(ranefnames,
            depvar)]
    }
    else {
        predictors = pred
    }
    plots = list()
    for (i in 1:length(predictors)) {
        cat("preparing panel for", predictors[i], "\n")
        if (is.na(xlabel) | length(predictors) > 1)
            xlabel = predictors[i]
        if ((length(predictors) == 1) & (!is.null(conditioningVals))) {
            if (is.null(conditioningColors)) {
                colors = rep(1, length(conditioningVals))
                lineTypes = rep(1, length(conditioningVals))
            }
            else {
                colors = conditioningColors
                lineTypes = conditioningLines
                if (length(colors) < length(conditioningVals)) {
                  nc = (length(conditioningVals)%%length(colors)) +
                    1
                  colors = rep(colors, nc)
                }
                if (length(lineTypes) < length(conditioningVals)) {
                  nc = (length(conditioningLines)%%length(lineTypes)) +
                    1
                  lineTypes = rep(lineTypes, nc)
                }
            }
            val = conditioningVals[1]
            cat(" preparing interaction for ", conditioningPred,
                " = ", val, "\n")
            m = makeDefaultMatrix.fnc(model, n, conditioningPred,
                val)
            subplots = list()
            dfr = preparePredictor.fnc(predictors[i], model,
                m, ylabel, fun, val, xlabel = xlabel, mcmc = mcmcMat,
                lty = 1, col = 0, ...)
            subplots[[1]] = dfr
            for (j in 2:length(conditioningVals)) {
                val = conditioningVals[j]
                cat(" preparing interaction for ", conditioningPred,
                  " = ", val, "\n")
                m = makeDefaultMatrix.fnc(model, n, conditioningPred,
                  val)
                dfr = preparePredictor.fnc(predictors[i], model,
                  m, ylabel, fun, val, mcmc = mcmcMat, lty = j,
                  xlabel = xlabel, ...)
                subplots[[j]] = dfr
            }
            plots[[i]] = subplots
        }
        else {
            lineTypes = 1
            m = makeDefaultMatrix.fnc(model, n, "", NULL)
            dfr = preparePredictor.fnc(predictors[i], model,
                m, ylabel, fun, val = NA, xlabel = xlabel, mcmc = mcmcMat,
                ...)
            plots[[i]] = dfr
        }
    }
    names(plots) = predictors
    plotAll.fnc(plots, sameYrange = lockYlim, ylabel, intrName = conditioningPred,
        pos = conditioningPos, ylimit = ylimit, addlines = addlines,
        cexsize = cexsize, conditioningVals = conditioningVals,
        conditioningColors = colors, conditioningLines = lineTypes)
    if (withList)
        return(plots)
}






## --------------------------------------------------
my.pvals <- function (object, nsim = 10000, ndigits = 4, withMCMC = FALSE,
    addPlot = TRUE, useOldVersion = TRUE, ...)
{
    require("lme4", quietly = TRUE, character = TRUE)
    if (useOldVersion) {
        coefs = summary(object)@coefs
        ncoef = length(coefs[, 1])
        if (nsim > 0) {
            mcmc = mcmcsamp(object, n = nsim)
            hpd = lme4::HPDinterval(mcmc)
#            sumry = summary(mcmc)$statistics
            nr <- nrow(mcmc)
            prop <- colSums(mcmc[, 1:ncoef] > 0)/nr
            ans <- 2 * pmax(0.5/nr, pmin(prop, 1 - prop))
            fixed = data.frame(Estimate = round(as.numeric(coefs[,
                1]), ndigits), MCMCmean = round(apply(t(mcmc@fixef),
                2, mean), ndigits), HPD95lower = round(hpd[1:ncoef, 1],
                ndigits), HPD95upper = round(hpd[1:ncoef, 2],
                ndigits), pMCMC = round(ans, ndigits), pT = round(2 *
                (1 - pt(abs(coefs[, 3]), nrow(object@frame) -
                  ncoef)), ndigits), row.names = names(coefs[,
                1]))
            if (class(object) == "lmer")
                colnames(fixed)[ncol(fixed)] = "Pr(>|t|)"
            else colnames(fixed)[ncol(fixed)] = "Pr(>|z|)"
            v = (ncoef + 1):nrow(hpd)
            random = data.frame(round(apply(t(mcmc@ranef),
                2, mean), ndigits), HPD95lower = hpd[v,
                1], HPD95upper = hpd[v, 2])
            nms = rownames(random)
            logs = substr(rownames(random), 1, 3)
            rows = logs == "log"
            random[rows, ] = sqrt(exp(random[rows, ]))
            nms[rows] = substr(nms[rows], 5, nchar(nms) - 1)
            nms[1] = substr(nms[1], 1, nchar(nms[1]) - 2)
            atanhs = substr(rownames(random), 1, 5)
            rows = atanhs == "atanh"
            if (sum(rows) > 0) {
                random[rows, ] = tanh(random[rows, ])
                nms[rows] = substr(nms[rows], 7, nchar(nms[rows]) -
                  1)
            }
            rownames(random) = nms
            if (addPlot)
                print(densityplot(mcmc, par.strip.text = list(cex = 0.75)))
            if (withMCMC) {
                mcmc = as.matrix(mcmc)
                return(list(fixed = format(fixed, digits = ndigits,
                  sci = FALSE), random = format(random, digits = ndigits,
                  sci = FALSE), mcmc = mcmc))
            }
            else return(list(fixed = format(fixed, digits = ndigits,
                sci = FALSE), random = format(random, digits = ndigits,
                sci = FALSE)))
        }
        else {
            fixed = data.frame(Estimate = as.numeric(coefs[,
                1]), pT = round(2 * (1 - pt(abs(coefs[, 3]),
                nrow(object@frame) - ncoef)), ndigits), row.names = names(coefs[,
                1]))
            if (class(object) == "lmer")
                colnames(fixed)[ncol(fixed)] = "Pr(>|t|)"
            else colnames(fixed)[ncol(fixed)] = "Pr(>|z|)"
            return(fixed)
        }
    }
    else {
        if (is(object, "mer")) {
		 mcmc = mcmcsamp(object, n = nsim)
          	 hpd = lme4::HPDinterval(mcmc)
            coefs = summary(object)@coefs
            ncoef = length(coefs[, 1])
            sgma = summary(object)@sigma
            if (nsim > 0) {
                if (colnames(coefs)[3] == "z value") {
                  stop("mcmc sampling is not yet implemented for generalized mixed models\n")
                }
                mcmcfixef = t(mcmc@fixef)
                nr <- nrow(mcmcfixef)
                prop <- colSums(mcmcfixef > 0)/nr
                ans <- 2 * pmax(0.5/nr, pmin(prop, 1 - prop))
                fixed = data.frame(Estimate = round(as.numeric(coefs[,
                  1]), ndigits), MCMCmean = round(apply(t(mcmc@fixef),
                  2, mean), ndigits), HPD95lower = round(hpd$fixef[,
                  1], ndigits), HPD95upper = round(hpd$fixef[,
                  2], ndigits), pMCMC = round(ans, ndigits),
                  pT = round(2 * (1 - pt(abs(coefs[, 3]), nrow(object@frame) -
                    ncoef)), ndigits), row.names = names(coefs[,
                    1]))
                colnames(fixed)[ncol(fixed)] = "Pr(>|t|)"
                isStdev = vector()
                theValues = vector()
                names(object@ST) = names(object@flist)
                for (i in 1:length(object@ST)) {
                  m = object@ST[[i]]
                  if (i == 1) {
                    nms = rownames(m)
                    nms = paste(rep(names(object@ST)[i], length(nms)),
                      nms, sep = " ")
                    for (j in 1:nrow(m)) {
                      isStdev = c(isStdev, TRUE)
                      theValues = c(theValues, m[j, 1])
                    }
                  }
                  else {
                    tmp = rownames(m)
                    tmp = paste(rep(names(object@ST)[i], length(tmp)),
                      tmp, sep = " ")
                    nms = c(nms, tmp)
                    for (j in 1:nrow(m)) {
                      isStdev = c(isStdev, TRUE)
                      theValues = c(theValues, m[j, 1])
                    }
                  }
                  if (nrow(m) > 1) {
                    for (j in 2:ncol(m)) {
                      for (k in 1:nrow(m)) {
                        if (m[k, j] != 0) {
                          nms = c(nms, paste(names(object@ST)[i],
                            colnames(m)[j - 1], rownames(m)[k],
                            sep = " "))
                          isStdev = c(isStdev, FALSE)
                          theValues = c(theValues, m[k, j])
                        }
                      }
                    }
                  }
                }
                rowSigma = round(c(mean(mcmc@sigma), median(mcmc@sigma),
                  as.vector(hpd$sigma)), ndigits)
                mcmcHPD = as.data.frame(lme4::HPDinterval(mcmc)$ST)
                rownames(mcmcHPD) = nms
                means = apply(mcmc@ST, 1, mean)
                medians = apply(mcmc@ST, 1, median)
                rowsST = round(cbind(means, medians, mcmcHPD),
                  ndigits)
                rownames(rowsST) = nms
                random = data.frame(rbind(rowsST, rowSigma))
                rownames(random)[nrow(random)] = "Sigma"
                colnames(random) = c("MCMCmean", "MCMCmedian",
                  "HPD95lower", "HPD95upper")
                mcmcM = as.matrix(mcmc)
                if (addPlot) {
                  beg = ncol(mcmcM) - length(nms)
                  end = ncol(mcmcM)
                  colnames(mcmcM)[beg:end] = c(nms, "sigma")
                  m = data.frame(Value = mcmcM[, 1], Predictor = rep(colnames(mcmcM)[1],
                    nrow(mcmcM)))
                  for (i in 2:ncol(mcmcM)) {
                    mtmp = data.frame(Value = mcmcM[, i], Predictor = rep(colnames(mcmcM)[i],
                      nrow(mcmcM)))
                    m = rbind(m, mtmp)
                  }
                  print(densityplot(~Value | Predictor, data = m,
                    scales = list(relation = "free"), par.strip.text = list(cex = 0.75),
                    xlab = "Posterior Values", ylab = "Density",
                    pch = "."))
                  if (withMCMC) {
                    return(list(fixed = format(fixed, digits = ndigits,
                      sci = FALSE), random = format(random, digits = ndigits,
                      sci = FALSE), mcmc = as.data.frame(mcmcM)))
                  }
                  else {
                    return(list(fixed = format(fixed, digits = ndigits,
                      sci = FALSE), random = format(random, digits = ndigits,
                      sci = FALSE)))
                  }
                }
                else {
                  return(list(fixed = format(fixed, digits = ndigits,
                    sci = FALSE), random = format(random, digits = ndigits,
                    sci = FALSE)))
                }
            }
            else {
                coefs = summary(object)@coefs
                ncoef = length(coefs[, 1])
                fixed = data.frame(Estimate = round(as.numeric(coefs[,
                  1]), ndigits), pT = round(2 * (1 - pt(abs(coefs[,
                  3]), nrow(object@frame) - ncoef)), ndigits),
                  row.names = names(coefs[, 1]))
                colnames(fixed)[ncol(fixed)] = "Pr(>|t|)"
                return(list(fixed = format(fixed, digits = ndigits,
                  sci = FALSE)))
            }
        }
        else {
            cat("the input model is not an lmer, glmer or mer object\n")
            return()
        }
    }
}
