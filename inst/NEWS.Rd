\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{broom.mixed News}
\encoding{UTF-8}

\section{CHANGES IN VERSION 0.2.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item more consistent term names in \code{brmsfit} tidiers
      \item improved \code{tidy.MCMCglmm}
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item all methods return tibbles (\code{tbl_df}) rather than data
      frames
      \item the value of the \code{group} variable for fixed-effect
      parameters has changed from \code{"fixed"} to \code{NA}
      \item \code{brmsfit} and \code{rstanarm} tidiers are more
      consistent with other tidiers (e.g. \code{conf.level} rather than
      \code{prob})
    }
  }
}
\section{CHANGES IN VERSION 0.0.1}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Sorted out some of the confusion over random effect naming:
      "ran_vals" extracts conditional modes/BLUPs/varying parameters
      (deviations from population-level estimates),
      while "ran_coefs" extracts group-level estimates
    }
  }
  \subsection{NEW FEATURES}{
     \itemize{
       \item improved \code{nlme} tidiers
       \item improved \code{glmmTMB} tidiers (can handle some
       zero-inflation parameters)
       \item \code{lme4} tidiers now optionally take a pre-computed
       \code{profile} argument when using \code{conf.method="profile"}

     }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item The default behaviour of most mixed-model tidiers has
      changed/will gradually be changed to the following
      (description modified from TJ Mahr at
      \url{https://github.com/tidyverse/broom/issues/96}):
      \itemize{
	\item Random effect variances and covariances can now be extracted. 
	\code{effects = "ran_pars"} returns the standard deviations/correlations
	of random effects (if \code{scales="sdcor"} [default])
	or their variances
	and covariances (if \code{scales = "varcov"})
	\item Random effects estimates are now extracted with 
	\code{effects = "ran_coefs"} for the group-level
	estimates (previously these effects were 
	extracted with \code{tidy(model, "random")}) or
	\code{effects = "ran_vals"} for the conditional modes
	(deviations of the group-level parameters from the
	population-level estimates)
	\item \code{effects} can take a vector of values (those listed
	above, plus "fixed" for fixed effects). The default value
	is \code{effects = c("ran_pars", "fixed")} which extracts random 
	effect variances/covariances and fixed effect estimates.
      }
      \item \code{term} names for random-effect parameters no longer
      contain a (redundant) group specifier (at least for \code{lme4}
      models); use something like
      \code{tidyr::unite(term,term,group,sep=".")} 
      to collapse the two columns
    }
  }
}
