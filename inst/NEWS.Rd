\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{broom.mixed News}
\encoding{UTF-8}

\section{CHANGES IN VERSION 0.0.1}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Sorted out some of the confusion over random effect naming:
      "ran_modes" extracts conditional modes (deviations from
      population-level estimates), while "ran_coefs" extracts
      group-level estimates
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
      \item \code{term} names for random-effect parameters no longer
      contain a (redundant) group specifier (at least for \code{lme4}
      models); use something like
      \code{tidyr::unite(term,term,group,sep=".")} 
      to collapse the two columns
    }
  }
}
