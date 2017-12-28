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
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{lme4} tidiers now optionally take a pre-computed
      \code{profile} argument when using \code{conf.method="profile"}
    }
  }
}