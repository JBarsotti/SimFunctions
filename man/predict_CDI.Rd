% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Base_Functions.R
\name{predict_CDI}
\alias{predict_CDI}
\title{Predict Probability of CDI}
\usage{
predict_CDI(
  df,
  cores = 1L,
  folds = 5L,
  model = "gbm",
  scoring = "ROC",
  repeats = 2L,
  trees = 1L,
  grid = 0L
)
}
\arguments{
\item{df}{Table of agent data to use in predictive model}

\item{cores}{Number of cores to use in parallel}

\item{folds}{Number of cross-validation folds}

\item{model}{Machine-learning model to use}

\item{scoring}{Method of scoring to judge model accuracy}

\item{repeats}{Number of repeats for cross-validation}

\item{trees}{Number of trees to use if applicable to machine-learning model}

\item{grid}{Tuning grid to use if applicable}
}
\value{
Machine-learning model
}
\description{
Use machine-learning methods to predict probability of agent having CDI
}
\examples{
model <- predict_CDI(df, cores = 1L, folds = 5L, model = "gbm", scoring = "ROC", repeats = 2L, trees = 1L, grid = 0L);
}
