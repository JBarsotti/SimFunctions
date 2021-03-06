% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Base_Functions.R
\name{manual_predict_CDI}
\alias{manual_predict_CDI}
\title{Performing Manual Cross Validation for Bayesian Logistic Regression}
\usage{
manual_predict_CDI(data, folds = 5L, repeats = 5L)
}
\arguments{
\item{data}{Table of agent data to use in predictive model}

\item{folds}{Number of cross-validation folds}

\item{repeats}{Number of repeats for cross-validation}
}
\value{
Model and statistics for accuracy of the model
}
\description{
Use manual cross-validation with Bayesian logistic regression for difficult to separate classes
}
\examples{
results <- manual_predict_CDI(data, folds = 5L, repeats = 5L)
}
