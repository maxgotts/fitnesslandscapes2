\name{2-dimensional frequency binning}
\alias{binCounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{binCounts}
\description{
  If we want to understand the frequencies of two variables, we can use individual binnings. However, this is not recommended becaues it ignores all epistasis. That is, swim speed and flipper length may be positively correlated, so we should not treat their histogram/density plots as independent. Hence, we create a function that bins two variables simultaneously. It creates bins according to the fraction by which we divide the standard deviation (i.e. bins are min+n*SD/k for k fixed and n a natural number) for convenience.
}
