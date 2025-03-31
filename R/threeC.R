#'  threeC
#'
#' A random subsample of 500 subjects,  from the French Three-Cities cohort aimed at assessing the relation between vascular factors and dementia in the elderly, for the example of the LSJM package.
#' The sample include participants without dementia at baseline and aged 65 years old or older.
#' Repeated measures of systolic blood pressure were collected over a maximum period of 20 years. At each visit, systolic blood pressure was measured two or three times.
#'
#' @format A data frame with 5248 rows and 11 variables:
#' \describe{
#'   \item{ID}{the ID of each subject}
#'   \item{SBP}{the value of measured sysotlic blood pressure in mmHg}
#'   \item{age.visit}{the age at measurement of SBP}
#'   \item{age.final}{the age of death or censoring (last news)}
#'   \item{age0}{the age at entry in the cohort}
#'   \item{age.last}{the age of last visit without dementia}
#'   \item{age.first}{the age of visit diagnosis for dementia}
#'   \item{sex}{the sex of the subject}
#'   \item{dem}{the indicator of dementia}
#'   \item{death}{the indicator of death}
#'   \item{num.visit}{the idenfication of the visit}
#' }
#'
"threeC"
