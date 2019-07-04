#' CLL data
#'
#' A dataset containing survival outcome and predictors on 694 patients who
#' received hematopoietic stem cell transplant.
#'
#' @format A data frame with 694 rows and 11 variables. Each row describes the data from a single patient. The below described variables are included in the data file. Missing observations are present in the variables performance status(9\%), remission status (6\%) and cytogenic abnormality (25\%).
#' \describe{
#'   \item{id}{record identification number}
#'   \item{age10}{age at transplantation}
#'   \item{perfstat}{performance status indicated by the Karnofsky Index (four categories)}
#'   \item{remstat}{remission status at transplantation (three categories)}
#'   \item{cyto}{cytogenetic abnormalities (four categories)}
#'   \item{asct}{previous autologous transplantation (two categories)}
#'   \item{donor}{donor type (three categories)}
#'   \item{sex_match}{patient-donor sex match (four categories)}
#'   \item{cond}{conditioning regimen (three categories)}
#'   \item{srv5y}{overall survival (OS) up to five years after first allogeneic stem cell transplantation}
#'   \item{srv5y_s}{censoring indicator (0=alive at end follow-up, 1=dead)}
#' }
#' @references Please reference the following papers when using this data. Schetelig, J. et al. (2017) Risk factors for treatment failure after allogeneic transplantation of patients with CLL: a report from the European Society for Blood and Marrow Transplantation. {Bone Marrow Transplantation},  52, 552-560. Schetelig, J. et al. (2017) Centre characteristics and procedure-related factors have an impact on outcomes of allogeneic transplantation for patients with CLL: a retrospective analysis from the European Society for Blood and Marrow Transplantation (EBMT). {British Journal of Haematology},  178, 521-533. Mertens, B.J.A. et al. (2019) Construction and assessment of prediction rules for binary outcome in the presence of missing predictor data using multiple imputation and cross-validation: theoretical perspective and data-based evaluation. {Biometrical Journal}. See ArXiv for an early version \url{https://arxiv.org/abs/1810.05099}.  We thank EBMT and DKMS for their work in collecting and preparing the CLL data and for approval to share the data.
#' @source European Society for Blood and Marrow Transplantation (EBMT). \url{https://www.ebmt.org}
#'
#'
#' @docType data
"cll"
