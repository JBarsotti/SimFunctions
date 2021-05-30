#' Installing Packages 
#'
#' Install packages needed to run agent-based CDI simulator
#' @return Nothing
#' @examples 
#' install_packages();
#' @export
install_packages <- function()
{

  # Data analysis packages
  packages <- c("gbm", "xgboost", "nortest", "RcppArmadillo",
                "arm", "janitor", "pROC", "caret", "randomForest",
                "reshape2", "ggpubr", "parallel", "doSNOW", "foreach", "Rmpi", "tidyverse")

  install.packages(packages)

}

#------------------------------------------------------------------------------



#' Adding packages
#'
#' Load the packages needed to run the agent-based CDI simulator
#' @return Nothing
#' @examples 
#' add_packages();
#' @export
add_packages <- function()
{

  library(caret)
  library(nortest)
  library(RcppArmadillo)
  library(doParallel)
  library(gbm)
  library(arm)
  library(janitor)
  library(pROC)
  library(randomForest)
  library(xgboost)
  library(reshape2)
  library(ggpubr)
  library(parallel)
  library(doSNOW)
  library(foreach)
  library(tidyverse)

}

#------------------------------------------------------------------------------



#' Access Patient-Level Data
#'
#' Load patient data to be used in machine learning and generating patients for the simulator
#' @param link The directory containing the patient data
#' @param core The table you would like to use from the patient data
#' @return The patient database
#' @examples 
#' db <- HCUP(link = "/home/ebarsotti/LSS/Statepi_Marketscan/databases/HCUP/IA.db", core = "IA_SID_2010_core");
#' @export
HCUP <- function(link = "/home/ebarsotti/LSS/Statepi_Marketscan/databases/HCUP/IA.db",
                 core = "IA_SID_2010_core")
{
  # Get all the tables into one database
  db <- src_sqlite(link)
  db

  # Pull just data from one year
  db_year <- db %>% tbl(core)
  print(colnames(db_year))

  return (db_year)
}

#------------------------------------------------------------------------------



#' Reformat Data for Machine Learning
#'
#' Remove some variables and convert some variable types to factors
#' @param master The data to be reformatted
#' @return A new, reformatted version of the data
#' @examples 
#' data <- CDI_predict_data(master);
#' @export
CDI_predict_data <- function(master)
{

  # Removing variables not currently being used in the model
  copy_master <- master %>% dplyr::select(-c(hcup_ed, key, proctype, dxccs, dshospid,
                                             visitlink, cdi_per_month, pstco2, medincstq,
                                             all_cases_per_month, pnum_r, died, los,
                                             amonth, drg_nopoa, pl_nchs2006, tran_out,
                                             daystoevent, zip3,
                                             pl_ruca4_2005, X, dispuniform, mdc, tran_in,
                                             pointoforiginub04, dx))


  # Cleaning variable names just in case R doesn't like them
  copy_master <- copy_master %>% clean_names()

  # Removing variable COMBINED_ID because not used in prediction
  # Also removing duplicate COMBINED_IDs because no longer needed
  copy_master <- copy_master %>% dplyr::select(-combined_id)

  # Subtraction 1 from each NDX for CDI patients
  # so that CDI is not counted as one of the number of diagnoses
  copy_master <- copy_master %>% dplyr::mutate(ndx = ifelse(cdi == 1, ndx - 1, ndx))

  # Factorizing categorical variables
  copy_master$atype          <- as.factor(copy_master$atype)
  copy_master$female         <- as.factor(copy_master$female)
  copy_master$hospitalunit   <- as.factor(copy_master$hospitalunit)
  copy_master$orproc         <- as.factor(copy_master$orproc)
  copy_master$race           <- as.factor(copy_master$race)
  copy_master$zipinc_qrtl    <- as.factor(copy_master$zipinc_qrtl)
  copy_master$pl_ruca10_2005 <- as.factor(copy_master$pl_ruca10_2005)
  copy_master$age            <- as.factor(copy_master$age)

  copy_master$cdi <- factor(copy_master$cdi,
                               levels = c(1, 0), labels = c("Yes", "No"))

  copy_master$cdi <- relevel(copy_master$cdi, ref = "No")

  return(copy_master)

}

#------------------------------------------------------------------------------



#' Convert Some Data to Quantiles
#'
#' Convert some categorical variables to quantiles to reduce number of categories
#' @param data Data containing variables to be converted
#' @param columns Variable columns to be converted
#' @param partition Number of quantiles to use
#' @return New version of the data
#' @examples 
#' data <- quantiles(data, cols, 8);
#' @export
quantiles <- function(data, columns, partition)
{

  for (i in 1:length(columns)) {
    data[[columns[i]]] <- with(data,
                               ntile(data[columns[i]], partition[i]))
  }

  return(data)

}

#------------------------------------------------------------------------------



#' Predict Probability of CDI
#'
#' Use machine-learning methods to predict probability of agent having CDI
#' @param df Table of agent data to use in predictive model
#' @param cores Number of cores to use in parallel
#' @param folds Number of cross-validation folds
#' @param model Machine-learning model to use
#' @param scoring Method of scoring to judge model accuracy
#' @param repeats Number of repeats for cross-validation
#' @param trees Number of trees to use if applicable to machine-learning model
#' @param grid Tuning grid to use if applicable
#' @return Machine-learning model
#' @examples 
#' model <- predict_CDI(df, cores = 1L, folds = 5L, model = "gbm", scoring = "ROC", repeats = 2L, trees = 1L, grid = 0L);
#' @export
predict_CDI <- function(df, cores = 1L, folds = 5L, model = "gbm",
                        scoring = "ROC", repeats = 2L, trees = 1L, grid = 0L)
{

  set.seed(1234)

  cl <- makePSOCKcluster(cores)
  registerDoParallel(cl)

  trControl <- trainControl(method  = "cv", number = folds, repeats = repeats,
                            savePredictions = TRUE, classProbs = TRUE,
                            summaryFunction = twoClassSummary, verboseIter = TRUE,
                            allowParallel = TRUE)

  if (model == "glmnet") {
    model <- caret::train(cdi ~ .,
                          data = df,
                          method = model,
                          metric = scoring,
                          trControl = trControl,
                          tuneGrid = grid)
  }
  else if (model != "rf") {
    model <- caret::train(cdi ~ .,
                        data = df,
                        method = model,
                        metric = scoring,
                        trControl = trControl)
  }
  else {
    model <- caret::train(cdi ~ .,
                          data = df,
                          method = model,
                          metric = scoring,
                          ntree = trees,
                          trControl = trControl)
  }

  stopCluster(cl)
  rm(cl)

  return(model)

}

#------------------------------------------------------------------------------



#' Machine-Learning Model Statistics
#'
#' Check model accuracy and statistics using ROC curves and variable importance
#' @param model Machine-learning model
#' @param model_type Model type for plotting titles
#' @param train_data Data trained on by the model
#' @param test_data Data to test on by the model
#' @return ROC curve and variable importance scores
#' @examples 
#' importance <- diagnose(gbm_model, "GBM", train, test)
#' @export
diagnose <- function(model, model_type, train_data, test_data)
{

  train <- train_data
  test <- test_data

  model.classTrain <- predict(model, newdata = train, type = "raw")
  print(confusionMatrix(train$cdi, model.classTrain))

  model.probs <- predict(model, newdata = test, type = "prob")
  rocCurve.model <- roc(test$cdi, model.probs[, "Yes"], direction = "<")

  importance_score <- varImp(model, scale = FALSE)
  plot(rocCurve.model, col = "red", main = model_type)

  print(auc(rocCurve.model))
  return(importance_score)

}

#------------------------------------------------------------------------------



#' XGBOOST Algorithm
#'
#' Use the R XGBOOST machinel-learning algorithm 
#' @param train Data to train the model on
#' @param test Data to test the model on
#' @param folds Number of cross-validation folds
#' @param scoring Method of scoring to judge model accuracy
#' @param scale Scaling of the importance of the positive class
#' @return Model and area under the curve of the model (AUC)
#' @examples 
#' results <- xgb_cdi(train, test, folds = 5L, scoring = "ROC", scale = 1L)
#' @export
xgb_cdi <- function(train, test, folds = 5L, scoring = "ROC", scale = 1L)
{

  train_xgb <- sparse.model.matrix(cdi ~., data = train)[, -1]
  outcome <- as.numeric(train$cdi)
  outcome[which(outcome == 1)] <- 0L
  outcome[which(outcome == 2)] <- 1L
  data_train <- xgb.DMatrix(data = train_xgb, label = outcome)

  set.seed(1234)

  xgb_model <- xgb.train (data = data_train, nrounds = 150, booster = "gbtree",
                          verbose = TRUE, eta = 0.3, objective = "binary:logistic",
                          gamma = 0, max_depth = 7, scale_pos_weight = scale,
                          eval_metric = scoring)


  test_xgb <- sparse.model.matrix(cdi ~., data = test)[, -1]
  outcome <- as.numeric(test$cdi)
  outcome[which(outcome == 1)] <- 0L
  outcome[which(outcome == 2)] <- 1L
  data_test <- xgb.DMatrix(data = test_xgb, label = outcome)
  preds_xgb <- predict(xgb_model, data_test, type = "prob")

  auc_xgb <- roc(test$cdi ~ preds_xgb, levels = c("No", "Yes"), auc = T, plot = T, direction = "<")

  return(list(xgb_model, auc_xgb))

}

#------------------------------------------------------------------------------



#' Performing Manual Cross Validation for Bayesian Logistic Regression
#'
#' Use manual cross-validation with Bayesian logistic regression for difficult to separate classes
#' @param data Table of agent data to use in predictive model
#' @param folds Number of cross-validation folds
#' @param repeats Number of repeats for cross-validation
#' @return Model and statistics for accuracy of the model
#' @examples 
#' results <- manual_predict_CDI(data, folds = 5L, repeats = 5L)
#' @export
manual_predict_CDI <- function(data, folds = 5L, repeats = 5L)
{

  aucs   <- 0L
  thresh <- 0L
  sens   <- 0L
  spec   <- 0L
  ppv    <- 0L
  npv    <- 0L

  for (i in 1:repeats) {

    df <- data

    # Create folds (some code from Dr. Miller's CrossValidation.R script)
    folded_df_index <- sample(rep(1:folds, length = nrow(df)))
    df$fold <- folded_df_index
    df <- df %>% tibble::add_column(preds = NA)

    for (j in 1:folds) {
      train <- df %>% dplyr::filter(fold != j) %>% dplyr::select(-c(fold, preds))
      test <- df %>% dplyr::filter(fold == j) %>% dplyr::select(-c(fold, preds))
      test_indeces <- which(df$fold == j)

      set.seed(1234)

      model_cdi <- bayesglm(cdi ~ ., data = train, family = "binomial")
      preds <- predict(model_cdi, test, type = "response")
      df$preds[test_indeces] <- preds

      auc_cdi <- roc(test$cdi ~ preds, levels = c("No", "Yes"), auc = T, direction = "<")

      best_coord <- coords(auc_cdi, x = "best",
                                 best.method  = "closest.topleft",
                                 ret = c("threshold", "sens", "spec", "ppv", "npv"))

      aucs   <- aucs + auc_cdi$auc
      thresh <- thresh + best_coord$threshold
      sens   <- sens + best_coord$sensitivity
      spec   <- spec + best_coord$specificity
      ppv    <- ppv + best_coord$ppv
      npv    <- npv + best_coord$npv

      if (j == 1) {
        plot(auc_cdi, col = rainbow(50)[j+5], main = "ROC Curves for Each Fold Using Bayes GLM")
        legend_nums <- paste("Fold", j)
        legend_cols <- rainbow(50)[j+5]
      }
      else {
        plot(auc_cdi, add = TRUE, col = rainbow(50)[j+5])
        legend_nums <- c(legend_nums, paste("Fold", j))
        legend_cols <- c(legend_cols, rainbow(50)[j+5])
      }
    }

    legend(0.0, 0.5, legend = legend_nums, col = legend_cols, lty = 1, cex = 0.8)
  }

  avg_auc    <- aucs/(folds*repeats)
  avg_thresh <- thresh/(folds*repeats)
  avg_sens   <- sens/(folds*repeats)
  avg_spec   <- spec/(folds*repeats)
  avg_ppv    <- ppv/(folds*repeats)
  avg_npv    <- npv/(folds*repeats)

  stats <- tibble(auc = avg_auc, threshold = avg_thresh, sensitivity = avg_sens,
                  specificity = avg_spec, ppv = avg_ppv, npv = avg_npv)

  return(list(model_cdi, stats))

}

#------------------------------------------------------------------------------



#' Basic Univariate Statistics of Data
#'
#' Generate basic statistics for numerical data, as well as random values to be used in the agent-based simulation 
#' @param data_in Table of data to be analyzed
#' @param key Name of the variable to be analzyed
#' @param samplesize Number of random values from the data to be taken
#' @param histo TRUE if you want to plot a histogram, FALSE if you want to plot a bar graph
#' @param plot Decides whether data will be plotted or not
#' @return A vector of random values taken from the data
#' @examples 
#' random_values <- basic_stats(data_in, key, samplesize, width, histo = TRUE, plot = FALSE)
#' @export
basic_stats <- function(data_in, key, samplesize, width, histo = TRUE, plot = FALSE)
{
  type <- class(data[[1]])
  data <- data_in[[1]]

  while (!is.null(dev.list()))  dev.off()

  # Plot histogram of data
  #------------------------------------------------------------------------------
  if (histo == FALSE) {
    xvals <- data_in[[1]]

    if (plot == TRUE) {
      plt <- ggplot(data_in, aes(x = xvals, fill = xvals)) +
        geom_bar() +
        theme(legend.position = "right") +
        coord_flip() +
        labs(title = paste("Barplot of", key), x = paste(key), y = "Count")
    }

    data_num <- as.numeric(data)
  }
  else {
    if (plot == TRUE) {
      plt <- qplot(data,
                   geom = "histogram",
                   binwidth = width,
                   main = paste("Histogram for", key),
                   xlab = paste(key),
                   ylab = "Frequency",
                   fill=I("blue"),
                   col=I("red"),
                   alpha=I(.2))
    }
  }

  if (plot == TRUE) {
    print(plt)
  }

  #------------------------------------------------------------------------------

  # Get summary statistics
  #------------------------------------------------------------------------------
  len <- length(unique(data))
  summary_vals <- summary(data)
  #------------------------------------------------------------------------------
  

  # Check normality of variable
  #------------------------------------------------------------------------------
  norm_test <- ad.test(data_num)
  #------------------------------------------------------------------------------

  # Sampling from LOS
  #------------------------------------------------------------------------------
  result <- sample(data, samplesize, replace = TRUE)
  #------------------------------------------------------------------------------

  return(result)

}

#------------------------------------------------------------------------------



#' Patient-Level Independent Variables
#'
#' Draw independent values for patients, age and sex, from patient data
#' @param data Data to be used
#' @return A list of two vectors: age of patients and sex of patients
#' @examples 
#' list_of_age_and_sex <- independent(data)
#' @export
independent <- function(data)
{
  # Age of Patients
  age      <- data %>% dplyr::select(age)
  key      <- "age"
  width    <- 1
  age_pat  <- basic_stats(age, key, npat, width, histo = FALSE, plot = FALSE)

  # Sex of Patients
  sex      <- data %>% dplyr::select(female)
  key      <- "sex"
  width    <- 1
  sex_pat  <- basic_stats(sex, key, npat, width, histo = FALSE, plot = FALSE)

  return(list(age_pat, sex_pat))

}

#------------------------------------------------------------------------------



#' Patient-Level Dependent Variables
#'
#' Draw dependent values for patients, ndx, elix, necode, orproc, npr, prior_visit, and hospitalunit, from patient data
#' @param data Data to be used
#' @param npat Number of patients
#' @param age_pat Ages of patients
#' @param sex_pat Sexes of patients
#' @return A tibble of patients with their individual characteristics
#' @examples 
#' tibble_patients <- dependent(data, npat, age_pat, sex_pat)
#' @export
dependent <- function(data, npat, age_pat, sex_pat)
{

  ndx          <- numeric(npat)
  elix         <- ndx
  necode       <- ndx
  npr          <- ndx
  prev_cdi     <- ndx
  prior_visit  <- ndx
  hospitalunit <- ndx
  orproc       <- ndx

  for (i in 1:npat) {

    age <- age_pat[i]
    sex <- sex_pat[i]

    subset <- data %>% dplyr::filter(age == age & female == sex)
    ndx[i] <- sample(subset$ndx, 1L)
    elix[i] <- sample(subset$elix, 1L)
    necode[i] <- sample(subset$necode, 1L)
    npr[i] <- sample(subset$npr, 1L)
    prev_cdi[i] <- sample(subset$prev_cdi, 1L)
    prior_visit[i] <- sample(subset$prior_visit, 1L)
    hospitalunit[i] <- sample(subset$hospitalunit, 1L)
    orproc[i] <- sample(subset$orproc, 1L)

  }

  tibble_pat <- tibble(ndx = ndx, elix = elix, npr = npr, prior_visit = prior_visit,
                       necode = necode, prev_cdi = prev_cdi, hospitalunit = hospitalunit,
                       orproc = orproc, age = age_pat, female = sex_pat)

  tibble_pat$female       <- as.factor(tibble_pat$female)
  tibble_pat$hospitalunit <- as.factor(tibble_pat$hospitalunit)
  tibble_pat$orproc       <- factor(tibble_pat$orproc, levels = c("1", "2"), labels = c("0", "1"))

  return(tibble_pat)

}

#------------------------------------------------------------------------------



#' Patient CDI Infection Probabilities
#'
#' Calculate probabilities of patients being infected based on machine-learning algorithm
#' @param model Machine-learning model to be used
#' @param data Tibble of patients
#' @param initial_infect_pat Initial number of patients to be set as infected
#' @return A new tibble of patients with infectivities and susceptibilities based on probabilities from the machine-learning model
#' @examples 
#' new_tibble_patients <- pat_preds(model, data, initial_infect_pat)
#' @export
pat_preds <- function(model, data, intial_infect_pat)
{

  preds_cdi <- predict(gbm_model, data, type = "prob")

  data <- data %>% dplyr::mutate(id = 1:npat) %>%
    dplyr::relocate(id) %>%
    dplyr::mutate(Yes = preds_cdi$Yes)

  data <- data %>% arrange(-Yes)
  infected_pat <- data$id[1:initial_infect_pat]

  data <- data %>% dplyr::mutate(disease_stat = ifelse(id %in% infected_pat, 1L, 0L))
  data <- data %>% dplyr::mutate(sus = Yes)

  return(data)

}

#------------------------------------------------------------------------------



#' Patient-Level Length of Stay Variable Generation
#'
#' Calculate length of stay for individual agents based on age, sex, OR procedure, and CDI status
#' @param data_one Data for all patients
#' @param data_two Tibble of data only for patients to be used in this model
#' @param npat Number of patients in the tibble of patients
#' @return An updated tibble of patients with a new length of stay column
#' @examples 
#' complete_tibble_patients <- pat_los(data_one, data_two, npat)
#' @export
pat_los <- function(data_one, data_two, npat)
{

  los <- numeric(npat)

  for (i in 1:npat) {

    subset <- data_one %>% dplyr::filter(age == data_two$age[i] &
                                       female == data_two$female[i] &
                                       cdi == data_two$disease_stat[i] &
                                       orproc == data_two$orproc[i])
    tryCatch(
      {
        los[i] <- sample(subset$los, 1L)
      },
      error = function(e)
      {
        print("No LOS for those parameters found. Switching to just age and CDI.")
        subset <- data_one %>% dplyr::filter(age == data_two$age[i] &
                                               cdi == data_two$disease_stat[i])
        los[i] <- sample(subset$los, 1L)
      }

    )
  }

  data_two <- data_two %>% dplyr::mutate(los = los)
  return(data_two)

}

#------------------------------------------------------------------------------



#' Constructing Tibble of Patients
#'
#' Generate patients based on previous patient-level data on characteristics such as age, sex, Elixhauser risk score, number of diagnoses, etc.
#' @param ml_data Data to use in machine-learning algorithm
#' @param all_data All patient-level data available
#' @param npat Number of patients in the simulation
#' @param initial_infected_pat Number of initial infected patients in the simulation
#' @param model Machine-learning model to use
#' @return A tibble of patients for use in the simuation; patients generated by multiple other functions in this library
#' @examples 
#' patients <- patients_table(ml_data, all_data, npat, initial_infect_pat, model)
#' @export
patients_table <- function(ml_data, all_data, npat, initial_infect_pat, model)
{

  ind_vars <- independent(ml_data)
  age_pat <- as.factor(ind_vars[[1]])
  sex_pat <- factor(ind_vars[[2]], levels = c(1, 0), labels = c(1, 0))

  tibble_pat <- dependent(ml_data, npat, age_pat, sex_pat)

  tibble_pat <- pat_preds(model, tibble_pat, inital_infect_pat)

  tibble_pat <- pat_los(all_data, tibble_pat, npat)

  tibble_pat <- tibble_pat %>%
    dplyr::rename(los_start = los) %>%
    dplyr::mutate(los_current = los_start) %>%
    dplyr::mutate(mobility = ifelse(orproc == 1L, 0L, 1L)) %>%
    dplyr::mutate(infect = sus) %>%
    dplyr::mutate(cdi_state = ifelse(disease_stat == 1L, 1L, 0L))

  return(tibble_pat)

}

#------------------------------------------------------------------------------



#' Constructing Tibble of Healthcare Workers
#'
#' Generate healthcare workers, of whom none are initially infected, with susceptibility equal to 0.05 quantile of patients and infectivity equal to 0.5 quantile of patients
#' @param tibble_pat Tibble of individual patients in the simulation
#' @param nhcw Number of healthcare workers in the simulation
#' @param data Patient-level data that is also used for finding out which hospital unit healthcare professionals work in
#' @return A tibble of individual healthcare workers
#' @examples 
#' tibble_hcw <- hcws_table(tibble_pat, nhcw, data)
#' @export
hcws_table <- function(tibble_pat, nhcw, data)
{

  mat_hcw <- matrix(-10L, nrow = nhcw, ncol = ncol(tibble_pat))
  names_col <- c(colnames(tibble_pat))
  colnames(mat_hcw) <- names_col

  tibble_hcw <- as_tibble(mat_hcw)
  tibble_hcw$id <- 1:nhcw

  hcw_age                           <- runif(nhcw, 0, 1)
  hcw_age[which(hcw_age <= 0.005)]  <- 1L
  hcw_age[which(hcw_age <= 0.0599)] <- 2L
  hcw_age[which(hcw_age <= 0.3269)] <- 3L
  hcw_age[which(hcw_age <= 0.5719)] <- 4L
  hcw_age[which(hcw_age <= 0.7774)] <- 5L
  hcw_age[which(hcw_age <= 0.9401)] <- 6L
  hcw_age[which(hcw_age < 1.0)]     <- 7L

  tibble_hcw$age <- hcw_age

  hcw_female <- runif(nhcw, 0, 1)
  hcw_female[which(hcw_female <= 0.7)]  <- 1L
  hcw_female[which(hcw_female < 1.0)] <- 0L

  tibble_hcw$female <- hcw_female
  tibble_hcw$hospitalunit <- sample(data$hospitalunit, nhcw)
  tibble_hcw$disease_stat <- 0L
  tibble_hcw$mobility <- 1L

  tibble_hcw <- tibble_hcw %>% dplyr::mutate(sus = quantile(tibble_pat$sus, 0.05))
  tibble_hcw$infect <- quantile(tibble_pat$infect, 0.5)

  return(tibble_hcw)

}

#------------------------------------------------------------------------------



#' Initialize Agents
#'
#' Make one list of agents that contain both healthcare workers and patients, which is needed for the C++ disease propogation function
#' @param nhcw Number of healthcare workers in the simulation
#' @param npat Number of patients in the simulation
#' @param hcw_feature_names Individual healthcare worker feature names
#' @param pat_feature_names Individual patient feature names
#' @param tibble_pat Tibble of patients
#' @param tibble_hcw Tibble of healthcare workers
#' @return A list containing two objects: the list of agents is the first; the initial set of infected agents is the second
#' @examples 
#' all_agents <- init_agents(nhcw = 800L, npat = 100L, hcw_feature_names, pat_feature_names, tibble_pat, tibble_hcw)
#' @export
init_agents <- function(nhcw = 800L, npat = 100L, hcw_feature_names, pat_feature_names,
                        tibble_pat, tibble_hcw)
{

  total_agents             <- nhcw + npat
  feature_names            <- c("type", union(hcw_feature_names, pat_feature_names))
  total_features           <- length(feature_names)

  # Initialize all values in list of agents to -10, which signifies that
  # the agent does not have the specified characteristic
  agents <- rep(list(numeric(total_agents) - 10L), total_features)

  # Name agent features
  names(agents) <- feature_names

  # Agent type 1 is patients and 2 is HCWs
  agents$type[1:npat]              <- 1L
  agents$type[npat+1:total_agents] <- 2L

  # Initialize initial infecteds vector
  initial_set <- c()

  # Initialize status of each individual patient feature
  for (i in 1:npat) {

    agents$id[i]           <- i
    agents$los_start[i]    <- as.numeric(tibble_pat$los_start[i])
    agents$sus[i]          <- as.numeric(tibble_pat$sus[i])
    agents$infect[i]       <- as.numeric(tibble_pat$infect[i])
    agents$disease_stat[i] <- as.numeric(tibble_pat$disease_stat[i])
    agents$necode[i]       <- as.numeric(tibble_pat$necode[i])
    agents$female[i]       <- as.numeric(tibble_pat$female[i])
    agents$los_current[i]  <- as.numeric(tibble_pat$los_current[i])
    agents$ndx[i]          <- as.numeric(tibble_pat$ndx[i])
    agents$prev_cdi[i]     <- as.numeric(tibble_pat$prev_cdi[i])
    agents$Yes[i]          <- as.numeric(tibble_pat$Yes[i])
    agents$mobility[i]     <- as.numeric(tibble_pat$mobility[i])
    agents$elix[i]         <- as.numeric(tibble_pat$elix[i])
    agents$hospitalunit[i] <- as.numeric(tibble_pat$hospitalunit[i])
    agents$npr[i]          <- as.numeric(tibble_pat$npr[i])
    agents$orproc[i]       <- as.numeric(tibble_pat$orproc[i])
    agents$prior_visit[i]  <- as.numeric(tibble_pat$prior_visit[i])
    agents$age[i]          <- as.numeric(tibble_pat$age[i])

    if (tibble_pat$disease_stat[i] == 1L) {
      initial_set <- c(initial_set, i)
    }

  }

  # Do the same for HCW features
  for (i in (npat+1):total_agents) {

    agents$id[i]           <- i
    agents$los_start[i]    <- as.numeric(tibble_hcw$los_start[i - npat])
    agents$sus[i]          <- as.numeric(tibble_hcw$sus[i - npat])
    agents$infect[i]       <- as.numeric(tibble_hcw$infect[i - npat])
    agents$disease_stat[i] <- as.numeric(tibble_hcw$disease_stat[i - npat])
    agents$necode[i]       <- as.numeric(tibble_hcw$necode[i - npat])
    agents$female[i]       <- as.numeric(tibble_hcw$female[i - npat])
    agents$los_current[i]  <- as.numeric(tibble_hcw$los_current[i - npat])
    agents$ndx[i]          <- as.numeric(tibble_hcw$ndx[i - npat])
    agents$prev_cdi[i]     <- as.numeric(tibble_hcw$prev_cdi[i - npat])
    agents$Yes[i]          <- as.numeric(tibble_hcw$Yes[i - npat])
    agents$mobility[i]     <- as.numeric(tibble_hcw$mobility[i - npat])
    agents$elix[i]         <- as.numeric(tibble_hcw$elix[i - npat])
    agents$hospitalunit[i] <- as.numeric(tibble_hcw$hospitalunit[i - npat])
    agents$npr[i]          <- as.numeric(tibble_hcw$npr[i - npat])
    agents$orproc[i]       <- as.numeric(tibble_hcw$orproc[i - npat])
    agents$prior_visit[i]  <- as.numeric(tibble_hcw$prior_visit[i - npat])
    agents$age[i]          <- as.numeric(tibble_hcw$age[i - npat])

    if (tibble_hcw$disease_stat[i - npat] == 1L) {
      initial_set <- c(initial_set, i)
    }

  }

  return (list(agents, initial_set))

}

#------------------------------------------------------------------------------



#' Contact Matrix
#'
#' Randomly generate a symmetric matrix of contacts between agents
#' @param agents The list of agents, both patients and healthcare workers, who will be in the simulation
#' @return A symmetric matrix of size NxN, where N is the number of agents, containing number of contacts between each individual agent for each timestep
#' @examples 
#' contact_matrix <- contacts(agents)
#' @export
contacts <- function(agents)
{

  total_agents   <- length(agents$id)
  contact_matrix <-  matrix(0L, nrow = total_agents, ncol = total_agents)

  for (i in 1:total_agents) {
    for (j in 1:total_agents) {

      if (i == j) {
        contact_matrix[i, j] <- 0L
      }
      else if (agents$type[i] == 1L & agents$type[j] == 1L) {
        contact_matrix[i, j] <- 0L
      }
      else if (agents$type[i] == 1L & agents$type[j] == 2L & (agents$hospitalunit[i] != agents$hospitalunit[j])) {
        contact_matrix[i, j] <- 0L
      }
      else if (agents$type[i] == 2L & agents$type[j] == 1L & (agents$hospitalunit[i] != agents$hospitalunit[j])) {
        contact_matrix[i, j] <- 0L
      }
      else if (agents$type[i] == 2L & agents$type[j] == 2L & (agents$hospitalunit[i] != agents$hospitalunit[j])) {
        contact_matrix[i, j] <- sample(1L:2L, 1L)
      }
      else if ((agents$mobility[i] == 0L & agents$type[j] == 2L)) {
        contact_matrix[i, j] <- sample(0L:1L, 1L)
      }
      else {
        contact_matrix[i, j] <- sample(1L:4L, 1L)
      }

    }
  }

  contact_matrix[lower.tri(contact_matrix)] <- t(contact_matrix)[lower.tri(contact_matrix)]
  return(contact_matrix)

}

#------------------------------------------------------------------------------



#' CDI Agent-Based Simulation Function
#'
#' Take in the list of all agents and the contact matrix, call the C++ function "propogate" that spreads the infection, and return statistics about how much the infection spreads and to whom it spreads
#' @param ml_data The data to be used in the machine learning algorithm to calculate risk of CDI
#' @param mod The machine-learning model to be used to calculate risk of CDI
#' @param timestep Number of iterations for the simulation to loop through
#' @param agents List of all agents in the simulation
#' @param contact_matrix Matrix of contacts between all individual agents
#' @param initial_set Initial set of infected agents
#' @param plot_now TRUE if you want to watch attack rate over time; FALSE if you don't
#' @param num Variable that keeps track of number of timesteps completed
#' @param secondary_infections Variable keeps track of number of secondary infections
#' @param hcw_to_pat_infections Variable keeps track of number of healthcare worker to patient transmissions
#' @param pat_to_hcw_infections Variable keeps track of number of patient to healthcare worker transmissions
#' @param pat_to_pat_infections Variable keeps track of number of patient to patient transmissions
#' @param hcw_to_hcw_infections Variable keeps track of number of healthcare worker to healthcare worker transmissions
#' @param secondary_patient_infected Variable keeps track of number of secondary infections spread to patients
#' @param secondary_hcw_infected Variable keeps track of number of secondary infections spread to healthcare workers
#' @param R0 Variable to store the basic reproduction number of the model
#' @param attack_curve Variable to store the attack curve, currently just the total number of secondary infections, in the model
#' @param pat_to_hcw Variable to store number of transmissions from patients to healthcare workers as a vector for later use
#' @param hcw_to_pat Variable to store number of transmissions from healthcare workers to patients as a vector for later use
#' @param hcw_to_hcw Variable to store number of transmissions from healthcare workers to other healthcare workers as a list for later use
#' @param iters Variable to store iteration number as a list for later use
#' @param discharged Variable to store the number of discharged patients
#' @param total_pat Variable to store the total number of patients, initial plus discharged
#' @param transmission_dir Variable containing directory name where the C++ transmission function is stored
#' @return Long list containing updated agents and simulation statistics
#' @examples 
#' sim_results <- simulate(ml_data, mod, timestep = 24L, agents, contact_matrix, initial_set, plot_now = FALSE, num, secondary_infections, hcw_to_pat_infections, pat_to_hcw_infections, pat_to_pat_infections, hcw_to_hcw_infections, secondary_patient_infected, secondary_hcw_infected, R0, attack_curve, pat_to_hcw, hcw_to_pat, hcw_to_hcw, iters, discharged, total_pat)
#' @export
simulate <- function(ml_data,
                     mod,
                     timestep = 24L,
                     agents,
                     contact_matrix,
                     initial_set,
                     plot_now = FALSE,
                     num,
                     secondary_infections,
                     hcw_to_pat_infections,
                     pat_to_hcw_infections,
                     pat_to_pat_infections,
                     hcw_to_hcw_infections,
                     secondary_patient_infected,
                     secondary_hcw_infected,
                     R0,
                     attack_curve,
                     pat_to_hcw,
                     hcw_to_pat,
                     hcw_to_hcw,
                     iters,
                     discharged,
                     total_pat,
                    transmission_dir)
{
  # For loop runs until out of timesteps.
  for (i in 1:timestep) {

    # Call Rcpp Function
    Rcpp::sourceCpp(paste0(transmission_dir, "Transmission_Function.cpp"))
    disease_stats       <- propogate(agents, contact_matrix, initial_set)

    infection_stats     <- unlist(tail(disease_stats, n = 8))
    new_disease_stats   <- unlist(disease_stats[-(length(disease_stats)-7):-length(disease_stats)])
    new_disease_stats[which(new_disease_stats == 2)] <- 1L
    agents$disease_stat <- new_disease_stats

    # Update current length of stay
    for(j in 1:length(agents$los_current)) {
      if(agents$los_current[j] > 0L) {
        agents$los_current[j] = agents$los_current[j] - 1L
      }
      else if(agents$los_current[j] == 0L & agents$type[j] == 1L) {
        discharged     <- discharged + 1L
        total_pat      <- total_pat + 1L
        tibble_new_pat <- patients_table(ml_data = ml_data, all_data = data, npat = 1L,
                                     initial_infect_pat = 0L, model = mod)

        agents$id[j]           <- j
        agents$los_start[j]    <- as.numeric(tibble_new_pat$los_start[1])
        agents$sus[j]          <- as.numeric(tibble_new_pat$sus[1])
        agents$infect[j]       <- as.numeric(tibble_new_pat$infect[1])
        agents$disease_stat[j] <- as.numeric(tibble_new_pat$disease_stat[1])
        agents$necode[j]       <- as.numeric(tibble_new_pat$necode[1])
        agents$female[j]       <- as.numeric(tibble_new_pat$female[1])
        agents$los_current[j]  <- as.numeric(tibble_new_pat$los_current[1])
        agents$ndx[j]          <- as.numeric(tibble_new_pat$ndx[1])
        agents$prev_cdi[j]     <- as.numeric(tibble_new_pat$prev_cdi[1])
        agents$Yes[j]          <- as.numeric(tibble_new_pat$Yes[1])
        agents$mobility[j]     <- as.numeric(tibble_new_pat$mobility[1])
        agents$elix[j]         <- as.numeric(tibble_new_pat$elix[1])
        agents$hospitalunit[j] <- as.numeric(tibble_new_pat$hospitalunit[1])
        agents$npr[j]          <- as.numeric(tibble_new_pat$npr[1])
        agents$orproc[j]       <- as.numeric(tibble_new_pat$orproc[1])
        agents$prior_visit[j]  <- as.numeric(tibble_new_pat$prior_visit[1])
        agents$age[j]          <- as.numeric(tibble_new_pat$age[1])
        
        for (k in 1:length(agents$age)) {
          
          if (j == k) {
            contact_matrix[j, k] <- 0L
          }
          else if (agents$type[j] == 1L & agents$type[k] == 1L) {
            contact_matrix[j, k] <- 0L
          }
          else if (agents$type[j] == 1L & agents$type[k] == 2L & (agents$hospitalunit[j] != agents$hospitalunit[k])) {
            contact_matrix[j, k] <- 0L
          }
          else if (agents$type[j] == 2L & agents$type[k] == 1L & (agents$hospitalunit[j] != agents$hospitalunit[k])) {
            contact_matrix[i, j] <- 0L
          }
          else if (agents$type[j] == 2L & agents$type[k] == 2L & (agents$hospitalunit[j] != agents$hospitalunit[k])) {
            contact_matrix[j, k] <- sample(1L:2L, 1L)
          }
          else if ((agents$mobility[j] == 0L & agents$type[k] == 2L)) {
            contact_matrix[j, k] <- sample(0L:1L, 1L)
          }
          else {
            contact_matrix[j, k] <- sample(1L:4L, 1L)
          }
          
        }

        contact_matrix[, j] <- contact_matrix[j, ]
      }
      else {
        gandalf <- 1L
      }
    }
    
    # Compile results
    secondary_infections       <- secondary_infections + infection_stats[1]
    hcw_to_pat_infections      <- hcw_to_pat_infections + infection_stats[2]
    pat_to_hcw_infections      <- pat_to_hcw_infections + infection_stats[3]
    pat_to_pat_infections      <- pat_to_pat_infections + infection_stats[4]
    hcw_to_hcw_infections      <- hcw_to_hcw_infections + infection_stats[5]
    secondary_patient_infected <- secondary_patient_infected + infection_stats[6]
    secondary_hcw_infected     <- secondary_hcw_infected + infection_stats[7]
    R0                         <- R0 + infection_stats[8]

    attack_curve <- c(attack_curve, secondary_infections)
    pat_to_hcw <- c(pat_to_hcw, pat_to_hcw_infections)
    hcw_to_pat <- c(hcw_to_pat, hcw_to_pat_infections)
    hcw_to_hcw <- c(hcw_to_hcw, hcw_to_hcw_infections)
    iters <- c(iters, i)

    num = num + 1L

    if (plot_now == TRUE) {
      Sys.sleep(0.5)
      par(bg = "black", fg = "white")
      plot(x = iters, y = attack_curve, type = "b", main = "%Attack Rate Over Time", xlab = "Iteration",
           ylab = "%Total", col = "red", panel.first = grid(),
           col.main = "red", col.lab = "red", col.axis = "red")

    }

  }

  return(list(num, secondary_infections, hcw_to_pat_infections, pat_to_hcw_infections,
              pat_to_pat_infections, hcw_to_hcw_infections, secondary_patient_infected,
              secondary_hcw_infected, R0, attack_curve, pat_to_hcw, hcw_to_pat,
              hcw_to_hcw, iters, discharged, total_pat))

}

#------------------------------------------------------------------------------

