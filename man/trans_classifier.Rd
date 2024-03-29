% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trans_classifier.R
\name{trans_classifier}
\alias{trans_classifier}
\title{Create trans_classifier object for machine-learning-based model prediction.}
\description{
This class is a wrapper for methods of machine-learning-based classification or regression models, including data pre-processing, feature selection, 
data split, model training, prediction, confusionMatrix and ROC (Receiver Operator Characteristic) or PR (Precision-Recall) curve.

Author(s): Felipe Mansoldo and Chi Liu
}
\examples{

## ------------------------------------------------
## Method `trans_classifier$new`
## ------------------------------------------------

\donttest{
data(dataset)
t1 <- trans_classifier$new(
		dataset = dataset, 
		x.predictors = "Genus",
		y.response = "Group")
}

## ------------------------------------------------
## Method `trans_classifier$cal_preProcess`
## ------------------------------------------------

\dontrun{
t1$cal_preProcess(method = c("center", "scale", "nzv"))
}

## ------------------------------------------------
## Method `trans_classifier$cal_feature_sel`
## ------------------------------------------------

\dontrun{
t1$cal_feature_sel(boruta.maxRuns = 300, boruta.pValue = 0.01)
}

## ------------------------------------------------
## Method `trans_classifier$cal_split`
## ------------------------------------------------

\dontrun{
t1$cal_split(prop.train = 3/4)
}

## ------------------------------------------------
## Method `trans_classifier$set_trainControl`
## ------------------------------------------------

\dontrun{
t1$set_trainControl(method = 'repeatedcv')
}

## ------------------------------------------------
## Method `trans_classifier$cal_train`
## ------------------------------------------------

\dontrun{
# random forest
t1$cal_train(method = "rf")
# Support Vector Machines with Radial Basis Function Kernel
t1$cal_train(method = "svmRadial", tuneLength = 15)
}

## ------------------------------------------------
## Method `trans_classifier$cal_feature_imp`
## ------------------------------------------------

\dontrun{
t1$cal_feature_imp()
}

## ------------------------------------------------
## Method `trans_classifier$plot_feature_imp`
## ------------------------------------------------

\dontrun{
t1$plot_feature_imp(use_number = 1:20, coord_flip = FALSE)
}

## ------------------------------------------------
## Method `trans_classifier$cal_predict`
## ------------------------------------------------

\dontrun{
t1$cal_predict()
}

## ------------------------------------------------
## Method `trans_classifier$plot_confusionMatrix`
## ------------------------------------------------

\dontrun{
t1$plot_confusionMatrix()
}

## ------------------------------------------------
## Method `trans_classifier$cal_ROC`
## ------------------------------------------------

\dontrun{
t1$cal_ROC()
}

## ------------------------------------------------
## Method `trans_classifier$plot_ROC`
## ------------------------------------------------

\dontrun{
t1$plot_ROC(size = 1, alpha = 0.7)
}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-trans_classifier-new}{\code{trans_classifier$new()}}
\item \href{#method-trans_classifier-cal_preProcess}{\code{trans_classifier$cal_preProcess()}}
\item \href{#method-trans_classifier-cal_feature_sel}{\code{trans_classifier$cal_feature_sel()}}
\item \href{#method-trans_classifier-cal_split}{\code{trans_classifier$cal_split()}}
\item \href{#method-trans_classifier-set_trainControl}{\code{trans_classifier$set_trainControl()}}
\item \href{#method-trans_classifier-cal_train}{\code{trans_classifier$cal_train()}}
\item \href{#method-trans_classifier-cal_feature_imp}{\code{trans_classifier$cal_feature_imp()}}
\item \href{#method-trans_classifier-plot_feature_imp}{\code{trans_classifier$plot_feature_imp()}}
\item \href{#method-trans_classifier-cal_predict}{\code{trans_classifier$cal_predict()}}
\item \href{#method-trans_classifier-plot_confusionMatrix}{\code{trans_classifier$plot_confusionMatrix()}}
\item \href{#method-trans_classifier-cal_ROC}{\code{trans_classifier$cal_ROC()}}
\item \href{#method-trans_classifier-plot_ROC}{\code{trans_classifier$plot_ROC()}}
\item \href{#method-trans_classifier-clone}{\code{trans_classifier$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-trans_classifier-new"></a>}}
\if{latex}{\out{\hypertarget{method-trans_classifier-new}{}}}
\subsection{Method \code{new()}}{
Create the trans_classifier object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{trans_classifier$new(
  dataset = NULL,
  x.predictors = "all",
  y.response = NULL,
  n.cores = 1
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{dataset}}{the object of \code{\link{microtable}} Class.}

\item{\code{x.predictors}}{default "all"; character string or data.frame; a character string represents selecting the corresponding data from microtable$taxa_abund; 
data.frame represents other customized input. See the following available options:
\describe{
  \item{\strong{'all'}}{use all the taxa stored in microtable$taxa_abund}
  \item{\strong{'Genus'}}{use Genus level table in microtable$taxa_abund, or other specific taxonomic rank, e.g. 'Phylum'}
  \item{\strong{other input}}{must be a data.frame; It should have the same format with the data.frame in microtable$taxa_abund, i.e. rows are features; 
    cols are samples with same names in sample_table}
}}

\item{\code{y.response}}{default NULL; the response variable in sample_table.}

\item{\code{n.cores}}{default 1; the CPU thread used.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
data_feature and data_response in the object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\donttest{
data(dataset)
t1 <- trans_classifier$new(
		dataset = dataset, 
		x.predictors = "Genus",
		y.response = "Group")
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-trans_classifier-cal_preProcess"></a>}}
\if{latex}{\out{\hypertarget{method-trans_classifier-cal_preProcess}{}}}
\subsection{Method \code{cal_preProcess()}}{
Pre-process (centering, scaling etc.) of the feature data based on the caret::preProcess function. 
	 See \href{https://topepo.github.io/caret/pre-processing.html}{https://topepo.github.io/caret/pre-processing.html} for more details.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{trans_classifier$cal_preProcess(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{parameters pass to preProcess function of caret package.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
converted data_feature in the object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
t1$cal_preProcess(method = c("center", "scale", "nzv"))
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-trans_classifier-cal_feature_sel"></a>}}
\if{latex}{\out{\hypertarget{method-trans_classifier-cal_feature_sel}{}}}
\subsection{Method \code{cal_feature_sel()}}{
Perform feature selection.
	 See \href{https://topepo.github.io/caret/feature-selection-overview.html}{https://topepo.github.io/caret/feature-selection-overview.html} for more details.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{trans_classifier$cal_feature_sel(
  boruta.maxRuns = 300,
  boruta.pValue = 0.01,
  boruta.repetitions = 4,
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{boruta.maxRuns}}{default 300; maximal number of importance source runs; passed to the maxRuns parameter in Boruta function of Boruta package.}

\item{\code{boruta.pValue}}{default 0.01; p value passed to the pValue parameter in Boruta function of Boruta package.}

\item{\code{boruta.repetitions}}{default 4; repetition runs for the feature selection.}

\item{\code{...}}{parameters pass to Boruta function of Boruta package.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
optimized data_feature in the object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
t1$cal_feature_sel(boruta.maxRuns = 300, boruta.pValue = 0.01)
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-trans_classifier-cal_split"></a>}}
\if{latex}{\out{\hypertarget{method-trans_classifier-cal_split}{}}}
\subsection{Method \code{cal_split()}}{
Split data for training and testing.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{trans_classifier$cal_split(prop.train = 3/4)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{prop.train}}{default 3/4; the ratio of the dataset used for the training.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
data_train and data_test in the object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
t1$cal_split(prop.train = 3/4)
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-trans_classifier-set_trainControl"></a>}}
\if{latex}{\out{\hypertarget{method-trans_classifier-set_trainControl}{}}}
\subsection{Method \code{set_trainControl()}}{
Control parameters for the following training. See trainControl function of caret package for details.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{trans_classifier$set_trainControl(
  method = "repeatedcv",
  classProbs = TRUE,
  savePredictions = TRUE,
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{method}}{default 'repeatedcv'; 'repeatedcv': Repeated k-Fold cross validation; 
see method parameter in \code{trainControl} function of \code{caret} package for available options.}

\item{\code{classProbs}}{default TRUE; should class probabilities be computed for classification models?;
see classProbs parameter in \code{caret::trainControl} function.}

\item{\code{savePredictions}}{default TRUE; see \code{savePredictions} parameter in \code{caret::trainControl} function.}

\item{\code{...}}{parameters pass to trainControl function of caret package.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
trainControl in the object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
t1$set_trainControl(method = 'repeatedcv')
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-trans_classifier-cal_train"></a>}}
\if{latex}{\out{\hypertarget{method-trans_classifier-cal_train}{}}}
\subsection{Method \code{cal_train()}}{
Run the model training.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{trans_classifier$cal_train(method = "rf", max.mtry = 2, max.ntree = 200, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{method}}{default "rf"; "rf": random forest; see method in caret::train function for other options.}

\item{\code{max.mtry}}{default 2; for method = "rf"; maximum mtry used for the tunegrid to do hyperparameter tuning to optimize the model.}

\item{\code{max.ntree}}{default 200; for method = "rf"; maximum number of trees used to optimize the model.}

\item{\code{...}}{parameters pass to \code{caret::train} function.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
res_train in the object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
# random forest
t1$cal_train(method = "rf")
# Support Vector Machines with Radial Basis Function Kernel
t1$cal_train(method = "svmRadial", tuneLength = 15)
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-trans_classifier-cal_feature_imp"></a>}}
\if{latex}{\out{\hypertarget{method-trans_classifier-cal_feature_imp}{}}}
\subsection{Method \code{cal_feature_imp()}}{
Get feature importance from the training model.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{trans_classifier$cal_feature_imp(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{parameters pass to varImp function of caret package.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
res_feature_imp in the object. One row for each predictor variable. The column(s) are different importance measures.
  For the method 'rf', it is MeanDecreaseGini (classification) or IncNodePurity (regression).
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
t1$cal_feature_imp()
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-trans_classifier-plot_feature_imp"></a>}}
\if{latex}{\out{\hypertarget{method-trans_classifier-plot_feature_imp}{}}}
\subsection{Method \code{plot_feature_imp()}}{
Bar plot for feature importance.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{trans_classifier$plot_feature_imp(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{parameters pass to \code{plot_diff_bar} function of \code{trans_diff} package.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
ggplot2 object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
t1$plot_feature_imp(use_number = 1:20, coord_flip = FALSE)
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-trans_classifier-cal_predict"></a>}}
\if{latex}{\out{\hypertarget{method-trans_classifier-cal_predict}{}}}
\subsection{Method \code{cal_predict()}}{
Run the prediction.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{trans_classifier$cal_predict(positive_class = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{positive_class}}{default NULL; see positive parameter in confusionMatrix function of caret package;
If positive_class is NULL, use the first group in data as the positive class automatically.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
res_predict, res_confusion_fit and res_confusion_stats stored in the object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
t1$cal_predict()
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-trans_classifier-plot_confusionMatrix"></a>}}
\if{latex}{\out{\hypertarget{method-trans_classifier-plot_confusionMatrix}{}}}
\subsection{Method \code{plot_confusionMatrix()}}{
Plot the cross-tabulation of observed and predicted classes with associated statistics.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{trans_classifier$plot_confusionMatrix(
  plot_confusion = TRUE,
  plot_statistics = TRUE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{plot_confusion}}{default TRUE; whether plot the confusion matrix.}

\item{\code{plot_statistics}}{default TRUE; whether plot the statistics.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
ggplot object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
t1$plot_confusionMatrix()
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-trans_classifier-cal_ROC"></a>}}
\if{latex}{\out{\hypertarget{method-trans_classifier-cal_ROC}{}}}
\subsection{Method \code{cal_ROC()}}{
Get ROC (Receiver Operator Characteristic) curve data and the performance data.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{trans_classifier$cal_ROC(input = "pred")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{input}}{default "pred"; 'pred' or 'train'; 'pred' represents using prediction results;
'train' represents using training results.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a list res_ROC stored in the object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
t1$cal_ROC()
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-trans_classifier-plot_ROC"></a>}}
\if{latex}{\out{\hypertarget{method-trans_classifier-plot_ROC}{}}}
\subsection{Method \code{plot_ROC()}}{
Plot ROC curve.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{trans_classifier$plot_ROC(
  plot_type = c("ROC", "PR")[1],
  plot_group = "all",
  color_values = RColorBrewer::brewer.pal(8, "Dark2"),
  add_AUC = TRUE,
  plot_method = FALSE,
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{plot_type}}{default c("ROC", "PR")[1]; 'ROC' represents ROC (Receiver Operator Characteristic) curve; 
'PR' represents PR (Precision-Recall) curve.}

\item{\code{plot_group}}{default "all"; 'all' represents all the classes in the model;
'add' represents all adding micro-average and macro-average results, see 
\href{https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc.html}{https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc.html};
other options should be one or more class names, same with the names in Group column of res_ROC$res_roc from cal_ROC function.}

\item{\code{color_values}}{default RColorBrewer::brewer.pal(8, "Dark2"); colors used in the plot.}

\item{\code{add_AUC}}{default TRUE; whether add AUC in the legend.}

\item{\code{plot_method}}{default FALSE; If TRUE, show the method in the legend though only one method is found.}

\item{\code{...}}{parameters pass to geom_path function of ggplot2 package.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
ggplot2 object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
t1$plot_ROC(size = 1, alpha = 0.7)
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-trans_classifier-clone"></a>}}
\if{latex}{\out{\hypertarget{method-trans_classifier-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{trans_classifier$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
