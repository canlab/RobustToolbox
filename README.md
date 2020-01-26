# RobustToolbox

Robust regression is an alternative to Ordinary Least Squares regression. It uses iterative algorithms to identify and down-weight potential outliers. In robust regression, these are defined as observations that lie far from the regression slope (or plane).

• Outliers can cause violations of statistical assumptions, which can increase false positives and false negatives

• Outliers have very large effects on regression coefficients (slopes). These slopes are interpreted as measures of activation or connectivity in neuroimaging studies. Even one or a few outliers can completely alter the results, particularly if they are high-leverage data points.

• Normally, a data analyst would check the assumptions and data distribution in any statistical model, and make adjustments as needed. But this is not possible with neuroimaging and other applications where hundreds or thousands of models are tested in parallel.

• When assumptions cannot be checked for every regression model, automatic procedures for weighting based on outlier status can be advantageous

• Robust regression is an automatic procedure for identifying cases that are potential outliers and down-weighting them.

The robust regression toolbox was created by Tor Wager. If you use it, please see (and consider citing) the accompanying paper:
Wager, T. D., Keller, M. C., Lacey, S. C., & Jonides, J. (2005). Increased sensitivity in neuroimaging analyses using robust regression. Neuroimage, 26(1), 99-113.

## Rationale

One way to think of outliers is as observations that come from a different *generative model* (i.e., different population), with a different distribution, from the rest of the observations. With outliers, this distribution is higher variance than the distribution for the main dataset, causing observations to have more spread. Since the "pull" of observations on the regression line is proportional to the square of the distance from the line, these observations from a high-variance distribution tend to dominate if their observed values are extreme. This "pull" also depends on the *leverage*, which is a function of how extreme an observations predicted value is. High-leverage outliers can completely change regression results even in large samples.  

If some observations come from a different generative model with a different distribution, this violates the IID assumptions underlying classical statistical inference -- i.e., that the observations conditional on the model (i.e., the residuals) are independent and come from an identical distribution. This situation also violates the equality of variance (homoscedasticity) assumption when one is using classical P-values for inference, and can also cause violations of the normality assumption.  

Robust regression uses an iterative algorithm to identify observations with large residuals and down-weight them.
The Robust Toolbox uses the Iteratively Reweighted Least Squares (IRLS) algorithm, with the following steps:

1. Fit the regression model using weighted least squares, with weights set to 1/leverage for each point

2. Normalize the residuals by their Median Absolute Deviation and apply a weight function based on normalized residuals.

3. Fit the regression model using weighted least squares, using weights from Step 2.
   Repeat Steps 2-3 until convergence.

4. Adjust variance, degrees of freedom, and P-values to account for reweighting

This is quite useful for second-level (group) analyses with one contrast image per participant entered as input data. But it can also be used for first-level (time series) analysis.

## Installation and setup

To use the CANlab Robust Regression toolbox, you'll need Matlab and three toolboxes on your Matlab path: SPM12, the CANlab Core Tools repository, and the CANlab Robust Regression toolbox. These tools include the sample dataset used in the robust regression help walkthrough, and other datasets as well.

For help installing CANlab tools, walkthroughs, and more [canlab.github.io](https://canlab.github.io)

CANlab code repositories are at [CANlab Github](https://github.com/canlab)

## Using the toolbox

After installation, there are two ways to use CANlab robust regression within Matlab. The first way uses this toolbox and operates on Nifti (.nii) or Analyze (.img) files, and writes output files to disk. The second way uses CANLab objects in an object-oriented interactive framework. Both methods rely on Matlab's **robustfit** function, which implements the core method. Thus, they will produce the same results. The CANLab Core tools can be used to visualize and make tables of results obtained using either method.

- The Robust toolbox includes a subfolder called **Robust_regression_walkthrough**, with a walkthrough including an example dataset for each use case.


### The classic Robust Regression toolbox

- Using this toolbox takes a set of image filenames and regressors (a GLM design matrix) as input, and runs robust regression at every voxel in the dataset.  
- The main function in the toolbox is called **robfit**. This creates a series of directories, one for each set of images you pass in, and writes images containing maps of T-values, P-values, observation weights, and brain coverage for each analysis. A SETUP.mat file includes the design matrix, image file names, and other meta-data to track what was done.
- **Output**: Each regressor is assigned a number. 0001 is the intercept, and 0002 and on are regressors you enter.  Regressors are mean-centered on average so that the intercept map can be interpreted as the group average. Here are some output files:

| File                  | Description                                                                 |
| --------              | ------------------------------------------------------------                |
| rob_beta_0001.nii     | Intercept activation values (group mean activation)                         |
| rob_tmap_0001.nii     | Intercept t-values                                                          |
| rob_p_0001.nii        | Intercept p-values                                                          |
| rob_beta_0002.nii     | Regression slopes (activation values) for first user-entered regressor      |
| rob_tmap_0002.nii     | t-values for first user-entered regressor                                   |
| rob_p_0002.nii        | p-values for first user-entered regressor                                   |
| weights.nii           | 4-D image of regression weights for each input image                        |

- **robust_results_batch** is a script that loads the files from disk into CANlab objects and generates a series of visualizations and tables.
- **publish_robust_regression_report** is a script that publishes an HTML report with the results of the analysis.
- **robust_regression_walkthrough_toolbox.mlx** is a Matlab live script that walks you through a sample analysis.

Both of these can be customized for your application, and the code contains more information about how to use CANlab tools and generate other kinds of output.

### Object-oriented tools

- The second way uses the regress() method for fmri_data in CANlab object-oriented toolbox. It does not use this toolbox, but uses the same robust regression algorithm. It returns statistic_image class objects, which can be visualized and written to disk (e.g., as Nifti files), but it does not write files to disk by default.

- **robust_regression_walkthrough_objectoriented.mlx** is a Matlab live script that walks you through a sample analysis.
