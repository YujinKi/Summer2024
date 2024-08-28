# Post-Selection Inference in Linear Regression: The Impact of Selection Methods in LASSO

<p xmlns:cc="http://creativecommons.org/ns#" xmlns:dct="http://purl.org/dc/terms/"><a property="dct:title" rel="cc:attributionURL" href="https://github.com/YujinKi/Summer2024">Post-Selection Inference in Linear Regression: The Impact of Selection Methods in LASSO</a> by <a rel="cc:attributionURL dct:creator" property="cc:attributionName" href="https://github.com/YujinKi">Yujin Ki</a> is licensed under <a href="https://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1" alt=""><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1" alt=""></a></p>

The project "Post-Selection Inference in Linear Regression: The Impact of Selection Methods in LASSO" is conducted as a summer project for MSc in Statistics at Imperial College London. In this project, five different scenarios have been considered - 4 standard cases and an extreme case. The r codes and figures can be found in this GitHub repository. 


## R codes 

There are five R files which are the codes for each corresponding scenario - all non-zero true coefficients, $p$ approaching $n$, $p$ exceeding $n$, 5-fold cross-validation, and $p$ exceeding $n$ with all non-zero true coefficients where $n$ is the number of observations and $p$ is the number of covariates. 


The following table presents brief information on each scenario but the detailed scenarios can be seen in the thesis. 

| Scenarios                                               | $n$ (Number of <br> observations) | $p$ (Number of <br> variables)| Number of <br> non-zero coefficients | Number of <br> cross-validation folds | 
| :-----------------------------------------------------: | :-: | :-: | :----------------------------:  | :--------------------------: | 
| all non-zero true coefficient                           | 60  | 18  | 18                              | 10 | 
| $p$ approaching $n$                                     | 60  | 45  | 10                              | 10 |
| $p$ exceeding $n$                                       | 60  | 80  | 10                              | 10 |
| 5-fold cross-validation                                 | 60  | 45  | 10                              | 5 |
| $p$ exceeding $n$ <br> with all non-zero true coefficients| 60  | 80  | 80                              | 10 |



## Figures 

Six folders in the main Figure folder, i.e. six types of figures are considered. Each folder contains figures for four or five different scenarios. The scenarios are stated as numbers. The numbers corresponding to each scenario are as follows:

1. all non-zero true coefficient
2. $p$ approaching $n$
3. $p$ exceeding $n$
4. 5-fold cross-validation
5. $p$ exceeding $n$ with all non-zero true coefficients


For reference, the histogram of selected models is only for the last scenario: $p$ exceeding $n$ with all non-zero true coefficients. 

