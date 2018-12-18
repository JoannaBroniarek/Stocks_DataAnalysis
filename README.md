# Stock, Dependency and Graphs
**December 3, 2018**

Project for the course Statistics in Data Science.
Authors: Joanna Broniarek and Guilherme Vescovi Nicchio

![sp500](https://www.avatrade.com/wp-content/uploads/2018/04/SP-500.jpg)

## HTML preview.
The whole analysis along with all generated graphs is available via: 

http://htmlpreview.github.io/?https://github.com/JoannaBroniarek/Stocks_DataAnalysis/blob/master/homework_2_realestate.html


## Study the dependency among stocks via marginal correlation graphs
The aim was to study the dependency among some standard measure of stock relative performance. Data was collected for the daily closing prices for six stocks, selected within those consistently in the S&P500 index from January 1, 2003 through today.

The stocks are categorized into 10 Global Industry Classification Standard sectors, including Consumer Discretionary,
Energy, Financials, Consumer Staples, Telecommunications Services, Health Care, Industrials, Information
Technology, Materials, and Utilities. It is expected that stocks from the same GICS sectors should tend to be clustered together, since stocks from the same gics sector tend to interact more with each other. This was the hypothesis to verify.

## Steps of Analysis
1. Selecting a sensible portfolio of stocks and taking data from January 1, 2003 through January 1, 2008, before the onset of the
“financial crisis”. 
2. Calculating the Pearson correlation coefficient between stocks, and implementing the bootstrap
procedure to build **marginal correlation graphs**
3. Building a marginal correlation graph based on the distance covariance (the multiple hypothesis testing - with and without Bonferroni correction)
4. For the same portfolio of stocks, grabing data from January 1, 2013 through January 1, 2018. Repeating the previous analysis.

### Technology
The whole analysis was provided in R (RStudio).
