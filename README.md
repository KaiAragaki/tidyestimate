# tidyestimate <img src='man/figures/logo.png' align="right" height="139" />

## Statement of Need

The ESTIMATE package has been fundamental for many analyses, but its documentation is lacking, and its functions sometimes overstep their bounds while not doing enough. This package is a refresh of ESTIMATE with the goal of maintaining the excellent backbone of the package while increasing its documentation and function scope.

Original paper [here](https://www.nature.com/articles/ncomms3612).

## Features

|            |          tidyestimate|   ESTIMATE|
|-----------:|:--------------------:|:---------:|
|       input|`data.frame`<br />`tibble`<br />`matrix`|`.GCT` file|
|      output|          `data.frame`|`.GCT` file|
|`%>%`/`\|>`?|                    ✔️|         ✖️|
|        size|                  <1MB|        ~7MB|

Additionally:

⚡ Faster. `tidyestimate` doesn't do any file conversion.

📝 Better documentation. Functions are more clear about input requirements and returns.

🕊️ Lighter. Less code, more readable (less to break, easier to fix).

💪 Robust. `tidyestimate` does conservative alias matching to allow compatibility with both old and new gene identifiers.
