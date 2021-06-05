# tidyestimate

## Statement of Need

The ESTIMATE package has been fundamental for many analyses, but its documentation is lacking, and its functions sometimes overstep their bounds while not doing enough. This package is a clone of ESTIMATE with an end goal of maintaining the excellent backbone of the package while increasing its documentation and function scope.

Original paper [here](https://www.nature.com/articles/ncomms3612)

## Features

|            |          tidyestimate|   ESTIMATE|
|-----------:|:--------------------:|:---------:|
|       input|`data.frame`, `tibble`|`.GCT` file|
|      output|          `data.frame`|`.GCT` file|
|`%>%`/`\|>`?|                    ✔️|         ✖️|

Additionally:
* Faster. `tidyestimate` doesn't do any file conversion.
* Better documentation. Functions are more clear about input requirements and returns.
* Lighter weight. Ditching `.GCT` files slims the package size by almost half.
* Less code, more readable (less to break, easier to fix).
