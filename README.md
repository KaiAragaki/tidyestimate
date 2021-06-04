# tidyestimate

## Statement of Need

The ESTIMATE package has been fundamental for many analyses, but its documentation is lacking, and its functions sometimes overstep their bounds while not doing enough. This package is a clone of ESTIMATE with an end goal of maintaining the excellent backbone of the package while increasing its documentation and function scope.

Original paper [here](https://www.nature.com/articles/ncomms3612)

## Differences

|            |          tidyestimate|   ESTIMATE|
|-----------:|:--------------------:|:---------:|
|       input|`data.frame`, `tibble`|`.GCT` file|
|      output|          `data.frame`|`.GCT` file|
|`%>%`/`\|>`?|                    ✔️|         ✖️|

Additionally:
* Functions are better documented (more clear about what they want as inputs)
* Reduced and more readable code (less to break, easier to fix)
