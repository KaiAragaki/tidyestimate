# tidy_estimate

## Statement of Need

The ESTIMATE package has been fundemental for many analyses, but its documentation is lacking, and its functions sometimes overstep their bounds while not doing enough. This package is a clone of ESTIMATE with an end goal of maintaining the excellent backbone of the package while increasing its documentation and function scope.

## Notes

Yep, I'm using a README for notes.

### outputGCT

camelCase to snake_case

The fact that it can take both a file or a dataframe seems strange.
* Will only allow DF

Should mention the following about input data:

Expects rownames (may change this down the line - perhaps like a DESeq-esque tidy flag)

rownames are put both as NAME and Description (capitalization all over the place) - we only need one. Maybe this is to fit in with a GSEA like format, but it's not doing us any services anymore.

Really confused as to what is going on in rows 10:14

A lot of this can be piped for clarity. Since R 4.1 is introducing a native pipe I think I'll do that and set a dependency to R 4.1. Don't expect this package to be used by really anyone, so it's not going to be an issue

Function should make an object, not a file.

No idea what's going on with lines 17 and 18. Are what appear to be colnames being added levels? That doesn't seem right.

The crux of the issue here is that it's trying to be compatible with GCT when that's not what my plan is. A lot of this can be stripped out.
