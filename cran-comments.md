## R CMD check results

0 errors | 0 warnings | 1 note

* This is a First release.

- Examples with CPU (user + system) or elapsed time > 5s
-> Answer: This is an intrinsic problem with metabarcoding analyses that deal with large datasets. I've often reduced the datasets for the examples, but sometimes the use of quite large dataset is necessary. The function ancombc_pq() and build_phytree_pq() relies on already time-consuming functions from other packages. The total examples run time is 166s on my computer which seems to be reasonable given the high number of examples in this package.

