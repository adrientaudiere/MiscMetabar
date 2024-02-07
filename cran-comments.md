## R CMD check results

0 errors | 0 warnings | 2 note

* This is a First release.


Notes: 

- Found the following (possibly) invalid URLs:
     URL: 
       From: inst/doc/MiscMetabar.html
       Message: Empty URL
       
-> Answer: I didn't find the place of this url.

- Examples with CPU (user + system) or elapsed time > 5s
-> Answer: This is an intrinsic problem with metabarcoding analyses that deal with large datasets. I've often reduced the datasets for the examples, but sometimes the use of full data is necessary. The total examples run time is 246s on my computer which seems to be reasonable given the high number of examples in this package.

