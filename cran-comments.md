# v0.1.0
### round 1

> Found the following (possibly) invalid URLs:
  URL: http://roliveros-ramos.github.io/ibm
    From: DESCRIPTION
    Status: 404
    Message: Not Found

The link is working now.

> The Description field should not start with the package name,
  'This package' or similar.

DESCRIPTION file has been modified following your suggestion.

> * checking top-level files ... NOTE
Non-standard file/directory found at top level:
  'cran-comments.md'

.Rbuildignore has been updated.

> * checking R code for possible problems ... NOTE
.summary: no visible binding for global variable 'median'
Undefined global functions or variables:
  median
Consider adding
  importFrom("stats", "median")
to your NAMESPACE file.

Done.

> Please fix and resubmit.
Best,
Uwe Ligges