# Man directory for package documentation

This directory contains .Rd (R documentation) files that document your package functions.

Each exported function should have a corresponding .Rd file.

These files are typically generated automatically from roxygen2 comments in your R code using:
\`\`\`r
roxygen2::roxygenise()
\`\`\`
