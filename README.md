# Optimal weighted Bonferroni tests and their graphical extension

This is the repository associated with manuscript "Optimal weighted Bonferroni tests and their graphical extension". All source code used in the manuscript are stored in this repo. All table and figures are reproducible and demonstrated in the Rmarkdown files.

## Source code
All source code are located in under folder `src`.

* `src/utils` stores R codes which define the utility functions, like initialization of starting points, etc;
* `src/objectives` stores R codes which define the objective functions, like disjunctive power, conjunctive power;
* `src/optimization` stores R codes which define the optimization algorithms proposed in the manuscript to yields optimal weights for Bonferroni test.

To find examples or use cases for the source code, please refer to the reproducible documents for the implementation.


## Reproducible documents
All table and figures demonstrated in the manuscript can be reproduced by Rmarkdown files located under folder `docs/markdowns`.

* `main_content.Rmd` yields the Figure 1 and Table 1-5 in the main content;
* `graph_optimization.Rmd` explains and demonstrates step-by-step procedure of the optimization of the graphical approach example in Section 5;
* `supplementary.Rmd` yields the tables and figures presented in the supplmentary materials.

All rendered pdf using above Rmarkdown files are also displayed in the folder `docs/markdowns`.
