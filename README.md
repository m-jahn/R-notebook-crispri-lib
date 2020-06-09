# R-notebooks

### Introduction | Disclaimer

R markdown notebooks to document my computational biology projects.
This repository contains data processing pipelines or data curation projects
that require a higher level of maintenance and documentation. 

In contrast to some of my other github projects, 
like [SysbioTreemaps](https://github.com/m-jahn/SysbioTreemaps),
[ShinyMC](https://github.com/m-jahn/ShinyMC), [ShinyProt](https://github.com/m-jahn/ShinyProt),
or [ShinyLib](https://github.com/m-jahn/ShinyLib),
the pipelines in this repository are not necessarily in the form of R packages or apps.
This means the pipelines were created for a particular purpose or data analysis task,
but are not further developed or maintained.

The pipelines and their documentation collected here therefore come with no warranty
regarding stability or runability in future R versions. However, all care was taken to
guarantee scientific accuracy and adhere to good scientific practice in terms of statistics,
reproducibility and code documentation. Please report any errors by filing a github issue.

### How to run the pipelines

The pipelines collected in this repository are self-contained and executable. 
The code _and_ the documentation are part of one and the same R markdown document
for each pipeline. The pipelines themselves can be downloaded and executed 
from the `pipeline` sub-folders. To simply view the rendered pipelines follow 
the links to the `*.html` reports under   [Contents](#Contents).

To download the repository on your local drive use `git clone` in a (linux) terminal:

``` bash
cd /your-target-folder
git clone https://github.com/m-jahn/R-notebooks.git
```

Open a pipeline with Rstudio and execute code (chunks) with the `Run` button.
Alternatively, open an interactive R session and render the R markdown pipeline:

``` bash
require(rmarkdown)
rmarkdown::render("pipeline.Rmd")
```


### Contents

#### Ralstonia resource allocation (in preparation)

- [_Ralstonia eutropha_ COG re-annotation](https://m-jahn.github.io/R-notebooks/Ralstonia_H16_genome_re_annotation.nb.html)
- [_Ralstonia eutropha_ gene essentiality analysis](https://m-jahn.github.io/R-notebooks/Ralstonia_H16_essentiality_analysis.nb.html)
- [_Ralstonia eutropha_ model constraints](https://m-jahn.github.io/R-notebooks/Ralstonia_model_constraints.nb.html)
- [_Ralstonia eutropha_ protein utilization](https://m-jahn.github.io/R-notebooks/Ralstonia_variability_analysis.nb.html)

#### CRISPRi library in _Synechocystis_ sp. PCC6803 (published in Yao *et al*., Nature Communications, 2020)

- [_Synechocystis_ CRISPRi library data processing](https://m-jahn.github.io/R-notebooks/CRISPRi_library_data_processing.nb.html)
- [_Synechocystis_ CRISPRi library enrichment analysis](https://m-jahn.github.io/R-notebooks/CRISPRi_library_enrichment_analysis.nb.html)
- [_Synechocystis_ CRISPRi library additional tests](https://m-jahn.github.io/R-notebooks/CRISPRi_library_additional_tests.nb.html)
