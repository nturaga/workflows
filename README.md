Errors and warnings
===

## rnaseqGene errors

```
~/D/workflows ❯❯❯ R CMD INSTALL rnaseqGene_1.0.0.tar.gz                           master ✖ ✱ ◼
Bioconductor version 3.5 (BiocInstaller 1.26.1), ?biocLite for help
* installing to library '/Library/Frameworks/R.framework/Versions/3.4/Resources/library'
* installing *source* package 'rnaseqGene' ...
** inst
** help
No man pages found in package  'rnaseqGene'
*** installing help indices
** building package indices
** installing vignettes
** testing if installed package can be loaded
Bioconductor version 3.5 (BiocInstaller 1.26.1), ?biocLite for help
No methods found in "RSQLite" for requests: dbGetQuery
Warning: replacing previous import 'GenomicAlignments::first' by 'dplyr::first' when loading 'rnaseqGene'
Warning: replacing previous import 'GenomicAlignments::last' by 'dplyr::last' when loading 'rnaseqGene'
Warning: replacing previous import 'dplyr::select' by 'AnnotationDbi::select' when loading 'rnaseqGene'
Warning: replacing previous import 'Rsamtools::path' by 'ReportingTools::path' when loading 'rnaseqGene'
* DONE (rnaseqGene)
```


## maEndToEnd

```
~/D/workflows ❯❯❯ R CMD INSTALL maEndToEnd_1.0.0.tar.gz                           master ✖ ✱ ◼
Bioconductor version 3.5 (BiocInstaller 1.26.1), ?biocLite for help
* installing to library '/Library/Frameworks/R.framework/Versions/3.4/Resources/library'
* installing *source* package 'maEndToEnd' ...
** inst
** help
No man pages found in package  'maEndToEnd'
*** installing help indices
** building package indices
** testing if installed package can be loaded
Bioconductor version 3.5 (BiocInstaller 1.26.1), ?biocLite for help
No methods found in "RSQLite" for requests: dbGetQuery

groupGOTerms: 	GOBPTerm, GOMFTerm, GOCCTerm environments built.
No methods found in "RSQLite" for requests: dbGetQuery
Warning: replacing previous import 'oligo::summarize' by 'dplyr::summarize' when loading 'maEndToEnd'
Warning: replacing previous import 'Biobase::combine' by 'dplyr::combine' when loading 'maEndToEnd'
Warning: replacing previous import 'Biobase::anyMissing' by 'matrixStats::anyMissing' when loading 'maEndToEnd'
Warning: replacing previous import 'dplyr::count' by 'matrixStats::count' when loading 'maEndToEnd'
Warning: replacing previous import 'Biobase::rowMedians' by 'matrixStats::rowMedians' when loading 'maEndToEnd'
Warning: replacing previous import 'matrixStats::rowSds' by 'genefilter::rowSds' when loading 'maEndToEnd'
Warning: replacing previous import 'matrixStats::rowVars' by 'genefilter::rowVars' when loading 'maEndToEnd'
Warning: replacing previous import 'dplyr::select' by 'biomaRt::select' when loading 'maEndToEnd'
Warning: replacing previous import 'oligo::backgroundCorrect' by 'limma::backgroundCorrect' when loading 'maEndToEnd'
Warning: replacing previous import 'geneplotter::plotMA' by 'limma::plotMA' when loading 'maEndToEnd'
Warning: replacing previous import 'oligo::normalize' by 'EnrichmentBrowser::normalize' when loading 'maEndToEnd'
* DONE (maEndToEnd)
```

## Proteomics

```
pool-1-thread-4: Computing spectral E-values finished (elapsed time: 58.00 sec)
pool-1-thread-3: Computing spectral E-values finished (elapsed time: 67.00 sec)
Computing q-values...
Computing q-values finished (elapsed time: 0.05 sec)
Writing results...
Writing results finished (elapsed time: 4.74 sec)
MS-GF+ complete (total elapsed time: 132.57 sec)
Assertion failure at kmp_runtime.cpp(6480): __kmp_thread_pool == __null.
Assertion failure at kmp_runtime.cpp(6480): __kmp_thread_pool == __null.
OMP: Error #13: Assertion failure at kmp_runtime.cpp(6480).
OMP: Hint: Please submit a bug report with this message, compile and run commands used, and machine configuration info including native compiler and operating system versions. Faster response will be obtained by including all program sources. For information on submitting this issue, please see http://www.intel.com/software/products/support/.
OMP: Error #13: Assertion failure at kmp_runtime.cpp(6480).
OMP: Hint: Please submit a bug report with this message, compile and run commands used, and machine configuration info including native compiler and operating system versions. Faster response will be obtained by including all program sources. For information on submitting this issue, please see http://www.intel.com/software/products/support/.
Quitting from lines 473-477 (proteomics.Rmd)
Error: processing vignette 'proteomics.Rmd' failed with diagnostics:
non-numeric argument to mathematical function
Execution halted
```

## SimpleSingleCell

```
Bioconductor version 3.5 (BiocInstaller 1.26.1), ?biocLite for help
* installing to library '/Library/Frameworks/R.framework/Versions/3.4/Resources/library'
* installing *source* package 'simpleSingleCell' ...
** help
No man pages found in package  'simpleSingleCell'
*** installing help indices
** building package indices
** testing if installed package can be loaded
Bioconductor version 3.5 (BiocInstaller 1.26.1), ?biocLite for help
Error: package or namespace load failed for 'GenomicFeatures' in dyn.load(file, DLLpath = DLLpath, ...):
 unable to load shared object '/Library/Frameworks/R.framework/Versions/3.4/Resources/library/Rsamtools/libs/Rsamtools.so':
  `maximal number of DLLs reached...
Error : package 'GenomicFeatures' could not be loaded
Error: loading failed
Execution halted
ERROR: loading failed
* removing '/Library/Frameworks/R.framework/Versions/3.4/Resources/library/simpleSingleCell'
```