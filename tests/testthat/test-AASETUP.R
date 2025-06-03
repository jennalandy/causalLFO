# NMF::nmf() internally uses setupLibPaths("NMF"), which calls path.package("NMF").
# This requires the NMF package to be attached, not just imported.
suppressPackageStartupMessages(library(NMF))
