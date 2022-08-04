# The script should be written here to benchmark the package on Github.
# NOTE: save_data argument must be set to TRUE for plot_metrics functions.

## TEST 1
Rperform::plot_metrics(
    test_path = "tests/testthat/test-TRAVIS-sequentialSearch.R",
    metric = "time", save_data = TRUE,
    save_plots = TRUE
)


## TEST 2
Rperform::plot_metrics(
    test_path = "tests/testthat/test-CRAN-PeakSegFPOP_dir.R",
    metric = "time", save_data = TRUE,
    save_plots = TRUE
)

## TEST 3
# Rperform::time_compare(
#    test_path = "inst/tests/test-dup.r",
#    num_commits = 2,
#    save_data = TRUE
# )

## TEST 4
# Rperform::mem_compare(
#    test_path = "inst/tests/test-dup.r",
#    num_commits = 2,
#    save_data = TRUE
# )

## TEST 5
#  Rperform::plot_branchmetrics(
#    test_path = "inst/tests/test-check.r",
#    metric = "time",
#    branch1 = "rperform_test",
#    branch2 = "master",
#    save_data = T,
#    save_plots = F
# )
