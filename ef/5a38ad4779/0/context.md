# Session Context

## User Prompts

### Prompt 1

The github action @.github/workflows/test-vsearch-windows.yaml failed with message : 

Run testthat::test_file(
✔ | F W  S  OK | Context

⠏ |          0 | vsearch                                                        
✖ | 1 1      0 | vsearch
────────────────────────────────────────────────────────────────────────────────
Warning ('test_vsearch.R:7:1'): (code run outside of `test_that()`)
data set 'data_fungi' not found
Backtrace:
    ▆
 1. └─utils::data("data_fungi") at test_vsearch.R:7:1

...

