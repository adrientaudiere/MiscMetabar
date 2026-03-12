# Session Context

## User Prompts

### Prompt 1

Bug in upset_pq(data_fungi_mini,
    fact = "Height", width_ratio = 0.2,
    taxa_fill = "Class"
  )


Error in `arrange()`:
ℹ In argument: `..1 = .`.
Caused by error:
! object '.' not found
Hide Traceback
Fix
Explain
     ▆
  1. ├─MiscMetabar::upset_pq(...)
  2. │ ├─dplyr::arrange(...)
  3. │ └─dplyr:::arrange.data.frame(...)
  4. │   └─dplyr:::arrange_rows(.data, dots = dots, locale = .locale)
  5. │     ├─dplyr::mutate(data, `:=`("{name}", !!dot), .keep = "none")
  6. │     └─dplyr:::mutat...

