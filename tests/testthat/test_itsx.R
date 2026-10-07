test_that("parse_itsx_headers() reads the id and the origin code", {
  headers <- c(
    "s1|F|ITS1 Extracted ITS1 sequence 46-191 (146 bp)",
    "s2|T|SSU Extracted SSU sequence (45 bp)",
    "id|with|pipes|A|ITS2 Extracted ITS2 sequence 292-489 (198 bp)",
    "s4",
    "s5|F fungi ITS sequence (503 bp) on main strand"
  )
  out <- parse_itsx_headers(headers)
  expect_equal(out$id, c("s1", "s2", "id|with|pipes", "s4", "s5"))
  expect_equal(out$code, c("F", "T", "A", NA, "F"))
})

test_that("find_itsx() honours the MiscMetabar.itsxpath option", {
  withr::local_options(MiscMetabar.itsxpath = "/some/where/ITSx")
  expect_equal(find_itsx(), "/some/where/ITSx")
})

# ITSx usually lives in a conda environment: set MISCMETABAR_ITSX_PRELUDE to
# the string returned by install_itsx() to run these tests.
itsx_prelude <- Sys.getenv("MISCMETABAR_ITSX_PRELUDE", unset = "")

test_that("itsx_pq() extracts the region, records the origin and filters", {
  skip_on_cran()
  skip_if_not(is_itsx_installed(args_before_itsx = itsx_prelude))
  pq <- phyloseq::prune_taxa(
    phyloseq::taxa_names(data_fungi_mini)[1:10],
    data_fungi_mini
  )
  res <- itsx_pq(
    pq,
    region = "ITS2",
    args_before_itsx = itsx_prelude,
    verbose = FALSE
  )
  expect_s4_class(res, "phyloseq")
  expect_true("ITSx_origin" %in% colnames(res@tax_table))
  expect_false(anyDuplicated(as.character(res@refseq)) > 0)
  shared <- intersect(phyloseq::taxa_names(res), phyloseq::taxa_names(pq))
  expect_true(all(
    Biostrings::width(res@refseq[shared]) <=
      Biostrings::width(pq@refseq[shared])
  ))
  origin <- as.character(res@tax_table[, "ITSx_origin"])
  expect_true(all(origin %in% c(unname(itsx_origin_codes), NA)))

  only_fungi <- itsx_pq(
    pq,
    region = "ITS2",
    remove_other_origin = TRUE,
    args_before_itsx = itsx_prelude,
    verbose = FALSE
  )
  expect_true(all(only_fungi@tax_table[, "ITSx_origin"] == "Fungi"))

  detected_only <- itsx_pq(
    pq,
    region = "ITS2",
    keep_undetected = FALSE,
    add_origin = FALSE,
    args_before_itsx = itsx_prelude,
    verbose = FALSE
  )
  expect_false("ITSx_origin" %in% colnames(detected_only@tax_table))
  expect_lte(phyloseq::ntaxa(detected_only), phyloseq::ntaxa(pq))

  flagged <- itsx_pq(
    pq,
    region = "ITS2",
    add_detected = TRUE,
    duplicated_seqs = "remove",
    args_before_itsx = itsx_prelude,
    verbose = FALSE
  )
  detected <- as.character(flagged@tax_table[, "ITSx_detected"])
  expect_true(all(detected %in% c("TRUE", "FALSE")))
  undetected <- phyloseq::taxa_names(flagged)[detected == "FALSE"]
  expect_equal(
    unname(as.character(flagged@refseq))[detected == "FALSE"],
    unname(as.character(pq@refseq))[match(undetected, phyloseq::taxa_names(pq))]
  )
  expect_equal(
    sum(detected == "TRUE"),
    phyloseq::ntaxa(detected_only)
  )
})
