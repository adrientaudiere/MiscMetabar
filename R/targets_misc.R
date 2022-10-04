
# MM_idtaxa(seqtab, train_maarjam, type = "collapsed")
MM_idtaxa <- function(test, trainingSet, ...) {
    idtaxa <- DECIPHER::IdTaxa(test = test, trainingSet = trainingSet, ...)
    col2add <- 8 - lengths(regmatches(idtaxa, gregexpr(";", idtaxa)))
    for (i in seq_along(idtaxa)) {
        idtaxa[i] <-
            paste(
                idtaxa[i],
                paste(as.character(rep(";", each = col2add[i])),
                    collapse = ""
                )
            )
    }
    t_idtaxa <- tibble::tibble(data.frame(
        stringr::str_split_fixed(idtaxa, ";", 9)
    ))
    colnames(t_idtaxa) <- c(
        "Kingdom_idtaxa",
        "Phyla_idtaxa",
        "Class_idtaxa",
        "Order_idtaxa",
        "Family_idtaxa",
        "Genus_idtaxa",
        "Species_idtaxa",
        "VT_idtaxa"
    )
    return(t_idtaxa)
}
