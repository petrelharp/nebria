
transl <- function (x) {
    tt <- c(
        "Army" = "Army Pass",
        "Conness" = "Conness Lake",
        "Crabtree" = "Crabtree Lakes",
        "Donohue" = "Donohue Pass",
        "Dusy" = "Dusy Pass",
        "FRecess" = "Recess Lakes",
        "HungryP" = "Hungry Packer Lake",
        "Hungry Pass" = "Hungry Packer Lake",
        "Italy" = "Italy Lake",
        "Kuna" = "Kuna Lake",
        "Lamarck" = "Lamarck Lakes",
        "Lyell" = "Lyell Peak",
        "Millys" = "Milly's Footpass",
        "Monarch" = "Monarch Lake",
        "NForester" = "North Forester",
        "Ottoway" = "Ottoway Lake",
        "Pear" = "Pear Lake",
        "Piute" = "Piute Pass",
        "Ritter" = "Ritter Range",
        "Ruby" = "Ruby Lake",
        "SamMack" = "Sam Mack Lake",
        "Selden" = "Selden Pass",
        "SForester" = "South Forester",
        "Sixty" = "Sixty Lakes",
        "Sphinx" = "Sphinx Lakes",
        "Sphinx  Lakes" = "Sphinx Lakes",
        "Taboose" = "Taboose Pass",
        "Treasure" = "Treasure Lake",
        "Wright" = "Wright Lakes"
    )
    repl <- (x %in% names(tt))
    x[repl] <- tt[match(x[repl], names(tt))]
    return(x)
}

make_names = function (a, b) {
    paste(
          ifelse(a < b, a, b),
          ifelse(a < b, b, a),
          sep="_"
    )
}

