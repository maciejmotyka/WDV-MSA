read_haps <- function(path) {
    files <-
        list.files(path = path,
                   full.names = T,
                   pattern = "*.fasta")
    haps <- Biostrings::readDNAStringSet(files)
    haps
}

read_haps_abqr <- function(path) {
    haps <- read_haps(path = path)
    names(haps) %<>%
        stringr::str_remove_all(" ") %>%
        stringr::str_remove("freq:") %>%
        stringr::str_remove("_bwa") %>%
        stringr::str_replace("strain", "h") %>%
        stringr::str_c("s", .) %>%
        tibble::enframe(name = NULL) %>%
        tidyr::separate(value,
                        into = c("sample", "hapl", "freq"),
                        sep = "_") %>%
        dplyr::mutate(freq = freq %>%
                          as.numeric() %>%
                          magrittr::multiply_by(100) %>%
                          round(1)) %>%
        tidyr::unite(col = ids) %>%
        dplyr::pull(ids)
    haps
}

read_haps_csnv <- function(path) {
    haps <- read_haps(path = path)
    names(haps) %<>%
        stringr::str_c("s", .) %>%
        stringr::str_remove("bwa") %>%
        stringr::str_remove("fr_") %>%
        tibble::enframe(name = NULL) %>%
        tidyr::separate(value,
                        into = c("sample", "hapl", "freq"),
                        sep = "_") %>%
        dplyr::mutate(freq = freq %>% as.numeric() %>% round(3) * 100) %>%
        dplyr::mutate(hapl = as.numeric(hapl) %>%
                          magrittr::add(1) %>%
                          as.character() %>%
                          stringr::str_c("h", .)) %>%
        tidyr::unite(col = ids) %>%
        dplyr::pull(ids)
    haps
}

align_convert <- function(inputSeqs, alignMethod, convertTo) {
    aligned_haps <-
        msa::msa(inputSeqs = inputSeqs, method = alignMethod) %>% 
        msa::msaConvert(type = convertTo)
    aligned_haps
}