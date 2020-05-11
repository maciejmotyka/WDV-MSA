plan <- drake::drake_plan(
    
    # Read haplotype seqs and clean their names ####
    haps = target(
        read_function(path = drake::file_in(path)),
        transform = map(
            path = c(
                "/home/lejno/Desktop/aBayesQR-nf/out_bwa_001/haplotypes",
                "/home/lejno/Desktop/cliqueSNV-nf/out_t5_tf001/haplotypes",
                "/home/lejno/Desktop/cliqueSNV-nf/out_t5_tf01/haplotypes",
                "/home/lejno/Desktop/cliqueSNV-nf/out_t10_tf01/haplotypes"
            ),
            read_function = !!rlang::syms(
                c(
                    "read_haps_abqr",
                    "read_haps_csnv",
                    "read_haps_csnv",
                    "read_haps_csnv"
                )
            ),
            .names = c("abqr_001",
                       "csnv_t5_tf001",
                       "csnv_t5_tf01",
                       "csnv_t10_tf01")
        )
    ), 
    
    # Align the haps using Muscle, ClustalO and ClustalW with default settings ####
    aln = target(
        msa::msa(inputSeqs = haps, method = methods),
        transform = cross(
            haps,
            methods = c("Muscle", "ClustalOmega", "ClustalW"),
            .names = outer(
                c("Muscle", "ClustO", "ClustW"),
                c("abqr_001",
                  "csnv_t5_tf001",
                  "csnv_t5_tf01",
                  "csnv_t10_tf01"),
                paste,
                sep = "_"
            ) %>% as.vector()
        )
    ),
    
    # Calculate distance matrices ####
    mlDist = target(
        msa::msaConvert(aln, type = "phangorn::phyDat") %>%
            phangorn::dist.ml(),
        transform = map(aln)
    ),
    # FIXME: returns a matrix of NaNs; ask at biostars 
    dnaDist = target(msa::msaConvert(aln, type = "ape::DNAbin") %>%
                         ape::dist.dna(),
                     transform = map(aln)
    ),

    # Construct NJ trees ####
    nj_mlDist = target(
        phangorn::NJ(mlDist) %>%
            ape::ladderize(),
        transform = map(
            mlDist,
            .names = outer(
                c("Muscle", "ClustO", "ClustW"),
                c("abqr_001",
                  "csnv_t5_tf001",
                  "csnv_t5_tf01",
                  "csnv_t10_tf01"),
                paste,
                sep = "_"
            ) %>%
                paste("nj_mlDist", ., sep = "_") %>%
                as.vector()
        )
    ),
    
    # FIXME: fix msaConvert(type = "ape::DNAbin") 
    # nj_dnaDist = target(phangorn::NJ(dnaDist) %>% ape::ladderize(),
    #                     transform = map(dnaDist,
    #                                     .id = c(ids2, ids)))
    # 
    
    
    # Make UPGMA trees ####
    # abqr001_Muscle_dml_upgma = abqr001_Muscle_dml %>% phangorn::upgma() %>% ape::ladderize(),
    # abqr001_ClustalO_dml_upgma = abqr001_ClustalO_dml %>% phangorn::upgma() %>% ape::ladderize(),
    # abqr001_ClustalW_dml_upgma = abqr001_ClustalW_dml %>% phangorn::upgma() %>% ape::ladderize(),
    # abqr001_Muscle_dml_nj = abqr001_Muscle_dml %>% phangorn::NJ() %>% ape::ladderize(),
    # abqr001_ClustalO_dml_nj = abqr001_ClustalO_dml %>% phangorn::NJ() %>% ape::ladderize(),
    # abqr001_ClustalW_dml_nj = abqr001_ClustalW_dml %>% phangorn::NJ() %>% ape::ladderize(),
    # 
    # # Compare the trees
    # trees = c(
    #     abqr001_Muscle_dml_upgma,
    #     abqr001_ClustalO_dml_upgma,
    #     abqr001_ClustalW_dml_upgma),
    
    # trees = {
    #   tree <- ape::rmtree(1, 10)
    #   c(tree, tree, tree)
    # },
    
    # abqr_001_trspc = treespace::treespace(trees, nf = 2),
    
    # abqr_001_trspc = {
    # trees <- list(
    #   abqr001_Muscle_dml_upgma,
    #   abqr001_ClustalO_dml_upgma,
    #   abqr001_ClustalW_dml_upgma
    # abqr001_Muscle_dml_nj,
    # abqr001_ClustalO_dml_nj,
    # abqr001_ClustalW_dml_nj
    # )
    # class(trees) <- "multiPhylo"
    #   abqr_001_trspc <- treespace::treespace(trees)
    #   abqr_001_trspc
    # },
    
    # Make report ####
    # report = rmarkdown::render(
    #   input = drake::knitr_in("report.Rmd"),
    #   output_file = drake::file_out("report.html"),
    #   quiet = T)
)