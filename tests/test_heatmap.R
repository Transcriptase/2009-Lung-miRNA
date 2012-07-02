context("Heatmap Prep Functions")
test_that("Name Character Sub Works",
          expect_that(name_prep("mmu-miR-33*"),
                      equals("mmu.miR.33.")))
t_samp <- c(3, 4, 5, 9, 10, 11, 12, 13)
mir_name <- "mmu-miR-32"
test_that("miRNA Lookup",{
          expect_that(length(miRNA_lookup(mir_name, t_samp)),
                       equals(length(t_samp))
                      )
          expect_that(miRNA_lookup(mir_name, t_samp),
                      is_a("numeric")
                      )
          })

mir_names <-c("mmu-miR-9",
          "mmu-miR-15b",
          "mmu-miR-29a",
          "mmu-miR-30a",
          "mmu-miR-103")
a <- make_table(mir_names, t_samp)

test_that("Table creation", {
    expect_is(a,"matrix")
    expect_equal(length(a[,1]), length(mir_names))
    expect_equal(length(a[1,]), length(t_samp))
    expect_equal(rownames(a), mir_names)
})