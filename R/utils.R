
# A function to create the color palette for the package
cite_colorPal <- function(n) {

    cols_10 <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
                 "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC")


    cols_20 <- c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F",
                 "#8CD17D", "#B6992D", "#F1CE63", "#499894", "#86BCB6",
                 "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", "#D37295",
                 "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6")

    if (n <= length(cols_10)) {
        return(cols_10[seq_len(n)])
    } else if (n <= length(cols_20)) {
        return(cols_20[seq_len(n)])
    } else {
        return(scales::gradient_n_pal(cols_20)(seq(0, 1, length.out = n)))
    }

}


cite_shapePal <- function(n) {
    return(c(16, 17, 15, 3, 6, 8, 1, 0, 5)[seq_len(n)])
}
