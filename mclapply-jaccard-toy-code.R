library(jaccard)
library(CellTagR)
library(parallel)

# presumes there is a celltag object with metric.filtered.count slot\
# grabs cells with celltags that passed filtering
test_mtx <- celltag@metric.filtered.count[complete.cases(as.matrix(celltag@metric.filtered.count)),]

# create function to pass into mclapply
# jaccard.test supports multiple types of p value generations
# jaccard.test.bootstrap/jaccard.test.mca seem to be the most stable and appropriate (not an expert though)
my_test <- function(pair, data){
  jaccard.test.mca(data[pair[1],],
                   data[pair[2],])
}

# function to run jaccard test between all rows in a celltag@metric.filtered.count
jaccard_test_pairwise <- function(dat){
  # returns all unique permutations, will result in calculating upper triangle only
  my_pair <- combn(1:nrow(dat), 2, simplify = FALSE)
  # for all pairs, run jaccard.test.mca
  out <- mclapply(my_pair, my_test, data = dat)
  # create long format table with statistics (jaccard coef) and p value (pvalue)
  table <- data.frame(i = sapply(my_pair, '[[', 2),
                      j = sapply(my_pair, '[[', 1),
                      statistic = sapply(out, '[[', "statistics"),
                      p.value = sapply(out, '[[', "pvalue"))
  return(table)  
}

# mclapply should automatically detect the number of cores available and use them all, watch out for memory issues
# you can run into cases where you select too many cores and not enough RAM, and you won't get full usage
test_run <- jaccard_test_pairwise(test_mtx)

# check that simil() returns roughly the same values as jaccard.test.mca()
test_join <- summary(celltag@jaccard.mtx) %>%
  left_join(test_run, by = c("i", "j"))
# perform p value adjustment
test_stats <- test_join %>%
  filter(!is.na(p.value)) %>%
  mutate(FDR = p.adjust(p.value))