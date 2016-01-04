#sparsify.R

###VECTORIZE MATRIX TO SPARSE MATRIX/DATAFRAME
#eg turn a matrix and row and column indices into a sparse representation
sparsify <- function(m, cols, rows)
{
  #as vector goes column-by-column
  data.frame( i = rep(cols,times=length(rows)),
              j = rep(rows,each=length(cols)),
              x = as.vector(m))
}
