

# A wrappper to a function that solves a system of linear equations using Gauss Seidel
fitLSYS <- function(C, rhs, b, active, RSS, maxIter, tol) {
    active <- active - 1L # for the 0-based index
    ans <- .Call("fitLSYS", C, rhs, b, active, RSS, maxIter, tol)
    return(list(b = ans[[1]], RSS = ans[[2]]))
}
