######################################################################

# TITLE: Descriptive Summary Statistics
# DATE MODIFIED: 03/21/2022
# AUTHOR: EDB

######################################################################


##### Eigenvalues and Eigenvectors #####
#   Eigenvalues and eigvectors sastify Ae(k) = Lam(k)e(k), k = 1...p
#   Eigenvectors have unit length and are pairwise orthogonal
#   The spectral decomposition of A can be written as A = Lam(1)e(1)e'(1)+...
#   The elements of the eigenvectors e(k) are coefficients of linear combinations of random
#       variables X(1), X(2), ..., X(p) that will allow use to reduce the dimensionality
#       p-dimensional data.

# define the matrix A
A <- matrix(c(9, -2, -2, 6),nrow=2,ncol=2,byrow=T)

# determine the eigenvalues and eigenvectors of a
eigA <- eigen(A)

# write out the spectral decomposition of A
P = eigA$vectors #define the eigenvectors
t(P)  # evaluate its transpose
solve(P) # evaluate its inverse

options(digits=2)
t(P)%*%P
Lam = diag(eigA$values) # create a matrix with diagonal elements equal to eigenvalues of A
P%*%Lam%*%t(P)  # P(Lam)P' = A?
A

##### Inverse of a Matrix #####

# find A-1
Ainv = solve(A)
eigAinv = eigen(Ainv)

# write out the spectral decomposition of A -1
Pinv = eigAinv$vectors 
t(Pinv)
solve(Pinv)
options(digits=2)
t(Pinv)%*%Pinv
Laminv = diag(eigAinv$values) # create a matrix with diagonal elements equal to eigenvalues of A
Pinv%*%Laminv%*%t(Pinv)  # P(Lam)P' = A?
Ainv

##### Square Root Matrix (& its inverse) #####

# find the matrices A 1/2 and A -1/2
eigA = eigen(A)
P = eigA$vectors
Lam = diag(eigA$values)
srLam = diag(sqrt(eigA$values))
srLaminv = diag(1/sqrt(eigA$values))

sqrtA = P%*%srLam%*%t(P)
sqrtA

# check square root matrix condition
sqrtA%*%sqrtA

sqrtAinv = P%*%srLaminv%*%t(P)
sqrtAinv

sqrtAinv%*%sqrtAinv # Does this equal the inverse of A?
solve(A)

sqrtA%*%sqrtAinv # identity condition check for the matrix & its inverse

##### Transpose of a Matrix #####

A =  matrix(c(4, 8, 8, 
              3, 6, 9),nrow=2,ncol=3,byrow=T)
A

# Calculate AA'
B = t(A) # calculate A'
AB = A%*%B # multiply by A
AB

# obtain its eigenvalues and eigenvectors
eigAB = eigen(AB)
eigAB
P = eigAB$vectors
P
Lam = eigAB$values
Lam

# Calculate A'A
BA = B%*%A
BA

# obtain its eigenvalues and eigenvectors
eigBA = eigen(BA)
eigBA
P = eigBA$vectors
P
Lam = eigBA$values
Lam

##### Singular Value Decomposition #####
#   Singular value decomposition will decompose an n×p matrix X into a
#   product of three matrices: X=UDV′ where,
#   U=n×n matrix containing the eigenvectors of XX′
#   D=n×pmatrix containing the singular values in a p×p submatrix. 
#     The singular values are the non-zero eigenvalues of XX′
#   V=p×p matrix containing the eigenvectors of X′X
#
#   Typically the first few singular values are large relative to the rest. 

Asvd = svd(A)

Asvd$d
Asvd$u
Asvd$v

U = Asvd$u
D = diag(Asvd$d)
V = Asvd$v

# Check that UDV' = X
U%*%D%*%t(V) 

##### Olive Oil Data #####
# There are two geographic classifications in this data: the first is nine
#   individual growing areas in Ital.  A broader classification is the growing 
#   region in Italy.

# Use visualization methods to identify the fatty acits that would be most
#   useful in discriminating between olive oils grown in the nine different growing
#   areas represented.



