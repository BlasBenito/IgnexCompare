#' Solver
#'
#' Search for optimal slotting.
#'
#' @usage
#'
#' @param
#' @param
#' @return
#' @author Blas Benito <blasbenito@gmail.com>
#' @examples
#' @export
Solver=function(distance.matrix){

  #example distance matrix
  distance.matrix=manhattan.matrix
  dm=distance.matrix

  #search area
  a=1
  b=1
  a.max=dim(distance.matrix)[1]+1
  b.max=dim(distance.matrix)[2]+1

  #matrix to store I don't know what yet
  f=matrix(nrow=100, ncol=100)
  f[0:100, 0:100]=0

  #I don't know what it does yet...
  #loop 1
  for (i in a:a.max){

    i1 = i + 1

    f[i1, 1] = f[i, 1] + 2.0 * dm[i, 1]

  }#end of loop 1

  #loop2
  for(k in b:b.max){

    k1 = k + 1

    f[1, k1] = f[1, k] + 2.0 * dm[1, k]

    for(j in a:a.max){

      j1 = j + 1

      x = min((f[j, k1] + dm[j, k1]), (f[j1, k] + dm[j1, k]))

      f[j1, k1] = x + dm[j, k]
    }
  }
  #end of loop2




}#end of function


#
# //C      DIMENSION D(110,110)
#
# //C     COMMON F(110,110),NS(2,210)
#
# //C New
#
# //DIMENSION D(310, 310)
#
# //COMMON F(310, 310), NS(2, 410)
#
# //C
#
# int A, AMAX, B, BMAX;
#
# F(1, 1) = 0.0;
#
# for(5 J = 1, AMAX)
# {
#
#   J1 = J + 1;
#
#   F(J1, 1) = F(J, 1) + 2.0*D[J, 1];
# }
#
# for(6 K = 1, BMAX)
# {
#
#   K1 = K + 1;
#
#   F(1, K1) = F(1, K) + 2.0*D[1, K];
#
#   for(6 J = 1, AMAX)
#   {
#
#     J1 = J + 1;
#
#     X = AMIN1((F(J, K1) + D[J, K1]), (F(J1, K) + D[J1, K]));
#
#     F(J1, K1) = X + D[J, K];
#   }
# }
#
# J = AMAX;
#
# K = BMAX;
#
# L = AMAX + BMAX;
#
# M = L;
#
# for(10 I = 1, L)
# {
#
#   if (K  <=  0.) goto g12;
#
#   g13:
#     if (J  <=  0.) goto g11;
#
#   g14:
#     J1 = J + 1;
#
#     K1 = K + 1;
#
#     if (ABS(F(J1, K) + D[J, K] + D[J1, K] - F(J1, K1))  >=  1E-05) goto g12;
#
#     //c    The following original 'IF statement' does not produce desired
#
#     //c           results in profort:
#
#       //c      IF (F(J1,K)+D(J,K)+D(J1,K)-F(J1,K1)) 12,11,12
#
#     g11:
#       NS(1, M) = B;
#
#     NS(2, M) = K;
#
#     K = K - 1;
#
#     goto g10;
#
#     g12:
#       NS(1, M) = A;
#
#     NS(2, M) = J;
#
#     J = J - 1;
#
#     M = M - 1;
# }
#
# return;
#
# };
