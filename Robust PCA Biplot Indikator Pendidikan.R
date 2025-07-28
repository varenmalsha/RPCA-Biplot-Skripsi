rm(list = ls())

#===== Input data =====
data<-read.csv(file.choose(),header=T,sep=";",dec=",")
head(data)
str(data)

#===== Membuat matriks Data =====
#--- Matriks Data Y ----
Y = as.matrix(data[,2:6])
Y

#===== STANDARISASI DATA =====
y1 = (Y[,1]-mean(Y[,1]))/(sqrt(sum((Y[,1]-
                                      mean(Y[,1]))^2)/(27)))
y2 = (Y[,2]-mean(Y[,2]))/(sqrt(sum((Y[,2]-
                                      mean(Y[,2]))^2)/(27)))
y3 = (Y[,3]-mean(Y[,3]))/(sqrt(sum((Y[,3]-
                                      mean(Y[,3]))^2)/(27)))
y4 = (Y[,4]-mean(Y[,4]))/(sqrt(sum((Y[,4]-
                                      mean(Y[,4]))^2)/(27)))
y5 = (Y[,5]-mean(Y[,5]))/(sqrt(sum((Y[,5]-
                                      mean(Y[,5]))^2)/(27)))
datastandar = data.frame(x1=c(y1), x2=c(y2),
                         x3=c(y3), x4=c(y4),
                         x5=c(y5))
X = as.matrix(datastandar) 
X
round(datastandar,4)

#===== MCD =====
library(robustbase)
mcd = covMcd(datastandar)
mcd$center    # Robust mean
mcd$cov       # Robust covariance matrix

#===== Robust Mahalanobis Distance =====
rd = mahalanobis(datastandar, center = mcd$center,
                 cov = mcd$cov)
rd

#---Deteksi Outlier
cutoff = qchisq(p = 0.95 , df =
                  ncol(datastandar)) ;cutoff
datastandar[rd > cutoff ,]

#===== RSVD =====
#---NILAI EIGEN DAN VEKTOR EIGEN DARI MATRIKS KOVARIANS ROBUST
CMCD.e = eigen(mcd$cov)
CMCD.e

#---Matriks L
Eigen = CMCD.e$values  ;Eigen
AkarEigen = sqrt(Eigen)   ;AkarEigen
Ldiag = AkarEigen   ;Ldiag
L = cbind(Ldiag,Ldiag,Ldiag,Ldiag,Ldiag);L
L[1,2:5] = 0
L[2,3:5] = 0
L[3,4:5] = 0
L[3,5] = 0
L[2,1] = 0
L[3,1:2] = 0
L[4,1:3] = 0
L[5,1:4] = 0
L[4,5] = 0
L
colnames(L) = c("PC1","PC2","PC3","PC4","PC5")
L
round(L,4)

#---Matriks A
A = CMCD.e$vectors
A
round(A,4)

#---Matriks U
#Mencari Matriks 1/L
Linv = 1/L
Linv[1,2:5] = 0
Linv[2,3:5] = 0
Linv[3,4:5] = 0
Linv[3,5] = 0
Linv[2,1] = 0
Linv[3,1:2] = 0
Linv[4,1:3] = 0
Linv[5,1:4] = 0
Linv[4,5] = 0
Linv

#Hasil Matriks U
U = X %*% A %*% Linv
U
round(U,4)

#===== Pembentukan Matriks G dan Matriks H =====
#---Matriks G
alpha = 0.5
G = U%*%(L^alpha)
G
round(G,4)

#---Matriks H
Ht = L^(1-alpha)%*%t(A)
H = t(Ht)
H
round(H,4)

#===== IDENTIFIKASI PERSENTASE KERAGAMAN DATA =====
#Variasi Data
V = Eigen
V

#Proporsi total varians populasi yang bisa dijelaskan oleh PCr
PV1 = V[1]/sum(V)
PV2 = V[2]/sum(V)
PV3 = V[3]/sum(V)
PV4 = V[4]/sum(V)
PV5 = V[5]/sum(V)
data.frame(PV1,PV2,PV3,PV4,PV5)
round(data.frame(PV1,PV2,PV3,PV4,PV5),4)

#Kumulatif Variasi Data
PC1 = PV1
PC2 = PC1+PV2
PC3 = PC2+PV3
PC4 = PC3+PV4
PC5 = PC4+PV5
data.frame(PC1,PC2,PC3,PC4,PC5)
round(data.frame(PC1,PC2,PC3,PC4,PC5),4)

#===== ANALISIS BIPLOT =====
#---Matriks G2 dan Matriks H2
G2 = G[,c(1,2)]
round(G2,4)
H2 = H[,c(1,2)]
round(H2,4)

#---Visualisasi Biplot 2D
G2 = G[,c(1,2)]
H2 = H[,c(1,2)]
biplot(G2,H2, cex=1,main="PC1 vs PC2",xlab = "PC 1 = 62.75%",ylab="PC 2 = 21.66%")
abline(h=0)
abline(v=0)

#===== Identifikasi Hasil PCA Biplot =====
#---Korelasi Kosinus antar Variabel
korelasi_var = function(x) {
  library(geometry)
  p = nrow(x)
  y = matrix(,nrow=p,ncol=p)
  for (i in 1:p) {
    for (j in 1:p) {
      y[i,j] = dot(x[i,],x[j,]) / (sqrt(sum(x[i,]^2))*
                                     sqrt(sum(x[j,]^2)))
    }
  }
  print(y) }
R = korelasi_var(H2)
round(R,4)

#---Nilai Variabel pada Suatu Objek
korelasi_obvar = function(X, Y) {
  p = nrow(X)
  q = nrow(Y)
  r = matrix(, nrow = p, ncol = q)
  for (i in 1:p) {
    for (j in 1:q) {
      r[i, j] = dot(X[i, ], Y[j, ]) /
        (sqrt(sum(X[i, ]^2)) * sqrt(sum(Y[j, ]^2)))
    }
  }
  print(r) }
K = korelasi_obvar(G2, H2)
round(K,4)
