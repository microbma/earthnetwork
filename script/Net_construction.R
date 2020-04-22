# FDR

p.fdr <- p.adjust(cor.p, method = "BH")

par(mfcol=c(1,2))
hist(c(cor.p),breaks = 10000,xlim=c(0, .01))
hist(c(p.fdr),breaks = 10000,xlim=c(0, .01))

# Network construction
ot <- data.frame(soil.matrix.rare.present)

## correlation matrix

## KL divergence 

klperm <- function(x=ot$otu1, y=ot$otu3, N=1000){
        reps <- replicate(N, KLdiv(sample(x), y))
        obs <- KLdiv(x,y)
        p <- mean(reps > c(obs)) # shortcut for sum(reps > obs)/N
        p
}

kl.res <- matrix(NA, nrow=NCOL(ot), ncol=NCOL(ot))

for(iii in 1:NCOL(ot)) {
        for(jjj in 1:NCOL(ot)) {
                kl.res[iii,jjj] <- klperm(ot[,iii], ot[,jjj])
                print(iii)
        }
}

rownames(kl.res)<-names(ot)
colnames(kl.res) <- names(ot)

## Merge P values


## Generate network



## Silencing method

G <- abs(mean)

I <- diag(2364)
D <- (G-I)%*%G
D1 <- D
D1[upper.tri(D)] <- 0
D1[lower.tri(D)] <- 0

J <- (G - I + D1) %*% solve(G)

diag(G) <- 0
diag(J) <- 0

plot(c(G[,1:10]),J[,1:10])
