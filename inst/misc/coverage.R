## code to check which methods are available

library(broom.mixed)
library(Matrix)

## https://stackoverflow.com/questions/51729207/how-to-list-all-s3-methods-defined-in-a-specific-package-namespace-for-a-parti

findFun <- function (Fun, pkg) {
  all_fun <- ls(getNamespace(pkg))
  all_fun[startsWith(all_fun, sprintf("%s.", Fun))]
}

ff <- function(m) {
    gsub(sprintf("^%s\\.",m),"",findFun(m,"broom.mixed"))
}

mvec <- c("tidy","glance","augment")
mm <- setNames(lapply(mvec,ff),mvec)
all_m <- sort(Reduce(union,mm))

res <- Matrix::Matrix(data=0,ncol=3,nrow=length(all_m),
                      dimnames=list(all_m,mvec))

for (i in mvec) {
    for (j in all_m) {
        if (j %in% mm[[i]]) res[j,i] <- 1
    }
}
print(res)
image(res)
par(mar=c(20,1,1,4))
heatmap(as.matrix(res),Rowv=NA,Colv=NA,scale="none",
        margins=c(10,5),
        col=c("white","gray"))
