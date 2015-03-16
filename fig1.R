fig1 <- function(){
     x <- read.table("fig1.txt")
     postscript(file="./fig1.eps",width=5,height=5,horizontal=FALSE,onefile=FALSE,paper="special",family="Helvetica")
     plot(x$V1, x$V2, xlab="r", ylab="Ball Probability")
     dev.off()
}

fig1()