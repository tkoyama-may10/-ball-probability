prob <- read.table("tmp", colClasses=c("character","character"))
time <- read.table("tmp2",colClasses=rep("character",4))

t1 <- data.frame(dim=seq(5,100,5),prob=prob$V2[1:20], time=time$V4[1:20])

sink("tab3.tex")

cat("\\begin{table}[htbp]\n")
cat("\\begin{center}\n")
cat("\\begin{tabular}{ccc}\n")
cat("dim& $p$ & time(s)\\\\\n")
cat("\\hline\n")

 for(i in c(1:20)){cat(sprintf("%d & %s & %s \\\\",t1[i,1],t1[i,2],t1[i,3]), "\n")}

cat("\\hline\n")
cat("\\end{tabular}\n")
cat("\\end{center}\n")
cat("\\caption{Computational times for Anderson-Darling statistic}\n")
cat("\\label{tab:AD}\n")
cat("\\end{table}\n")


sink()

time <- read.table("tmp2")

postscript(file="./tab3.eps",width=5,height=5,horizontal=FALSE,onefile=FALSE,paper="special",family="Helvetica")
plot(seq(5,100,5), time$V4[1:20], xlab="Dimension", ylab="Computational time")
dev.off()
