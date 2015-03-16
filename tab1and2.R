prob <- read.table("tmp", colClasses=c("character","character"))
time <- read.table("tmp2",colClasses=rep("character",4))

 t1 <- data.frame(dim=c(10:20),prob1=prob$V2[1:11], time1=time$V4[1:11],prob2=prob$V2[12:22], time2=time$V4[12:22])
 t2 <- data.frame(dim=c(10:20),prob1=prob$V2[23:33], time1=time$V4[23:33],prob2=prob$V2[34:44], time2=time$V4[34:44])

sink("tab1and2.tex")

cat("\\begin{table}[htbp]\n")
cat("\\begin{center}\n")
cat("\\begin{tabular}{c|cccc}\n")
cat("\\hline\n")
cat("\\multirow{2}{*}{dim}&\n")
cat("\\multicolumn{2}{c}{$\\mu=\\mu^{(1)}=0$}&\n")
cat("\\multicolumn{2}{c}{$\\mu=\\mu^{(2)}\\neq 0$}\\\\\n")
cat("& $p$ & time(s) & $p$ & time(s)  \\\\\n")
cat("\\hline\n")
 for(i in c(1:11)){cat(sprintf("%d & %s & %s & %s & %s \\\\",t1[i,1],t1[i,2],t1[i,3],t1[i,4],t1[i,5]), "\n")}
cat("\\hline\n")
cat("\\end{tabular}\n")
cat("\\end{center}\n")
cat("\\caption{Probabilities  and times for $\\Sigma^{(1)}$}\n")
cat("\\label{tab:hirotsu-1}\n")
cat("\\end{table}\n")

cat("\n")

cat("\\begin{table}[htbp]\n")
cat("\\begin{center}\n")
cat("\\begin{tabular}{c|cccc}\n")
cat("\\hline\n")
cat("\\multirow{2}{*}{dim}&\n")
cat("\\multicolumn{2}{c}{$\\mu=\\mu^{(1)}=0$}&\n")
cat("\\multicolumn{2}{c}{$\\mu=\\mu^{(2)}\\neq 0$}\\\\\n")
cat("& $p$ & time(s) & $p$ & time(s)  \\\\\n")
cat("\\hline\n")


 for(i in c(1:11)){cat(sprintf("%d & %s & %s & %s & %s \\\\",t2[i,1],t2[i,2],t2[i,3],t2[i,4],t2[i,5]), "\n")}

cat("\\hline\n")
cat("\\end{tabular}\n")
cat("\\end{center}\n")
cat("\\caption{Probabilities  and times for $\\Sigma^{(2)}$}\n")
cat("\\label{tab:hirotsu-2}\n")
cat("\\end{table}\n")

sink()