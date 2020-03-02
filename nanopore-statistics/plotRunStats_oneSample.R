args = commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
print(args)
input_fn <- args[[1]]
output_fn <- args[[2]]

pdf(output_fn)

histD <- read.delim(input_fn, header = F)

readLengths <- rep(histD[[1]], histD[[2]])

print(readLengths)

np_megabases <- sum(as.numeric(readLengths))
rL_min <- min(readLengths)
rL_median <- median(readLengths)
rL_max <- max(readLengths)
readLengths <- sort(readLengths)
rL_N50 <- -1

sum_readLengths <- 0
for(i in 1:length(readLengths))
{
	sum_readLengths <- sum_readLengths + readLengths[[i]]
	if(sum_readLengths >= 0.5 * np_megabases)
	{
		rL_N50 <- readLengths[[i]]
		break
	}
}
H <- hist(readLengths, nclass = 100, plot = F)

truncateX <- 60000
np_gigabases <- (np_megabases / 1024**3)
if(np_gigabases < 1)
{
	plot(H,  xlab = "Read length", main = paste(args[[2]], "\n", sprintf("%.1f",(np_megabases / 1024**2)), "Mb data, median read length ", sprintf("%.1fkb", rL_median/1000), ", ", sep = ""), xlim = c(0, truncateX))
} else
{
	plot(H,  xlab = "Read length", main = paste(args[[2]], "\n", sprintf("%.1f",(np_megabases / 1024**3)), "Gb data, median read length ", sprintf("%.1fkb", rL_median/1000), ", ", sep = ""), xlim = c(0, truncateX))
}
y_span <- max(H$counts) - min(H$counts)
y_oneStep <- y_span / 20

text(truncateX, max(H$counts) - 0 * y_oneStep, paste("Min. read: ", rL_min, sep = ""), adj = 1, cex = 0.8)
text(truncateX, max(H$counts) - 1 * y_oneStep, paste("Max. read: ", rL_max, sep = ""), adj = 1, cex = 0.8)
text(truncateX, max(H$counts) - 2 * y_oneStep, paste("Read N50: ", sprintf("%.1f", rL_N50/1000), "kb", sep = ""), adj = 1, cex = 0.8)
text(truncateX, max(H$counts) - 3 * y_oneStep, paste("Total reads: ", length(readLengths), sep = ""), adj = 1, cex = 0.8)

dev.off()

