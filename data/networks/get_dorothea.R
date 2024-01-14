

library(dorothea)


net <- decoupleR::get_dorothea(levels = c('A', 'B', 'C', 'D'))
head(net)
dim(net)

file_ab = "dorothea_AB.txt"
df = as.data.frame(net)
head(df)
df2 = df[df$confidence %in% c("A","B"),]
dim(df2) # 15113     4
df2 = df2[,c("source","target")]
colnames(df2) = c("node1","node2")
write.table(df2, file_ab, sep="\t", row.names = FALSE, quote=FALSE)