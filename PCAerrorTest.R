
library(magrittr)

a = c(0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, -4, -1, 2)
b = c(1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, -1, -5, 3)
c = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  1,  0,  0)
d = c()

df = data.frame(a, b, c)

prin_comp <- prcomp(df, scale = F)
prin_comp2 <- prcomp(df, scale = T)

sapply(1:ncol(t(df)), function(x){
  length = unique(t(df)[, x]) %>% length
}) %>% table



a = c(1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, -1, -1, 1)
b = c(1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, -1, -1, 2)
c = c(1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  1,  1, 1)
d = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  1,  1, 1)

df = data.frame(a, b, c,d)
df
prin_comp2 <- prcomp(t(df), center = TRUE, scale = F)

#get rownames for constant variables:
print("These variables get constant values, thus got droped during PCA scale:")
row.names(df[apply(df, MARGIN = 1, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE)),])


#exclude columns with constant value:
newdf <- df[,!apply(df, MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))]
newdf <- df[!apply(df, MARGIN = 1, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE)),]

colnames(df[!(colnames(df) %in% colnames(newdf))])
constantRow <- row.names(df[!(row.names(df) %in% row.names(newdf)),])

#exclude columns with constant value
df[,sapply(df, function(v) var(v, na.rm=TRUE)!=0)]

colnames(df[, sapply(df, function(v) var(v, na.rm=TRUE)==0)])
