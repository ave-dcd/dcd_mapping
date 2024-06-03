library(dplyr)
library(ggplot2)
library(magrittr)

mave <- read.csv('mave_concordance.csv')
ex <- c()
for (i in 1:nrow(mave)) {
  temp = strsplit(mave$Scoreset[i], '')[[1]]
  vec = c()
  for (j in 1:length(temp)) {
     if (temp[j] == '-') {
      vec = c(vec, j)
    }
  }
  ex <- c(ex, substr(mave$Scoreset[i], 1, vec[2] -1))
}
mave$experiment <- ex

mave_assay <- data.frame(matrix(ncol = 2, nrow = 0))
x <- c('experiment', 'context')
colnames(mave_assay) <- x
for (i in 1:nrow(mave)) {
  if (mave$experiment[i] %in% mave_assay$experiment) {
    next
  }
  else {
    vec = c(mave$experiment[i], mave$Cellular.Context[i])
    mave_assay[i,] <- vec
  }
}
mave_assay <- mave_assay %>%
  na.omit(experiment)

context_counts <- c(table(mave_assay$context))
names = names(table(mave_assay$context))
df <- data.frame('Experiment Cellular Context' = names, 'value' = context_counts)

ggplot(df, aes(x = factor(Experiment.Cellular.Context, levels = c('Human', 'Yeast', 'Bacteria', 'Mouse', 'Bacteriophage', 'N/A')), y = value, fill = rownames(df))) +
  geom_bar(stat = 'identity', fill = c("#F8766D","#B79F00","#90ee90","#00BFC4","#619CFF","#F564E3")) +
  geom_text(aes(label = value), vjust = ifelse(df$value != 86, -1, 3), size = 10, colour = ifelse(df$value == 97, 'white', 'black')) + 
  xlab('MAVE Experiment Cellular Context') +
  ylab('Number of Experiments') +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(legend.position = 'none') +
  theme(axis.text=element_text(size=12, face = 'bold'),axis.title=element_text(size=14,face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))
