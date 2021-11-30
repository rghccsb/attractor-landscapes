library(readr)
library(tibble)
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(scales)
library(factoextra)
library(ggrepel)
library(PCDimension)

df <- read_tsv('./results/objective_table.tsv')
df['Or'] <- df['Or'] / 54
df['Op'] <- 1 / df['Op']
df['Of'] <- 1 / (df['Of'] + 1)

pca <- df %>%
    select('Od', 'Oc', 'Or', 'Op', 'Of') %>%
    prcomp(center = TRUE, scale. = TRUE)

eig <- get_eig(pca)

df['Objective'] <- pca$x[,1]
df['Model ID'] <- factor(
    t(df['Model ID']),
    levels = sapply(df %>%
                  arrange(Objective) %>%
                select(`Model ID`), as.character))

shapiro_p.value <- shapiro.test(pca$x[,1])$p
print(shapiro_p.value)

proportional.variance <- pca$sdev^2 / sum(pca$sdev^2)
brokenstick <- brokenStick(1:5, 5)

df['p.value'] <- sapply(
    1:length(pca$x[,1]),
    function(i) t.test(x = pca$x[-i,1], mu = pca$x[i,1], alternative = 'greater')$p.value
)
df <- df %>%
    arrange(p.value)
sig_level <- cut(df$p.value, breaks = c(-Inf, 0.01, Inf), labels = c('*', ''))
df['pretty.p.value'] <- paste0(signif(df$p.value, digits = 2), sig_level)

x_range <- 1:length(t(df[,1]))
# Plot the objective values as a line and the p-values as annotations
objective_line <- ggplot(df, aes(x = `Model ID`, y = Objective, group = 1)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = -0.6, linetype = 'dashed', color = 'red') +
    geom_label_repel(
        aes(x = `Model ID`, label = pretty.p.value),
        force = 10,
        nudge_x = 0.1,
        nudge_y = -0.1,
        min.segment.length = 0
    ) +
    labs(subtitle = 'PC1 Objective Values and P-values') +
    ylab('PC1') +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = 'black'))


qqplot <- ggplot(df, aes(sample = Objective)) +
    stat_qq() +
    labs(subtitle = 'Q-Q Plot of PC1 Values') +
    stat_qq_line() +
    xlab('Normal Theoretical Quantiles') +
    ylab('Observed PC1 Quantiles') +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = 'black'))


## p_value_line <- ggplot(df, aes(x = `Model ID`, y = p.value, group = 1)) +
##     geom_line() +
##     geom_hline(yintercept = 0.01, linetype = 'dashed', color = 'red') +
##     scale_y_continuous(
##         trans = 'log2',
##         breaks = trans_breaks('log2', function(x) 2^x),
##         labels = trans_format('log2', math_format(2^.x))
##     ) +
##     labs(subtitle = 'T-test P-values')

pca_weights_df <- as.data.frame(pca$rotation[,1])
pca_weights_bar <- ggplot(pca_weights_df, aes(x = reorder(rownames(pca_weights_df), -`pca$rotation[, 1]`), y = `pca$rotation[, 1]`)) +
    geom_bar(stat='identity') +
    labs(subtitle = 'PCA Weights for PC1') +
    xlab('Objective') +
    ylab('PCA Weight') +
    scale_x_discrete(
        labels = c(
           'Od' = parse(text = TeX('$O_d$')),
           'Oc' = parse(text = TeX('$O_c$')),
           'Or' = parse(text = TeX('$O_r$')),
           'Op' = parse(text = TeX('$O_p$')),
           'Of' = parse(text = TeX('$O_f$'))
        )
    ) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = 'black'))


fig <- ggarrange(
    qqplot,
    objective_line,
    pca_weights_bar,
    labels = c('a.', 'b.', 'c.'),
    ncol = 1,
    nrow = 3
)
fig

ggsave('./figures/objectives.pdf', fig)
