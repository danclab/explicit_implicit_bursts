
library(lme4)
library(emmeans)

df_rt <- read.csv('./tmp_data.csv')
df_rt$block_label <- factor(df_rt$block_label,
                             levels = c('Baseline', 'First Adaptation Block', 'Last Adaptation Block', 'Washout'))

model <- glmer(reach_rt ~ block_label * group + (1 | subject), 
               data = df_rt, 
               family = Gamma(link = "log"),
               control = glmerControl(optimizer = "bobyqa"))

print(summary(model))

emm_block <- emmeans(model, ~ block_label | group)
block_vs_baseline <- contrast(emm_block, method = "trt.vs.ctrl", ref = "Baseline", adjust = "holm")
summary(block_vs_baseline)

emm_group <- emmeans(model, ~ group | block_label)
group_diff <- contrast(emm_group, method = "revpairwise", adjust = "holm")
summary(group_diff)

write.csv(as.data.frame(block_vs_baseline), './posthoc_block_vs_baseline.csv', row.names = FALSE)
write.csv(as.data.frame(group_diff), './posthoc_group_diff.csv', row.names = FALSE)
