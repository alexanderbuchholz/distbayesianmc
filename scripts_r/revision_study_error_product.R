# study error normal approximation

vec_datasets <- c("hla_ultra_small1")#, "hla_ultra_small2")
#vec_splits <- c(1,2,4,5,10,20,50)
#vec_splits <- c(2,4,8,15,25,30,35,40,45)
vec_splits <- c(1,2,4,5,8,10,15,20,25,30,35,40,45,50)

df <- f_combine_const_data_in_frame(vec_splits, vec_datasets, vec_types_splits, 6)
dfinter <- df %>% filter(dataset == "hla_ultra_small1")

true_norm_const <-  mean(pull(dfinter %>% filter(splits == "1") %>% select(normconstcombined)))

integral_true_post <- true_norm_const-as.numeric(pull(dfinter %>% filter(splits != "1") %>% select(logpost))) - as.numeric(pull(dfinter %>% filter(splits != "1") %>% select(alphasub)))
integral_normal_approx <- as.numeric(pull(dfinter %>% filter(splits != "1") %>% select(logsubposteriorapprox)))
error = integral_true_post-integral_normal_approx
splits_covar <- as.numeric(as.character(pull(dfinter %>% filter(splits != "1") %>% select(splits))))


# plot just the error
reg <- lm(error ~ splits_covar)
summary(reg)
reg$coefficients
f_interpol <- function(s){
  reg$coefficients[2]*s
}
plot(splits_covar, error)

y = sapply(1:50, f_interpol)
x = 1:50
lines(x,y)

df_data <- data.frame(error, splits_covar, logerror)
plot_log_diff <- ggplot(df_data, aes(splits_covar, error)) + geom_point(alpha=0.7, size=5) + labs(x = "# Splits of data", y = "Error", title="Relative log error") + scale_size_continuous(range = c(3,10)) +
  scale_colour_Publication()+ theme_Publication() + geom_smooth(method='lm')

plot_diff_log <- ggplot(df_data, aes(splits_covar, logerror)) + geom_point(alpha=0.7, size=5) + labs(x = "# Splits of data", y = "log of error", title="Log error difference") + scale_size_continuous(range = c(3,10)) +
  scale_colour_Publication()+ theme_Publication() + 
  geom_smooth(method='lm') 

plot_diff_log_log <- ggplot(df_data, aes(splits_covar, logerror)) + geom_point(alpha=0.7, size=5) + labs(x = "log # Splits of data", y = "log of error", title="Log error difference") + scale_size_continuous(range = c(3,10)) +
  scale_colour_Publication()+ theme_Publication() + 
  geom_smooth(method='lm') + scale_x_log10()

library(ggpubr)
figure_error <- ggarrange(plot_log_diff, plot_diff_log,plot_diff_log_log, nrow=1)
ggsave("error_approximation.png", plot = figure_error, width = 12, height = 6)

library(matrixStats)
logSumExp
matdiff <- cbind(-integral_true_post, integral_normal_approx)
logerror <- apply(matdiff, 1, logSumExp)
logsplits <- log(splits_covar)
moderror <- logerror-3/2*logsplits


reg <- lm(logerror ~ splits_covar+logsplits)
summary(reg)
reg$coefficients
f_interpol <- function(s){
  reg$coefficients[1]+reg$coefficients[2]*s+reg$coefficients[3]*log(s)
}
plot(splits_covar, moderror)

y = sapply(1:50, f_interpol)
x = 1:50
plot(log(splits_covar), logerror)
lines(x,y)
