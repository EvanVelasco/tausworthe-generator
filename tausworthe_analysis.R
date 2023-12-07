# load functions
source("tausworthe_functions.R")

# Generate PRNs
tausworthe("111",n = 5, r = 1, q = 3, l = 3) # code for example in theory subsection

# Generate two vectors of PRNs
t = tausworthe("0101101001001010010110101011111",
               n = 10000, r = 4, q = 31, l = 20)
t2 = tausworthe("1101111001111110011110111011000",
                n = 10000, r = 4, q = 31, l = 20)

# Place PRNs into dataframe
prn = data.frame(u1 = t, u2=t2) # store PRNs in dataframe


# graph results

# standard plot
plot = ggplot(prn, aes(x=1:length(u1),y=u1)) +
  geom_point(shape=1, size=2) +
  ggtitle("Plot of Tausworthe PRNs") +
  labs(x = "Index", y = "Value") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    title = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16), 
    legend.text = element_text(size = 12) 
  )
plot


# histogram
custom_breaks <- seq(0, 1, by = 0.1)
hist_plot = ggplot(prn, aes(x=u1)) +
  geom_histogram(col="black", fill = "grey", alpha=0.8, breaks = custom_breaks) +
  ggtitle("Histogram of Tausworthe PRNs") +
  labs(x = "Value", y = "Frequency") +
  xlim(0,1) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    title = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16), 
    legend.text = element_text(size = 12) 
  )
hist_plot


# unit square
unit_square = adjacentPrns(t, dimensions = 2)
ggplot(unit_square, aes(x=u, y=u1)) +
  geom_point(shape=1, size=3) +
  ggtitle("Plot of Adjacent Tausworthe PRNs \n on Unit Square") +
  labs(x = "Ui", y = "Ui+1") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    title = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16), 
    legend.text = element_text(size = 12) 
  )


# unit cube
unit_cube = adjacentPrns(t, dimensions = 3)
plot_ly(unit_cube, x=~u, y=~u1, z=~u2, type="scatter3d", mode="markers",
        marker = list(size = 2, opacity = 0.5, color="black")) %>% 
  layout(title = list(text = "Plot of Adjacent Tausworthe PRNs<br>on Unit Cube",
                      font = list(size = 20, color='black')),
         scene = list(
           xaxis = list(title = "Ui"),
           yaxis = list(title = "Ui+1"),
           zaxis = list(title = "Ui+2")
         ),
         margin = list(t = 60))


# statistical test for GoF
uniformChiSq(prn$u1, alpha=0.05, k=5)


# statistical tests for Independence
runsUpDown(prn$u1, alpha=0.05)
runsAboveBelow(prn$u1, alpha=0.05)
correlationTest(prn$u1, alpha=0.05)



# Normal variate generation
norm = data.frame(normals = boxMullerNorm(prn$u1, prn$u2))
ggplot(norm, aes(x=normals)) +
  geom_histogram(col="black", fill = "grey", alpha=0.8, bins=15) +
  ggtitle("Histogram of Box Muller Normals") +
  labs(x = "Value", y = "Frequency") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    title = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16), 
    legend.text = element_text(size = 12) 
  )

# test for normal distribution
set.seed(4123123)
ks.test(norm, rnorm(10000))

ggplot(norm, aes(sample=normals)) +
  stat_qq(shape=1, size=3) + 
  stat_qq_line() +
  ggtitle("Normal QQ Plot") +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    title = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16), 
    legend.text = element_text(size = 12) 
  )

