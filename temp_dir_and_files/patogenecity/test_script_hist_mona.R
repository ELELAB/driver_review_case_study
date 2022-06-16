## Test script to test geom_histogram to plot revel scores for DEG, Moonlight and anti-moonlight sets


# plot histogram using position = identity for df_all
ggplot(data = df_all, 
       mapping = aes(x = scores, 
                     color = name)) +
  geom_histogram(fill = "white", 
                 binwidth = 0.01, position = "identity", alpha = 0.5)



# subset df_all to contain only anti moonlight set
df_all_am <- df_all[which(df_all$name == "Anti_moonlight"),]

# plot histogram of anti moonligt set
ggplot(data = df_all_am, 
       mapping = aes(x = scores, 
                     color = name)) + 
  geom_histogram(binwidth = 0.01)

# subset df_all to contain only DEG set
df_all_DEG <- df_all[which(df_all$name == "DEG"),]

# plot histogram of DEG set
ggplot(data = df_all_DEG, 
       mapping = aes(x = scores, 
                     color = name)) + 
  geom_histogram(binwidth = 0.01)

# subset df_all to contain only Moonlight set
df_all_moonlight <- df_all[which(df_all$name == "Moonlight"),]

# plot histogram of Moonlight set
ggplot(data = df_all_moonlight, 
       mapping = aes(x = scores, 
                     color = name)) + 
  geom_histogram(binwidth = 0.01)
