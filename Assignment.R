# Install the required package
install.packages("readr")
install.packages("styler")
install.packages("iNEXT")
# Load the installed Package
#Install the library
library(readr)
library(tidyverse)
library(ggplot2)
library(vegan)
library(lubridate)
library(styler)
##Reading in data i.e Brachiopoda location from the PC

df_bold <- read_tsv(file = "../R/result.tsv")

#Exploring data with various functions

# Check the class of my data frame (df_bold)
class(df_bold)

# Display the column names in my data frame above
names(df_bold)

# Show the dimensions
dim(df_bold)
 #take note of a shorter view of some important variables
df_bold.subset <- df_bold[, c("processid", "bin_uri", "order", "species", "country/ocean", "realm", "coord")]
df_bold.subset

#Understanding the dataset with some functions
# Count the number of unique countries in the country/ocean column
number_countries <- length(unique(df_bold$`country/ocean`))
number_countries

#Display the list of unique countries found in 'country/ocean'
unique(df_bold$`country/ocean`)
#Shows the total number
length(unique(df_bold$`country/ocean`))
unique(df_bold$order)
length(df_bold$genus)
unique(df_bold$bin_uri)
length(df_bold$bin_uri)

# Now, I am trying to create a frequency table of the 'order' variable
table(df_bold$order)

# Also create a frequency table of the 'country/ocean' variable
table(df_bold$`country/ocean`)

# Sort the frequency table of 'country/ocean' in descending order
sort(table(df_bold$`country/ocean`), decreasing = TRUE)

# Plot the top 20 most frequent 'country/ocean' values
plot(sort(table(df_bold$`country/ocean`), decreasing = TRUE)[1:20])

# Sort the frequency table of 'bin_uri' in descending order
sort(table(df_bold$bin_uri), decreasing = TRUE)

#Count occurrences of each 'country/ocean' value, sorted by frequency
df_bold %>%
  count(`country/ocean`, sort = TRUE)
#so here, i removed the NAs
df_bold %>%
  count(`country/ocean`, sort = TRUE, na.rm = TRUE)

# Also count occurrences of each 'bin_uri', sorted by frequency and removing NAs
df_bold %>%
  count(bin_uri, sort = TRUE, na.rm = TRUE)



# Group data by 'bin_uri' and count the number of records per BIN
df_CountbyBIN <- df_bold %>%
  group_by(bin_uri) %>%
  count(bin_uri)


#separate the coordinates into longitude and latitude
df_boldd <- df_bold %>%
separate(coord, into = c("Lat","Long"), sep = ",", convert = FALSE) %>%
  mutate(Lat  = as.numeric(str_replace_all(Lat, "[^0-9.-]", "")), Long = as.numeric(str_replace_all(Long, "[^0-9.-]", "")))

df_bold <- df_bold.subset %>%
  separate(coord, into = c("lat","long"), sep = ",", convert = FALSE) %>%
  mutate(lat  = as.numeric(str_replace_all(lat, "[^0-9.-]", "")), long = as.numeric(str_replace_all(long, "[^0-9.-]", "")))


#so I can have a new subset to be;
df_bold.subset2 <- df_bold[, c("processid", "bin_uri", "order", "species", "country/ocean", "realm", "long", "lat")]
df_bold.subset
#let us get things started
#Taking a look at the unique bins per the country
# Filter out rows with missing 'country/region' or 'bin_uri' values, then get distinct combinations of 'country/region' and 'bin_uri' and  now count the number of unique BINs per country/region
bins_country <- df_bold %>%
  filter(!is.na('country/region') & !is.na(bin_uri)) %>%
  distinct('country/region', bin_uri) %>%
  count('country/region', name = "n_bins", sort = TRUE)

# Sort and group the dataset by 'country/ocean', count unique BINs per country,
# and arrange results from highest to lowest BIN count
country_bins <- df_bold %>%
group_by(`country/ocean`) %>%
summarise(num_bins = n_distinct(bin_uri)) %>%
arrange(desc(num_bins))            


# Create a bar plot showing the number of unique BINs per country.
ggplot(country_bins, aes(x = reorder(`country/ocean`, -num_bins), y = num_bins)) +
geom_bar(stat = "identity", fill = "darkblue") +
theme_minimal() +
labs(title = "Number of BINs per Country", x = "Country", y = "Number of Unique BINs") +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Summarize BIN diversity by taxonomic order and visualize the results.
bins_order <- df_bold %>%
  group_by(order) %>%
  summarise(num_bins = n_distinct(bin_uri)) %>%
  arrange(desc(num_bins))

# Create a bar plot showing BIN diversity per order
ggplot(bins_order, aes(y = reorder(order, num_bins), x = num_bins)) +
geom_bar(stat = "identity", fill = "steelblue") +
coord_flip() +
theme_minimal() +
labs(y = "Order", x = "Number of Unique BINs", title = "BIN Diversity per Taxonomic Order")

#An important question to ask is what is the relationship between latitude and bin diveristy like
bin_richness <- df_bold %>%
  group_by(`country/ocean`) %>%
  summarise(mean_lat = mean(lat, na.rm = TRUE), BIN_richness = n_distinct(bin_uri)) %>%
  arrange(desc(BIN_richness))

head(bin_richness)

#Visualize the relationship between latitutde and Bin diversity
ggplot(bin_richness, aes(x = mean_lat, y = BIN_richness)) +
  geom_point(size = 3, color = "steelblue", alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "darkred", linewidth = 1) +
  labs(
    title = "Relationship Between Latitude and BIN Diversity",
    x = "Mean Latitude (°)",
    y = "BIN Richness (Number of Unique BINs)"
  ) +
  theme_minimal(base_size = 13)

# Summarize the number of unique BINs per country, arrange in descending order,
# select the top 10 countries, and visualize using a bar chart.
df_bold %>%
group_by(`country/ocean`) %>%
summarise(unique_BINs = n_distinct(bin_uri)) %>%
arrange(desc(unique_BINs)) %>%
head(10) %>%  #Here, I need top 10 countries

ggplot(aes(y = reorder(`country/ocean`, unique_BINs), x = unique_BINs, fill = `country/ocean`)) +
geom_col(show.legend = FALSE) +
coord_flip() +
theme_minimal() +
labs(y = "Country", x = "Number of Unique BINs", title = "Top 10 Countries by Number of Unique BINs DNA Barcoded")

#Create a species accumulation curve
# Generate a presence–absence matrix of BINs (columns) across countries (rows)
# Compute a normalized index showing how widespread each taxon is
taxa_country_count <- taxa_country_count %>%
  mutate(global_widespread_index = n_countries / max(n_countries))


# Display the 10 most globally widespread taxonomic orders
head(taxa_country_count[order(-taxa_country_count$global_widespread_index), ], 10)

# Calculate a species accumulation curve using random sampling order
species_matrix <- table(df_bold$`country/ocean`, df_bold$bin_uri)
sac <- specaccum(species_matrix, method = "random")

# Plot the species accumulation curve
plot(sac, col = "steelblue", lwd = 1.2, ylab = "Number of Samples", xlab = "Accumulated BIN Richness", main = "Species Accumulation Curve")


#ANother interesting discovery is right here
# Filter the dataset to include only rows with non-missing values for 'order' and 'bin_uri'.
df_boldsub <- df_bold %>%
  filter(!is.na(order), !is.na(bin_uri))

# Group data by 'order' and count the number of unique BINs for each order.
# Then arrange in descending order of BIN diversity
order_diversity <- df_bold %>%
  group_by(order) %>%
  summarise(Unique_BINs = n_distinct(bin_uri)) %>%
  arrange(desc(Unique_BINs))

# this is to show the top 10 most diverse orders based on number of unique BINs.
head(order_diversity, 10)

# Draw a  vertical bar chart showing barcode diversity across taxonomic orders
ggplot(order_diversity, aes(y = reorder(order, Unique_BINs), x = Unique_BINs, fill = order)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  theme_minimal(base_size = 13) +
  labs(title = "Barcode Diversity Across Orders", subtitle = "Number of unique BINs per order", y = "Order", x = "Unique BINs")

# Create a new data frame including the country and order and filtering the missing values

df_boldnew <- df_bold %>%
    filter(!is.na(`country/ocean`), !is.na(order)) 


#Count how many countries each taxon is found in

taxa_country_count <- df_bold %>%
  group_by(order) %>%                     
  summarise(n_countries = n_distinct(`country/ocean`), n_records = n()) %>%
  arrange(desc(n_countries))

# Display the first few rows (top orders by number of countries)
head(taxa_country_count)

#VIsualize the global widespread of the taxa by number of countries
ggplot(taxa_country_count, aes(y = fct_reorder(order, n_countries), x = n_countries, fill = n_countries)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  labs(title = "Global Widespreadness of Taxa (by Number of Countries)", y = "Order", x = "Number of countries") +
  scale_fill_gradient(low = "blue", high = "lightblue") +
  theme_minimal(base_size = 13)





