
# It is the first script for introduction in R
# Author: Aleksei Zamalutdinov
# PhD student, Project Center for Agrotechnologies, Skoltech

#-------------------------------------------------------------------------

# Let's start with style of R scripts
# Script consists of two parts: executable (active part) and comments
# All comments starts from # symbol

It is the executable part # select this line and click Run
# it raise an error, if we try to run it

# Try again with it
"It is the executable part" # it works now and print our text in Console

# Let's try with numbers
7+3
# We can do arithmetic inside Console!
# Try with other arithmetic operations and explain result
# -
# *
# /
# %%
# %/%
# **
# ^
# So we can use R as calculator

#-------------------------------------------------------------------------

# We made many calculations, but all our results were just printed in console
# Let's save them

a <- 7+3 # this arrow (<-) means that we assign value to the variable a. Alt+- to draw it
a # we can see what is inside variable a

b <- 9
a + b

#-------------------------------------------------------------------------

# We have already worked with numbers and text, so let's discuss other format types
# There are 4 main data types in R

1      # integer - whole numbers, use less space
1.01   # numeric - all numbers
"Text" # character - texts
TRUE   # logical - stores TRUE or FALSE states

# It is rather obvious why we need and how to save first three of them
# So let's try to perform this command

1>2      # FALSE, are we sure that it is logical?
str(1>2) # There is special command to check data type
# Try to use it with other values. Type your code below

#-------------------------------------------------------------------------

# We faced with a new type of things like str(1>2). They are called functions.
# Functions may combine several operations inside and produce some result.
# So, str - is the name of function and everything in brackets is called arguments.
# Arguments may contain input data and parameters that affect function performance
# However, some of the functions may work without any arguments

getwd() # shows the path to this script

#-------------------------------------------------------------------------

# There are several formats how to organize our data, when we work with more than one value
# Basic and easy format for storing data is vectors
a <- c(1,2,3)
b <- c(1,"2",3)
c <- c() # Let's create vector with logical, characters and numbers

# Vectors can store values of the one type: numeric, characters or logical
# That is why, if we try to add different types, it will convert them to one

# We can add values to vector and concatenate them
d <- c(a,5,a)

# We can perform all arithmetic not only on a single value, but also on vector

# Write below calculations with different operators
# +
# -
# *
# /
# %%
# %/%
# **
# ^

# We can also perform some basic statistical operations like:

min(d)    # minimum value
max(d)    # maximum value
range(d)  # minimum and maximum value as a new vector
mean(d)   # mean of values
length(d) # number of elements in vector

# And even compare vector with some value

d>2

# Basing on these functions let's do small task
# We will use the value of Lake Huron level from LakeHuron dataset
?LakeHuron # question mark open the help page of function and datasets
Lake <- as.vector(LakeHuron) # we will save it in a new variable

# Calculate the difference between lowest and highest level of the lake 
# and save it in my_answer variable


my_answer <- 
print(paste("Difference of the lake levels is",my_answer,"feets. Not so big difference"))

#-------------------------------------------------------------------------

# Also we can choose one or several particular numbers from vector 
# There are several ways to do that:
# 1) position of element in vector (index). It starts from 1!
# Check how it works
d[4]
d[c(1,3)]
d[1:3]

# 2) Logical vectors
e <- d>2
d[e] # or we can do the same in one line
d[d>2]

# 3) Positions to exclude
d[-4]
d[-(1:3)]

# 4) Selection using names
fruit_and_vegetables <- c(150, 70, 80, 200, 99)
names(fruit_and_vegetables) <- c("cucumber", "banana", "apple", "tomato", "mandarin")
fruit_salad <- fruit_and_vegetables[c("banana","apple","mandarin")]
sum(fruit_salad)

# It is time for a new small task
islands <- islands # it is one of the default datasets in R

# Using indexing, remove continents, select 10 smallest islands and 
# save one that you would like to visit to the new variable, 
# then convert the area from square miles to square km

# You may need to use sort(your_vector) to arrange values from smallest to largest

my_answer <- 
  
print(paste("My favorite island is",names(my_answer),", its area is",as.numeric(my_answer), "sq km"))

#-------------------------------------------------------------------------
# Another way to store data is matrices. 
# In R they are the concatenation of vectors of the same size and data type
# let's add information to our fruit example
fruit_and_vegetables_price <- c(150, 70, 80, 200, 99)
fruit_and_vegetables_weight <- c(0.4,0.3,0.25,0.2,0.1)

# Now we can merge two vectors into one 2D matrix and receive table 
fruit_and_vegetables <- array(c(fruit_and_vegetables_price,fruit_and_vegetables_weight), dim=c(5,2))

# We can change values inside matrix using indexing [row,column]
fruit_and_vegetables[2,1] <- 80

# We can also use indexing to select columns and use them in calculations
fruit_and_vegetables_price_per_fruit <- fruit_and_vegetables[,1]*fruit_and_vegetables[,2]

# Then we can add result into matrix as one more column
fruit_and_vegetables <- cbind(fruit_and_vegetables, fruit_and_vegetables_price_per_fruit)
colnames(fruit_and_vegetables) <- c("price", "weight", "price per fruit")
fruit_and_vegetables
row.names(fruit_and_vegetables) <-c("cucumber", "banana", "apple", "tomato", "mandarin")
fruit_and_vegetables

# Sometimes it is also useful to create an empty matrix with 0
empty_matrix <- array(0, dim = c(4,2))

# Create an empty matrix with dimensions 3,6 and fill it with some values using indexing



# We can also perform some operations with matrices if they are of the same size

matrix_1 <- array(1, dim = c(4,2))
full_matrix <- array(seq(1,6), dim = c(4,2))

matrix_1 + full_matrix
matrix_1 - full_matrix
matrix_1 * full_matrix
matrix_1 / full_matrix

# And we can transpose it
t(full_matrix)

#-------------------------------------------------------------------------
# However, the most useful way to organize data in R is data frames
# They are set of vectors, like arrays and matrices, 
# but vectors may be of different data types

fruit_and_vegetables_names <- c("cucumber", "banana", "apple", "tomato", "mandarin")
fruit_and_vegetables_price <- c(150, 70, 80, 200, 99)
fruit_and_vegetables_weight <- c(0.4,0.3,0.25,0.2,0.1)

fruit_df <- data.frame(fruit=fruit_and_vegetables_names, price=fruit_and_vegetables_price, weight=fruit_and_vegetables_weight)
fruit_df
# We can also convert matrix into dataframe
df_from_matrix <- as.data.frame(fruit_and_vegetables)
df_from_matrix

# We can select values from data frame in a same manner as in matrix
fruit_df[3,1]
# However, we can also use column names to select one column (actually, one vector)
fruit_df$fruit


# Let's create a data frame with your favorite books

book_title <- c()
author <- c()
my_books <- data.frame()
my_books

# We can add columns to dataframe in a same way as selecting columns. 
# Let's add your recommendations column to your data frame. 
# Will you recommend this book to read?

my_books$recommendation <- c()
my_books

# Let's add one more book
rbind()
#-------------------------------------------------------------------------

# We learnt how to perform some basic operations in R and now we are ready to move further

# Now we will see how to manipulate data more efficiently

# For that we will use family of tidy packages package. What is package?
install.packages("tidyverse")
library(tidyverse)
# Package is a set of functions that extend the functionality of R

table4a <- table4a
table4a # cases of tuberculosis

# In order to convert it to tidy data we can use pivot_longer function
table4a_long <- pivot_longer(table4a, cols = 2:3, names_to ="year", values_to = "cases")
table4a_long
# In such format it is easier to manipulate with data

# However, it is not so easy to look at it as usual table
pivot_wider(table4a_long, names_from = "year", values_from = "cases")

# Sometimes we may need to divide some value into several columns
table3 <- table3
table3_sep <- separate(table3, rate, sep = "/", into = c("cases", "pop"))
table3_sep

# In contrast, sometimes we need to combine 2 columns
unite(table3_sep, country, year, col = "country_year", sep = "_")

# Let's imagine that we want to make our experiment in several combinations of conditions

experiment <- tibble(
  temperature   = c(25, 25, 37, 42),
  time  = c(10,20,10,20),
  result = c(456, 849, 67, 47)
) 

# Looks nice, however, we have "hidden" missing data
# We made 4 tests, but actual number of combinations is 6
expand(experiment, temperature, time)

# We can expand our table and introduce NA
complete(experiment, temperature, time)

# Small task
field <- tibble(
  variety   = c("Bright", "Red", "First", "Leader"),
  yield_2021  = c(4.18, 2.54, 1.92, 0.17),
  yield_2022 = c(1.88, 0.59, 0.92, 1.54)
)
field
# Transform this dataset to tidy format. It should be like this

# variety parameter year  yield
# <chr>   <chr>     <chr> <dbl>
# 1 Bright  yield     2021   4.18
# 2 Bright  yield     2022   1.88
# 3 Red     yield     2021   2.54
#-------------------------------------------------------------------------

# Now, we will learn how to upload data from files

# We can use full paths to select files, however sometimes it is not convenient and easy to read
# This is a trick to detect which folder contains the R script and associated data
main_dir <- dirname(rstudioapi::getSourceEditorContext()$path) 
setwd(main_dir)

# As for me, I prefer to use button "Import Dataset",  
# because it also allows to preview data before import

# However, we can just use simple command like this
Life_expectancy <- read.csv("Life expectancy.csv")
# And we can save it 
write.csv(Life_expectancy, file="Life expectancy_new.csv")

# Let's have a look at our data
str(Life_expectancy) # We have 1 column with text and 4 columns with numbers
# Sometimes it is also helpful to open the table and look at it from other point of view
View(Life_expectancy)

#-------------------------------------------------------------------------

# Let's look through its functions
filter(Life_expectancy, Happiness.Score >7) # It helps to select rows that match some condition
# Select countries that have men life expectancy higher than 81
filter()

# We can also use filter with such functions as is.na() or %in%
filter(Life_expectancy, Country %in% c("Australia", "Austria"))

# However, filter works with some logical statements. 
# Sometimes we need to select rows per position. Select rows from 10 to 15 using indexing
Life_expectancy[,]

# dplyr has more understandable analogue - slice
slice(Life_expectancy, 10:15)

# slice is like a family of functions and 
# they provide additional functionality in comparison with usual indexing

slice_sample(Life_expectancy, n=10) # select random 10 rows
slice_min(Life_expectancy, Life.Expectancy..years....Men, prop=0.2) # select 20% min values
slice_max(Life_expectancy, Life.Expectancy..years....Men, prop=0.2) # select 20% max values
slice_head(Life_expectancy, n=5) # select first 5 rows
slice_tail(Life_expectancy, n=5) # select last 5 rows

#-------------------------------------------------------------------------

# In some cases we would like to add observations to our dataset
# Let's try to make subsets
Life_expectancy_bot <- slice_min(Life_expectancy, Life.Expectancy..years....Men, prop=0.2)
Life_expectancy_top <- slice_max(Life_expectancy, Life.Expectancy..years....Men, prop=0.2)
# And then merge them in one
bind_rows(Life_expectancy_bot, Life_expectancy_top)

#-------------------------------------------------------------------------

# Another cool option with dplyr is the piping. Ctrl+Shift+M will draw %>% 
# It means that dataset in output will be used as input in next command
# Let's try with summarise function, it allows to calculate statistics for the whole column
# We will calculate mean of men life expectancy
filter(Life_expectancy, Happiness.Score >7) %>% 
  summarise(avg = mean(Life.Expectancy..years....Men))

#-------------------------------------------------------------------------

# We can make a subset of our dataset by selecting several columns
# As a reminder, select columns with country and women life expectancy using indexing
Life_expectancy[,]

# dplyr offers us one more way to select columns just using their names.
# It is much more useful when we work with large datasets and we don't want to count columns

country_women <- select(Life_expectancy, Country, Life.Expectancy..years....Women) %>%
  slice_max(Life.Expectancy..years....Women, prop=0.2)
#-------------------------------------------------------------------------

# Next we will learn how to modify and create new columns
# Just for comparison, calculate the difference between women and men life expectancy
Life_expectancy$diff <- 
  
Life_expectancy <- Life_expectancy %>% mutate(diff_dplyr=Life.Expectancy..years....Women - Life.Expectancy..years....Men)
# So, now we can check the result
View(Life_expectancy)

#-------------------------------------------------------------------------

# Sometimes we need to sort our table by some of the variables

Life_expectancy <- Life_expectancy %>% 
  mutate(diff_dplyr=Life.Expectancy..years....Women - Life.Expectancy..years....Men) %>% 
  arrange(diff_dplyr)
View(Life_expectancy)

#-------------------------------------------------------------------------

# We can make groups based on our variables. 
# For example, let's use Happiness.Score_level as grouping factor 
Life_expectancy_grouped <- Life_expectancy %>% group_by(Happiness.Score_level)
Life_expectancy_grouped
# We see such text Groups:   Happiness.Score_level [3]
# It means that this variable is grouping factor
# What opportunities does it provide us?
# Now we can calculate statistics per group
Life_expectancy_grouped %>% 
  summarise(avg=mean(Life.Expectancy..years....Men))

# We can also count the number of observations in each group
Life_expectancy_grouped %>% tally()

#-------------------------------------------------------------------------

# Let's practice our skills in mtcars dataset. It is some historical data about US cars
mtcars_df <- mtcars
# We would like to choose the most economical and practical car
# I suggest to do several steps

# 1) convert mpg to km per liter value using mutate function and then to liters per 100 km
# 2) select 25% most economic cars based on this parameter
# 3) filter cars with automatic transmission from these 25%

# And we want to compare, how weight is important for the fuel efficiency
# 4) Let's return back to the full dataset.
#    Divide cars into 3 groups by weight and calculate mean liter per 100 km
# 5) Add column with money expense on fuel per 100 km (in rubles). 
#    Let's assume that fuel price is 48 rubles per liter

#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
# We discussed how we can manipulate our data. However, it is also 
# important to present data in understandable manner
# We will discuss how to use ggplot2 package and customize our plots
install.packages("ggplot2")
library(dplyr)
library(ggplot2)
main_dir <- dirname(rstudioapi::getSourceEditorContext()$path) 
setwd(main_dir)

diamonds_df <- diamonds # It is rather big
diamonds_subset <- slice_sample(diamonds_df, prop=0.05) # Not so big
#-------------------------------------------------------------------------
# First, let's discuss the structure of the command
# ggplot(our_dataset, aes(which columns we would like to draw)) +
# geom_bar() type of the plot

# Now we can try this very basic command
ggplot(data=diamonds_subset, aes(x=cut))+ 
  geom_bar()

#-------------------------------------------------------------------------
## One continuous variable

# There are a lot of plot types that we can use for visualization.
# We will start from the cases when we have one continuous variable.

# This plots help us to demonstrate the distribution of our variable
# There are several ways to draw such plots

# This one shows us the smoothed version of our real distribution 
ggplot(data=diamonds_subset, aes(x=price))+ 
  geom_density()

# Below we will see several plot types that merge nearest points into bins,
# On the one hand, it reduce detalization, on the other hand it helps us to smooth distribution 

# More common analogue with the same result is histogram. 
# Height of the bar depends on the number of dots
ggplot(data=diamonds_subset, aes(x=price))+ 
  geom_histogram(binwidth=100)

# Try to visualize some other variable from our diamonds_subset dataframe


#-------------------------------------------------------------------------
## One discrete variable

# It is also possible to visualize distribution of discrete variables but 
# the choice is limited
# The only way to do it is geom_bar()
ggplot(data=diamonds_subset, aes(x=cut))+ 
  geom_bar()

# let's see what will change if we use price instead of cut. Please, write the command below


#-------------------------------------------------------------------------
## Two continuous variables

# Another idea of using plots is to show relationship between different variables
# Again we have different plots for different variable types
# Let's start from both continuous variables

# geom_points shows us plot with points and each dot represent one observation
ggplot(data=diamonds_subset, aes(x=carat, y=price))+ 
  geom_point()

# However, trend may be not so obvious. In this case we can use geom_smooth.
# It tries to fit some model (linear regression in our case)
ggplot(data=diamonds_subset, aes(x=carat, y=price))+ 
  geom_smooth(method = "lm")

# Each geom is just a layer, we can combine them to produce more informative plots
ggplot(data=diamonds_subset, aes(x=carat, y=price))+ 
  geom_point()+
  geom_smooth(method = "lm")

# Try these plots with some other pair of variables

#-------------------------------------------------------------------------
## One continuous variable and one dicrete

# When we have one continuous and one discrete variable we should use other plots
# For sure, we can use previous set of plots but it will be rather useless
ggplot(data=diamonds_subset, aes(x=cut, y=price))+ 
  geom_point()

# The most common but the most useless way is to use columns
# It shows the mean for each group, but we don't have any idea about distribution and range of observations
ggplot(data=diamonds_subset, aes(x=cut, y=price)) + 
  geom_col()

# More applicable for scientific research are boxplots. 
# They show the median, first and third quartiles, 
# whiskers up to 1.5 inter-quartile range and all outliers as dots
ggplot(data=diamonds_subset, aes(x=cut, y=price)) + 
  geom_boxplot()

# We see distribution, but we are missed some valuable statistics
# It is possible to combine these 2 plots in violin plot
ggplot(data=diamonds_subset, aes(x=cut, y=price)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))

#-------------------------------------------------------------------------
## Two discrete variables

# Okay, that's great but what if we have only discrete variables

# There are several options to create plots 
# We can represent each pair of our discrete values as a dot, and 
# its size will represent the number of observations with such parameters
ggplot(data=diamonds_subset, aes(x=cut, y=color)) + 
  geom_count()

#-------------------------------------------------------------------------
# Let's make some practice

# I suggest to play with Barcelona cars dataset
# 790 cars are included in the dataset, being the cars that were for sale on the Flexicar website 
# for the province of Barcelona at the time of extraction, described by the following attributes:
#   
# brand - car brand
# model - car model
# price - car price
# engine - car engine
# year - production year
# mileage - car mileage
# fuel - fuel type
# gearbox - gearbox type
# location - car location

# Import the dataset, explore the distribution of variables and build some of
# the pairwise plots in order to, for example, 
# check dependency of price on different parameters graphically. 

#-------------------------------------------------------------------------

# ggplot2 is really flexible tool to make plots
# Let's discuss how to customize your plot
# the first way to make your plot more informative is aesthetics
?geom_point()
# We can add one more variable to the plot as color
ggplot(data=diamonds_subset, aes(x=carat, y=price, colour=cut)) + 
  geom_point()
# Add one more variable with another aesthetics



# Let's start from axis names, it can be done inside labs()
ggplot(data=diamonds_subset, aes(x=carat, y=price, colour=cut)) + 
  geom_point() +
  labs(x="Carats", 
       y="Price, $", 
       title = "Price of diamonds with different weights and quality",
       subtitle = "2697 diamonds in total",
       colour = "Quality of cut")

# Let's add such labs to our plot with more variables

#-------------------------------------------------------------------------
# Okay, we made our plots rather informative, but still we can improve its design
# We can do it with using theme() family of function
# As an easy option, there are several ready themes

ggplot(data=diamonds_subset, aes(x=carat, y=price, colour=cut)) + 
  geom_point() +
  theme_light()

# Try other themes and see the effect

# theme_bw()
# theme_gray()
# theme_dark()
# theme_classic()
# theme_linedraw()
# theme_minimal()
# theme_void()

#-------------------------------------------------------------------------
# It is all really cool, but how to save picture

# We can store not only some data in variables but also plots
a <- ggplot(data=diamonds_subset, aes(x=carat, y=price, colour=cut)) + 
  geom_point()
a

# Then we can save it as file
ggsave("plot.png") # by default it is the last plot

# Using variables we can save any plot. We can also specify parameters of the picture quality
ggsave("plot.jpg",plot=a,dpi = 100,width = 10, height = 5, units = "cm")

#-------------------------------------------------------------------------
# Faceting
# We can multiple plots at the same time in order to show the effect of several variables
# Here we visualized 5 variables at the same time
ggplot(data=diamonds_subset, aes(x=carat, y=price, colour=clarity)) + 
  geom_point() + 
  facet_grid(rows = vars(cut), cols = vars(color))

# Try to visualize 5 variables from barselona cars dataset with and without facets
# Let's compare these 2 plots

#-------------------------------------------------------------------------

# Let's discuss using this simple example how important to use graphs in your analysis
head(anscombe)

# Rename columns in this data
anscombe.1 <- data.frame(x = anscombe[["x1"]], y = anscombe[["y1"]], Set = "Anscombe Set 1")
anscombe.2 <- data.frame(x = anscombe[["x2"]], y = anscombe[["y2"]], Set = "Anscombe Set 2")
anscombe.3 <- data.frame(x = anscombe[["x3"]], y = anscombe[["y3"]], Set = "Anscombe Set 3")
anscombe.4 <- data.frame(x = anscombe[["x4"]], y = anscombe[["y4"]], Set = "Anscombe Set 4")

# Combine dataframes together 
anscombe.data <- rbind(anscombe.1,anscombe.2 ,anscombe.3 ,anscombe.4 ) # add sets 
head(anscombe.data)

# However, we can also do like this
x_data <- pivot_longer(anscombe, cols = 1:4, names_to ="Label_x", values_to = "x") %>% separate(col = Label_x, sep = 1, into = c("Var","Set"))
y_data <- pivot_longer(anscombe, cols = 5:8, names_to ="Label_y", values_to = "y") %>% separate(col = Label_y, sep = 1, into = c("Var","Set"))
anscombe.data <- cbind(x_data[,6:7],y_data[,7])

###Calculate mean
anscombe.data %>% group_by(Set) %>% summarize(mean(x),mean(y),sd(x),sd(y),cor(x,y))

###Save the linear regression between x and y for each Set
model1 <- lm(x ~ y, filter(anscombe.data, Set == 1))
model2
model3
model4

##Get a SUMMARY statistics of obtained correlations
summary(model1)
model2
model3
model4

##Visualize data
ggplot(anscombe.data, aes(x = x, y = y)) + 
  geom_point(color = "black") + 
  facet_grid(cols = vars(Set)) + 
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE, data = anscombe.data)

#-------------------------------------------------------------------------

# Some practice
# https://www.kaggle.com/datasets/rajkumarpandey02/wheat-production-statistics

# 1) Import International_wheat_production_statistics.csv file
# 2) Check its structure and column types
# 3) Convert it to tidy format
# 4) Draw plot year*production
# 5) Choose 8 countries and make a subset with them
# 6) Draw plot year*production with selected countries in color



