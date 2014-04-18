#+ license, echo=FALSE
# This is an introduction to data manipulation with dplyr and data plotting with
# ggplot2.
# 
# Copyright (C) 2014 Simon Garnier
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#+ libraries, echo=FALSE, message=FALSE
require("data.table")
require("dplyr")
require("ggplot2")

#+ data, echo=FALSE
data <- fread("movies.csv")

data.simplified <- select(data, title, year, length, budget, rating, Short)

data.shorts <- filter(data.simplified, Short == 1)

data.features <- filter(data.simplified, Short != 1)

#' Plot 1
#' ---
#' Relationship between year of production and budget size for feature movies.
#+ plot1, echo = FALSE, message=FALSE, warning=FALSE
ggplot(data.features,
       aes(x = year,
           y = budget)) +
  geom_point() + 
  geom_smooth()

#' Plot 2
#' ---
#' Relationship between year of production and duration for feature movies.
#+ plot2, echo = FALSE, message=FALSE, warning=FALSE
ggplot(data.features,
       aes(x = year,
           y = length)) +
  geom_point() + 
  geom_smooth() + 
  ylim(0, 1000)

#' Plot 3
#' ---
#' Relationship between year of production and duration for feature AND short
#' movies.
#+ plot3, echo = FALSE, message=FALSE, warning=FALSE
ggplot(data.simplified,
       aes(x = year,
           y = length,
           color = factor(Short))) +
  geom_point() + 
  geom_smooth(color = "black") +
  xlab("Year of production") +
  ylab("Duration (min)") +
  scale_color_discrete(name = "Movie\ntype",
                       labels = c("Feature", "Short")) +
  facet_wrap(~ Short, nrow = 2)

#' Plot 4
#' ---
#' Relationship between year of production and budget size for feature movies as boxplots
#+ plot4, echo = FALSE, message=FALSE, warning=FALSE
ggplot(data.features,
       aes(x = factor(year),
           y = budget,
           fill = year,
           color = year)) +
  geom_boxplot()


