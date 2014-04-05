#+ license, echo=FALSE
# This script compute and display a sigmoidal function
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

#+ libraries, echo=FALSE
source("functions.R")

#' Plot 1
#' ---
#' This plot represents the logistic function for b = .1
#+ plot1, echo=FALSE
y <- GenLogistic(time = -20:20, A = 0, K = 1000,
                 b = .1, v = .5, Q = .5, M = 0)
plot(x = -20:20, y = y)

#' Plot 2
#' ---
#' This plot represents the logistic function for b = 1
#+ plot2, echo=FALSE
y <- GenLogistic(time = -20:20, A = 0, K = 1000,
                 b = 1, v = .5, Q = .5, M = 0)
plot(x = -20:20, y = y)






