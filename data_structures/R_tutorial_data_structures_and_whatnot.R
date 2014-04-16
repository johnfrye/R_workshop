# Script to introduce students to basic data structures and data exploration

# intstalling R

#  go to www.r-project.org, and download a 'precompiled binary' 
#	for your system from the 'cran' site nearest you
#  follow links to 'download'
#  follow the directions to install

# starting R: double click the R icon, or double click a .R file (a script)

#<><>R as a calculator and simple operations<><><><>
# try these simple commands by typing them into the 'console'
2 + 2
c(1, 2, 3)
x <- c(1, 2, 3)

x
x <- c(1:10, 13)
y <- x^2
y

#  Scripting
#  The command line is fine, but it is not very repeatable.  Scripts allow us to 
#  1) document our analyses
#  2) repeat our analyses
#  2) share our methods


ls()  # this lists all objects you have created that are still in memory
rm(x) # Use the rm command to remove objects after you're finished.
rm(ls())  # removes everything!!!

# Now, open this script and we can get rolling!

# To execute a command from the script, simply highlight the part you want to run,
# 	and press command-enter (Mac) or ctrl-r (PC).
#	On the PC you can also press ctrl-r anytime and the entire line the cursor is
#	on will be executed.

# assignments
# everything in R is an 'object'.  That means that R knows something 
# 	about each object you create, and that there are 'methods' 
# 	available to work with objects of various classes.

# scalars
a <- 10
b <- 10
# print to console
b

# vectors
x <- c(3, 5, 2, 8, 5, 6, 8, 10)
y <- c(1, 3, 4, 5, 3, 9, 8, 3)
x
y

# other ways to create vectors

# from:to colon creates vector of all integers or factors from from to to
1:10

# seq creates a sequence from, to, by some step or to some length
seq(1, 10, .5)
seq(1, 10, length=5)

# rep just repeats some value or set of values
rep('a', 10)
rep(c('a', 'b'), 5)

# lists
# lists are arbitrary collections of objects
mylist <- list(aa = "foo", bb = 66)

mylist

mylist$aa

mylist[[2]]

str(mylist)
summary(mylist)
structure(mylist)

# character
char1 <-  "this is some text"
char1

# factors
# factors are like character vectors but only contain values in predefined 'levels'
fac1 <- factor(rep(c("a", "b"), 4))
fac1

levels(fac1)


# matrices (2 dimensional arrays)
xx <- matrix(c(1:6), 2, 3) # 
xx

# arrays
yy <- array(1:10, c(2,3,4))

yy

# dataframes
# dataframes are special objects in R that we will use a lot.  
# They are special because each column of a dataframe can 
#	itself be of a different class.
df <- data.frame(x, y, fac1)
# df <- data.frame(x, y, treatment = fac1)
df
dim(df)
summary(df)
str(df)
head(df)

# Transforming variables

df$lny <- log10(df$y)
df

plot(df$y ~ df$x)
plot(y ~ x, data = df)




# indexing 
# indexing dataframes
df[1,2]  #  returns row 1, column 2

df[1, 'y']  # same

df[,2]  # returns column y  as a vector

df[,'y'] # same

df['y']

df[1,]

df$y  # same

df[2:4,]  # rows 1-4

df[df$fac1=="a",]  # rows where fac1 is 'a'

df[df$x>mean(df$x),]  # rows where x is above the mean for x

mean(df$x)

df[order(df$x), ]  # returns whole dataframe sorted by x

df[df$y%in%c(3:7),]  # rows where 'y' matches a list of values

names(df)

df[2]

# indexing lists
mylist

mylist[1] #  item 1 in list, with label

mylist[[1]] # the value of item 1

myotherlist <- c(mylist, fac1)
myotherlist

myotherlist[[3]] [1]

# [ can return multiple values
# [[ can return only one value

plot(y ~ x, data = df)
mylm <- lm(y ~ x, data = df)
mylm
str(mylm)
mylm$coefficients
mylm[1]

mylm[[2]][3]

mylm[1]

abline(mylm)
abline(mylm$coefficients)

summary(mylm)
anova(mylm)

hist(mylm$residuals)

#################################
# did thru here in R workshop.
################################










# <><><> Getting Help <><><>

# what to do when you get an error
# 1. check braces (), commas, and other syntax for errors
# 2. check help page 
?data.frame  # searches local help for 'data.frame'
example(data.frame)  # gives example uses 
data.frame #  shows the function's internals, not always helpful
# 3. read the FAQ !!!
# 4. search the r-help archives
# 5. send a message to the r-help list.  Only as a last resort.  
# 	Formulate your question carefully and be sure to include a repeatable example

# finding methods, eg, for plotting, stats, etc
# search the R-help archives
# Rnews can be useful as a resource
# r-help list - only as a last resort; search the archives first!

# Excercise:  Now let's do some real work.  We will take three datasets, load them,
#	merge them, then do some exploratory data analysis to see what the data may 
#	tell us.  
# The datasets include:
#	1. Community data for a series of plots in 10 old fields at Cedar Creek
#	2. Some environmental data for the fields and plots
#	3. Trait data for some of the more common species in these fields
# We will use these data sets to see if trait abundances change with succession

# getting data into R

# Use getwd() to find out the current directory, and setwd() to change it.
# you can also do this through the menus (e.g., Misc > Set working directory... on a Mac). 
getwd()

# First, see what files are available here:
dir()

# now, let's load the plot data into R
# take a look at the data file "CCSppAbun.txt" in a text editor such as 
#	notepad or textedit ...   hmm, looks like it is comma delimited with a header row

sppAbun=read.table("CCSppAbunAnon.csv", header=TRUE, sep=",")
head(sppAbun) # look at the header and make sure it looks right

dim(sppAbun)  # the dimensions of the dataframe
str(sppAbun)  # lots of info about the df and its elements


# now, let's repeat that with our other two data files  
# Note that there are several read.XXX() commands.  We can use read.csv() as it
#	assumes a header line and comma delimiters in the default arguments.
#	The help page, ?read.table, will give you the defaults for these commands
envData=read.csv("CCEnvDataAnon.csv")
head(envData)
traitsRaw=read.csv("CCSppTraitsAnon.csv")
head(traitsRaw)

# whoa!  The trait file has all sorts of stuff in it that we do not want

# first, the columns that have 'var' in the name are the variance for each trait,
#	as the data came from 6 reps per species.
# So, lets drop all voriables (or columns) with 'var' in the name

traitsRaw <- traitsRaw[, -grep("var", names(traitsRaw))]  # the "-" means not those

# A WORD OF CAUTION: recycling an object name as I just did ("traits=traits[...") is 
#	a bit dangerous as you have written over the old traits, and it no longer exists.
#	In this case you can simply scroll up and reload the original, but that may not
#	always be the case
head(traitsRaw)

# let's now choose a few traits we would like to focus on, and name it 'traits'
#	so we can get back to the original if we want
traits <- traitsRaw[, c("sp","species","family","longevity","group","origin","C3C4",
	"seedsize","leafCN","height","leafarea","lma")]
head(traits)

# let's take a minute and look at these data
# one easy way to eyeball some data is with the pairs() command
pairs(traits[8:12])

# Graphics devices
# Let's stop for a minute and talk about graphics output.  If you simply issue 
#	a plotting cammand, a graphics window will appear and your figure will plot in it.
#	However, you can control the size of the window.
quartz(width=6, height=4) # on PC use windows(width=6, height=4, record=TRUE)

# Finally, you will want to output your figure to a platform-indepentent format
#	when you get it looking all pretty, and the best option for this is pdf().
#	after you finish the figure, use dev.off() to close the connection
# Two notes on devices and output.
#	1. On the mac (or PC) you can run a series of figures to pdf() so you can then 
#		scroll through them and avoid opening dozens of windows on your screen.


# ...mean while, back to our analyses

# wow, it looks like we should glance at some histograms, as some of these look highly skewed
#	Also, we have one species with gigundus leaves - I wonder if that is a typo?

head(traits[order(traits$leafarea, decreasing=TRUE),])
# Nope, its just Frances with its huge leaves.

# We also have something with rather large seeds...
head(traits[order(traits$seedsize, decreasing=TRUE),])
# Right, Myrtle, a legume, would have large seeds

# So, let's look at some histgrams to see if we want to log transform some of these traits

hist(traits$seedsize)  #  very skewed
hist(log10(traits$seedsize)) # that's better

# We can even write a little loop to look at all of the continuous traits. The way a loop works is to
# 	walk through a set of values, performing some function. Here we want to go through each
# 	trait, making a new histogram each time. Here, the trait values are in columns 8 - 12 in the
#	"traits" dataframe, so we'll set up a loop with the placeholder "i" (although we could use
#	any character here) as the column of the data frame "traits", and tell R to replace this 
#	placeholder with each number from 8 to 12, making a histogram each time.
for (i in (8:12)) {
	hist(traits[,i], main=paste(names(traits[i])))
	}

# It looks like just leaf area and seedsize are highly skewed, so let's make new transformed
#	variables for them

traits$l10seedsize <- log10(traits$seedsize)
hist(traits$l10seedsize)

traits$l10leafarea=log10(traits$leafarea)
hist(traits$l10leafarea)

# Much better.

# Now, one of the things you will notice about R is that there are a thousand ways to do any
#	one thing.  Sometimes that is good, sometimes you just are just getting familiar
#	with a command or package and find that there is a better way to do things

# R has a set of plotting functions in the base package (including polt() and hist(), but it 
#	turns out there is an completely different set of plotting functions based on a package
#	called Grid. Grid graphics tend to be much more powerful, but customizing them can be 
#	a little more complicated.

# One example of the power is an alternative to pairs() called splom(), or Scatter PLOt Matrix,
#	in package "lattice"

# So, let's download our first package, lattice
# On Mac, go to the menu 'packages & data > package intstaller', choose 'get list'.  It will 
#	show a slew of available packages for download. Let's get 'lattice' and 'Hmisc' for now.
#	Just select them and then click 'intstall selected'

# When they are done downloading, you can load lattice into memory with:
library(lattice)

# Now we can really look at our data
head(traits)
splom(~traits[c(9,10,12:14)])  
# That looks about the same as pairs()
#	But here is where the fun is...
splom(~traits[c(9,10,12:14)] | origin, data = traits)  
# 	or better yet
splom(~traits[c(9,10,12:14)], groups= longevity, data = traits, auto.key = TRUE)  

# 	or even combine the two
splom(~traits[c(9,10,12:14)] | origin, groups = longevity, data = traits, auto.key = TRUE)  

# Well, that's all well and good, but what is up with this data?
#	Looks like we have some issues here -- Like three species which are neither 
#	native nor introduced

# Now lets go ahead and drop all species where we lack complete data
# We can select rows that have NA's for a given trait

traits[is.na(traits$leafCN),]

# or those that have leaf CN

traits[!is.na(traits$leafCN),]  # the ! means not

# better yet, we can select all rows that lack any NA's
traits[complete.cases(traits),]

# but that really tosses a lot of species, mostly because we lack seedsize and CN,
# 	so lets drop those traits and see if we get a few more species
traits=traits[,-(c(8:9, 13))]
traits=traits[complete.cases(traits),]
traits
# Good, that dropped fungi, and moss, and all the species where we lack any traits
# lets drop Florence as well
traits[traits$sp=="Florence",]

traits=traits[traits$sp!="Florence",]

# Now let's straighten out 'longevity
levels(traits$longevity)

# First, let's change all the 'biannual' to 'biennial' 

traits[traits$longevity=='biannual', 'longevity']
traits[traits$longevity=='biannual', 'longevity'] = 'biennial'
traits$longevity=factor(traits$longevity) # drops unused levels
levels(traits$longevity)

# ok, let's eyeball the data again
head(traits)
splom(~traits[c(8,10,11)] | origin, groups= longevity, data=traits, auto.key=TRUE)  



# stopping here for today












# Next, let's loox at the plot data
head(sppAbun)

# we have a lot of data here, and we are really interested in values at the field level
#	so let's summarize the abundances at the field level

# There are lots of tools for summarizing ( eg, apply(), tapply(), aggregate(), 
#	or even writing your own loop to page through fields and species).
# I have found that summarize() from the Hmisc package works great for this 
#	sort of thing.

library(Hmisc)
fieldAbun <- with(sppAbun, summarize(X = cover, by = llist(sp, field), 
	stat.name = "totalcover", FUN = sum))
head(fieldAbun)

# ok, that gave us total cover by field, but we need to know how much area 
#	was sampled in each field.  Each plot is a 100cm line transect, and the 
#	cover values are in cm on that transect, so to get relative cover,
# 	we simply need to sum the total number of plots*100

head(sppAbun)

plotAbun <- with(unique(sppAbun[,1:3]), summarize(X = plot, by = llist(field), 
	stat.name = "plotsPerField", FUN = length))
	
# Now let's merge those data with the field abun data
fieldAbun <- merge(fieldAbun, plotAbun, All=TRUE)
# And now calculate cover at the field level
fieldAbun$cover <- fieldAbun$totalcover/(fieldAbun$plotsPerField*100)
fieldAbun <- fieldAbun[,-c(3,4)]
head(fieldAbun)

# Check to see if total cover is reasonable
with(fieldAbun, summarize(X = cover, by = llist(field), stat.name = "totalcover", FUN = sum))
#  Odd, cover seems not to be numeric...
str(fieldAbun$cover)
fieldAbun$cover <- as.numeric(fieldAbun$cover)
with(fieldAbun, summarize(X = cover, by = llist(field), stat.name = "totalcover", FUN = sum))
# That's better.  Cover ranges from 21% to 54%, OK.

# Now, Let's merge these data with the env data
head(fieldAbun)
head(envData)

# first we will aggregate the env data by taking the mean soil N for each field
fieldEnv <- with(envData, summarize(X = soilPctN, by = llist(field), 
	stat.name = "soilN", FUN = mean, na.rm = TRUE))  # na.rm=TRUE drops the NA's
# Now merge that result with envData
fieldEnv <- merge(fieldEnv, unique(envData[, c(1,6,7)]))
# Now merge that with fieldAbun
fieldAbun <- merge(fieldEnv, fieldAbun)

# now lets add the trait data

fieldSpp <- merge(fieldAbun, traits)  
# note that all species for which we lack trait data have been dropped
head(fieldSpp)

# OK, now that are data are together, let's test a few hypotheses about 
#	succesion and traits

# First, let's plot weighted leaf area versus field age

LeafArea <- with(fieldSpp, summarize(X=cbind(l10leafarea, cover), 
  by=llist(field, soilN, lastCrop, age07), 
	stat.name="wLogLeafArea", FUN=weighted.mean, na.rm=TRUE))  # na.rm=TRUE drops the NA's

# Let's plot that result
plot(wLogLeafArea ~ age07, data=LeafArea)
# Cool, so leaf size increases with field age
# We can run a simple stats test using lm()
leafAreaMod=lm(wLogLeafArea ~ age07, data=LeafArea)
summary(leafAreaMod)
anova(leafAreaMod)
leafAreaMod # what's inside the model results
str(leafAreaMod)  # even more inside the model output

#  We can use the model output to add a line
abline(leafAreaMod$coef)

# An alternative plotting method is to use xyplot
xyplot(wLogLeafArea ~ age07, data=LeafArea)
# We can dress it up a bit
xyplot(wLogLeafArea ~ age07, data=LeafArea, col="black",
	ylim=c(0.25, 0.37),
	xlab="Field Age (years)", ylab=expression("Weighted mean leaf size (cm"^-2*")"),
	scales=list(
		x=list(tck=c(1,0)),
		y=list(tck=c(1,0), at=log10(seq(1.8, 2.3, 0.1)),  labels=seq(1.8, 2.3, 0.1))))

# Now, does height vary within longevity classes as succession proceeds?

longXht=with(fieldSpp, summarize(X=cbind(height, cover), 
  by=llist(field, longevity, soilN, lastCrop, age07), 
	stat.name="wHt", FUN=weighted.mean, na.rm=TRUE))  # na.rm=TRUE drops the NA's

# Let's plot that result
xyplot(wHt ~ age07, groups=longevity, data= longXht,
	pch=1, col=c("red", "blue", "green"),
	ylim=c(8, 23),
#	abline=list(leafAreaMod$coef),
	xlab="Field Age (years)", ylab=expression("Weighted mean height (cm)"),
	scales=list(
		x=list(tck=c(1,0)),
		y=list(tck=c(1,0))),
	key=list(cex=0.8, points=list(pch=1, col=c("red", "blue", "green")),
		text=list(levels(longXht$longevity),
		x = 0.58, y = 0.72, corner = c(0, 0))))
		
# we can dress it up a bit more 
plotcol=c("red", "blue", "green")
xyplot(wHt ~ age07, groups=longevity, data= longXht,
	panel=panel.superpose, col = plotcol,
	panel.groups=function(x, y, subscripts, groups, ... ) {
		panel.xyplot(x, y, ...)
		p = unique(longXht$longevity[subscripts]) # p='annual'
		fit=lm(longXht[longXht$longevity==p,"wHt"] ~ longXht[longXht$longevity==p,"age07"])
		panel.curve(coef(fit)[1] + (coef(fit)[2] * x),			
    lwd = 2, type = "l", col = plotcol[match(p,levels(longXht$longevity))])
		},
	ylim=c(8, 23),
	xlab="Field Age (years)", ylab=expression("Weighted mean height (cm)"),
	scales=list(
		x=list(tck=c(1,0)),
		y=list(tck=c(1,0))),
	key=list(x = 0.3, y = 0.92, cex=0.8, points=list(pch=1, col=c("red", "blue", "green")),
		text=list(levels(longXht$longevity))))

# Lattice graphics are really powerful, but the learing curve is tough
# You can make the same plot in base graphics by plotting just one longevity
#	first, and then add in the other series
plot(wHt ~ age07, data= longXht[longXht$longevity=="annual",], col= "red")
points(wHt ~ age07, data= longXht[longXht$longevity=="biennial",], col= "blue")
points(wHt ~ age07, data= longXht[longXht$longevity=="perennial",], col= "green")
# and then add in the ablines....

# Now, is this relation real?
longXhtFit=lm(wHt ~ age07*longevity, data=longXht)
anova(longXhtFit)

# OK, so height differs by longevity, but mean height does not vary with field age or 
#	the interaction between field age and longevity

# Finally, lets look at mean soil N by field age
# Here we will go back to the original envData 
head(envData)

# But first, a bit about functions
# In R you can write your own functions for stuff that you do often

# This code make a new funtcion called meanstderr, which calculates the mean,
#	standard error, and replication for groups of data.
# By itself this is pretty trivial, because you can calculate any of these 
#	easily.  I use this function inside summarize(), because summarize()
#	takes only one function as an argument, and I often want to get means
#	and errors for exploratory data analysis and plotting
meanstderr <- function(x) c(Mean=mean(x,na.rm=TRUE),
	se=sqrt(var(x,na.rm=TRUE)/length(na.omit(x))), n=length(na.omit(x)))

# let's try it!
meanstderr(c(2,4,5,2,3))

# Now, if you forget what it does, you can always look inside:
meanstderr

# ok, so now let's see haw soilN varies with age
fieldN=with(envData, summarize(X=soilPctN, by=llist(field, lastCrop, age07), 
	stat.name="meanN", FUN=meanstderr))  # na.rm=TRUE drops the NA's

# Now we can plot mean N, with confidence intervals, across the successional gradient

# Note we are using xYplot from Hmisc, not xyplot from lattice  
#	(and Cbind from Hmisc, not cbind from base)
# Frank Harrell likes to improve on things, and also likes upper case!

trellis.par.set("grid.pars", list(lex = 2, cex=1.2))
xYplot(Cbind(meanN, se)~age07, groups=lastCrop, data=fieldN, ylim=c(.02,.08), 
	pch=20,
	xlab="Years since abandonment",
	ylab="Soil Nitrogen (%)")
Key(.7,.4)

# To wrap up, we can print this figure to pdf()

pdf(file = "MyPlot.pdf", width = 6, height = 6) 
trellis.par.set("grid.pars", list(lex = 2, cex=1.2))

xYplot(Cbind(meanN, se)~age07, groups=lastCrop, data=fieldN, ylim=c(.02,.08), 
	pch=20,
	xlab="Years since abandonment",
	ylab="Soil Nitrogen (%)")
Key(.8,.2)

dev.off()

# you might want to save your data objects for later use
save(list=c("fieldAbun", "fieldSpp"), file="mydata.Rdata")
# you can load it later with
load(file="mydata.Rdata")


# I hope you have enjoyed the fun!!!

# Now, on your own load your own dataset into R and do some exploratory data analysis!

