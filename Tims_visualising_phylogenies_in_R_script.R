####################################
### VISUALISING PHYLOGENIES IN R ###
####################################

# first, clear your workspace so things are less cluttered
rm(list=ls()) 
# ls() returns a list of all objects stored in your current R workspace
# rm() removes them

.pardefault <- par() # save the default graphical paramters so it's easy to reset them later

# install the packaes which we'll need to run this script, if you don't already have them. 
# We'll use the CSIRO cran mirror
install.packages(c("ape","phytools","paleotree","geoscale","strap"), repos="https://cran.csiro.au/") 

# load in the packages used in this workshop
library(ape) # contains most of the core functions for manipulating trees in R
library(phytools) # some other useful tools and functions
library(paleotree) # used for timescaling trees
library(geoscale) # incroporating geological scales
library(strap) # for plotting geological timescales on phylogenies

# you can usually get the citations for any packages used like this:
citation("ape")


###################################
### tree objects and structures ###
###################################

tree <- rtree(n=10) # this will generate us a random tree with 10 tips

# you can find out more about a function by prepending a question mark to the function name
# ?rtree

# quick tree plot
plot(tree)

# let's look at the structure of our tree
str(tree) # it's an object of class "phylo" - it's basically a list of stuff

# within the list masquerading as a "phylo" object we have 4 things:

# 1:
# a vector containing the tip label names:
tree$tip.label

# you can check how many tips you have using
Ntip(tree)


# 2:
# the number of internal nodes in the tree
tree$Nnode # there are m = n - 1 nodes for a fully bifurcating tree

# let's plot the node labels onto the tree and take a look
 nodelabels() # internal nodes

# the node number of the root node is (number of taxa + 1)
# this is because the tips also get node numbers
tiplabels()


# 3:
# a matrix defining the edges (branches):
tree$edge

edgelabels() # row numbers are edge numbers
# the numbers in the tree$edge matrix tell you which two nodes are connected by each edge


# 4: 
# a vector detailing the length of each branch
tree$edge.length 
axisPhylo() # adds a scale so you can see the branch lengths


### we can also generate more than one tree at once

trees <- rmtree(N=2,n=30) # N is the number of trees we want
str(trees) # object is of class multiPhylo - a list of multiple "phylo" objects

# or bind together existing phylo objects
trees2 <- c(tree,tree)
trees2


### some arguments can be added or removed

# we can remove some elements from "phylo" objects

# our tree object has branch lengths
tree # "rooted; includes branch lengths"

# if you set edge.length to NULL, the tree won't include any branch lengths
tree2 <- tree
tree2$edge.length <- NULL

tree2 # "rooted; no branch lengths"
plot(tree2) # we can see the branch lengths have been removed and the tree is now somewhat arbitrary


# we can add a root.time element to our phylo object
# this tells us the age of our root node

# by default, the most recent taxon is assumed to be age 0 when plotting
plot(tree); axisPhylo() # draw tree with scale bar

# the root node age is inferred from that, but not included as a default argument

tree$root.time # doesn't specify the age of the root node

# we can calculate the root node age by looking at the variance-covariance matrix of the tree
# a variance-covariance matrix is a matrix representation of a phylogeny

round(vcv(tree),2) # print the vcv matrix with each element rounded to 2dp
# this shows the amount of branch length shared between each pair of taxa
# i.e. how much evolutionary history they share

# the diagonal of the matrix shows how much branch length a taxon shares with itself
# i.e. it shows how far each taxon is from the root node

# if we know the age of the most recent taxon, we can calculate the root node age with
max(vcv(tree)) # this is how far away the root node is from the most recent taxon

#if the most recent taxon is not at 0, we can calculate the root node age from
# max(vcv(tree)) + [age of most recent taxon]

# if our most recent tip is 5Ma old
tree$root.time <- max(vcv(tree)) + 5

plot(tree); axisPhylo()

# the most recent taxon now appears at 5, not 0
# and the root time is longer ago
tree$root.time

# you can set any root time you want
tree$root.time <- 25

plot(tree); axisPhylo()



#####################################
### importing and exporting trees ###
#####################################

# you can import trees with one of two functions, depending on your file type
# tree <- read.nexus("myNexFile.nex")
# tree <- read.tree("myTreFile.tre")
# tree <- read.tree("myNewickString.txt")

# or by writing in a newick string in plain text
newick.tree <- read.tree(text = "(((A,B,C),D,E),F);")

plot(newick.tree)

# trees can be saved in a similar way

# write.nexus(tree, file = "myNexFile.nex") # make sure to include "file = " here
write.nexus(newick.tree) # this is the format of a nexus file

# write.tree(tree, "myTreFile.tre")

# write.tree(tree, "myNewickString.txt") 
write.tree(newick.tree) # this is the format of a newick string

# if the trees have branch lengths, the respective files look like:
write.nexus(tree)
write.tree(tree)


#############################
### plotting a basic tree ###
#############################

# thus far, we've just used the most basic plot command
plot(tree) #R knows to use plot.phylo() on an object of class 'phylo' when you type plot()

# this function contains a number of different arguments, including:
plot(tree, edge.width = 4)  # edge.width alters the branch thickness
plot(tree, label.offset = 3) # label offset alters the distance between the branch tip and the tip label
plot(tree, tip.color = "cornflowerblue") # changes the colours of the tip labels
plot(tree, direction = "up") # rotates the way your tree is pointing
plot(tree, direction = "up", srt = -45) # srt (string rotation) rotates your tip labels
plot(tree, direction = "up", cex = 0.7) # cex (character expansion) changes your text size

# the default tree type is "phylogram", but there are a number of others you can use
# for example 

plot(tree, type = "cladogram", edge.width = 2, label.offset = 0.1)
plot(tree, type = "fan", edge.width = 2, label.offset = 0.1)
plot(tree, type = "unrooted", edge.width = 2, cex = 0.6, font = 1)

# there are a number of other arguments which can be included in plot()
# ?plot.phylo # they can be seen here


###################################
### manipulating tree structure ###
###################################


### binary and non-binary trees ###

#we can read in a simple tree with a polytomy from a newick string:
t1 <- read.tree(text = "(((A,B,C),D,E),F);")

#and plot it:
plot(t1, type = "cladogram",tip.color=c("black","cornflowerblue"),
	edge.width = 2.5,label.offset = 0.1, cex = 1.5)
# tip colours are included here to represent an arbitrary character state for each taxon

#is it fully dichotomous?
is.binary.tree(t1)

#multi2di randomly resolves polytomies in a multichotomous tree to produce a binary tree:
t2 <- multi2di(t1)

#now it will be binary:
is.binary.tree(t2)

plot(t2, type = "cladogram", tip.color = c("black","cornflowerblue"),
	edge.width = 2.5, label.offset = 0.1, cex = 1.5)

# note: randomly resolving polytomies [e.g. in your consensus tree] is not necessarily a good idea!)
# can be useful for certain purposes, though
# e.g. phylogenetic clustering of traits across trees
# resolve the tree multiple times and see how results vary with topology

par(mfrow=c(2,2),mar=c(0.2,0.2,0.2,0.2))
plot(multi2di(t1), edge.width = 2.5, label.offset = 0.1, 
	type = "cladogram",tip.color=c("black","cornflowerblue"),cex=1.5)
plot(multi2di(t1), edge.width = 2.5,  label.offset = 0.1,
	type = "cladogram",tip.color=c("black","cornflowerblue"),cex=1.5)
plot(multi2di(t1), edge.width = 2.5,  label.offset = 0.1,
	type = "cladogram",tip.color=c("black","cornflowerblue"),cex=1.5)
plot(multi2di(t1), edge.width = 2.5,  label.offset = 0.1,
	type = "cladogram",tip.color=c("black","cornflowerblue"),cex=1.5)

# the apparent phylogenetic clumping of the colours ("character states") depends on the resolved topology



### ultrametric and non-ultrametric trees ###

# rtree generates non-ultrametric trees - the tips are all different distances from the root
t3 <- rtree(n=30) # generate a random tree with 30 tips
par(.pardefault)
plot(t3, show.tip.label = FALSE) # plot it without the tip labels
axisPhylo() # adds a basic axis to the tree
is.ultrametric(t3) # is the tree ultrametric?

# we can use compute.brlen() to scale the branch lengths such that the tree is ultrametric
# this is a crude branch-length scaling method that's really only suitable for visualisation

t4 <- compute.brlen(t3, method = "Grafen", power = 0.5)
is.ultrametric(t4)
plot(t4, type = "fan", edge.width = 2, label.offset = 0.01)


# more info on this method can be found with
# ?compute.brlen()

# this function can also be used to standardise branch lengths

plot(compute.brlen(t3,1),type = "cladogram") # set ll branch lengths in the tree to 1
axisPhylo()

# this can be useful as a crude visualisation for things such as speciation rates


### adding, removing and rearranging taxa ###

# Let's read in a tree in newick string format with actual names, this time
treeA <- read.tree(text = "((Thelodonti,(Antiarchi,(Arthrodira,(((Actinistia,(Tetrapoda,Dipnoi)),Actinpoterygii),(Chondrichthyes,Acanthodii))))),Cyclostomata);")

par(mfrow=c(1,2))
plot(treeA, no.margin=T, label.offset = 0.5); tiplabels()

# it looks pretty messy...

# the ladderize() function reorganises the internal structure of the tree to neaten it up
tree <- ladderize(treeA)

plot(tree,no.margin=T,label.offset=0.5); tiplabels() # looks a bit neater now
# branches have been rearranged and node labels aren't quite so squished into the branches

# when you ladderize the the tree, it doesn't change the tip label order in the "phylo" object
tree$tip.label # now doesn't match the order they plot in

# you can neaten this up in the tree by using:
tree <- read.tree(text = write.tree(tree)) 
# this rewrites the tree into its ladderized form and reads the phylo object back in

tree$tip.label # now they're saved in they same order as they appear in the tree

# plot the  tree and you'll see that the tip numbers are now in order
par(.pardefault)
plot(tree,label.offset = 0.4); tiplabels()

# this is useful for matching up rows in data tables with tip labels in your tree


### graft a taxon in just below a particular node ###

# adding in a jawless group, Osteostraci, which is just more basal than antiarchs

# plot the node numbers so we can see what's going on
nodelabels()

# find the node which represents the most recent common ancestor of antiarchs and tetrapods
node <- getMRCA(tree,c("Antiarchi","Tetrapoda")) 
node

# now, we add the taxon in at this node
tree <- bind.tip(tree, tip.label = "Osteostraci", where = node, position = 0.1)
# "position" is the distance below the node you want to graft in the branch
# it doesn't mean much here because our tree has no branch lengths
# but if your tree has branch lengths you can specify exactly where you want your new node to be

# let's take a look at the new tree
plot(tree, label.offset = 0.1); tiplabels()


### taking taxa out ###
# for example, if we wanted to remove the thelodonts...

tree2 <- drop.tip(phy = tree, tip = "Thelodonti")
plot(tree2, label.offset = 0.4); tiplabels()

# we can also extract clades from a bigger tree
subtree <- extract.clade(phy = tree, node = getMRCA(phy = tree, tip = c("Antiarchi","Tetrapoda")))
plot(subtree, label.offset = 0.1)

# for timescaled trees, we can also choose to drop all the extinct or extant taxa  using a function in the paleotree library
t3$root.time <- max(vcv(t3)) # set our tree's root time, assuming the most recent taxon is extant

t5 <- dropExtant(t3, tol = 1.5) # tol is where we want to draw the line at "extant" - if the taxon dies out before then, it's considered extinct. If it occurs after tol, it's considered extant

t6 <- dropExtinct(t3, tol = 1.5)

par(mfrow=c(1,3)) # plot 3 images side-by-side

plot(t3, label.offset = 0.1, x.lim=c(0,5),cex = 1.5); title("all taxa") # plot the whole tree
axisPhylo() # add an axis
abline(v=t3$root.time-1.5,lty = 2) # add a vertical line at 1.5 - our extinct/extant cutoff

plot(t5, label.offset = 0.1, x.lim=c(0,5),cex = 1.5); title("extinct taxa") # plot the extinct taxa
axisPhylo() # add an axis - you'll see it ends at 1.5

plot(t6, label.offset = 0.1, x.lim=c(0,5),cex = 1.5); title("extant taxa") # plot the extant taxa
axisPhylo() # add an axis


# you can also reroot trees using a different outgroup
rr.tree <- root(tree, outgroup = "Tetrapoda")

par(.pardefault)
plot(rr.tree, label.offset = 0.1)


### additional note: fixing that typo

# it's easy to spot typos when you've only got a small tree like this
# but if you've got a big tree, spotting typos can be hard

# if you're comparing your phylogeny to a dataset/species list,
# you can compare the taxa to make sure they're all the same 

# let's set up a comparison dataset

timeData <- matrix(data = c(415,0,367,0,380,0,425,0,422,0,445,253,419,359,419,359,433,359,458,359,419,0), 
	 byrow=T,nrow = 11, ncol = 2, dimnames = list(c("Tetrapoda","Dipnoi",
	"Actinistia","Actinopterygii","Chondrichthyes","Acanthodii","Arthrodira","Antiarchi",
	"Osteostraci","Thelodonti","Cyclostomata"),c("FAD","LAD")))
timeData

setdiff(rownames(timeData),tree$tip.label) # which taxa are in the dataset but not the tree?

setdiff(tree$tip.label,rownames(timeData)) # which taxa are in the tree but not the dataset?

# find the misspelled tip label and replace it
tree$tip.label[which(tree$tip.label=="Actinpoterygii")] <- "Actinopterygii"
plot(tree,label.offset=0.1)

setdiff(rownames(timeData),tree$tip.label) # are there any mismatches now?


#######################################################
### timescaling a tree with a posteriori algorithms ###
#######################################################

# the timeData dataset we just imported can be used to estimate branch lengths in a tree
# we just need a list of taxa, and their first appearance dates (FADs)
# and last appearance dates (LADs)
timeData

# we can optionally have a vector of minimum ages for the internal nodes
# this is for an optional variable (node.mins) in the timescaling function we're about to use
# it places a constraint on the minimum divergence time of a node
# e.g. from a fossil calibration

nodeMins <- rep(0,tree$Nnode) # set a minimum age for each internal node in the tree to 0

# in this case, we know that the root node (the first internal node)  is at least early Cambrian
# so we'll say that our root node (the split between cyclostomes and gnathostomes) must be at least 530Ma

nodeMins[1] <- 530

# we won't constrain the other nodes - we'll leave their minimum age as 0

#this function call timescales our tree using tip-ages and the Minimum Branch Length (MBL) algorithm
timescaled.tree <- timePaleoPhy(tree = tree, timeData = timeData, node.mins = nodeMins,
	type = "equal", vartime = 10, plot = T, add.term = T) 
# type: specifies the algorithm used for timescaling. There are multiple options to choose from
# vartime: minimum length a branch can be.
# plot: do we want the function to plot its result? It will plot both a timescaled and non-timescaled tree
# add.term: by default, this function plots the tree tips to coincide with the estimated FADs
# if we set add.term = T, (add terminal branches) it adds on the range info so the tips correspond to the LADs

#we can now plot our timescaled tree
plot(timescaled.tree, edge.width = 2, label.offset = 0.5)
axisPhylo()

# if we want to add a fancy geological timescale, we can use:

par(oma=c(1,1,1,1)) # make sure to add some margins to your plot - the function doesn't do this automatically
geoscalePhylo(timescaled.tree, ages = timeData, tick.scale = "Period", 
	label.offset = 5, units = "Period", boxes = "Period", quat.rm=T,
	cex.ts = 0.7, cex.tip = 1, cex.age = 0.7, width = 2)

# for further information on what all these variables mean, check:
# ?geoscalePhylo()


###################################################
### matching data in your dataset to your trees ###
###################################################

# the timescaling performed above gives us an example of how we might relate data to a tree
# in that case, it helped us estimate branch lengths, but it can be used for a number of things

# what's important is that the rows in our dataset match up to the taxa in our tree

# let's build another dummy dataset
x <- fastBM(t3) # simulate the evolution of a continuous trait on the tree (returns values at tips)
x # look at the trait values

#are there any mismatches between the names of your trait and the tips of the tree?
all.equal(names(x), t3$tip.label) 

# all.equal() is a bit like setdiff, which we used above
# however, all.equal() relies on the two variables being in the same order as each other

# if we shuffle the order of the dataset...
x <- sample(x)

all.equal(names(x), t3$tip.label) # there's a bunch of mismatches!

x <- x[t3$tip.label] #reorder the data

all.equal(names(x), t3$tip.label) # fixed it!


# sometimes, we only want to use subsets of our data

# what if your tree contains fewer taxa than your trait data?
# for example, if we have a big dataset, but we only have a phylogeny of one family...

t3a <- drop.tip(t3, sample(t3$tip, size = 4)) # randomly drop four taxa from the tree
all.equal(names(x), t3a$tip.label) # now our dataset doesn't all match up

setdiff(names(x), t3a$tip.label) # which taxa are in our dataset but not in the tree?

#we can fix that
xa <- x[t3a$tip.label] # only take the elements of x which match tip labels in our tree

# do they match now?
all.equal(names(xa), t3a$tip.label) # yes :)


# what if your tree contains more taxa than your trait data?
# e.g., if we have a phylogeny of lots of species, but we only have data for ones from a certain region?

xb <- x[!(names(x) %in% names(sample(x, 4)))] #randomly drop four taxa from trait vector
all.equal(names(xb), t3$tip.label) # now we've got some mismatches

dropme <- setdiff(t3$tip.label, names(xb)) #find which taxa are in our tree but not our dataset
dropme

t3b <- drop.tip(t3, dropme) # drop those tips from the tree
all.equal(names(xb), t3b$tip.label) #now our tree and our data match up again



#######################################
### adding phenotypic data to trees ###
#######################################

# by now, we've learned how to get your tree into R
# and how to match your tree up to your dataset

# now, we're going to look at how to visualise some of this data on your tree

# we've already generated some trait values for our tree tips using fastBM()
x

# as an example, let's pretend these values represent body lengths of each taxon
# ignore the fact that some values are negative!


# we can visualise the distribution of body size evolution across the tree using:
par(.pardefault)
phenogram(t3, x, col = "cornflowerblue"); title(main = "BM trait values") 

# this plots the variation in a continuous trait across the phylogeny

# however, this isn't very easy to read!

# we can make a "heat map" along the branches of our phylogeny to visualise 
# body size changes in a slightly neater fashion

par(mfrow = c(1,2))
phenogram(t3, x, col = "cornflowerblue"); title(main = "BM trait values")
plotBranchbyTrait(ladderize(t3), x, mode = "tips") # colour edges of tree according to ancestral states reconstructed from tip values

# this shows us an evolutionary history of the trait


# if we want to also add discrete states to the tip labels, we need more data
# for example, let's say that any taxa with a BM trait value - our body length proxy - of over 1.2 are "big"
# and any taxa under 1.2 are "small"

# if we wanted to plot "big" taxa with black names and "small" ones with blue names
# we'd do something like this

# set up a vector to put our tip colours into
tip.colors <- rep("black",length(x)) # make sure it's the same size as our tree and dataset
tip.colors[which(x < 0.5)] <- "cornflowerblue" # any taxa smaller than 0.5 can be blue instead


# double check this looks about right
cbind(x,tip.colors) 

# all the taxa under 1.2 are "cornflowerblue"
# all the taxa over 1.2 are "black"

# let's plot this
par(.pardefault)
plotBranchbyTrait(ladderize(t3), x, mode = "tips",tip.color=tip.colors) 



###############################
### making your tree pretty ###
###############################

### rearranging clades ###

# start with the vertebrate tree
par(mfrow=c(1,3))
plot(treeA, edge.width = 2, label.offset = 0.1, cex=1.2, no.margin = T); nodelabels()

# ladderizing the tree orders clades according to their size
# we already did this earlier in the script with the ladderize() function:
# tree <- ladderize(tree)

# plot the ladderized tree
plot(tree, edge.width = 2, label.offset = 0.1,cex = 1.2, no.margin = T); nodelabels()

# you can use rotate() to switch clades at a node (can be useful for producing figures)
rt.tree <- rotate(tree, node = 17) # rotate node 17
plot(rt.tree, edge.width = 2, label.offset = 0.1,cex = 1.2, no.margin = T); nodelabels()


### changing colours/line types of branches ###

# this can be done with the edge.color argument

# if, for example, I wanted to highlight the branches relating to crown-group gnathostomes...

crown.gnathostomes <- c("Acanthodii","Chondrichthyes","Actinopterygii","Actinistia","Dipnoi","Tetrapoda") # list of taxa in your clade of interest
clade.edges <- which.edge(tree,crown.gnathostomes) # find which edges are relevant to this clade
clade.edges # these are the branch numbers we want to colour

edge.cols <- rep("black",nrow(tree$edge)) # set up a vector of edge (branch) colours with all the colours defaulted to black
edge.cols[clade.edges] <- "cornflowerblue" # turn the ones which are relevant to our clade of interest blue
edge.cols # we've now got a vector of the colours each branch should be

par(.pardefault)
plot(tree, edge.width = 3, edge.color = edge.cols, label.offset = 0.2) # plot the tree with our coloured branches

# you can also change the line type in a similar manner using the edge.lty argument

edge.lty <- rep(1, nrow(tree$edge)) # vector of clade line types (1 is a standard solid line)
edge.lty[clade.edges] <- 3 # change the relevant ones to line type 3 (dashed)
plot(tree, edge.width = 3, edge.lty = edge.lty, edge.color = edge.cols, label.offset = 0.2)

# now you've got blue dotted lines


### adding coloured shapes to tree tips ###

par(mfrow = c(1,2))
plot(tree, label.offset = 0.5, no.margin = T) # plot a boring tree

rainbow.tips <- rainbow(length(tree$tip.label)) # make a vector of colours
# in this example, they're just random colours
# but you can add colours to represent specific things
# e.g. threat level, plumage colour... whatever you want, really

names(rainbow.tips) <- tree$tip.label # match the colours to the tip names
rainbow.tips # take a look at the colours

# plot them on the tree
tiplabels(pch = 23, cex = 2, bg = rainbow.tips)


# if we want to change a tip colour...
rainbow.tips["Tetrapoda"] <- "cornflowerblue"
plot(tree, label.offset = 0.5, no.margin = T) # plot the tree again
tiplabels(pch = 23, cex = 2, bg = rainbow.tips) # add the tips - tetrapods are now blue



### higlighting a clade ###

# this is a slightly long-winded way of doing this on such a small tree
# however, it works well for drawing multiple boxes quickly
# because you can do everything with for-loops

# the package "ggtree" - an extension of ggplot2 - has ways of doing this sort of thing
# but it's still in early stages of development and is quite buggy
# so this way is safer!
# (also, lots of people don't like ggplot)

# first, we neeed some branch lengths so we can calculate where we want our box
# we're going to use the vertebrate tree as an example

tree <- compute.brlen(tree, method = "Grafen", power = 0.5)

par(.pardefault)
plot(tree,label.offset = 0.01)
axisPhylo(); nodelabels()

# to define a box, we need four values:
# xmin, xmax, ymin and ymax
# these tell us where the edges go

# we can calculate these from our tree and use them to draw boxes
# we use the axis for our x-values
# and our tip labels for our y-values
# the bottom tip has a y coordinate of 1 - each tip label increments by 1
points(1,1,pch=19)
points(1,2,col="red",pch=19)


# In this example, I want to draw a box around my crown group gnathostomes
# everything from Acanthodii up to Tetrapoda 

# which node defines the most recent common ancestor of all the crown gnathostomes?
node <- getMRCA(tree,c("Tetrapoda","Acanthodii"))
node

# how many taxa are there in this clade?
clade.size <- Ntip(extract.clade(tree,node)) # extract the subtree from the node of interest.
clade.size #  How many taxa are in this clade?

# because this clade is at the bottom of the tree, our box can start at 0.5
ybottom <- 0.6 # if we make it 0.6, it leaves a bit of a gap above the axis and looks a bit neater

# we can calculate how high the box needs to be from the number of taxa in it
ytop <- ybottom + clade.size - 0.2 # take a bit off again
# if we make the values 0.5, rather than 0.4 and 0.6, there's no gap between adjacent boxes
# and it looks messy

# the leftmost part of the box should be just below node 17
# dist.nodes() tells you how far any two nodes apart are

# remember, node (Ntip(tree) + 1) is our root node number

# find the distance from the root to node 17
# take a little bit off so the box starts just below the node

xleft <- dist.nodes(tree)[Ntip(tree)+1,node]-0.04

# the right side of the box should be level with our tips
# find how far the tips (nodes 1 - 11) are from the root
xright <- dist.nodes(tree)[Ntip(tree)+1,1]

# we can also do this with 
max(vcv(tree))


# we've got our coordinates - now we can plot our box!
rect(xleft, ybottom, xright, ytop, col = "cornflowerblue", border=NA) 

# we now have a box. Unfortunately, we've covered up our tree

# while we can set the box to appear transparent by changing the alpha value,
# this still means our tree branches will no longer be black
# they'll go a bit blue because of the colour overlay

# instead, we will re-plot our tree over the box

par(new=TRUE) # this tells R that we're still drawwing on the same plot
# if we don't do this, it will clear the window before plotting the tree again

plot(tree,label.offset = 0.01); nodelabels() # plot the tree again, over the top of our box this time!

# you  might have noticed our tip labels look aa bit funky because they've been drawn on twice

# we can overcome this by telling R to plot them in white the first time round
plot(tree,label.offset = 0.01,tip.color = "white") # tree with invisible tip labels
rect(xleft, ybottom, xright, ytop, col="cornflowerblue",border = NA) # box which covers up the tree
par(new=TRUE) # keep plotting in the same window
plot(tree,label.offset = 0.01) # add the tree back on the top

### if you're drawing a fan-type tree, you can use plot.sector from the Circlize package
### to plot your background blocks

### this requires r and theta values rather than x and y


### adding images to tree ###

### this borrows heavily from a script at:
# https://gist.github.com/scrogster/7fc5b7597b63585a00b6

install.packages(c("png","RCurl"), repos="https://cran.csiro.au/") # we need these two packages. Download them and read them in
library(png)
library(RCurl)


# we're going to use this little function taken from the URL cited above to add our pictures

logoing_func <- function(logo, x, y, size) {
	dims <- dim(logo)[1:2] #number of x-y pixels for the logo (aspect ratio)
	AR <- dims[1]/dims[2]
	par(usr=c(0, 1, 0, 1))
	rasterImage(logo, x-(size/2), y-(AR*size/2), x+(size/2), y+(AR*size/2), interpolate=TRUE)
}


# I want to put a picture of a Dunkleosteus on my phylogeny

# find the URL of a phylopic silhouette I like
arthrodira.url <- "http://phylopic.org/assets/images/submissions/8c479a10-d9c1-4324-86cf-60e04e4008e1.512.png"

arthrodira.logo <-  readPNG(getURLContent(arthrodira.url)) # and get the file itself from the URL


par(mar = c(1,1,1,3), oma = c(1,1,1,1), xpd = T) # make sure I have enough margin space to add in my picture - enable xpd to draw outside plot area
plot(tree) # plot the tree
logoing_func(arthrodira.logo, x=0.95, y=0.59, size=0.15) # add in my picture
# x and y are the the proportion of how far into your plot you want the image to be
# size changes the scale of the picture

# yay, we've got a picture on our plot!


# now lets add some pictures for the rest of the groups in with some URLs I found earlier
cyc.logo <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/5a5d9328-7438-40f8-b538-13edeedc0f0c.512.png"))
the.logo <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/68a14c4d-30ae-45be-b7a3-1b1c837fc7c6.512.png"))
ost.logo <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/279a1d82-1283-4d0a-9d83-4139de8cc416.512.png"))
ant.logo <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/2d73e1fe-2e14-46bf-8d44-2ef9414167be.512.png"))
aca.logo <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/efd34171-3ad3-4de4-a044-490ae3904266.512.png"))
cho.logo <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/7f26031d-169a-4438-8f17-7d5e418f5964.512.png"))
act.logo <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/884d77c5-932b-4331-a3c7-e9df01d805e6.original.png"))
coe.logo <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/ff131fa9-0ac4-4e6d-8df9-c44226ea4e7c.512.png"))
dip.logo <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/de33e9bd-643d-4a68-a3c4-c3206f9b508e.512.png"))
tet.logo <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/f5024c43-2614-44d9-b66f-03f0269996ef.512.png"))


# and add them to our plot
logoing_func(cyc.logo, x=0.95, y=0.96, size=0.15)
logoing_func(the.logo, x=0.95, y=0.87, size=0.15)
logoing_func(ost.logo, x=0.95, y=0.78, size=0.15)
logoing_func(ant.logo, x=0.95, y=0.69, size=0.15)
logoing_func(aca.logo, x=0.95, y=0.50, size=0.15)
logoing_func(cho.logo, x=0.95, y=0.41, size=0.15)
logoing_func(act.logo, x=0.95, y=0.32, size=0.15)
logoing_func(coe.logo, x=0.95, y=0.22, size=0.15)
logoing_func(dip.logo, x=0.95, y=0.12, size=0.15)
logoing_func(tet.logo, x=0.95, y=0.03, size=0.15)

# Et voila! An illustrated phylogeny of early vertebrates! :)

