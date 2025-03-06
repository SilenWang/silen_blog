---
title: Notes on Drawing Maps with R
date: 2018-07-29 00:38:09
tags: [R]
categories: Script
---

During my graduate studies, I had the need to display results using maps and used R to create them. I recorded the code from that time below.

<!-- more -->

# Epidemiological Data Map Visualization

## Code Source

The script code is mainly adapted from an article on Statistical Analysis:
[Drawing Maps of China and Displaying Epidemiological Data with R](http://cos.name/2014/08/r-maps-for-china/)

## Code Body (Includes Comments)

```r
####loading packages Load necessary libraries####
library(maptools) # Used for reading and manipulating map data
library(ggplot2) # Used for plotting maps
####set working directory####
getwd()
setwd("/home/silen/R _Data/Map") # Set the current working directory to the location of the map data
####loading map data(shape data)####
map <- readShapePoly("countries_shp/countries.shp") # Read *.shp format map data, note that *.shp files must be placed together with other files like *.shx to be readable; otherwise, errors will occur.
####check loaded data Debugging part, used to check the loaded map data####
names(map)
str(map$NAME)
table(map$EU)

####Normal method, using geom_polygon / geom_path Using basic polygon and path methods for map plotting (unclear how to display epidemiological data)####
DB <- fortify(map) # Important statement: Through this step, the imported map data is converted into a format that ggplot can recognize.
p <- ggplot(DB, aes(x=long,y=lat, group = group)) # Create a plot, initialize the graph, determine the use of the database and map longitude and latitude to the x and y axes; group=group ensures that the plotted graph is not chaotic.
# World map, using geom_path instead of geom_polygon (from ggplot help documentation, unclear why)
World <- p + geom_polygon(color="white", fill = "grey") + # Use polygon (polygon) to draw the map, set the fill color to grey and line color to white
  theme(panel.grid.major=element_blank(),# Theme() within are all style effect settings, mainly hiding x.y axis and corresponding text, background removal, grid removal,
        panel.background=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank())+
  ggtitle("World")+
  scale_y_continuous(breaks = (-2:2) * 30) +# Forgot what it does
  scale_x_continuous(breaks = (-4:4) * 45) +
  #coord_map("orthographic") # Orthographic view
  #coord_map("gilbert") # Ball shape view
  coord_map("lagrange") # Flat
# coord_map() within are special coordinate axis mapping, see ggplot package documentation for details

####Epi method using geom_map Using the dedicated geom_map method to draw maps and display epidemiological data (heat map)####
# Load data and prepare for plotting Generate the required map data, load part at the beginning
DB <- fortify(map, region = "NAME")# Transform map shape data into an object that ggplot can read
# DB <- transform(DB, id = iconv(id, from = 'GBK'), code = iconv(code, from = 'GBK'))# Transform coding from GBK
head(DB) # Check the data situation
names(DB)[c(1, 2)] = c("x", "y")# Change the header of the data, for the purpose of using expand_limits() Change the names of variables in the data table as required by geom_map method
# Prepare epidemiological data -> crude death rate Prepare the epidemiological data, note the correspondence issue
mor <- read.csv("mapRate.csv",head =T, sep =",") # Read epidemiological data from *.csv file
epi <- data.frame(id = unique(sort(DB$id)))# Use epi data frame to save epidemiological data, obtain all region (or location) names from DB (map data)
epi <- merge(epi, mor, by.x= "id", all.x = T)# Merge region names and corresponding epidemiological data based on id variable
# write.csv(epi, "1234.csv") # Debugging use
# epi
p1 <- ggplot(epi, fill = "#595959") + 
              geom_map(aes(map_id = id, fill = rate ),# Use geom_map to draw the map and display data
                       colour = "white",
                           map = DB) +
                    expand_limits(DB) +# Important statement: Without this statement, the map cannot be drawn. The meaning of the statement is pending investigation
                    # Changing theme, make background blank/transparent
                    theme(panel.grid.major = element_blank(),# Theme settings, mainly hiding x.y axis and removing background, grid to make background transparent
                          panel.background = element_blank(),
                          plot.background = element_rect(fill = "transparent",colour = NA),
                          legend.background = element_rect(fill = "transparent",colour = NA),
                          legend.title = element_blank(),
                          axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.x = element_blank(),
                          axis.text.y = element_blank(),
                          axis.line = element_blank(),
                          axis.ticks = element_blank())+
                    # ggtitle("World") + # Title of plot
                    scale_fill_gradient(high = "#F70909", low = "#E99799") + # Set heat map color
                    scale_y_continuous(breaks = (-2:2) * 30) +# Meaning unclear
                    scale_x_continuous(breaks = (-4:4) * 45) +
                    # coord_map("orthographic")
                    # coord_map("gilbert") # Ball
                    coord_map("lagrange") # Flat# Coordinate setting
####print plot Output image, only used when a transparent background image is needed####
png('world.png',width=600,height=600,units="px",bg = "transparent")# Set name, size, size unit and set background to transparent
print (p1)# Output p1
dev.off()# Indicates completion of output

