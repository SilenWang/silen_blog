---
title: 用R绘制地图的笔记
date: 2018-07-29 00:38:09
tags: [r]
categories: Script
---

研究生期间曾有作地图展示结果的需要, 使用R绘制过地图, 将当时的代码记录如下

<!-- more -->

#流行病学数据地图展示


## 代码来源

脚本代码主要改写自统计之都文章:
[R绘制中国地图，并展示流行病学数据](http://cos.name/2014/08/r-maps-for-china/)


## 代码本体(已包含备注)

```r
####loading packages 载入需要的程序包####
library(maptools) #用于读取并操作地图数据
library(ggplot2) #用于绘制地图
####set working directory####
getwd()
setwd("/home/silen/R _Data/Map") #设定当前工作目录,将目录设置为地图数据所在位置即可
####loading map data(shape data)####
map <- readShapePoly("countries_shp/countries.shp") #读取*.shp格式的地图数据,注意*.shp数据须和其他文件放置在一起如*.shx才可读取,否则会出现错误
####check loaded data 调试部分,用于查看已读取地图数据情况####
names(map)
str(map$NAME)
table(map$EU)

####Normal method, using geom_polygon / geom_path 使用基本的polygon方法与path方法进行地图绘制(暂不清楚如何展示流行病数据)####
DB <- fortify(map) #重要语句,通过此步骤将读入地图数据转化为ggplot可以识别的格式
p <- ggplot(DB, aes(x=long,y=lat, group = group)) #creat a plot 图形初始化,确定使用数据库并将经纬度映射到x,y轴;group=group保证绘制图形无错乱
# World map, using geom_path instead of geom_polygon(来自ggplot帮助文档,暂不明为何)
World <- p + geom_polygon(color="white", fill = "grey") + #使用polygon(多边形)绘制地图,并设置填充色为灰,线条色为白
  theme(panel.grid.major=element_blank(),#theme()内都为样式效果设定,主要为隐藏x.y轴以及相应文字,背景去除,网格去除,
        panel.background=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank())+
  ggtitle("World")+
  scale_y_continuous(breaks = (-2:2) * 30) +#忘了干嘛的
  scale_x_continuous(breaks = (-4:4) * 45) +
  #coord_map("orthographic") #极视图
  #coord_map("gilbert") #ball 球形视图
  coord_map("lagrange") #flat
#coord_map()内为特殊坐标轴映射,具体见ggplot包说明

####Epi method using geom_map 使用专有的geom_map方法绘制地图,可展示流行病学数据(热图)####
#load data and get prepared for plotting 生成需要的地图数据,载入部分在开头
DB <- fortify(map, region = "NAME")#transform map shape data to object can be read by ggplot
#DB <- transform(DB, id = iconv(id, from = 'GBK'), code = iconv(code, from = 'GBK'))#transform coding from GBK
head(DB) #查看数据情况
names(DB)[c(1, 2)] = c("x", "y")#change head of data, for purpose of using expand_limits() 改变数据表中变量的名称,此处为geom_map方法所需
#prepare Epidemiology data -> crude death rate 准备好流行病学数据,注意对应问题
mor <- read.csv("mapRate.csv",head =T, sep =",") #从*.csv文件中读入流行病学数据
epi <- data.frame(id = unique(sort(DB$id)))# 使用名为epi的数据框保存流行病学数据,从DB(地图数据)中获得所有区域(或地点)的名称(对应用)
epi <- merge(epi, mor, by.x= "id", all.x = T)# 以id变量为对应标准,合并区域名称及相应流行病数据
#write.csv(epi, "1234.csv") #调试用
#epi
p1 <- ggplot(epi, fill = "#595959") + 
              geom_map(aes(map_id = id, fill = rate ),#使用geom_map绘制地图并展示数据
                       colour = "white",
                           map = DB) +
                    expand_limits(DB) +#重要语句,无此语句无法绘制出地图,语句意义待查
                    #changing theme, make backgroud blank/transparent
                    theme(panel.grid.major = element_blank(),#主题设置,大致为隐藏x.y轴并去掉背景、网格,使背景透明
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
                    #ggtitle("World") + #title of plot
                    scale_fill_gradient(high = "#F70909", low = "#E99799") + #设定热图颜色
                    scale_y_continuous(breaks = (-2:2) * 30) +#意义不明
                    scale_x_continuous(breaks = (-4:4) * 45) +
                    #coord_map("orthographic")
                    #coord_map("gilbert") #ball
                    coord_map("lagrange") #flat#坐标设定
####print plot 输出图像,只有在需要输出透明背景图像时使用####
png('world.png',width=600,height=600,units="px",bg = "transparent")#设定名称、大小、大小单位以及设背景为透明
print (p1)#输出p1
dev.off()#表示完成输出

```