# Cleaning up
> objects()
> rm(list=objects())
> prodir=getwd()
> dir.create("Lab05")
> setwd("Lab05")
> url="https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/Lab04.R"
> download.file(url,"Lab04.R")
> file.edit("Lab04.R")
