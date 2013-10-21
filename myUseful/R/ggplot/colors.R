df <- read.table(header=T, con <- textConnection('
con yval
A	2
B	2.5
C   1.6
'))
close(con)
df2 <- read.table(header=T, con <- textConnection('
 cond1 cond2 yval
    A      I 2
    A      J 2.5
    A      K 1.6
    B      I 2.2
    B      J 2.4
    B      K 1.2
    C      I 1.7
    C      J 2.3
    C      K 1.9
'))
close(con)

ggplot(df,aes(x=con,y=yval)) + geom_bar()

## using hexadecimal
ggplot(df, aes(x=con,y=yval)) + geom_bar(colour = "#DCFA43")

ggplot(df, aes(x=con,y=yval)) + geom_bar(fill = "#DCFA43",colour= "blue")

ggplot(df,aes(x=con,y=yval)) + geom_line(aes(group=1)) + geom_point(size=3)

ggplot(df,aes(x=con,y=yval)) + geom_line(aes(group=1),colour="#000099") + geom_point(size=3,colour="#CC0000")


ggplot(df,aes(x=con,y=yval,fill=con)) + geom_bar()

ggplot(df2,aes(x=cond1,y=yval)) + geom_bar(aes(fill=cond2),colour="blue", position=position_dodge())

ggplot(df2,aes(x=cond1,y=yval)) + geom_line(aes(colour=cond2,group=cond2)) + geom_point(aes(colour=cond2),size=3)

ggplot(df2,aes(x=cond1,y=yval,colour=cond2)) + geom_line(aes(group=cond2)) + geom_point(size=3)
### using RColorBrewer
ggplot(df,aes(x=con,y=yval,fill=con)) + geom_bar() + scale_fill_brewer(palette="Set1")

ggplot(df2,aes(x=cond1,y=yval,colour=cond2)) + geom_line(aes(group=cond2)) + geom_point(size=3) + scale_colour_brewer()


