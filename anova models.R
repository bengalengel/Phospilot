# To begin with, we will use the example I had in class. There are three schools, with two students nested in each school.
scores <- c(25, 29, 14, 11, 11, 6, 22, 18, 17, 20, 5, 2)
school <- factor(c("A", "A", "A", "A", "B", "B", "B", "B", "C", "C", "C", "C"))
teacher <- factor(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6))
teacher2 <- factor(c(1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2))  # This is the way the data is coded for problems in the book (student!)

boxplot(scores ~ school)
boxplot(scores ~ teacher)
boxplot(scores ~ school:teacher2)

library(lattice)
dotplot(scores ~ teacher2 | school)

# What if we ignore the fact that the design is nested? 
anova(lm(scores ~ school))
anova(lm(scores ~ school + teacher))
anova(lm(scores ~ school * teacher))
anova(lm(scores ~ school + teacher2))
anova(lm(scores ~ school * teacher2))

# What if we do it correctly?
res1 <- lm(scores ~ school + school/teacher)
anova(res1)

res2 <- lm(scores ~ school + school/teacher2)
anova(res2)

# Which schools are different?
TukeyHSD(aov(res1), "school")
TukeyHSD(aov(res2), "school")

TukeyHSD(aov(res1), "school:teacher")  #Yikes
TukeyHSD(aov(res2), "school:teacher2")

contrast(res2, list(school = "A", teacher = "1"), list(school = "A", teacher = "2"))  #Bummer







