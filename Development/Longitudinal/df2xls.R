# Converting body fat data to an Excel file

# Load data:
load("C:/usr/rick/doc/Committees/PIN/PIN Director/Courses/Stats/TAD/TAD Code/Development/Longitudinal/bodyfat.RData")

# get library
install.packages('writexl')
library(writexl)

# write to file
write_xlsx(bodyfat, "C:/usr/rick/doc/Committees/PIN/PIN Director/Courses/Stats/TAD/TAD Code/Development/Longitudinal/bodyfat.xlsx")
