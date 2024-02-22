any(config$total_hosts > config$total_populations)
any(config$total_hosts > config$total_populations, na.rm=T)
any(config$total_exposed > config$total_populations)
any(config$total_exposed > config$total_populations, na.rm=T)
any(config$total_infecteds > config$total_populations)
any(config$total_infecteds > config$total_populations, na.rm=T)
# "Total hosts"
m2 <- matrix(c(2,2,2,2), ncol=2)
# "Total infecteds"
m3 <- matrix(c(3,3,3,3), ncol=2)
# "Total exposed"
m4 <- matrix(c(4,4,4,4), ncol = 2)
# "Total populations"
m5 <- matrix(c(5,5,5,5), ncol=2)
# TOtal hosts with an NA
m2na <- matrix(c(2,2,2,NA), ncol=2)
# Total hosts with a value greater than total populations
m26 <- matrix(c(2,6,2,2), ncol=2)
# Total hosts with more values than total population
m2x6 <- matrix(c(2,2,2,2,2,2), ncol=2)

if( any(m2 > m5) || any(m3 > m5) || any(m4 > m5)) { print("none greater than m5")}
# returns ""

if( any(m2 > m5) || any(m3 > m5) || any(m4 > m5)) { print("none greater than m5")} else {print("something else")}
# returns "something else"

if( any(m2na > m5) || any(m3 > m5) || any(m4 > m5)) { print("none greater than m5")} else {print("something else")}
# returns 
'Error in if (any(m2na > m5) || any(m3 > m5) || any(m4 > m5)) { : 
    missing value where TRUE/FALSE needed'

if( any(m26 > m5) || any(m3 > m5) || any(m4 > m5)) { print("none greater than m5")} else {print("something else")}
# returns "none greater than m5"

if( any(m2x6 > m5) || any(m3 > m5) || any(m4 > m5)) { print("none greater than m5")} else {print("something else")}
# returns
'Error in m2x6 > m5 : non-conformable arrays'

while( any(m2 > m5) || any(m3 > m5) || any(m4 > m5)) { print("none greater than m5")}

while( any(m2 > m5) || any(m3 > m5) || any(m26 > m5)) { print(m5)}
# Keeps printing m5
