library(maps)
data(world.cities)
samples <- read.table("../medicago/data/sample_info.tsv",sep="\t",header=TRUE)
samples$country <- samples$Country.of.Origin
levels(samples$country)[ levels(samples$country) ==  "Greece, Crete" ] <- "Greece_Crete"
levels(samples$country)[ levels(samples$country) ==  "France, Corsica" ] <- "France_Corsica"
levels(samples$country) <- gsub(",.*","",levels(samples$country))
levels(samples$country)[ levels(samples$country)=='' ] <- "blank"

countries <- setdiff( sort( levels( samples$country ) ), c("Bulgaria", "China", "Russia", "Uzbekistan") )
country.cols <- rainbow(length(countries))

other.cities <- c(
                  "France_Corsica" = "Ajaccio",
                  "Greece_Crete" = "Malia",
                  "U.S." = "Los Angeles"
                  )
cities <- subset( do.call( rbind, tapply( 1:nrow(world.cities), world.cities$country.etc,
                                 function (kk) { 
                                     thisone <- which.max( world.cities$pop[kk] )
                                     world.cities[kk[thisone],]
                                 } ) ), tolower(country.etc) %in% tolower(countries) )
for (k in seq_along(other.cities)) {
    thisone <- subset( world.cities, name == other.cities[k] )
    thisone$country.etc <- names(other.cities[k])
    cities <- rbind( cities, thisone[nrow(thisone),] )
}
cities <- cities[ match(countries,cities$country.etc), ]

lims <- structure(c(29.025, 62.085, -12.099, 41.649), .Dim = c(2L, 2L), .Dimnames = list(NULL, c("lat", "long")))

map("world", xlim=lims[,'long'], ylim=lims[,'lat'] )
points( cities[,"long"], cities[,"lat"], col=country.cols, pch=seq_along(countries), cex=3, lwd=2 )
text( cities[,"long"], cities[,"lat"], labels=cities$country.etc, col=country.cols, 
    pos=3, adj=1 )


