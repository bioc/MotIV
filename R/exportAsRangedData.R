
#####RANGED DATA#####
exportAsRangedData <- function (x, y, correction=TRUE)
{
        .Deprecated("exportAsGRanges")
        as(exportAsGRanges(x, y, correction=correction), "RangedData")
}
