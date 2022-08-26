library(AnnotationHub)
package = "AnnotationHub"

oldcache = path.expand(rappdirs::user_cache_dir(appname=package))
setAnnotationHubOption("CACHE", oldcache)
ah = AnnotationHub(localHub=TRUE)
## removes old location and all resources
removeCache(ah, ask=FALSE)

## create the new default caching location
newcache = tools::R_user_dir(package, which="cache")
setAnnotationHubOption("CACHE", newcache)
ah = AnnotationHub()
