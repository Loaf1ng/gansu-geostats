# 甘肃气候与县域分析图
invisible(Sys.setlocale("LC_CTYPE", ".UTF-8"))
suppressPackageStartupMessages({library(grid); library(sf); library(terra)})
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
root <- dirname(dirname(normalizePath(script)))
work <- file.path(root, "run")
out <- file.path(work, "figures")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
tables <- readRDS(file.path(root, "results", "figures.rds"))
read_result <- function(name) read.csv(file.path(root, "results", name), fileEncoding="UTF-8", stringsAsFactors=FALSE)
stations <- read_result("stations.csv")
counties_all <- read_result("counties.csv")
included <- subset(counties_all, included == 1)

# 地图
local({
sta <- stations
edges <- tables$edges
models <- c("precip_OK_log_spherical_T0_N16", "precip_CK_log_spherical_T0_N16",
            "evap_OK_log_gaussian_T0_N16", "evap_CK_log_exponential_T0_N16")
read_grid <- function(i, field) {
  rast(file.path(work, "rasters", paste0(models[i], "_", field, ".tif")))
}
pred <- lapply(seq_along(models), read_grid, field = "prediction")
se <- lapply(c(1, 3), read_grid, field = "se")
map_crs <- st_crs(pred[[1]])
province <- st_transform(st_read(file.path(work, "source/province/Gansu.shp"), quiet = TRUE), map_crs)
sites <- st_as_sf(sta, coords = c("x_m","y_m"), crs = map_crs)
stopifnot(nrow(sites) == 53L, all(is.finite(sta$x_m)), all(is.finite(sta$y_m)))
for (r in c(pred,se)) {
  stopifnot(compareGeom(r, pred[[1]], stopOnError = FALSE),
            sum(!is.na(values(r))) == 16356L)
}
poly <- st_collection_extract(st_voronoi(st_combine(st_geometry(sites)),
                            envelope = st_as_sfc(st_bbox(province))), "POLYGON")
site_id <- st_intersects(poly, sites)
stopifnot(all(lengths(site_id) == 1L), length(unique(unlist(site_id))) == 53L)
vor <- suppressWarnings(st_intersection(st_sf(evap = sta$evap[unlist(site_id)], geometry = poly),
                                        st_geometry(province)))
counties <- st_transform(st_read(file.path(work, "source/topic/Gansu-County.shp"), quiet = TRUE), map_crs)
obs <- included
ix <- match(seq_len(nrow(counties))-1L, obs$ZoneFID)
stopifnot(nrow(counties) == 89, sum(!is.na(ix)) == 88L, !anyDuplicated(obs$ZoneFID),
          all(as.character(counties$ADCODE99[!is.na(ix)]) == as.character(obs$ADCODE99[ix[!is.na(ix)]])))
bb <- unname(st_bbox(province))
ink <- "#263337"; muted <- "#536568"; edge <- "#687977"
rain <- c("#E4F1EF","#B6D9D8","#79B5BD","#418C9D","#1A6179","#123E58")
evap <- c("#FAEFE0","#F0D4AC","#E5B37F","#CF8C5B","#AA6446","#794335")
grass <- c("#F1F1DC","#D4DDB4","#A8C599","#75A78A","#478674","#245D55")
uncert <- c("#F1EFEB","#D9D8E3","#AAAAC5","#7B82AB","#525E8F","#333E6C")
resid <- c("#175973","#438DA4","#9AC7D0","#F3F1E9","#E7B99A","#C27B53","#934D37")
sans <- "Microsoft YaHei"
W <- 166
H <- 88

# Page coordinates are millimetres measured down from the upper-left corner.
start_fig <- function(n,h) {
  cairo_pdf(file.path(out, sprintf("fig%02d.pdf", n)), width = W/25.4,
            height = h/25.4, family = sans, pointsize = 10, bg = "white",
            fallback_resolution = 600)
  par(mar = rep(0,4), xaxs = "i", yaxs = "i", lend = "round", ljoin = "round",
      family = sans, fg = ink, xpd = FALSE)
  plot.new()
  plot.window(c(0,W), c(h,0), xaxs="i", yaxs="i")
}
label <- function(x,y,s,sz=7.6,adj=0,col=ink,bold=FALSE) {
  text(x,y,s,cex=sz/10,adj=c(adj,.5),col=col,font=if(bold) 2 else 1)
}
panel <- function(x,y,k,title) {
  label(x,y,k,9.2,bold=TRUE)
  label(x+5.4,y,title,8.6)
}
locator <- function(x,y,w,h) {
  s <- min(w/(bb[3]-bb[1]),h/(bb[4]-bb[2]))
  cx <- x+w/2; cy <- y+h/2
  list(x=function(v) cx+(v-(bb[1]+bb[3])/2)*s,
       y=function(v) cy-(v-(bb[2]+bb[4])/2)*s,
       east=function(v) (v-cx)/s+(bb[1]+bb[3])/2,
       north=function(v) (cy-v)/s+(bb[2]+bb[4])/2,
       s=s,left=x,top=y,width=w,height=h)
}
map_pair <- list(locator(3,10.8,77,62.2), locator(86,10.8,77,62.2))
shade <- function(x,pal,lim,na_col="#D9DADA") {
  lut <- colorRampPalette(pal,space="Lab")(256)
  q <- pmin(1,pmax(0,(x-lim[1])/diff(lim)))
  ifelse(is.na(q),na_col,lut[1+round(255*q)])
}
draw_polygons <- function(g, m, fill, stroke=edge, lwd=.42) {
  geoms <- st_geometry(g)
  fill <- rep_len(fill,length(geoms))
  for (i in seq_along(geoms)) {
    polys <- geoms[[i]]
    if (inherits(polys,"POLYGON")) polys <- list(polys)
    for (p in polys) {
      rings <- do.call(rbind,lapply(p,function(z) rbind(z[,1:2,drop=FALSE],NA_real_)))
      rings <- rings[-nrow(rings),,drop=FALSE]
      polypath(m$x(rings[,1]),m$y(rings[,2]),col=fill[i],border=stroke,lwd=lwd,rule="evenodd")
    }
  }
}
draw_raster <- function(r,m,pal,lim) {
  mat <- as.matrix(r,wide=TRUE)
  # Flip row order as the page's y axis is reversed for millimetre layout.
  cmat <- matrix(shade(as.vector(mat),pal,lim,"#FFFFFF00"),nrow=nrow(mat))
  rasterImage(as.raster(cmat), m$x(xmin(r)), m$y(ymin(r)), m$x(xmax(r)), m$y(ymax(r)),
              interpolate=FALSE)
}
circle <- function(x,y,r,fill,stroke=NA,lwd=.5) {
  # MS Word imposes a minimum stroke weight on imported SVGs. Concentric
  # filled discs keep thin station halos at their actual printed size.
  if (!is.na(stroke)) {
    half_width <- lwd * 25.4 / 96 / 2
    symbols(x,y,circles=r+half_width,inches=FALSE,add=TRUE,bg=stroke,fg=NA)
    r <- pmax(r-half_width,0.05)
  }
  symbols(x,y,circles=r,inches=FALSE,add=TRUE,bg=fill,fg=NA)
}
site_marks <- function(m) {
  circle(m$x(sta$x_m),m$y(sta$y_m),rep(.34,nrow(sta)),"#526674","#FFFFFF",.30)
}
scale_bar <- function(m,x,y) {
  # Correct the projected 200 km bar for the local ground scale at its midpoint.
  rx <- m$east(x); ry <- m$north(y)
  ends <- st_sfc(st_point(c(rx,ry)),st_point(c(rx+2e5,ry)),crs=map_crs)
  ground <- as.numeric(st_distance(st_transform(ends,4326))[1,2])
  len <- 2e5*m$s*(2e5/ground)
  segments(x,y,x+len,y,col=muted,lwd=.68)
  segments(c(x,x+len),y-.65,c(x,x+len),y+.25,col=muted,lwd=.68)
  label(x+len/2,y-2.6,"200 km",6.7,.5,col=muted)
}
true_north <- function(m,x,y) {
  # The projected meridian is inclined across western Gansu.
  base <- st_sfc(st_point(c(m$east(x),m$north(y))),crs=map_crs)
  ll <- st_coordinates(st_transform(base,4326))[1,]
  upper <- st_transform(st_sfc(st_point(ll+c(0,.1)),crs=4326),map_crs)
  tip <- st_coordinates(upper)[1,]
  ux <- m$x(tip[1])-x; uy <- m$y(tip[2])-y
  length <- sqrt(ux^2+uy^2)
  ux <- ux/length; uy <- uy/length
  ex <- x+4.2*ux; ey <- y+4.2*uy
  segments(x,y,ex,ey,lwd=.68,col=muted)
  polygon(c(ex,ex-1.6*ux+.56*uy,ex-1.6*ux-.56*uy),
          c(ey,ey-1.6*uy-.56*ux,ey-1.6*uy+.56*ux),col=muted,border=NA)
  label(ex+.2,ey-2.0,"N",6.7,.5,col=muted)
}
map_bar <- function(x,y,w,pal,limits,ticks,title,suffix="") {
  label(x,y,title,7.4)
  gradient <- matrix(colorRampPalette(pal,space="Lab")(256),nrow=1)
  rasterImage(as.raster(gradient),x,y+4.7,x+w,y+2.45,interpolate=TRUE)
  tx <- x+w*(ticks-limits[1])/diff(limits)
  segments(tx,y+4.7,tx,y+5.5,col=muted,lwd=.44)
  label(tx,y+7.75,paste0(ticks,suffix),6.7,.5,col=muted)
}
pred_limits <- list(c(0,800),c(1000,3600))
p_ticks <- c(0,200,400,600,800); e_ticks <- seq(1000,3500,500)

start_fig(2,100)
panel(3,4.2,"a","台站降水与高程");panel(86,4.2,"b","泰森多边形蒸发分区")
a <- map_pair[[1]]; b <- map_pair[[2]]
draw_polygons(province,a,"#F6F7F3")
radius_h <- function(h) .63*sqrt(h/1000)
ord <- order(sta$elev,decreasing=TRUE)
circle(a$x(sta$x_m[ord]),a$y(sta$y_m[ord]),radius_h(sta$elev[ord]),
       shade(sta$precip[ord],rain,pred_limits[[1]]),"#425E66",.45)
draw_polygons(vor,b,shade(vor$evap,evap,pred_limits[[2]]),stroke="#EEE9DF",lwd=.38)
draw_polygons(province,b,NA,lwd=.52)
site_marks(b)
# Two fixed text lanes in the whitespace east of the northern part of the province.
for (i in seq_along(c("马鬃山","梧桐沟"))) {
  name <- c("马鬃山","梧桐沟")[i]
  row <- sta[sta$name==name,,drop=FALSE]
  x <- b$x(row$x_m); y <- b$y(row$y_m); ly <- c(16.7,23.0)[i]
  lines(c(x,x+3.8,135.8),c(y,ly,ly),col="#675C51",lwd=.48)
  circle(x,y,.66,"white","#795540",.65)
  label(137,ly,name,7.2)
}
scale_bar(a,7.3,65.4)
true_north(a,26.4,65.4)
map_bar(6.6,79.3,62,rain,pred_limits[[1]],p_ticks,"降水量（mm）")
map_bar(89.6,79.3,62,evap,pred_limits[[2]],e_ticks,"蒸发量（mm）")
label(6.6,94.5,"高程（m）",7.3)
for (i in 1:3) {
  x <- c(26,43.4,60.8)[i]
  circle(x,94.5,radius_h(i*1000),"#FFFFFF","#425E66",.50)
  label(x+2.0,94.5,format(i*1000),6.7)
}
circle(90.2,94.5,.34,"#526674","#FFFFFF",.30); label(92.5,94.5,"气象台站",6.9)
dev.off()

for (var in 1:2) {
  n <- c(4,5)[var]; pal <- list(rain,evap)[[var]]; lim <- pred_limits[[var]]
  start_fig(n,H)
  for (j in 1:2) {
    m <- map_pair[[j]]
    panel(c(3,86)[j],4.2,letters[j],c("普通克里格","高程协同克里格")[j])
    draw_raster(pred[[(var-1)*2+j]],m,pal,lim)
    draw_polygons(province,m,NA,lwd=.48)
    site_marks(m)
  }
  map_bar(6.6,76.5,73,pal,lim,list(p_ticks,e_ticks)[[var]],
          c("降水量（mm）","蒸发量（mm）")[var])
  scale_bar(map_pair[[1]],137,80.4)
  circle(98,79.3,.34,"#526674","#FFFFFF",.30)
  label(100.5,79.3,"气象台站",7.2)
  dev.off()
}
start_fig(6,90)
for (j in 1:2) {
  m <- map_pair[[j]]; lim <- list(c(0,300),c(150,550))[[j]]
  panel(c(3,86)[j],4.2,letters[j],c("降水量","蒸发量")[j])
  draw_raster(se[[j]],m,uncert,lim);draw_polygons(province,m,NA,lwd=.48);site_marks(m)
  map_bar(c(6.6,89.6)[j],77.0,62,uncert,lim,
          list(c(0,100,200,300),c(150,250,350,450,550))[[j]],"预测标准误（mm）")
}
scale_bar(a,7.3,65.4)
dev.off()

start_fig(9,93)
m <- locator(7,5,95.5,80)
draw_polygons(province,m,"#F7F8F4")
segments(m$x(edges$x1),m$y(edges$y1),m$x(edges$x2),m$y(edges$y2),
         col="#91A5A5",lwd=.49)
radius_r <- function(r) .49+.70*sqrt(abs(r)/100)
ord <- order(abs(sta$ols_resid_p),decreasing=TRUE)
circle(m$x(sta$x_m[ord]),m$y(sta$y_m[ord]),radius_r(sta$ols_resid_p[ord]),
       shade(sta$ols_resid_p[ord],resid,c(-220,220)),"#465855",.48)
scale_bar(m,8,80.4)
true_north(m,29.8,80.4)
map_bar(111,22,48,resid,c(-220,220),c(-200,-100,0,100,200),"降水残差（mm）")
label(111,43,"残差绝对值（mm）",7.4)
for (i in 1:3) {
  x <- c(113.0,129.4,145.8)[i]; v <- i*50
  circle(x,49,radius_r(v),"#F8F8F4","#465855",.52)
  label(x+2.4,49,as.character(v),7)
}
segments(111,62.5,118,62.5,col="#91A5A5",lwd=.62)
label(120.4,62.5,"四近邻连线",7.4)
dev.off()

start_fig(11,95)
county_lim <- list(c(0,100),c(0,800))
county_values <- list(100*obs$percent[ix],obs$jiangshui[ix])
for (j in 1:2) {
  m <- map_pair[[j]]; pal <- list(grass,rain)[[j]]
  panel(c(3,86)[j],4.2,letters[j],c("县域林草覆盖率","县域平均降水量")[j])
  draw_polygons(counties,m,shade(county_values[[j]],pal,county_lim[[j]]),stroke="#798A83",lwd=.28)
  map_bar(c(6.6,89.6)[j],76.6,62,pal,county_lim[[j]],list(c(0,25,50,75,100),p_ticks)[[j]],
          c("林草覆盖率（%）","降水量（mm）")[j])
}
scale_bar(a,7.3,65.4)
rect(6.6,90.1,9.2,92,col="#D9DADA",border=NA);label(11,91,"无有效观测",6.8)
dev.off()




})

# 统计诊断
local({
station <- stations
county <- included
pairs <- tables$pairs
scores <- tables$scores
loads <- tables$loads
stopifnot(nrow(station) == 53, nrow(pairs) == choose(nrow(station), 2),
          nrow(county) == 88, nrow(scores) == nrow(county),
          max(abs(cor(scores[c("z_P", "z_E")],scores[c("PC1", "PC2")]) -
                  as.matrix(loads[c("loading_pc1","loading_pc2")]))) < 1e-8)

ink <- "#273437"
muted <- "#5C686B"
hair <- "#B2BBBC"
guide <- "#E1E6E6"
precip <- "#236176"
evap <- "#A55736"
family <- "Microsoft YaHei"
mm <- function(x) unit(x, "mm")
register_plot <- function(number, height) {
  cairo_pdf(file.path(out, sprintf("fig%02d.pdf", number)), width = 166/25.4,
            height = height/25.4, family = family, bg = "white", fallback_resolution = 600)
  grid.newpage()
  pushViewport(viewport(xscale=c(0,166), yscale=c(height,0)))
}
uv <- function(x) unit(x, "native")
label <- function(x,y,s, size=7.5, col=ink, just="centre", bold=FALSE, rot=0) {
  grid.text(s,uv(x),uv(y),just=just,rot=rot,gp=gpar(fontfamily=family,
    fontface=if (bold) "bold" else "plain", fontsize=size, col=col))
}
segment <- function(x1,y1,x2,y2,col=hair,lwd=.55,lty=1) {
  grid.segments(uv(x1),uv(y1),uv(x2),uv(y2),gp=gpar(col=col,lwd=lwd,lty=lty,lineend="butt"))
}
polyline <- function(x,y,col,lwd=.8,lty=1) {
  grid.lines(uv(x),uv(y),gp=gpar(col=col,lwd=lwd,lty=lty,linejoin="round",lineend="butt"))
}
points_mm <- function(x,y,diam=1,col=ink,rim=NA,lwd=.42) {
  grid.circle(uv(x),uv(y),mm(diam/2),gp=gpar(fill=col,col=rim,lwd=lwd))
}
frame <- function(left,top,width,height,xlim,ylim) {
  list(left=left,right=left+width,top=top,bottom=top+height,
       x=function(v) left+(v-xlim[1])*width/diff(xlim),
       y=function(v) top+height-(v-ylim[1])*height/diff(ylim))
}
tick_labels <- function(x) chartr("-", "−", format(x, trim=TRUE, scientific=FALSE))
axes <- function(p,xt,yt,xlab=TRUE,ylab=TRUE,guide_y=FALSE) {
  if (guide_y) segment(p$left,p$y(yt),p$right,p$y(yt),col=guide,lwd=.38)
  segment(p$left,p$top,p$left,p$bottom,col=muted,lwd=.62)
  segment(p$left,p$bottom,p$right,p$bottom,col=muted,lwd=.62)
  segment(p$x(xt),p$bottom,p$x(xt),p$bottom+1.05,col=muted,lwd=.5)
  segment(p$left-1.05,p$y(yt),p$left,p$y(yt),col=muted,lwd=.5)
  if (xlab) label(p$x(xt),p$bottom+3.3,tick_labels(xt),7.1,muted)
  if (ylab) label(p$left-1.9,p$y(yt),tick_labels(yt),7.1,muted,"right")
}
heading <- function(x,y,letter,title) {
  label(x,y,letter,9.3,ink,"left",TRUE)
  label(x+5,y,title,8.25,ink,"left")
}

# Each rank and the same 1:1 reference are shown separately for each transformation.
register_plot(1, 99)
left <- c(22,70,118)
for (j in 1:3) label(left[j]+20,4.0,c("原值","自然对数","Box–Cox")[j],8.25)
for (row in 1:2) {
  top <- c(9,51)[row]
  variable <- c("precip","evap")[row]
  tone <- c(precip,evap)[row]
  label(3,top+13,letters[row],9.3,ink,"left",TRUE)
  label(3,top+19,c("降水量","蒸发量")[row],7.9,ink,"left")
  for (j in 1:3) {
    p <- frame(left[j],top,41,37.8,c(-3.35,3.35),c(-2.75,3.43))
    axes(p,c(-2,0,2),c(-2,0,2),xlab=row==2,ylab=j==1)
    ref <- c(-2.75,3.35)
    segment(p$x(ref[1]),p$y(ref[1]),p$x(ref[2]),p$y(ref[2]),col="#98A6A6",lwd=.62,lty="32")
    obs <- sort(as.numeric(scale(station[[paste0(variable,c("","_log","_bc")[j])]])))
    points_mm(p$x(qnorm(ppoints(length(obs)))),p$y(obs),1.06,tone)
  }
}
label(2.6,3.8,"样本分位数",7.1,muted,"left")
label(90,96,"理论标准正态分位数",8)
dev.off()

# All station pairs are displayed. An overlay is the mean of the two coordinates
# within ten equal-width distance classes over the full observed lag range.
register_plot(3,77)
fine <- "#A6B8BC"
for (j in 1:2) {
  l <- c(17,99)[j]
  pr <- j==1
  g <- pairs[[if(pr) "gamma_p" else "gamma_e"]]/1000
  maxy <- if(pr) 300 else 3300
  p <- frame(l,12,63.3,47,c(0,1340),c(0,maxy))
  heading(l-1,4.2,letters[j],c("降水量","蒸发量")[j])
  label(l+62,4.2,sprintf("Spearman ρ = %.3f",cor(pairs$dist_km,g,method="spearman")),7.25,muted,"right")
  axes(p,c(0,400,800,1200),if(pr) c(0,100,200,300) else c(0,1000,2000,3000),guide_y=TRUE)
  points_mm(p$x(pairs$dist_km),p$y(g),.63,adjustcolor(fine,alpha.f=.48))
  group <- cut(pairs$dist_km,breaks=seq(0,max(pairs$dist_km),length.out=11),include.lowest=TRUE)
  d <- as.numeric(tapply(pairs$dist_km,group,mean))
  av <- as.numeric(tapply(g,group,mean))
  tone <- if(pr) precip else evap
  polyline(p$x(d),p$y(av),tone,1.02)
  points_mm(p$x(d),p$y(av),1.50,"white",tone,.88)
}
label(2.7,36,"半方差（10³ mm²）",7.8,rot=90)
label(89.5,68,"台站间距（km）",8)
points_mm(62,74,.95,fine)
label(64,74,"站点对",7.1,muted,"left")
segment(82,74,87,74,muted,.75)
points_mm(84.5,74,1.3,"white",muted,.7)
label(89,74,"距离分箱均值",7.1,muted,"left")
dev.off()

# OLS and the confidence band for the mean response use the original 88 rows.
fits <- lapply(c("jiangshui","zhengfa"), function(var) {
  d <- data.frame(x=county[[var]],y=county$percent)
  model <- lm(y~x,data=d)
  grid <- data.frame(x=seq(min(d$x),max(d$x),length.out=200))
  list(obs=d, xs=grid$x, band=100*predict(model,grid,interval="confidence"), model=model)
})
fitted_subei <- fitted(fits[[2]]$model)[county$ZoneFID==24]
stopifnot(round(fitted_subei,4)==-.0188)
register_plot(12,82)
all_y <- range(sapply(fits, function(f) range(f$band)),100*county$percent)
yrange <- c(min(-10,all_y[1]-1), max(105,all_y[2]+1))
for (j in 1:2) {
  f <- fits[[j]]
  l <- c(18,96)[j]
  tone <- c(precip,evap)[j]
  p <- frame(l,12,65.0,45.5,c(if(j==1) c(15,745) else c(1100,3080)),yrange)
  heading(l-1,4.1,letters[j],c("降水量","蒸发量")[j])
  label(l+64.8,4.1,sprintf("r = %.3f",cor(f$obs$x,f$obs$y)),7.25,muted,"right")
  xt <- if(j==1) c(200,400,600) else c(1500,2000,2500,3000)
  axes(p,xt,c(0,25,50,75,100),ylab=j==1,guide_y=TRUE)
  gc <- if(j==1) "#DAE7EA" else "#F0DDD2"
  grid.polygon(uv(p$x(c(f$xs,rev(f$xs)))),uv(p$y(c(f$band[,"lwr"],rev(f$band[,"upr"])))),
               gp=gpar(col=NA,fill=gc))
  segment(p$left,p$y(0),p$right,p$y(0),"#9BA5A5",.57,"33")
  polyline(p$x(f$xs),p$y(f$band[,"fit"]),tone,1.05)
  points_mm(p$x(f$obs$x),p$y(100*f$obs$y),1.22,adjustcolor(tone,alpha.f=.67))
  if(j==2){
    xx <- p$x(county$zhengfa[county$ZoneFID==24])
    yy <- p$y(fitted_subei*100)
    segment(xx-.78,yy-.78,xx+.78,yy+.78,ink,.92)
    segment(xx+.78,yy-.78,xx-.78,yy+.78,ink,.92)
  }
  label(l+32.5,67,c("多年平均降水量（mm）","多年平均蒸发量（mm）")[j],8)
}
label(3,35,"林草覆盖率（%）",8,rot=90)
points_mm(26,77,1.3,muted)
label(28.4,77,"县域观测",7.25,muted,"left")
grid.rect(uv(60),uv(77),mm(8),mm(2.1),gp=gpar(col=NA,fill="#D7DDDF"))
segment(56,77,64,77,muted,.92)
label(66,77,"OLS 拟合与 95% 均值置信带",7.25,muted,"left")
segment(116.2,76.25,117.7,77.75,ink,.86)
segment(117.7,76.25,116.2,77.75,ink,.86)
label(120,77,"肃北北片的越界拟合值",7.25,muted,"left")
dev.off()

# Left panel retains the entire standardized observation cloud; the correlation
# circle in the right panel uses measured correlations with both component scores.
register_plot(13,93)
heading(16,4.5,"a","县域水热分布")
heading(98,4.5,"b","主成分载荷")
p <- frame(18,11.6,52.0,59.9,c(-2.3,2.3),c(-1.65,3.65))
axes(p,c(-2,-1,0,1,2),c(-1,0,1,2,3))
segment(p$x(0),p$top,p$x(0),p$bottom,guide,.55,"33")
segment(p$left,p$y(0),p$right,p$y(0),guide,.55,"33")
ramp <- colorRampPalette(c("#E3EAD4","#ADC5A3","#63977D","#24674F","#164B3E"),space="Lab")(128)
index <- pmin(128,floor(scores$percent*128)+1)
points_mm(p$x(scores$z_P),p$y(scores$z_E),1.75,ramp[index],"#517565",.40)
label(44,80.3,"标准化降水量",8)
label(5,42,"标准化蒸发量",8,rot=90)
label(6.6,87.5,"林草覆盖率",7.0,muted,"left")
for (i in 1:64) {
  grid.rect(uv(24+(i-.5)*.72),uv(87.4),mm(.74),mm(1.7),gp=gpar(col=NA,fill=ramp[2*i]))
}
label(c(24,47.04,70.08),90.4,c("0%","50%","100%"),6.65,muted)
p <- frame(100,16,53.5,53.5,c(-1.15,1.15),c(-1.15,1.15))
axes(p,c(-1,0,1),c(-1,0,1))
ang <- seq(0,2*pi,length.out=361)
polyline(p$x(cos(ang)),p$y(sin(ang)),"#A6B2B0",.72)
segment(p$x(-1),p$y(0),p$x(1),p$y(0),guide,.63)
segment(p$x(0),p$y(-1),p$x(0),p$y(1),guide,.63)
for (j in 1:2) {
  tone <- c(precip,evap)[j]
  xx <- loads$loading_pc1[j]
  yy <- loads$loading_pc2[j]
  grid.segments(uv(p$x(0)),uv(p$y(0)),uv(p$x(xx)),uv(p$y(yy)),
                arrow=arrow(length=mm(1.6),type="closed",angle=22),
                gp=gpar(col=tone,fill=tone,lwd=.90))
  label(p$x(.57*sign(xx)),p$y(yy)-4.3,c("降水量","蒸发量")[j],7.4,tone)
}
vp <- colSums(as.matrix(loads[c("loading_pc1","loading_pc2")])^2)/2*100
label(126.8,80.3,sprintf("PC1（%.2f%%）",vp[1]),8)
label(87.3,42,sprintf("PC2（%.2f%%）",vp[2]),8,rot=90)
dev.off()

})

# 模型比较
local({
target <- out
ink <- "#263338"
muted <- "#556368"
light <- "#D8DFDF"
faint <- "#E9EEEE"
rain <- "#216A78"
evap <- "#A55438"
font <- "Microsoft YaHei"
u <- function(value) unit(value, "native")
txt <- function(x, y, label, size=8.2, colour=ink, align="left", face="plain") {
  grid.text(label, x=u(x), y=u(y), just=c(align, "centre"),
            gp=gpar(fontfamily=font, fontsize=size, col=colour, fontface=face))
}
seg <- function(x1,y1,x2,y2,colour=light,weight=0.68,lty=1) {
  grid.segments(u(x1),u(y1),u(x2),u(y2),gp=gpar(col=colour,lwd=weight,lty=lty,lineend="butt"))
}
dot <- function(x,y,diameter=1.75,fill=ink,stroke=fill,form=21,weight=0.8) {
  grid.points(u(x),u(y),pch=form,size=unit(diameter*1.45,"mm"),
              gp=gpar(col=stroke,fill=fill,lwd=weight))
}
tile <- function(x,y,w,h,fill) {
  grid.rect(u(x),u(y),unit(w,"mm"),unit(h,"mm"),gp=gpar(fill=fill,col=NA))
}
draw_pdf <- function(filename, height, fn) {
  cairo_pdf(file.path(target,filename),width=166/25.4,height=height/25.4,family=font,bg="white")
  on.exit(dev.off())
  grid.newpage()
  pushViewport(viewport(xscale=c(0,166),yscale=c(height,0)))
  fn()
  popViewport()
}
map <- function(value, bounds, left, right) {
  left + (value-bounds[1]) / diff(bounds) * (right-left)
}
panel <- function(x,y,letter,label) {
  txt(x,y,letter,10.1,ink,face="bold")
  txt(x+5.6,y,label,9.1,ink,face="bold")
}
axis_bottom <- function(left,right,y,bounds,ticks,labels=ticks) {
  pos <- map(ticks,bounds,left,right)
  seg(left,y,right,y,muted,0.65)
  seg(pos,y,pos,y+1.1,muted,0.58)
  txt(pos,y+4.5,labels,7.4,muted,"centre")
}

models <- c("precip_OK_log_spherical_T0_N16","precip_CK_log_spherical_T0_N16",
            "evap_OK_log_gaussian_T0_N16","evap_CK_log_exponential_T0_N16")
ck <- tables$model_comparison
cv <- tables$validation_metrics
nested <- tables$nested_metrics
ck <- subset(ck,model %in% models)
ck$scheme <- "cond_LOO"
cols <- c("scheme","variable","method","RMSE")
cv <- rbind(ck[,cols],cv[cv$model %in% models,cols],nested[,cols])
schemes <- c("cond_LOO","refit_LOO","spatial8","spatial5","spatial5_buffer100km","nested_spatial5")
stopifnot(nrow(cv)==24,!anyDuplicated(cv[c("scheme","variable","method")]))
draw_pdf("fig07.pdf",77,function() {
  ys <- c(18,25,35,42,49,56)
  txt(2,6.5,"验证方案",8.4,face="bold")
  txt(2,ys,c("条件留一","重拟合留一","空间 8 区","空间 5 区","5 区 + 100 km 缓冲","外 5 内 3 嵌套"),8.1)
  seg(2,30,164,30,light,0.6)
  for (j in 1:2) {
    precip <- j==1
    left <- c(40,108)[j]; right <- c(99,163)[j]
    bounds <- if(precip) c(55,157) else c(238,337)
    ticks <- if(precip) c(60,90,120,150) else c(240,270,300,330)
    col <- if(precip) rain else evap
    di <- cv[cv$variable == c("precip","evap")[j],]
    panel(left,6.5,letters[j],c("降水量","蒸发量")[j])
    x <- map(ticks,bounds,left,right)
    seg(x,13,x,60,faint,0.55)
    for (i in seq_along(schemes)) {
      values <- di[di$scheme==schemes[i],]
      xs <- map(values$RMSE[match(c("OK","CK"),values$method)],bounds,left,right)
      stopifnot(length(xs)==2, all(is.finite(xs)), all(xs>=left), all(xs<=right))
      seg(xs[1],ys[i]-1.25,xs[2],ys[i]+1.25,col,0.9)
      dot(xs[1],ys[i]-1.25,1.85,col)
      dot(xs[2],ys[i]+1.25,2,"white",col,23,1.12)
    }
    axis_bottom(left,right,60.5,bounds,ticks)
  }
  dot(3,73.3,1.85,ink); txt(6.3,73.2,"普通克里格",7.8)
  dot(28,73.3,2,"white",ink,23,1.1); txt(31.3,73.2,"高程协同克里格",7.8)
  txt(163,73.2,"预测均方根误差（mm）",8.2,ink,"right")
})

range_fit <- tables$variogram_sensitivity
range_cv <- tables$sensitivity_spatial5_metrics
ss <- merge(range_fit,range_cv,by=c("model","range_multiplier","bins"))
stopifnot(nrow(ss)==36,all(ss$model %in% models),all(as.logical(ss$complete)))
ss$ref <- NA_real_
for (m in models) {
  ref <- ss$RMSE[ss$model==m & ss$range_multiplier==2 & ss$bins==10]
  stopifnot(length(ref)==1)
  ss$ref[ss$model==m] <- ref
}
ss$change <- ss$RMSE-ss$ref
max_range <- max(ss$upper_km)
col_range <- grDevices::colorRamp(c("#EDF3F0","#AACABF","#509383","#215E53"),space="Lab")
col_change <- grDevices::colorRamp(c("#397994","#ACCDDB","#F4F3EC","#E0B59A","#A35940"),space="Lab")
colour <- function(ramp,fraction) rgb(ramp(pmax(0,pmin(1,fraction))),maxColorValue=255)
draw_pdf("fig08.pdf",98,function() {
  step <- 11.35; gap <- 2.7
  xs <- c(40+step*(0:2),40+3*step+gap+step*(0:2),40+6*step+2*gap+step*(0:2))
  names <- c("降水  OK","降水  CK","蒸发  OK","蒸发  CK")
  panel(2,5,"a","拟合变程")
  for(k in 1:41) tile(94+k*0.67,4.8,.70,2.1,colour(col_range,(k-1)/40))
  txt(92.6,4.9,"0",7,muted,"right")
  txt(123.5,4.9,paste(format(round(max_range,1),big.mark=" ",nsmall=1),"km"),7.1,muted)
  dot(146.4,4.9,1.05,"white",ink,21,.7); txt(149,4.9,"触及上限",7.1,muted)
  txt(2,13.4,"搜索上限",7.7,muted)
  for (g in 1:3) {
    subx <- xs[(g-1)*3+1:3]
    txt(mean(subx),13.4,c("1 倍","2 倍","4 倍")[g],8.1,ink,"centre")
    seg(min(subx)-5.2,16.2,max(subx)+5.2,16.2,light,.65)
  }
  txt(2,20.4,"滞后分箱数",7.7,muted)
  txt(xs,20.4,rep(c("8","10","14"),3),7.2,muted,"centre")
  for (i in 1:4) {
    m <- ss[ss$model==models[i],]
    m <- m[order(m$range_multiplier,m$bins),]
    ya <- 26.8+(i-1)*5.4
    yb <- 64.8+(i-1)*5.4
    txt(2,ya,names[i],8)
    txt(2,yb,names[i],8)
    for (k in 1:9) {
      tile(xs[k],ya,step-.45,5.0,colour(col_range,m$range_km[k]/max_range))
      if(isTRUE(as.logical(m$upper_hit[k]))) dot(xs[k],ya,1.06,"white",ink,21,.70)
      tile(xs[k],yb,step-.45,5.0,colour(col_change,(m$change[k]+8)/16))
    }
    txt(160,yb,sprintf("%.2f",m$ref[1]),7.6,ink,"right")
  }
  panel(2,51.7,"b","空间 5 区验证误差变化")
  for(k in 1:61) tile(101+k*.49,51.7,.53,2.1,colour(col_change,(k-1)/60))
  txt(99.7,51.8,"−8",7.1,muted,"right")
  txt(116.1,55.1,"0",6.8,muted,"centre")
  txt(132.6,51.8,"+8 mm",7.1,muted)
  txt(xs,58.9,rep(c("8","10","14"),3),7.2,muted,"centre")
  txt(160,58.9,"基准 / mm",7.2,muted,"right")
  txt(2,93,"误差变化为各方案的 RMSE 减去同一模型在 2 倍上限、10 箱下的 RMSE。",7.7,muted)
})

coef <- tables$coef_forest
stopifnot(nrow(coef)==24,all(coef$ci_low<=coef$estimate),all(coef$estimate<=coef$ci_high))
methods <- c("经典 OLS","HC3 稳健 SE","SEM (4近邻)","SEM (8近邻)")
draw_pdf("fig10.pdf",96,function() {
  left <- c(33,79,124); right <- c(70,115,162)
  bounds <- list(c(-105,250),c(-83,52),c(-0.31,0.16))
  ticks <- list(c(-100,0,100,200),c(-60,-30,0,30),c(-.3,-.15,0,.15))
  xtick <- list(c("−100","0","100","200"),c("−60","−30","0","30"),c("−0.30","−0.15","0","0.15"))
  txt(left,5,c("纬度","经度","高程"),8.8,ink,face="bold")
  txt(left,10.2,c("mm / °","mm / °","mm / m"),7.6,muted)
  for (k in 1:3) {
    zp <- map(0,bounds[[k]],left[k],right[k])
    for (j in 1:2) {
      y0 <- c(23,60.5)[j]; ys <- y0 + (0:3)*6.1
      colour <- c(rain,evap)[j]
      if(k==1) {
        panel(2,y0-7,letters[j],c("降水量","蒸发量")[j])
        txt(2,ys,c("OLS","OLS · HC3","SEM · 4 近邻","SEM · 8 近邻"),7.8)
      }
      seg(zp,y0-3.1,zp,max(ys)+3.1,"#B4BEBF",.72,"22")
      d <- coef[coef$term==c("Lat","Lon","Elev")[k] & coef$dependent==c("降水量","蒸发量")[j],]
      d <- d[match(methods,d$method),]
      xs <- function(v) map(v,bounds[[k]],left[k],right[k])
      stopifnot(all(xs(d$ci_low)>=left[k]),all(xs(d$ci_high)<=right[k]))
      seg(xs(d$ci_low),ys,xs(d$ci_high),ys,colour,1.08)
      seg(xs(d$ci_low),ys-.8,xs(d$ci_low),ys+.8,colour,.77)
      seg(xs(d$ci_high),ys-.8,xs(d$ci_high),ys+.8,colour,.77)
      dot(xs(d$estimate),ys,1.63,colour)
    }
    axis_bottom(left[k],right[k],84,bounds[[k]],ticks[[k]],xtick[[k]])
  }
  seg(2,47.6,163,47.6,light,.62)
  txt(163,93.6,"偏回归系数及 95% 置信区间",8.1,ink,"right")
})

cv_pref <- tables$prediction_validation_metrics
ridge <- tables$strict_pc1_validation_metrics
logit <- tables$strict_pc1_one_round_metrics
schemes <- c("spatial4","spatial5","spatial8")
cv_county <- rbind(subset(cv_pref,scheme %in% schemes & model %in% c("P_only","E_only","PCR_85","OLS_two"))[c("scheme","model","RMSE")],
                   subset(ridge,scheme %in% schemes & model=="PC1_nested_shrinkage")[c("scheme","model","RMSE")],
                   subset(logit,scheme %in% schemes & model=="PC1_fractional_logit")[c("scheme","model","RMSE")])
order_m <- c("E_only","P_only","OLS_two","PCR_85","PC1_nested_shrinkage","PC1_fractional_logit")
stopifnot(nrow(cv_county)==18,!anyDuplicated(cv_county[c("scheme","model")]))
draw_pdf("fig14.pdf",77,function() {
  ys <- c(17,24.2,31.4,42,49.2,56.4)
  txt(2,6,"县域模型",8.4,face="bold")
  txt(2,ys,c("蒸发一元 OLS","降水一元 OLS","水热双变量 OLS","主成分回归","PC1 岭回归","PC1 分数阶 logit"),8.0)
  seg(2,36.7,164,36.7,light,.6)
  left <- c(36,81,126); right <- left+36.5
  bounds <- c(.155,.206)
  ticks <- c(.16,.18,.20)
  for (j in 1:3) {
    data <- cv_county[cv_county$scheme==schemes[j],]
    vals <- data$RMSE[match(order_m,data$model)]
    stopifnot(length(vals)==6,all(is.finite(vals)),all(vals>=bounds[1] & vals<=bounds[2]))
    panel(left[j],6,letters[j],c("空间 4 折","空间 5 折","空间 8 折")[j])
    xp <- function(v) map(v,bounds,left[j],right[j])
    seg(xp(ticks),12,xp(ticks),60.5,faint,.56)
    seg(xp(vals[3]),12,xp(vals[3]),60.5,muted,.74,"22")
    for(i in 1:6) {
      hi <- i==4
      dot(xp(vals[i]),ys[i],if(hi) 2.05 else 1.65,
          if(i==3) "white" else if(hi) rain else "#68797B",
          if(hi) rain else "#617275",if(hi) 23 else 21,.92)
    }
    axis_bottom(left[j],right[j],61,bounds,ticks,format(ticks,nsmall=2))
  }
  seg(2,73,8,73,muted,.74,"22");txt(10,73,"双变量 OLS 基准",7.6,muted)
  txt(163,73,"林草覆盖率预测 RMSE（比例）",8.0,ink,"right")
})
})
