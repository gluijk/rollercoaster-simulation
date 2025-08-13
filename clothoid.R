# Plotting clothoids in R
# www.overfitting.net
# https://www.overfitting.net/2025/08/clotoides-las-curvas-que-cuidan-de-tu.html

library(pracma)  # fresnelC() and fresnelS()
library(Cairo)
library(png)


# Function to generate a standard clothoid (Euler spiral)
# of length L, reaching a min radius of Rmin, with n points of resolution
# (L and Rmin share units)
# Function outputs: x,y coords, s space, radius and centres of circumferences
build_clothoid <- function(L, Rmin, n=1000) {
    # Scale factor so curvature radius at s=L equals Rmin
    A <- sqrt(pi * Rmin * L)
    
    # Arc length sampling
    s <- seq(0, L, length.out = n)
    
    # Normalized arc length
    u <- s / A
    
    # Calculate clothoid (x,y) positions
    x <- A * fresnelC(u)
    y <- A * fresnelS(u)
    
    # Derivatives dx/ds and dy/ds (tangent vector components)
    dx_ds <- cos(pi * u^2 / 2)
    dy_ds <- sin(pi * u^2 / 2)
    
    # Radius of curvature R(s) = Rmin * L / s
    # Handle s=0 separately (radius = Inf)
    radius <- ifelse(s == 0, NA, Rmin * L / s)
    
    # Normal unit vector (perpendicular to tangent vector rotated 90 degrees CCW)
    nx <- -dy_ds
    ny <- dx_ds
    
    # Center of curvature coordinates
    center_x <- x + radius * nx
    center_y <- y + radius * ny
    
    return(list(
        x = x,  # coords of the clothoid
        y = y,
        s = s,  # length position of the clothoid
        radius = radius,  # radius at each point of the clothoid
        center_x = center_x,  # centres of the tangent circumferences
        center_y = center_y
    ))
}



##############################

# 1. BASIC EXAMPLE

NPOINTS=1000
L=100   # clothoid total length (m)
Rmin=3  # min curvature radius (m)

clothoid=build_clothoid(L, Rmin, NPOINTS)
x=clothoid$x
y=clothoid$y
s=clothoid$s
R=clothoid$radius
xc=clothoid$center_x
yc=clothoid$center_y

CairoPNG("clothoid.png", width=512, height=512)
    # Plot clothoid
    plot(x, y, type="l", asp=1, col="red", lwd=3,
         xlim=c(0,max(x)), ylim=c(0,max(y)),
         xlab="X (m)", ylab="Y (m)",
         main=paste0("Clothoid of L=", L, "m and Rmin=", Rmin, "m"))
    abline(h=0, col='gray')
    symbols(x=xc[NPOINTS], y=yc[NPOINTS], circles=R[NPOINTS],
            inches=FALSE, add=TRUE, fg='gray')
    segments(x[NPOINTS], y[NPOINTS], xc[NPOINTS], yc[NPOINTS],
             col='gray', lty='dotted')
dev.off()

CairoPNG("clothoid_params.png", width=512, height=800)
    par(mfrow = c(2, 1))
    
    # Plot radius vs arc length
    Rs=(Rmin * L) / s
    plot(s, Rs, type="l", col="red", lwd=3,
         ylim=c(0, 100),
         xlab="Arc length (m)", ylab="Radius (m)",
         main="Radius of curvature along the clothoid")
    abline(h=Rmin, col='gray')
    
    # Plot curvature (kappa=1/r) vs arc length
    kappa=1/Rs  # curvature is proportional to centrifugal force
    plot(s, kappa, type="l", col="red", lwd=3,
         xlab="Arc length (m)", ylab="Curvature kappa = 1/Radius (m^-1)",
         main="Curvature along the clothoid")
    #grid()
    
    par(mfrow = c(1, 1))
dev.off()


##############################

# 2. ANIMATION

NPOINTS=1000
L=50   # clothoid total length (m)
RMIN=0.1  # min curvature radius (m)
RMAX=200
NCURVES=100

NFRAMES=1452  # 60,5s at 24fps
for (i in 0:(NFRAMES-1)) {
    name=paste0("clothoidanim_", ifelse(i<10, "000", ifelse(i<100, "00",
                ifelse(i<1000, "0", ""))), i, ".png")
    print(name)
    CairoPNG(name, width=1080, height=1080)
        plot(1, type="l", asp=1, xlim=c(-30,30), ylim=c(0,50),  # NULL plot
             axes=FALSE, xlab="", ylab="")
        for (r in seq(0, 1, length.out=NCURVES)) {
            Rmin=r^1.5 * (((RMAX-RMIN)/(NFRAMES-1)*i+RMIN)-RMIN)+RMIN
            clothoid=build_clothoid(L, Rmin, NPOINTS)

            # Plot clothoid
            colour=rgb(Rmin/RMAX*0.7, 0.5, max((1-Rmin/RMAX)*0.6,0), 0.5)
            lines( clothoid$y, clothoid$x, col=colour, lwd=1)
            lines(-clothoid$y, clothoid$x, col=colour, lwd=1)
        }
    dev.off()
}

# Add background to each frame
bg=readPNG("background.png")
for (i in 0:(NFRAMES-1)) {
    name=paste0("clothoidanim_", ifelse(i<10, "000", ifelse(i<100, "00",
                ifelse(i<1000, "0", ""))), i, ".png")
    img=readPNG(name)
    imgout=bg
    imgout[,1:1080,]=img
    imgout[bg<1]=bg[bg<1]
    
    name=paste0("clothoidanimfinal_", ifelse(i<10, "000", ifelse(i<100, "00",
            ifelse(i<1000, "0", ""))), i, ".png")
    writePNG(imgout, name)
}
    

# MP4 Video (MPEG-4 AVC/H.264):
# ffmpeg -framerate 24 -i clothoidanim_%04d.png -i theexpanse.wav -c:v libx264
# -crf 18 -pix_fmt yuv420p clothoids.mp4
