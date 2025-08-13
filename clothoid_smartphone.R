# Plotting clothoids in R for smartphone soft borders
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
        x = x,
        y = y,
        s = s,
        radius = radius,
        center_x = center_x,
        center_y = center_y
    ))
}


##############################

# 1. FIND L FOR WHICH WE HAVE 45º (Rmin=0.5)

difangulo45 <- function(L) {
    Rmin=0.5
    n=2
    
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
    
    # Radius of curvature R(s) = Rmin * L / s; handle s=0 separately (radius = Inf)
    radius <- ifelse(s == 0, NA, Rmin * L / s)
    
    # Normal unit vector (perpendicular to tangent vector, rotated 90° CCW)
    nx <- -dy_ds
    ny <- dx_ds
    
    # Center of curvature coordinates
    center_x <- x + radius * nx
    center_y <- y + radius * ny

    return (atan((center_y[n]-y[n])/(center_x[n]-x[n]))*180/pi + 45)
}

sol <- uniroot(difangulo45, interval = c(0.5, 0.8))

# The value of L where difangulo45(L) = 0
sol$root  # L=0.785398163397

# The function value at that root (should be 0 or very close)
sol$f.root


##############################

# 2. PLOT SEMI-COMBINED CLOTOID (MIRRORED IN PHOTOSHOP)

NPOINTS=1000
L=0.785398163397   # clothoid total length (cm)
Rmin=0.5  # min curvature radius (cm)

clothoid=build_clothoid(L, Rmin, NPOINTS)
x=clothoid$x
y=clothoid$y
s=clothoid$s
R=clothoid$radius
xc=clothoid$center_x
yc=clothoid$center_y

angulo=atan((yc[NPOINTS]-y[NPOINTS])/(xc[NPOINTS]-x[NPOINTS]))*180/pi
print(paste0("Angle: ", angulo, "º"))

CairoPNG("clothoid_smartphone.png", width=1920, height=1920)
    # Plot clothoid
    plot(x, y, type="l", asp=1, col="red", lwd=6,
         xlim=c(0,max(x)*1.1), ylim=c(0,1),
         xlab="X (cm)", ylab="Y (cm)",
         main=paste0("Clothoid of L=", L, "m and Rmin=", Rmin, "m"))
    # symbols(x=xc[NPOINTS], y=yc[NPOINTS], circles=R[NPOINTS],
    #         inches=FALSE, add=TRUE, fg='gray')
    segments(x[NPOINTS], y[NPOINTS], xc[NPOINTS], yc[NPOINTS], col='gray', lty='dotted')
    
    # abline(h=c(0, yc[NPOINTS], Rcirc), v=c(0, x0), col='gray')
    abline(h=0, v=0, col='gray')
    
    x0=x[NPOINTS]/3
    Rcirc=x[NPOINTS]+y[NPOINTS]-x0
    symbols(x=x0, y=Rcirc, circles=Rcirc,
            inches=FALSE, add=TRUE, fg='black')
dev.off()

    

