# Designing a Constant G rollercoaster loop
# www.overfitting.net
# https://www.overfitting.net/2025/08/simulando-loops-de-una-montana-rusa-con_69.html


# library(pracma)  # fresnelC() and fresnelS() for clothoid
# library(splines)
library(ggplot2)
library(Cairo)


# Design a constant-G loop with a clothoid entry (centrifugal OUTWARD convention)
design_constantG_loop <- function(v0, Gt, g = 9.80665,
                                  Lc = 5,        # clothoid length (m)
                                  ds_cl = 0.01,  # clothoid arc step (m)
                                  dphi = 0.002,  # loop integration step in radians
                                  y0 = 0.0,
                                  nloop = 1) {  # number of G constant loops
    stopifnot(Gt > 1, v0 > 0, g > 0, Lc >= 0, ds_cl > 0, dphi > 0)
    
    # Storage
    xs <- numeric()
    ys <- numeric()
    phis <- numeric()
    ss <- numeric()
    phases <- character()
    Ginst <- numeric()
    
    # Helpers ---------------------------------------------------------------
    a_req <- function(phi) {
        # required inward acceleration to keep |g + a_out| = Gt*g
        root <- Gt^2 - sin(phi)^2
        if (root < 0) root <- 0 # numeric guard
        g * (-cos(phi) + sqrt(root))
    }
    totalG_from <- function(phi, v2, kappa) {
        # magnitude of g + a_out (outward) in units of g
        a_in <- v2 * kappa
        # cos(theta) between down and inward is -cos(phi)
        # |g e_down - a_in e_in| = sqrt(g^2 + a_in^2 - 2 g a_in cos(theta))
        # with cos(theta) = -cos(phi)
        num <- sqrt(g^2 + a_in^2 - 2*g*a_in*(-cos(phi)))
        num / g
    }
    
    # Phase 0: initial conditions ------------------------------------------
    x <- 0; y <- y0; phi <- 0; s <- 0
    v2 <- v0^2
    
    # Phase 1: clothoid entry (linear curvature from 0 to bottom curvature)
    # Bottom curvature that achieves target G at phi ~ 0:
    a_bottom <- g * (Gt - 1)                  # inward accel needed at the bottom
    kappa_target_bottom <- a_bottom / max(v2, 1e-9)
    
    if (Lc > 0) {
        alpha <- kappa_target_bottom / Lc       # kappa(s) = alpha * s
        s_loc <- 0
        while (s_loc < Lc) {
            kappa <- alpha * s_loc
            dphi_cl <- kappa * ds_cl
            phi <- phi + dphi_cl
            # advance position
            x <- x + cos(phi) * ds_cl
            y <- y + sin(phi) * ds_cl
            s <- s + ds_cl; s_loc <- s_loc + ds_cl
            # update speed via energy
            v2 <- v0^2 - 2*g*(y - y0)
            if (v2 <= 0) stop("Insufficient v0: train stalls during clothoid.")
            # log
            xs <- c(xs, x); ys <- c(ys, y); phis <- c(phis, phi); ss <- c(ss, s); phases <- c(phases, "clothoid")
            Ginst <- c(Ginst, totalG_from(phi, v2, kappa))
        }
    }
    
    # Phase 2: constant-G loop (integrate in φ with RK4)
    # ODEs: dy/dφ = sinφ * v^2 / a_req(φ) ; dx/dφ = cosφ * v^2 / a_req(φ)
    fy <- function(phi, y) {
        v2loc <- v0^2 - 2*g*(y - y0)
        if (v2loc <= 0) stop("Insufficient v0: train stalls in loop.")
        v2loc * sin(phi) / max(a_req(phi), 1e-9)
    }
    fx <- function(phi, y) {
        v2loc <- v0^2 - 2*g*(y - y0)
        if (v2loc <= 0) stop("Insufficient v0: train stalls in loop.")
        v2loc * cos(phi) / max(a_req(phi), 1e-9)
    }
    
    phi0 <- phi
    phi_end <- phi0 + 2*pi*nloop
    n_steps <- ceiling((phi_end - phi0) / dphi)
    for (i in 1:n_steps) {
        # RK4 for y(φ)
        k1y <- fy(phi, y)
        k1x <- fx(phi, y)
        
        k2y <- fy(phi + dphi/2, y + k1y*dphi/2)
        k2x <- fx(phi + dphi/2, y + k1y*dphi/2)
        
        k3y <- fy(phi + dphi/2, y + k2y*dphi/2)
        k3x <- fx(phi + dphi/2, y + k2y*dphi/2)
        
        k4y <- fy(phi + dphi, y + k3y*dphi)
        k4x <- fx(phi + dphi, y + k3y*dphi)
        
        y <- y + (dphi/6)*(k1y + 2*k2y + 2*k3y + k4y)
        x <- x + (dphi/6)*(k1x + 2*k2x + 2*k3x + k4x)
        phi <- phi + dphi
        
        # arc-length increment ds = dφ / κ = v^2/a_req * dφ
        v2 <- v0^2 - 2*g*(y - y0)
        a_in <- a_req(phi - dphi/2)  # mid-step approx
        ds_inc <- (v2 / max(a_in, 1e-9)) * dphi
        s <- s + ds_inc
        
        # log
        kappa_here <- a_in / max(v2, 1e-9)
        xs <- c(xs, x); ys <- c(ys, y); phis <- c(phis, phi); ss <- c(ss, s); phases <- c(phases, "loop")
        Ginst <- c(Ginst, totalG_from(phi, v2, kappa_here))
    }
    
    out <- data.frame(
        x = xs, y = ys, phi = phis, s = ss, phase = phases, G_total = Ginst
    )
    attr(out, "params") <- list(v0=v0, Gt=Gt, g=g, Lc=Lc, ds_cl=ds_cl, dphi=dphi)
    out
}


################################################

# 1. SIMULATION PARAMETERS

# Example: 28m/s entry speed, target 4.2G total, 10m clothoid

# Parameters needed by our simulator from rollercoaster topology:
g=9.80665
v0=28
t_max=999
dt=0.05      # time step (s)

# Parameters for current new design:
Gt=4.2  # number of constant G forces (Gt times g)
Lc=10

simu='3/4 loop'
#simu='3 loops'  # special run for 3 loops G constant rollercoaster
if (simu=='3/4 loop') {
    # Clothoidal section + constant G section (1 loop)
    track <- design_constantG_loop(v0 = v0, Gt = Gt, g = g, Lc = Lc, nloop = 0.73)
} else if (simu=='3 loops') {
    # Constant G section (3 loops)
    track <- design_constantG_loop(v0 = v0, Gt = Gt, g = g, Lc = 0, nloop = 3)
}

# Quick check of G constancy in the loop:
with(track[track$phase=="loop",], range(G_total))
# Should be ~ c(4.2, 4.2) apart from tiny numerical noise

# Plot
plot(track$x, track$y, type="l", asp=1, xlab="x (m)", ylab="y (m)")
grid()


################################################

# 2. DEFINE CONSTANT G ROLLERCOASTER

if (simu=='3/4 loop') {
    # Horizontal section + clothoidal section + constant G section (3/4 loop)
    x_points=c(-15, 0, track$x)+15
    y_points=c(0, 0, track$y)
} else if (simu=='3 loops') {
    # Constant G section (3 loops)
    x_points=track$x
    y_points=track$y
}

H=max(y_points)


# Check loop height vs initial speed compatibility

Hmax=v0^2/(2*g)
v0min=(2*g*H)^0.5
if (H>Hmax) stop("Loop height / initial speed incompatibility")


################################################

# 3. SIMULATION LOOP: CALCULATE ROLLERCOASTER POSITION, SPEED AND G FORCES


# Parameterize path by cumulative distance 's'
dx <- diff(x_points)
dy <- diff(y_points)
ds <- sqrt(dx^2 + dy^2)
s_points <- c(0, cumsum(ds))

# Total track length (m)
s_total <- max(s_points)

# Create spline interpolation functions for x(s) and y(s)
x_spline <- splinefun(s_points, x_points, method = "natural")
y_spline <- splinefun(s_points, y_points, method = "natural")

# Derivative functions (1st and 2nd derivatives)
x_prime  <- function(s) x_spline(s, deriv = 1)
x_dprime <- function(s) x_spline(s, deriv = 2)
y_prime  <- function(s) y_spline(s, deriv = 1)
y_dprime <- function(s) y_spline(s, deriv = 2)

# Curvature function kappa(s) = |x'y'' - y'x''| / (x'^2 + y'^2)^(3/2)
curvature <- function(s) {
    num <- x_prime(s) * y_dprime(s) - y_prime(s) * x_dprime(s)
    den <- (x_prime(s)^2 + y_prime(s)^2)^(3/2)
    abs(num) / den
}

# Energy conservation: total_energy = kinetic + potential at start
y0 <- y_points[1]
total_energy <- 0.5 * v0^2 + g * y0  # mass cancels out

# Initialize vectors to store results
time_vec <- numeric()
s_vec <- numeric()
x_vec <- numeric()
y_vec <- numeric()
v_vec <- numeric()
Fg_vec <- numeric()
Fc_vec <- numeric()

# Simulation loop
s <- 0       # initial path length coordinate
t <- 0

while (t < t_max && s <= s_total) {
    # Position
    x <- x_spline(s)
    y <- y_spline(s)
    
    # Calculate speed from energy conservation:
    # 0.5 v^2 + g*y = total_energy  => v = sqrt(2*(total_energy - g*y))
    v_sq <- 2 * (total_energy - g * y)
    # set velocity to zero if train “stops” climbing
    if (v_sq < 0) {
        v <- 0
        # optionally break to stop simulation here
        print(paste0("Simulation stopped at t=", t, " for losing all speed"))
        break
    } else {
        v <- sqrt(v_sq)
    }
    
    # Calculate radius of curvature r = 1/kappa
    k <- curvature(s)
    if (length(k) == 0 || is.na(k) || k == 0) {
        r <- Inf
    } else {
        r <- 1 / k
    }
    
    # Forces magnitudes (mass cancels)
    Fg <- g          # constant gravity force magnitude
    Fc <- v^2 / r    # centrifugal force magnitude
    
    # Store results
    time_vec <- c(time_vec, t)
    s_vec <- c(s_vec, s)
    x_vec <- c(x_vec, x)
    y_vec <- c(y_vec, y)
    v_vec <- c(v_vec, v)
    Fg_vec <- c(Fg_vec, Fg)
    Fc_vec <- c(Fc_vec, Fc)
    
    # Advance s by ds = v * dt (distance along path)
    s <- s + v * dt
    
    # Advance time
    t <- t + dt
}

# Combine results into a data frame
results <- data.frame(time = time_vec,
                      s = s_vec,
                      x = x_vec,
                      y = y_vec,
                      v = v_vec,
                      Fg = Fg_vec,
                      Fc = Fc_vec)


################################################

# 4. CALCULATE VECTOR COMPONENTS OF Fg, Fc -> Fsum

# Calculate unit tangent vector (dx/ds, dy/ds)
tangent_x <- x_prime(results$s)
tangent_y <- y_prime(results$s)
tangent_length <- sqrt(tangent_x^2 + tangent_y^2)
tangent_x <- tangent_x / tangent_length
tangent_y <- tangent_y / tangent_length

# Calculate unit normal vector by rotating tangent 90 degrees CCW:
normal_x <-  tangent_y
normal_y <- -tangent_x

# Force vectors components

# Gravity force (Fg)
Fg_x <- rep(0, length(results$x))
Fg_y <- -results$Fg

# Centrifugal force (Fc)
Fc_x <- normal_x * results$Fc
Fc_y <- normal_y * results$Fc

# Total G force (Fsum)
Fsum_x <- Fg_x + Fc_x
Fsum_y <- Fg_y + Fc_y
results$Fsum <- sqrt(Fsum_x^2 + Fsum_y^2)  # total G force experienced


################################################

# 5. DISPLAY AND PLOT RESULTS

# Show simulation parameters and results
text=paste0("SIMULATION PARAMETERS AND RESULTS:\n\n",
            "Loop height: ", round(H,1), "m (min initial speed: ", round(v0min, 1), "m/s = ", round(v0min/1000*3600,1),"km/h)\n",
            "Initial speed: ", round(v0, 1), "m/s = ", round(v0/1000*3600,1),"km/h (max loop height: ", round(Hmax, 1), "m)\n",
            "Total track length: ", round(s_total, 1), "m\n",
            "Total simulation time: ", round(max(results$time), 1), "s\n",
            "[min, max] speed reached: [", round(min(results$v), 1), ", ",
            round(max(results$v), 1), "] m/s\n",
            "[min, max] acceleration experienced: [", round(min(results$Fsum)/g, 1), ", ",
            round(max(results$Fsum)/g, 1), "] G's\n")
cat(text)


# Plot speed and G force vs time
CairoPNG(paste0("ConstantG_", round(H), "m_", round(v0), "ms_Plots.png"),
         width=512, height=800)
    par(mfrow = c(2, 1))
    
    # Plot speed
    plot(results$time, results$v, ylim=c(0,max(results$v)), type='l', col='red',
         xlab="Time (s)", ylab="Speed (m/s)")
    
    # Plot G force
    plot(results$time, results$Fsum/g, ylim=c(0,7), type='l', col='red',
         xlab="Time (s)", ylab="Total experienced force (G's)")
    abline(h=c(1,6), col="gray", lty="dotted")
    
    par(mfrow = c(1, 1))
dev.off()


# Plot graph of rollercoaster with all vector G forces

# Scale factor for plotting
scale_factor <- 0.2
Fg_x_plot <- Fg_x * scale_factor
Fg_y_plot <- Fg_y * scale_factor
Fc_x_plot <- Fc_x * scale_factor
Fc_y_plot <- Fc_y * scale_factor
Fsum_x_plot <- Fsum_x * scale_factor
Fsum_y_plot <- Fsum_y * scale_factor

# Data frames for arrows
arrows_Fg <- data.frame(
    x = results$x, y = results$y,
    xend = results$x + Fg_x_plot, yend = results$y + Fg_y_plot,
    Force = "Gravity"
)

arrows_Fc <- data.frame(
    x = results$x, y = results$y,
    xend = results$x + Fc_x_plot, yend = results$y + Fc_y_plot,
    Force = "Centrifugal"
)

arrows_Fsum <- data.frame(
    x = results$x, y = results$y,
    xend = results$x + Fsum_x_plot, yend = results$y + Fsum_y_plot,
    Force = "Total"
)

# Combine Fg and Fc into one for dim gray arrows
arrows_individual <- rbind(arrows_Fg, arrows_Fc)

x_mid=11
y_mid=5
CairoPNG(paste0("ConstantG_", round(H), "m_", round(v0), "ms_Graph.png"),
         width=1080, height=1080)
    ggplot() +
        geom_path(data = results, aes(x = x, y = y), color = "black", linewidth = 1) +
        geom_segment(data = arrows_individual,
                     aes(x = x, y = y, xend = xend, yend = yend),
                     arrow = arrow(length = unit(0.15, "cm")),
                     color = alpha("gray70", 0.6),
                     size = 0.3) +
        geom_segment(data = arrows_Fsum,
                     aes(x = x, y = y, xend = xend, yend = yend),
                     arrow = arrow(length = unit(0.3, "cm")),
                     color = "red",
                     size = 0.4) +
        annotate("text", x = x_mid, y = y_mid, label = text, color = "blue", size = 5, fontface = "bold") +
        coord_equal() +
        labs(title = "'constant G' rollercoaster simulation",
             x = "X (m)", y = "Y (m)") +
        theme_minimal() +
        xlim(0, 50) + ylim(-15, 35)
dev.off()




# Creative blueprint plotting
CairoPNG("beautiful_rollercoaster.png", width=1080*3, height=1080*3)
    ggplot() +
        geom_path(
            data = results,
            aes(x = x, y = y),
            color = "black",
            linewidth = 4
        ) +
        geom_segment(
            data = arrows_Fsum,
            aes(x = x, y = y, xend = xend, yend = yend),
            arrow = arrow(length = unit(1.5, "cm")),
            color = "gray",
            linewidth = 2
        ) +
        coord_equal() +
        theme_void()  # removes axes, gridlines, and background
dev.off()





