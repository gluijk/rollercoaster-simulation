# Simulating G forces on different rollercoaster loop topologies
# www.overfitting.net
# https://www.overfitting.net/2025/08/simulando-loops-de-una-montana-rusa-con_69.html


library(pracma)  # fresnelC() and fresnelS() for clothoid
# library(splines)
library(ggplot2)
library(Cairo)


# ChatGPT 5 prompt:
# -----------------
# Write R code to simulate the forces affecting a train on a rollercoaster:
# * The rollercoaster is defined by discrete points on two vectors
#   x (ground axis) and y (height axis), in the way: x=c(0, 1, 2, 3, 4, 5,...)
#   and y=c(0, 0, 0.5, 1, 2, 3.5, 5,...)
# * The simulation must be done time dependant in dt constant time steps.
#   This means you'll have to linearly interpolate along the definitory points
#   of the rollercoaster, not just give a force value for each rollercoaster
#   point, calculating the advance of the train along the path
# * The train is defined by a punctual mass m (but this should be irrelevant
#   since it cancels on all equations) starting from point (0,0) at initial speed
#   v0 left to right
# * You must calculate the sum of two forces: Fg=m*g and Fc=m*v^2/r, where Fg
#   is the gravity force and always points down and Fc is the centrifugal force
#   trying to evade the rollercoaster in a perpendicular direction.
# * Apply energy conservation: 1/2*m*v0^2=m*g*H to obtain each v speed
# 
# Do you need something else before beginning?



################################################

# 1. SIMULATION PARAMETERS

g <- 9.80665    # gravity (m/s^2)
v0 <- 28        # initial speed (m/s)
t_max <- 999    # max simulation time (s)
dt <- 0.05      # time step (s)
m <- 1          # mass (kg), cancels out



################################################

# 2. DEFINE CIRCULAR / CLOTHOID ROLLERCOASTER


NRESOL=5000
for (type in c('circular', 'clothoid')) {
    # Define discrete track points
    if (type=='circular') {
        R=15
        H=2*R
        
        # Add initial flat track
        Offset=25
        par=seq(0, 3/4*2*pi, length.out=NRESOL)
        x_points=c(seq(0, Offset*0.9, length.out=10), R*sin(par)+Offset)
        y_points=c(seq(0, 0, length.out=10), R*(1-cos(par)))
        
        plot(x_points, y_points, asp=1, type='l')
        abline(v=c(0,Offset), col='gray')
    } else if (type=='clothoid') {
        L=90   # clothoid total length (m)
        s=seq(0, L, length.out=NRESOL)
        Rmin=6.244394152  # Rmin for max height=30m
        # Scale factor so curvature radius at s = L equals Rmin
        A=(pi * Rmin * L)^0.5
        
        # Calculate clothoid (x,y) positions
        u=s/A
        x=A*fresnelC(u)
        y=A*fresnelS(u)
        
        # Truncate and simmetrize clothoid
        H=max(y)
        pos=which(y==H)
        x=x[1:pos]
        y=y[1:pos]
        
        # Plot clothoid
        plot(x, y, col='black', lwd=2, type="l", asp=1, xlim=c(-30,45), ylim=c(0,50),  # NULL plot
             xlab="", ylab="")
        abline(h=c(0,30), v=x[pos], col='gray',lty='dotted')
        
        SIME=1000
        x=c(x, -x[(pos-1):(pos-SIME)]+x[pos]*2)
        y=c(y,  y[(pos-1):(pos-SIME)])
        
        # Add initial flat track
        Offset=10
        x_points=c(seq(0, Offset*0.9, length.out=10), x+Offset)
        y_points=c(seq(0, 0, length.out=10), y)
        
        
        plot(x_points, y_points, col='black', lwd=2, type="l",
             asp=1, xlim=c(-30,45), ylim=c(0,50),  # NULL plot
             xlab="", ylab="")
        abline(h=c(0,30), v=x[pos]+Offset, col='gray',lty='dotted')
    }
    
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
    CairoPNG(paste0(type, "_", round(H), "m_", round(v0), "ms_Plots.png"),
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
    scale_factor <- 0.2  # 0.1
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
    
    x_mid=11  # <- mean(range(results$x))*1.2
    y_mid=5  # <- mean(range(results$y))
    CairoPNG(paste0(type, "_", round(H), "m_", round(v0), "ms_Graph.png"),
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
            labs(title = paste0("'", type, "' rollercoaster simulation"),
                 x = "X (m)", y = "Y (m)") +
            theme_minimal() +
            xlim(0, 50) + ylim(-15, 35)
    dev.off()

}




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



