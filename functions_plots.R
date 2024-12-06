############################## Plots functions ######################################
# Densities plots


# Plot of lines with diferent colours -----------------------------------------
# x: Name of variable. No characters # Plots lines
# y: Name of variable. No characters
# colour: Name of variable. No character
# x_lab: As character. Title of x axis.
# y_lab: As character. Title of y axis.
# colour_lab: As character. Name of variable for title legend.
# breaks: Vector with breaks for x axis. This vector contains the values of x axis.

plots <- function(filtered_data, x, y, colour, x_lab, y_lab, colour_lab, breaks){
    library(ggplot2)
    ggplot(filtered_data, aes(x = x, y = y, color = as.factor(colour), shape = as.factor(colour))) +
        geom_line() + geom_point() + scale_color_brewer(palette = "Dark2") +
        scale_fill_brewer(palette = "Dark2") + theme_classic2() + xlab(x_lab) + ylab(y_lab) +
        labs(color = colour_lab, shape = colour_lab) +
        scale_x_continuous(breaks = breaks)
}

# Plots side by side ----------------------------------------------------------
# With grid and labels

