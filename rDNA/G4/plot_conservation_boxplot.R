#!/usr/bin/env Rscript

# Plot the relationship between predicted G-quadruplex (G4) stability 
# and sequence conservation across maternal and paternal rDNA assemblies
#
# Produces three boxplots:
# 1) all strands
# 2) + strand only
# 3) - strand only
#
# Conservation score is shown on the y-axis.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
    stop("Usage: Rscript plot_conservation_boxplot.R input.bed")
}

input_file <- args[1]


# Derive output name
output_file <- sub(
    "\\.bed$",
    ".conservation_boxplot.png",
    input_file
)

figure_title <- sub(
    "\\.consensus.G4_vs_conservation.bed$",
    "",
    input_file
)


if (output_file == input_file) {
    output_file <- paste0(
        input_file,
        ".conservation_boxplot.png"
    )
}


# Read BED file
dat <- read.table(
    input_file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
)


colnames(dat) <- c(
    "seq",
    "start",
    "end",
    "quartile",
    "score",
    "strand",
    "conservation"
)


# Order quartiles
dat$quartile <- factor(
    dat$quartile,
    levels = c(
        "Q1",
        "Q2",
        "Q3",
        "Q4"
    )
)


# Remove missing conservation values
dat <- dat[
    !is.na(dat$conservation),
]


# Function to create boxplots
make_plot <- function(plot_dat, title) {


    # Sample sizes
    counts <- table(plot_dat$quartile)


    q_labels <- paste0(
        names(counts),
        "\n(n=",
        counts,
        ")"
    )


    # Initial boxplot
    boxplot(
        conservation ~ quartile,
        data = plot_dat,

        xaxt = "n",

        xlab = "G4 score quartile",
        ylab = "Conservation score",

        main = title,

        las = 1,

        outline = FALSE,

        col = "grey85",
        border = "#2171b5",

        lwd = 2,

        cex.main = 1.1,
        cex.lab = 1.1,
        cex.axis = 1.1
    )


    # Add background grid
    grid(
        nx = NA,
        ny = NULL,
        col = "grey85",
        lty = 1
    )


    # Redraw boxes over grid
    boxplot(
        conservation ~ quartile,
        data = plot_dat,

        add = TRUE,

        axes = FALSE,

        outline = FALSE,

        col = "grey85",
        border = "#2171b5",

        lwd = 2
    )


    # Custom x-axis
    axis(
        side = 1,

        at = seq_along(counts),

        labels = FALSE,

        las = 1,

        cex.axis = 1.1,
    )

    #plot the quartile labels
    mtext(
    q_labels,

    side = 1,

    at = seq_along(counts),

    line = 1.3,

    cex = 1.1
)


    # Add jittered individual observations
    set.seed(123)

    x <- as.numeric(plot_dat$quartile)

    x_jitter <- jitter(
        x,
        amount = 0.15
    )


    points(
        x_jitter,

        plot_dat$conservation,

        pch = 16,

        cex = 0.8,

        col = rgb(
            0,
            0,
            0,
            0.25
        )
    )
}



# Create high-resolution PNG
png(
    output_file,

    width = 2400,

    height = 3600,

    res = 300
)


# Plot layout
par(
    mfrow = c(3, 1),

    mar = c(
        6,
        6,
        5,
        2
    ),

    oma = c(
        1,
        1,
        1,
        1
    ),

    cex = 1.1
)



# 1. All strands
make_plot(
    dat,

    paste0("Conservation score by G4 score quartile\n(all strands)\n",figure_title)
)



# 2. + strand
plus_dat <- dat[
    dat$strand == "+",
]


if (nrow(plus_dat) > 0) {

    make_plot(
        plus_dat,

        paste0("Conservation score by G4 score quartile\n(+ strand)\n",figure_title)
    )
}



# 3. - strand
minus_dat <- dat[
    dat$strand == "-",
]


if (nrow(minus_dat) > 0) {

    make_plot(
        minus_dat,

        paste0("Conservation score by G4 score quartile\n(- strand)\n",figure_title)
    )
}



dev.off()


cat(
    "Output written to:",
    output_file,
    "\n"
)