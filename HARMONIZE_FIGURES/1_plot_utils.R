#MANY THANKS TO ALL AUTHORS SHARING R-CODE ONLINE
library(ggplot2)
library(scales)
library(cowplot)
library(wesanderson)

############
#PLOTS
############
#plot default aes
meg_stef_theme <- function () { 
  theme_linedraw(base_size=10) %+replace% 
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA), 
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA),
      legend.title=element_text(size=10, hjust=0),
      #legend.title=element_blank(),
     # axis.ticks = element_blank(),
      panel.grid.major = element_line(color = "grey90", size = 0.3), 
      panel.grid.minor = element_blank(),
      plot.title.position = "panel",
      plot.title = element_text(size = 10, hjust = 0, vjust = 0.5, 
                                margin = margin(b = 0.2, unit = "cm")),
      plot.subtitle = element_text(size = 10, hjust = 0, vjust = 0.5, 
                                   margin = margin(b = 0.4, unit = "cm")),
      plot.caption = element_text(size = 7, hjust = 1, face = "italic", 
                                  margin = margin(t = 0.1, unit = "cm")),
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9,hjust = 1)
    )
}

# labelled axis ticks for discrete values (boxplot)
every_nth = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}


heat_theme <- function () { 
  theme_minimal(base_size=10) %+replace% 
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA), 
      axis.ticks=element_blank(),
      axis.text=element_blank(),
      axis.line=element_blank(),
      axis.title.y=element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      plot.title.position = "panel",
      plot.title = element_text(size = 10, hjust = 0, vjust = 0.5, 
                                margin = margin(b = 0.2, unit = "cm")),
      plot.subtitle = element_text(size = 10, hjust = 0, vjust = 0.5, 
                                   margin = margin(b = 0.4, unit = "cm")),
      plot.caption = element_text(size = 7, hjust = 1, face = "italic", 
                                  margin = margin(t = 0.1, unit = "cm")),
      plot.margin = margin(t = 0,r = 0,b = 0, l = 20) 
   
    )
}


########
#COLOR PALETTES
########
#extracts hex
mycols <- function(...) {
  cols <- c(...)
  if (is.null(cols))
  return (mycolors)
  mycolors[cols]
}


#colorRamp calculates continuous colors
paper_pal <- function(palette = "Paper", reverse = FALSE, ...) {
  pal <- mypalettes[[palette]]
  if (reverse) pal <- rev(pal)
  colorRampPalette(pal, ...)
}

#customize color theme for ggplot2
scale_color_paper <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- paper_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("colour", paste0("paper_", palette), palette = pal, ...)
  } else {
    scale_color_gradientn(colours = pal(256), ...)
  }
}

scale_fill_paper <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- paper_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("fill", paste0("paper_", palette), palette = pal, ...)
  } else {
    scale_fill_gradientn(colours = pal(256), ...)
  }
}

#paper colors
mycolors<-c(`Byellow`="#F1BB7B", `Bpink`="#FD6467", `Bplum`="#5B1A18",  `Mlblue`="#C6CDF7",`Mblue`="#7294D4", `Cgreen`="#81A88D",`Mgrey`="#798E87",`Mwhite`="#D5D5D3",`Boldpink`="#D8A499", `Tgreen`="#02401B", `white`="#FFFFFF", `grey`="#cccccc")

#paper palettes
mypalettes <- list(
  `Paper` = mycols("Byellow", "Bpink", "Bplum",  "Mlblue","Mblue", "Cgreen","Tgreen","Mgrey","Mwhite","Boldpink"),
  `CF` = mycols("Bpink","Bplum","Byellow"),
  `Noise` = mycols("Mwhite","Boldpink"),
  `Time` = mycols("Mblue","Tgreen"),
  `Age` = mycols("Mblue","grey","Tgreen"),
  `DMM` = mycols("Byellow","Bplum", "Mlblue","Mblue"),
  `Bi2` = mycols("Bpink", "Mblue")
)


#display colors
#show_col(paper_pal("DMM")(6))
#show_col(mypalettes$DMM)
#get hex
# paper_pal("DMM")(6) 

#outliers
out<-c("1141","1149","23024","1156","1131","1105")