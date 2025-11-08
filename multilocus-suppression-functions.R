# Load in (and install, if necessary) all necessary packages
packages = c("tidyverse", "cowplot", "viridis","latex2exp", "here", "openxlsx", "readr", "purrr", "RColorBrewer", "bit64",
             "colorspace") 
installed = packages %in% installed.packages()[, "Package"]
if (any(!installed)) {
  install.packages(packages[!installed])
}
lapply(packages, library, character.only = TRUE)

#######################################################

### Useful helper functions ###

calculate_text_size = function(artboard_width, overleaf_scaling,desired_text_size = 10,latex_width = 500.484){
  # Calculates the necessary text size for the plot
  # Arguments:
  #     artboard_width: Adobe Illustrator artboard width
  #     overleaf_scaling: the scaling parameter on the figure in Overleaf
  #     desired_text_size: what size you want the text to be in the paper
  #     latex_width: Overleaf text width
  # Returns:
  #     text_size: size of text for everything in the plot
  #     subplot_text_size: size of the subplot labels (A, B, C, D)
  
  if (latex_width*overleaf_scaling > artboard_width){
    p = 1
  } else {
    p = (latex_width*overleaf_scaling)/artboard_width  
  }
  
  t = desired_text_size/p
  subplot.t = (desired_text_size + 2)/p
  return(list(text_size = t, subplot_text_size = subplot.t))
}

st_for_desired_load = function(desired_load, s_drive, n_targets){
  # Finds the necessary fitness cost at the individual targets
  # Arguments:
  #   desired_load: the lambda^* you're targeting
  #   s_drive: fitness cost at the drive sites
  #   n_targets: the number of targets
  # Returns:
  #   the disrupted fitness cost, s
  
  # desired load = 1 - (1-s_drive)(1-st)^n_targets
  return(1 - ((1-desired_load)/(1-s_drive))^(1/n_targets))
  
}

#######################################################

### Function for heatmaps ### 

heatmap_detailed = function(csv, 
                            xindex, yindex, zindex,
                            x_tick_marks, y_tick_marks,z_tick_marks,
                            zlow, zhigh,
                            color_scale=c("white","blue"),
                            xl = "", yl = "", zl = "",title = "", 
                            text_size = 12,
                            legend.pos = "bottom",
                            margin = F){
  
  # Arguments:
  #   csv: a tibble containing the x, y, and z (color) variable
  #   xindex: the column index of the x-axis variable in the csv
  #   yindex: the column index of the y-axis variable in the csv
  #   zindex: the column index of the z-axis (color axis) variable in the csv
  #   x_tick_marks: vector of tick-mark values for the x-axis
  #   y_tick_marks: vector of tick-mark values for the y-axis
  #   z_tick_marks: vector of tick-mark values for the colorbar
  #   zlow: minimum color value
  #   zhigh: maximum color value
  #   color_scale: vector of colors to use for the min and max of the z-variable range
  #   xl: x-axis label
  #   yl: y-axis label
  #   zl: z-axis label
  #   title: plot title
  #   text_size: size of all text
  #   legend.pos: position of the colorbar ("bottom", "right", or "none")
  #   margin: whether to include padding around the plot 
  
  indices = c(xindex, yindex, zindex)
  
  # focus on data of interest
  df = as.data.frame(csv[,indices])
  rng = range(c(0,1))
  class(df[,1]) = "double"
  class(df[,2]) = "double"
  class(df[,3]) = "double"
  names(df) = c('x','y','v')
  df = aggregate(v~x+y, data=df, mean)
  xs = sort(unique(df$x))
  ys = sort(unique(df$y))
  x_range = xs[2]-xs[1]
  y_range = ys[2]-ys[1]
  
  
  first = ggplot(df, aes(x, y)) +
    geom_tile(aes(fill = v)) +
    coord_fixed(ratio=1*x_range/y_range) +
    labs(x=xl, y=yl,title=title, fill = zl) 
  
  
  
  if (legend.pos %in% c("bottom", "none")){
    
    plot = first + scale_x_continuous(breaks = x_tick_marks) +
      scale_y_continuous(breaks = y_tick_marks)+
      scale_fill_gradientn(colours = color_scale,
                           breaks = z_tick_marks,
                           limits = c(zlow, zhigh),
                           guide = guide_colourbar(direction="horizontal",
                                                   label.position = "bottom",
                                                   barheight = unit(0.03, "npc"),  # sets colorbar height as fraction of plot's
                                                   barwidth = unit(0.45, "npc"))) + # sets colorbar width as fraction of plot's
      theme(legend.position=legend.pos,
            text = element_text(size=text_size),
            plot.title = element_text(size=text_size, face = "bold"),
            axis.text = element_text(size=text_size),
            legend.text = element_text(size=text_size))
    
  } else if (legend.pos == "right") {
    
    plot = first + scale_x_continuous(breaks = x_tick_marks) +
      scale_y_continuous(breaks = y_tick_marks)+
      scale_fill_gradientn(colours = color_scale,
                           breaks = z_tick_marks,
                           limits = c(zlow, zhigh),
                           guide = guide_colourbar(direction="vertical",
                                                   label.position = "left",
                                                   barwidth  = unit(0.03, "npc"),
                                                   barheight = unit(0.45, "npc"))) +
      theme(legend.position=legend.pos,
            text = element_text(size=text_size),
            plot.title = element_text(size=text_size, face = "bold"),
            axis.text = element_text(color="black",size=text_size),
            legend.text = element_text(size=text_size))
    
  }
  
  if (zl!=""){
    if (legend.pos == "right"){
      plot = plot + theme(legend.title = element_text(hjust = 1))
    } else if (legend.pos == "bottom"){
      plot = plot + theme(legend.title = element_text(vjust = 1))
    }
  }
  
  if (!margin){
    plot = plot + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
  }
  return(plot)
} 


#############################################

###### Single locus drive plotting functions ####

summarize_single_locus_directory = function(directory,all_gens=NULL, summarize = TRUE){
  # Reads in all of the replicate csv files in a directory, fills them (down from the final absorbing state) to a max generation,
  # and averages across generations
  
  # Arguments:
  #   directory: path to the files (don't end it in /)
  #           note -- ensure that all of these files have the same ending state (ie all ended in fixation or all ended in lost) 
  #           such that averaging is not misleading
  #   all_gens: vector of generations to ensure each replicate has (ex: 1:100)
  #   summarize: whether to average the trajectories
  # Returns:
  #   if !summarize, just returns the final_states. This shows the last row in each replicate.
  #   else returns the summary tibble (averaged values across generations) and the final states
  
  # --- Read all CSV files ---
  files = list.files(path = directory, pattern = "\\.csv$", full.names = TRUE)
  
  data_list <- map2(
    files,
    seq_along(files),
    ~ read_csv(.x, show_col_types = FALSE) %>%
      mutate(
        replicate = .y,
        filename = basename(.x)      # Save the file name
      )
  )
  
  all_data = bind_rows(data_list)
  
  # need a point for each
  if (is.null(all_gens)){
    all_gens = sort(unique(all_data$gen_num))  
  }
  
  final_states = all_data %>%
    group_by(replicate,filename) %>%
    filter(gen_num == max(gen_num, na.rm = TRUE)) %>%
    summarize(
      final_genetic_load = last(genetic_load),
      final_drive_rate = last(drive_rate),
      final_r1_drive_rate = last(resist_rate),
      .groups = "drop"
    )
  
  n = sum(final_states$final_drive_rate)
  all_ended_at_fixation = near(n, nrow(final_states))
  all_ended_at_loss = near(n, 0)
  if (!all_ended_at_fixation & !all_ended_at_loss){
    print("SEPARATE FINAL STATES")
  }
  
  if (!summarize){
    return(final_states)
  }
  
  
  # Extend for absorbing states
  filled_data = all_data %>%
    group_by(replicate) %>%
    arrange(gen_num) %>%
    # filter out generations above max
    filter(gen_num <= max(all_gens)) %>%
    # extend to all_gens, filling missing gens
    complete(gen_num = all_gens) %>%
    # carry forward the last observed values (which should be absorbing states)
    fill(genetic_load, drive_rate, resist_rate, mean_fit, .direction = "down") %>%
    mutate(
      across(
        c(genetic_load, drive_rate, resist_rate),
        ~ ifelse((near(lag(.x), 0) | near(lag(.x), 1)) & is.na(.x), lag(.x), .x)
      )
    ) %>%
    ungroup()
  
  summary = filled_data %>%
    group_by(gen_num) %>%
    summarise(
      avg_genetic_load = mean(genetic_load, na.rm = TRUE),
      min_genetic_load = min(genetic_load, na.rm = TRUE),
      max_genetic_load = max(genetic_load, na.rm = TRUE),
      avg_drive_rate = mean(drive_rate, na.rm = TRUE),
      min_drive_rate = min(drive_rate, na.rm = TRUE),
      max_drive_rate = max(drive_rate, na.rm = TRUE),
      avg_r1_drive_rate = mean(resist_rate, na.rm = TRUE),
      min_r1_drive_rate = min(resist_rate, na.rm = TRUE),
      max_r1_drive_rate = max(resist_rate, na.rm = TRUE),
      n_replicates = n()
    ) %>%
    ungroup() 
  
  return(list(summary = summary,
              final_states = final_states)) 
}


plot_single_locus_joint_gl_joint_res = function(resistance_rate,
                                                max_gen = NULL,
                                                legend.pos = "right",
                                                text_size = 14, 
                                                linewidth = 1,
                                                dsx_color = "mediumslateblue",
                                                s0.2_color = "green4") {
  
  # Makes the single locus main plots -- the averaged genetic loads for each drive and the averaged frequencies
  # Arguments:
  #     resistance_rate: the resistance rate at the drive site (1 of 1e-5, 1e-4, or 1e-3)
  #     max_gen: the max x-value of the plot
  #     legend.pos: 'right' for legend or 'none' for no legend
  #     text_size: size of all text in plot
  #     linewidth: width of trajectory and horizontal lines
  #     dsx_color: the color of the s_d = 1 drive
  #     s0.2_color: the color of the modification drive, with lower s
  # Returns:
  #     List of the genetic load averaged plot, resistance frequency averaged plot, and drive frequency averaged plot
  #     the transparency of a line reflects the % of replicates where it ended in that state
  
  
  # Set up vector for reading in the replicate data
  if (!is.null(max_gen)){
    all_gens = 1:max_gen
  } else {
    all_gens = NULL
  }
  
  
  dsx.dir.base = here(paste0("data/single-locus/dsx/"))
  s0.2.dir.base = here(paste0("data/single-locus/s_d0.2/"))
  if (near(resistance_rate, 1e-3)){
    # both dsx and s = 0.5 were always lost
    # 2 categories: s = c(1, 0.2)
    dsx.dir = paste0(dsx.dir.base, "res1e-3")
    s0.2.dir = paste0(s0.2.dir.base, "res1e-3")
    
    dsx.results = summarize_single_locus_directory(dsx.dir,all_gens)$summary %>%
      mutate(s_d = "1",
             alpha_val = n_replicates/20,
             group = 1)
    
    s0.2.results = summarize_single_locus_directory(s0.2.dir,all_gens)$summary %>%
      mutate(s_d = "0.2",
             alpha_val = n_replicates/20,
             group = 2)
    
    all_results = bind_rows(
      dsx.results,
      s0.2.results
    )
    
  } else if (near(resistance_rate, 1e-4)){
    # dsx always lost
    # s = 0.2 fixed in 12/20 replicates; lost in 8/20
    # 3 categories: s = c(1, 0.2_fixed, 0.2_lost)
    dsx.dir = paste0(dsx.dir.base, "res1e-4")
    
    dsx.results = summarize_single_locus_directory(dsx.dir,all_gens)$summary %>%
      mutate(s_d = "1",
             alpha_val = n_replicates/20,
             group = 1) # drive loss
    
    s0.2.dir.frequent = paste0(s0.2.dir.base, "res1e-4/fixed")
    s0.2.frequent.results = summarize_single_locus_directory(s0.2.dir.frequent,all_gens)$summary %>% 
      mutate(s_d = "0.2",
             alpha_val = n_replicates/20,
             group = 2)
    
    s0.2.dir.less.frequent = paste0(s0.2.dir.base, "res1e-4/lost")
    s0.2.less.frequent.results = summarize_single_locus_directory(s0.2.dir.less.frequent,all_gens)$summary %>% 
      mutate(s_d = "0.2",
             alpha_val = n_replicates/20,
             group = 3)
    
    all_results = bind_rows(
      dsx.results,
      s0.2.frequent.results,
      s0.2.less.frequent.results
    )
    
    
  } else if (near(resistance_rate, 1e-5)) {
    # dsx fixed in 10/20
    # dsx lost in 10/20
    # s = 0.2 fixed in 18/20
    # s = 0.2 lost in 2/20
    # 4 categories: s = c(1_fixed, 1_lost, 0.2_fixed, 0.2_lost)
    
    dsx.dir.frequent = paste0(dsx.dir.base, "res1e-5/lost")
    dsx.frequent.results = summarize_single_locus_directory(dsx.dir.frequent, all_gens)$summary %>%
      mutate(s_d = "1",
             alpha_val = n_replicates/20,
             group = 1)
    
    dsx.dir.less.frequent = paste0(dsx.dir.base, "res1e-5/fixed")
    dsx.less.frequent.results = summarize_single_locus_directory(dsx.dir.less.frequent,all_gens)$summary %>%
      mutate(s_d = "1",
             alpha_val = n_replicates/20,
             group = 2) 
    
    s0.2.dir.frequent = paste0(s0.2.dir.base, "res1e-5/fixed")
    s0.2.frequent.results = summarize_single_locus_directory(s0.2.dir.frequent,all_gens)$summary %>% 
      mutate(s_d = "0.2",
             alpha_val = n_replicates/20,
             group = 3) 
    
    s0.2.dir.less.frequent = paste0(s0.2.dir.base, "res1e-5/lost")
    s0.2.less.frequent.results = summarize_single_locus_directory(s0.2.dir.less.frequent,all_gens)$summary %>%
      mutate(s_d = "0.2",
             alpha_val = n_replicates/20, 
             group = 4) 
    
    all_results = bind_rows(
      dsx.frequent.results,
      dsx.less.frequent.results,
      s0.2.frequent.results,
      s0.2.less.frequent.results
    )
    
  }
  
  
  plt.gl = ggplot(all_results, aes(x = gen_num, y = avg_genetic_load, color = s_d, group = interaction(s_d, alpha_val, group)))  +
    geom_line(aes(alpha = alpha_val), linewidth = linewidth) +      # alpha varies by frequent/less frequent
    scale_color_manual(values = c("1" = dsx_color, "0.2" = s0.2_color)) +
    scale_alpha_identity()  + 
    geom_hline(yintercept = 1, linetype = "dashed", color = dsx_color, linewidth = linewidth*1.05)+
    geom_hline(yintercept = 0.2, linetype = "dashed", color = s0.2_color, linewidth = linewidth*1.05)+
    labs(y = TeX(r"(Genetic load ($\lambda$))"),
         x = TeX(r"(Generation)"),
         color = TeX(r"($s_d$)")) +
    theme(
      legend.position = legend.pos,
      text = element_text(size = text_size),
      plot.title = element_text(size = text_size, face = "bold"),
      axis.text = element_text(size = text_size),
      legend.text = element_text(size = text_size),
      aspect.ratio = 1
    ) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0,1))
  
  plt.res = ggplot(all_results, aes(x = gen_num, y = avg_r1_drive_rate, color = s_d, group = interaction(s_d, alpha_val, group))) +
    geom_line(aes(alpha = alpha_val), linewidth = linewidth) +      # alpha varies by frequent/less frequent
    scale_color_manual(values = c("1" = dsx_color, "0.2" = s0.2_color)) +
    scale_alpha_identity()  + 
    labs(y = TeX(r"(Resistance frequency)"),
         x = TeX(r"(Generation)"),
         color = TeX(r"($s_d$)")) +
    theme(
      legend.position = legend.pos,
      text = element_text(size = text_size),
      plot.title = element_text(size = text_size, face = "bold"),
      axis.text = element_text(size = text_size),
      legend.text = element_text(size = text_size),
      aspect.ratio = 1
    ) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0,1))
  
  
  plt.dr = ggplot(all_results, aes(x = gen_num, y = avg_drive_rate, color = s_d, group = interaction(s_d, alpha_val, group))) +
    geom_line(aes(alpha = alpha_val), linewidth = linewidth) +      # alpha varies by frequent/less frequent
    scale_color_manual(values = c("1" = dsx_color, "0.2" = s0.2_color)) +
    scale_alpha_identity()  + 
    labs(y = TeX(r"(Drive frequency)"),
         x = TeX(r"(Generation)"),
         color = TeX(r"($s_d$)")) +
    theme(
      legend.position = legend.pos,
      text = element_text(size = text_size),
      plot.title = element_text(size = text_size, face = "bold"),
      axis.text = element_text(size = text_size),
      legend.text = element_text(size = text_size),
      aspect.ratio = 1
    ) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0,1))
  
  
  plot.list = list(gl = plt.gl, res = plt.res, drive = plt.dr)
  return(plot.list)
}

#############################################

### Functions for multilocus drive plotting ### 

summarize_multilocus_directory = function(directory, all_gens = NULL, summarize = T) {
  
  # Arguments:
  #   directory: path to replicate data (not ending in /)
  #         make sure these are separated by absorbing state (ie all lost or all fixed)
  #   all_gens: vector of generations to summarize (1:max_gen)
  #   summarize: whether to create an averaged trajectory
  # Returns:
  #   if !summarize, just returns the final states (last rows of each replicate)
  #   else returns the final states and the replicate-averaged summary

  files = list.files(path = directory, pattern = "\\.csv$", full.names = TRUE)
  
  data_list = map2(
    files,
    seq_along(files),
    ~ {
      df = read_csv(.x, show_col_types = FALSE)
      
      # Generation where genetic load < 20%
      gen_marker_gen = df %>%
        filter(gen_marker == 1) %>%
        pull(gen_num) %>%
        first()
      
      # Generation where drive reaches fixation or lost
      drive_absorbed_gen = df %>%
        filter(near(drive_rate,0) | near(drive_rate, 1)) %>%
        pull(gen_num) %>%
        first()
      
      # Generation where targets reach fixation or lost
      disrupted_absorbed_gen = df %>%
        filter((near(broken_rate,0) | near(broken_rate, 1)) & gen_num > 1) %>%
        pull(gen_num) %>%
        first()
      
      df %>%
        mutate(
          replicate = .y,
          filename = basename(.x),
          gen_genetic_load_below20 = gen_marker_gen,
          gen_drive_fixed_or_lost = drive_absorbed_gen,
          gen_disrupted_fixed_or_lost = disrupted_absorbed_gen
        )
    }
  )
  
  all_data = bind_rows(data_list)
  
  # need a point for each
  if (is.null(all_gens)){
    all_gens = sort(unique(all_data$gen_num))  
  }
  
  final_states = all_data %>%
    group_by(replicate,filename) %>%
    filter(gen_num == max(gen_num, na.rm = TRUE)) %>%
    summarize(
      final_genetic_load = last(genetic_load),
      final_drive_rate = last(drive_rate),
      final_r1_drive_rate = last(resist_rate),
      final_disrupted_rate = last(broken_rate),
      final_r1_target_rate = last(func_rate),
      gen_genetic_load_below20 = last(gen_genetic_load_below20),
      gen_drive_fixed_or_lost = last(gen_drive_fixed_or_lost),
      gen_disrupted_fixed_or_lost = last(gen_disrupted_fixed_or_lost),
      .groups = "drop"
    )
  
  n = sum(final_states$final_disrupted_rate)
  all_ended_at_fixation = near(n, nrow(final_states))
  all_ended_at_loss = near(n, 0)
  if (!all_ended_at_fixation & !all_ended_at_loss){
    print("SEPARATE TARGET FINAL STATES")
    print(paste0("target lost in ", n,"/20"))
  }
  
  n = sum(final_states$final_drive_rate)
  all_ended_at_fixation = near(n, nrow(final_states))
  all_ended_at_loss = near(n, 0)
  if (!all_ended_at_fixation & !all_ended_at_loss){
    print("SEPARATE DRIVE FINAL STATES")
    print(paste0("drive fixed in ", n,"/20"))
  }
  
  if (!summarize){
    return(final_states)
  }
  
  
  # Extend for absorbing states
  filled_data = all_data %>%
    group_by(replicate) %>%
    arrange(gen_num) %>%
    # filter out generations above max
    filter(gen_num <= max(all_gens)) %>%
    # extend to all_gens, filling missing gens
    complete(gen_num = all_gens) %>%
    # carry forward the last observed values (which should be absorbing states)
    fill(genetic_load, drive_rate, resist_rate, mean_fit, broken_rate, func_rate,
         gen_genetic_load_below20, gen_drive_fixed_or_lost, gen_disrupted_fixed_or_lost,
         .direction = "down") %>%
    mutate(
      across(
        c(genetic_load, drive_rate, resist_rate, mean_fit, broken_rate, func_rate),
        ~ ifelse((near(lag(.x), 0) | near(lag(.x), 1)) & is.na(.x), lag(.x), .x)
      )
    ) %>%
    ungroup()
  
  summary = filled_data %>%
    group_by(gen_num) %>%
    summarise(
      avg_genetic_load = mean(genetic_load, na.rm = TRUE),
      min_genetic_load = min(genetic_load, na.rm = TRUE),
      max_genetic_load = max(genetic_load, na.rm = TRUE),
      avg_drive_rate = mean(drive_rate, na.rm = TRUE),
      min_drive_rate = min(drive_rate, na.rm = TRUE),
      max_drive_rate = max(drive_rate, na.rm = TRUE),
      avg_r1_drive_rate = mean(resist_rate, na.rm = TRUE),
      avg_disrupted_rate = mean(broken_rate, na.rm = TRUE),
      min_disrupted_rate = min(broken_rate, na.rm = TRUE),
      max_disrupted_rate = max(broken_rate, na.rm = TRUE),
      avg_r1_target_rate = mean(func_rate, na.rm = TRUE),
      n_replicates = n()
    ) %>%
    ungroup() 
  
  
  return(list(final_states = final_states,
              summary = summary))
}


plot_multilocus_genetic_load_summaries = function(base.directory = here("data/n-lambda/cas9-drive1-target0.8/lambda0.95/"),
                                                  max_gen = NULL,
                                                  lambda = 0.95, 
                                                  col.palette = c("saddlebrown", "green4", "mediumslateblue"),
                                                  text_size = 14, legend.pos = "right", linewidth = 1,
                                                  alpha.booster = 1.1, title = "") {
  
  # Creates the replicate-averaged plots for the (n, lambda^*) unlinked multilocus jobs
  # Arguments:
  #     base.directory: the directory where the num_targets<t>/drive-<fixed or lost> files are
  #     max_gen: the max x-value of the plot
  #     lambba: the common max genetic load 
  #     col.palette: the corresponding colors for 1-target, 5-targets, 10-targets
  #     text_size: size of all text in plot
  #     legend.pos: either "right" for legend or "none" for no legend
  #     linewidth: width of all lines
  #     alpha.booster: multiplier on the transparency of the "drive fixed" outcomes 
  #     title: optional title
  # Returns:
  #     3 plots:
  #         combined genetic load trajectories
  #         combined drive and disrupted trajectories
  #         combined disrupted and target-resistant trajectories
                                          
  all_gens = NULL
  if (!is.null(max_gen)){
    all_gens= 1:max_gen
  }
  
  # 1 target
  one.target_drive.fixed = summarize_multilocus_directory(directory = paste0(base.directory, "num_targets1/drive-fixed"), all_gens = all_gens)$summary %>%
    mutate(num_targets = "1",
           alpha_val = (n_replicates/20)*alpha.booster,
           group = 1)
  
  one.target_drive.lost = summarize_multilocus_directory(directory = paste0(base.directory, "num_targets1/drive-lost"), all_gens = all_gens)$summary %>%
    mutate(num_targets = "1",
           alpha_val = n_replicates/20,
           group = 2)
    
  # 5 targets
  five.target_drive.fixed = summarize_multilocus_directory(directory = paste0(base.directory, "num_targets5/drive-fixed"), all_gens = all_gens)$summary %>%
    mutate(num_targets = "5",
           alpha_val = (n_replicates/20)*alpha.booster,
           group = 3)
  
  five.target_drive.lost = summarize_multilocus_directory(directory = paste0(base.directory, "num_targets5/drive-lost"), all_gens = all_gens)$summary %>%
    mutate(num_targets = "5",
           alpha_val = n_replicates/20,
           group = 4)
  
  # 10 targets
  ten.target_drive.fixed = summarize_multilocus_directory(directory = paste0(base.directory, "num_targets10/drive-fixed"), all_gens = all_gens)$summary %>%
    mutate(num_targets = "10",
           alpha_val = (n_replicates/20)*alpha.booster,
           group = 5)
  
  ten.target_drive.lost = summarize_multilocus_directory(directory = paste0(base.directory, "num_targets10/drive-lost"), all_gens = all_gens)$summary %>%
    mutate(num_targets = "10",
           alpha_val = n_replicates/20,
           group = 6)
  
  values.palette = c("1" = col.palette[1],
                     "5" = col.palette[2],
                     "10" = col.palette[3])
  
  
  all_results = bind_rows(one.target_drive.fixed,
                               one.target_drive.lost,
                               five.target_drive.fixed,
                               five.target_drive.lost,
                               ten.target_drive.fixed,
                               ten.target_drive.lost)
  
  p.gl = ggplot(all_results, aes(x = gen_num, y = avg_genetic_load, color = num_targets, group = interaction(num_targets, alpha_val, group)))  +
    geom_line(aes(alpha = alpha_val), linewidth = linewidth) +      # alpha varies by frequent/less frequent
    scale_alpha_identity()  + 
    geom_hline(yintercept = lambda, linetype = "dashed", size = linewidth) +
    labs(
      x = TeX(r"(Generation)"),
      y = TeX(r"(Genetic load ($\lambda$))"),
      color = TeX(r"($n$)"),
      title = title
    ) +
    scale_color_manual(values = values.palette,
                       breaks = c("1","5", "10")) +
    theme(
    legend.position = legend.pos,
    text = element_text(size = text_size),
    plot.title = element_text(size = text_size, hjust = 0.5),
    axis.text = element_text(size = text_size),
    legend.text = element_text(size = text_size),
    aspect.ratio = 1
  ) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0,1)) 

  

  # Target site alleles
  targets.freq = all_results %>%
    pivot_longer(cols = c(avg_r1_target_rate, avg_disrupted_rate),
                 names_to = "allele", values_to = "freq")

  
  p.res = ggplot(targets.freq, aes(x = gen_num, y = freq,
                                  color = num_targets, group = interaction(num_targets, alpha_val, group, allele),
                                  linetype = allele))  +
    geom_line(aes(alpha = alpha_val), linewidth = linewidth) +      # alpha varies by frequent/less frequent
    scale_alpha_identity()  + 
    labs(
      x = TeX(r"(Generation)"),
      y = "Frequency",
      color = TeX(r"($n$)"),
      linetype = "target site\nallele"
    ) +
    scale_linetype_manual(
      values = c(avg_r1_target_rate = "longdash",
                 avg_disrupted_rate = "solid"),
      labels = c(avg_r1_target_rate = "resistance",
                 avg_disrupted_rate = "disrupted")
    ) +
    scale_color_manual(values = values.palette,
                       breaks = c("1","5", "10")) +
    theme(
      legend.position = legend.pos,
      text = element_text(size = text_size),
      axis.text = element_text(size = text_size),
      legend.text = element_text(size = text_size),
      aspect.ratio = 1
    ) +
    guides(linetype = guide_legend(keywidth = 2.2)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0,1))
    
  
  # Drive and disrupted target alleles
  allele.freq = all_results %>%
    pivot_longer(cols = c(avg_drive_rate, avg_disrupted_rate),
                 names_to = "allele", values_to = "freq")
  
  
  p.drive.target = ggplot(allele.freq, aes(x = gen_num, y = freq,
                                   color = num_targets, group = interaction(num_targets, alpha_val, group, allele),
                                   linetype = allele))  +
    geom_line(aes(alpha = alpha_val), linewidth = linewidth) +      # alpha varies by frequent/less frequent
    scale_alpha_identity()  + 
    labs(
      x = TeX(r"(Generation)"),
      y = "Frequency",
      color = TeX(r"($n$)"),
      linetype = "allele"
    ) +
    scale_linetype_manual(
      values = c(avg_disrupted_rate = "longdash",
                 avg_drive_rate = "solid"),
      labels = c(avg_disrupted_rate = "disrupted",
                 avg_drive_rate = "drive")
    ) +
    ylim(0, 1) +
    scale_color_manual(values = values.palette,
                       breaks = c("1","5", "10")) +
    theme(
      legend.position = legend.pos,
      text = element_text(size = text_size),
      axis.text = element_text(size = text_size),
      legend.text = element_text(size = text_size),
      aspect.ratio = 1
    ) +
    guides(linetype = guide_legend(keywidth = 2.25)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0,1)) 
  
  return(list(gl.plot = p.gl,
              target.freq.plot = p.res,
              drive.and.target.plot = p.drive.target))
}


#### Function for plotting log jam
log_jam_summary_plots = function(data, 
                          text_size = 14, 
                          legend.pos = "right",
                          rho_colors = c("#BEEBD3", "#7ED3B7", "#43BDAA", "#2F86A2", "#33679C", "#3E488E", "black"),
                          linewidth = 1, point_size = 1.5){
  # Arguments
  #   data: dataset with each of these rhos and each num_targets between 1 and 10
  #   text_size: size of all text on plot
  #   legend.pos: either "right" or "none" for no legend
  #   rho_colors: color of each recombination rate (least to greatest)
  #   linewidth: width of lines
  #   point_size: size of points
  # Returns:
  #   List of plots summarizing the data:
  #     1. rate.gl.below.0.2: fraction of simulations where the genetic load fell below 0.2
  #     2. rate.gl.above.0.2: fraction of simulations where the genetic load did NOT fall below 0.2
  #     3. time.gl.below.0.2: avg number of generations until the genetic load fell below 0.2
  #     4. time.gl.below.0.2.given.rate.above.10: avg number of generations until the genetic load fell below 0.2 | at least 2 replicates
  #     5. rate.disrupted.lost: fraction of simulations where disrupted alleles were eventually lost (after 10,000 gens)
  #     6. time.disrupted.lost: avg number of generations until disrupted alleles were lost (given that they were)
  #     7. max.gl: the avg maximum genetic load attained by the drive
  #     8. time.max.gl: avg number of generations until the max genetic load
  #     9. ending.gl: the avg genetic load after 10,000 generations
  
  labs = c(TeX(r"($10^{-7}$)"), TeX(r"($10^{-6}$)"), TeX(r"($10^{-5}$)"), 
           TeX(r"($10^{-4}$)"), TeX(r"($10^{-3}$)"), TeX(r"($10^{-2}$)"), TeX(r"($0.5$)"))
  
  
  plt.below.0.2 =  ggplot(data, aes(x = num_targets, 
                                    y = avg_rate_genetic_load_fell_below_marker, color = as.factor(recomb_rate))) +
    geom_point(size = point_size) +
    geom_line(size = linewidth) +
    labs(
      x = TeX(r"(Number of target sites ($n$))"),
      y = TeX(r"($P(\lambda(t) < 0.2)$)"),
      color = TeX(r"($\rho$)"),
    ) + ylim(0,1) +
    scale_color_manual(values = rho_colors,
                       labels = labs) +
    scale_x_continuous(breaks = seq(1, 10, by = 2)) +
    theme(
      legend.position = legend.pos,
      text = element_text(size = text_size),
      plot.title = element_text(size = text_size, face = "bold"),
      axis.text = element_text(size = text_size),
      legend.text = element_text(size = text_size),
      aspect.ratio = 1
    ) +
    scale_y_continuous(breaks = seq(0,1, by = 0.2), limits = c(0,1))
  
  plt.above.0.2 =  ggplot(data, aes(x = num_targets, 
                                    y = 1 - avg_rate_genetic_load_fell_below_marker, 
                                    color = as.factor(recomb_rate))) +
    geom_point(size = point_size) +
    geom_line(size = linewidth) +
    labs(
      x = TeX(r"(Number of target sites ($n$))"),
      y = TeX(r"($P(\lambda(t = 10000) > 0.2$))"),
      color = TeX(r"($\rho$)"),
    ) + ylim(0,1) +
    scale_color_manual(values = rho_colors,
                       labels = labs) +
    scale_x_continuous(breaks = seq(1, 10, by = 2)) +
    theme(
      legend.position = legend.pos,
      text = element_text(size = text_size),
      plot.title = element_text(size = text_size, face = "bold"),
      axis.text = element_text(size = text_size),
      legend.text = element_text(size = text_size),
      aspect.ratio = 1
    ) +
    scale_y_continuous(breaks = seq(0,1, by = 0.2), limits = c(0,1))
  
  plt.ending.gl =  ggplot(data, aes(x = num_targets, 
                                    y = avg_ending_genetic_load, color = as.factor(recomb_rate))) +
    geom_point(size = point_size) +
    geom_line(size = linewidth) +
    labs(
      x = TeX(r"(Number of target sites ($n$))"),
      y = TeX(r"($\lambda(t = 10000)$)"),
      color = TeX(r"($\rho$)"),
    ) +
    scale_color_manual(values = rho_colors,
                       labels = labs) +
    scale_x_continuous(breaks = seq(1, 10, by = 2)) +
    theme(
      legend.position = legend.pos,
      text = element_text(size = text_size),
      plot.title = element_text(size = text_size, face = "bold"),
      axis.text = element_text(size = text_size),
      legend.text = element_text(size = text_size),
      aspect.ratio = 1
    ) +
    scale_y_continuous(breaks = seq(0, 0.6, by = 0.2), limits = c(0,0.6))
  
  plt.max.gl =  ggplot(data, aes(x = num_targets, 
                                    y = avg_max_load, color = as.factor(recomb_rate))) +
    geom_point(size = point_size) +
    geom_line(size = linewidth) +
    geom_hline(yintercept = 0.95, linewidth = linewidth, linetype = "dashed") +
    labs(
      x = TeX(r"(Number of target sites ($n$))"),
      y = TeX(r"($\max(\lambda(t))$)"),
      color = TeX(r"($\rho$)"),
    ) +
    scale_color_manual(values = rho_colors,
                       labels = labs) +
    scale_x_continuous(breaks = seq(1, 10, by = 2)) +
    theme(
      legend.position = legend.pos,
      text = element_text(size = text_size),
      plot.title = element_text(size = text_size, face = "bold"),
      axis.text = element_text(size = text_size),
      legend.text = element_text(size = text_size),
      aspect.ratio = 1
    ) +
    scale_y_continuous(breaks = seq(0,1,by = 0.2), limits = c(0,1))
  
  
  plt.time.max.gl = ggplot(data, aes(x = num_targets, 
                                     y = avg_time_to_maximum_genetic_load, color = as.factor(recomb_rate))) +
    geom_point(size = point_size) +
    geom_line(size = linewidth) +
    labs(
      x = TeX(r"(Number of target sites ($n$))"),
      y = TeX(r"(Time to $\max(\lambda(t))$)"),
      color = TeX(r"($\rho$)"),
    ) +
    scale_color_manual(values = rho_colors,
                       labels = labs) +
    scale_x_continuous(breaks = seq(1, 10, by = 2)) +
    theme(
      legend.position = legend.pos,
      text = element_text(size = text_size),
      plot.title = element_text(size = text_size, face = "bold"),
      axis.text = element_text(size = text_size),
      legend.text = element_text(size = text_size),
      aspect.ratio = 1
    ) 
  
  plt.tau = ggplot(data, aes(x = num_targets, y = avg_time_to_genetic_load_below_marker, color = as.factor(recomb_rate))) +
    geom_point(size = point_size) +
    geom_line(size = linewidth) +
    labs(
      x = TeX(r"(Number of target sites ($n$))"),
      y = TeX(r"($\tau$)"),
      color = TeX(r"($\rho$)"),
    ) +
    scale_color_manual(values = rho_colors,
                       labels = labs) +
    scale_x_continuous(breaks = seq(1, 10, by = 2)) +
    theme(
      legend.position = legend.pos,
      text = element_text(size = text_size),
      plot.title = element_text(size = text_size, face = "bold"),
      axis.text.x = element_text(size = text_size),
      axis.title.x = element_text(size = text_size),
      axis.title.y = element_text(size = 1.3*text_size),
      legend.text = element_text(size = text_size),
      axis.text.y = element_text(size = text_size),
      aspect.ratio = 1
    ) 
  
  
  plt.tau.conditional = data %>% filter(avg_rate_genetic_load_fell_below_marker >= 2/20) %>%
    ggplot(aes(x = num_targets, y = avg_time_to_genetic_load_below_marker, color = as.factor(recomb_rate))) +
    geom_point(size = point_size) +
    geom_line(size = linewidth) +
    labs(
      x = TeX(r"(Number of target sites ($n$))"),
      y = TeX(r"($\tau$)"),
      color = TeX(r"($\rho$)"),
    ) +
    scale_color_manual(values = rho_colors,
                       labels = labs) +
    scale_x_continuous(breaks = seq(1, 10, by = 2)) +
    theme(
      legend.position = legend.pos,
      text = element_text(size = text_size),
      plot.title = element_text(size = text_size, face = "bold"),
      axis.text.x = element_text(size = text_size),
      axis.title.x = element_text(size = text_size),
      axis.title.y = element_text(size = 1.3*text_size),
      axis.text.y = element_text(size = text_size),
      legend.text = element_text(size = text_size),
      aspect.ratio = 1
    )
  
  
  plt.total.time = ggplot(data, aes(x = num_targets, 
                                    y = avg_time_to_disrupted_lost, 
                                    color = as.factor(recomb_rate))) +
    geom_point(size = point_size) +
    geom_line(size = linewidth) +
    labs(
      x = TeX(r"(Number of target sites ($n$))"),
      y = TeX(r"(Time until all disrupted sites lost)"),
      color = TeX(r"($\rho$)"),
    ) +
    scale_color_manual(values = rho_colors,
                       labels = labs) +
    scale_x_continuous(breaks = seq(1, 10, by = 2)) +
    theme(
      legend.position = legend.pos,
      text = element_text(size = text_size),
      plot.title = element_text(size = text_size, face = "bold"),
      axis.text = element_text(size = text_size),
      legend.text = element_text(size = text_size),
      aspect.ratio = 1
    )
  
  
  rate.disrupted.lost = ggplot(data, aes(x = num_targets, 
                                         y = disrupted_lost_rate, 
                                         color = as.factor(recomb_rate))) +
    geom_point(size = point_size) +
    geom_line(size = linewidth) +
    labs(
      x = TeX(r"(Number of target sites ($n$))"),
      y = TeX(r"(Rate of disrupted sites loss)"),
      color = TeX(r"($\rho$)"),
    ) +
    scale_color_manual(values = rho_colors,
                       labels = labs) +
    scale_x_continuous(breaks = seq(1, 10, by = 2)) +
    theme(
      legend.position = legend.pos,
      text = element_text(size = text_size),
      plot.title = element_text(size = text_size, face = "bold"),
      axis.text = element_text(size = text_size),
      legend.text = element_text(size = text_size),
      aspect.ratio = 1
    ) +
    scale_y_continuous(breaks = seq(0,1, by = 0.2), limits = c(0,1))
  
  
  plot.list = list(rate.gl.below.0.2 = plt.below.0.2,
                   rate.gl.above.0.2 = plt.above.0.2,
                    time.gl.below.0.2 = plt.tau,
                    time.gl.below.0.2.given.rate.above.10 = plt.tau.conditional,
                    rate.disrupted.lost = rate.disrupted.lost,
                    time.disrupted.lost = plt.total.time,
                    max.gl = plt.max.gl,
                    time.max.gl = plt.time.max.gl,
                    ending.gl = plt.ending.gl)
  
  return(plot.list)
}

#### Plotting a single multilocus replicate's genetic load, disrupted, and r1 frequencies
plot_multilocus_single_replicate = function(data, lambda = NULL,
                                            text_size = 14, 
                                            linewidth = 1,
                                            add_frequencies = T,
                                            add_genetic_load = F,
                                            add_lambda_line = F,
                                            legend.pos = "right",
                                            alpha = 1,
                                            col.palette = c("red", "green4","black")){
  # Arguments:
  #     data: SLiM csv output. Should have columns for gen_num, broken_rate, and func_rate
  #     lambda: the max potential genetic load (can leave NULL if !add_lambda_line)
  #     text_size: size of all text in plot
  #     linewidth: width of each line
  #     add_frequencies: whether to also plot the disrupted rate and functional target rate 
  #                     (otherwise just plots the genetic load)
  #     add_lambda_line: whether to add the lambda dashed line
  #     legend.pos: "right" for legend or "none" for no legend
  #     alpha: transparency
  #     col.palette: c(disrupted color, resistance color,  genetic load color)
  #
  # Returns:
  #   plot of the genetic load over time, as well as the disrupted and functional-target frequencies if (add_frequencies)
  
  if (add_frequencies){
    
    if (add_genetic_load){
      data.long = data %>% pivot_longer(cols = c(genetic_load, broken_rate, func_rate), names_to = "type",
                                        values_to = "value")
      
      data.long$type = factor(
        data.long$type,
        levels = c("broken_rate", "genetic_load", "func_rate")  # desired plotting order
      )
      
      base = ggplot(data.long, aes(x = gen_num, y = value, color = type))
      
      if (add_lambda_line){
        base = base + geom_hline(yintercept = lambda, linetype = "dashed", linewidth = linewidth,
                                 color = col.palette[3], alpha = alpha)
      }
      
      p = base +
        geom_line(linewidth = linewidth, alpha = alpha) +
        scale_color_manual(values = c(broken_rate = col.palette[1], func_rate = col.palette[2],
                                      genetic_load = col.palette[3]),
                           breaks = c( "broken_rate","func_rate", "genetic_load"),  
                           labels = c("disrupted frequency", "resistance frequency", 
                                      TeX(r"(genetic load ($\lambda$))"))) +
        labs(y = "Value", color = "Value") 
    } else {
      # just disrupted and resistance frequencies
      data.long = data %>% pivot_longer(cols = c(broken_rate, func_rate), names_to = "allele",
                                        values_to = "frequency")
      
      data.long$allele = factor(
        data.long$allele,
        levels = c("broken_rate","func_rate")  # desired plotting order
      )
      
      base = ggplot(data.long, aes(x = gen_num, y = frequency, color = allele))
      
      p = base +
        geom_line(linewidth = linewidth, alpha = alpha) +
        scale_color_manual(values = c(broken_rate = col.palette[1], func_rate = col.palette[2]),
                           breaks = c( "broken_rate","func_rate"),  
                           labels = c("disrupted", 
                                      "resistance")) +
        labs(y = "Frequency", color = "allele") 
      
      
    }
  
    # just genetic load  
  } else {
    base = ggplot(data, aes(x = gen_num, y = genetic_load))
    if (add_lambda_line){
      base = base + geom_hline(yintercept = lambda, 
                               linetype = "dashed", li, alpha = alpha)
    }
    
    p = base +
      geom_line(linewidth = linewidth, color = col.palette[3], alpha = alpha) +
      labs(y =TeX(r"(Genetic load ($\lambda$))"))
  }
  
  p2 = p + labs(x = TeX(r"(Generation ($t$))")) +
    scale_y_continuous(breaks = seq(0,1, by = 0.2), limits = c(0,1)) +
    theme(
      legend.position = legend.pos,
      text = element_text(size = text_size),
      plot.title = element_text(size = text_size, hjust = 0.5),
      axis.text = element_text(size = text_size),
      legend.text = element_text(size = text_size),
      aspect.ratio = 1
    )
    
  return(p2)
    
}
