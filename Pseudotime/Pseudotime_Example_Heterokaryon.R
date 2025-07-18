library(dyno)
library(dynwrap)
library(dyneval)
library(dynmethods)
library(tidyverse)
library(Matrix)
library(scales)


## function to save trajectory plot
# save plots easily
saveTrajectoryPlot <- function(h, title, save_path) {
  h <- h+labs(title=title)
  ggsave(save_path, plot = h, width = 8, height = 6, dpi = 300)
}

## Parameters
root_dir = "Heterokaryon_Example_Data"
csv_dir = file.path(root_dir,"Input_O-SNAP_CSV")
markers = c("H3K27me3","H2B")
countsValue = 0.5
start_group = "CtrlhFb"
end_group = "HK48h"

## Setup
babelwhale::set_default_config(babelwhale::create_docker_config()) # ensure that the code runs with docker (as opposed to singularity)
set.seed(0) # for reproducibility

## Open connection to find optimal timeline inference methods
guidelines <- guidelines_shiny()

# ## Below is the hard-coded version of responses to the guidelines
# # If desired, comment out line 29 and uncomment below (33-43) to hardcode guideline responses
# answers <- dynguidelines::answer_questions(
#   multiple_disconnected = FALSE,
#   expect_topology = TRUE,
#   expected_topology = "linear",
#   n_cells = 350,
#   n_features = 144,
#   memory = "2GB",
#   prior_information = c("start_id","start_n","end_n","groups_id","groups_network"),
#   docker = TRUE
# )
# guidelines <- guidelines(answers = answers)

## Get methods based on results from guidelines
methods_selected <- guidelines$methods_selected
print(guidelines$methods_selected)
print("----")

for (m in 1:length(markers)){
  
  marker <- markers[[m]]
  analysisName <- paste("ANALYSIS_Heterokaryon",marker,sep="_")
  
  ## Directory management
  save_dir = file.path(root_dir,analysisName)
  csv_path = file.path(csv_dir,paste0(analysisName,".csv"))
  
  ## Prepare data (filter, normalize)
  data <- read_csv(csv_path,TRUE,show_col_types=FALSE) %>%
    filter(group != "CtrlmESC")%>%
    unite("cell_id",1:3,remove=FALSE) %>%
    subset(select = -c(biological_replicate,name)) %>%
    mutate(across(where(is.numeric), scale)) %>%
    select_if(colSums(!is.na(.)) > 0) %>%
    filter(complete.cases(.)) 
  
  ## Create expression data
  expression <-  data  %>% subset(select = -c(cell_id, group)) %>% as("sparseMatrix")
  expression@Dimnames[[1]] = data$cell_id
  
  ## Fill count data with dummy value (countsValue)
  counts <- expression
  counts@x = rep(countsValue, length(expression@x));
  
  ## Create dataset variable to pass to dyno methods
  dataset <- wrap_expression(expression=expression, counts=counts)
  ## OPTIONAL (highly recommended): add prior information
  # PRIORS
  #   start_id:       The start cells
  #   start_n:        Number of start states
  #   end_n:          Number of end states
  #   group_id:       The grouping of cells, a dataframe with cell_id and group_id
  #   groups_network: The network between groups, a dataframe with from and to
  groups_id <- tibble(.rows = NULL,
                      cell_id = data$cell_id,
                      group_id = data$group
  )
  groups_network <- tibble(.rows = NULL,
                           from = c("CtrlhFb","HK06h","HK24h"),
                           to = c("HK06h","HK24h","HK48h")
  )
  dataset <- add_prior_information(dataset,
                                   start_id = data %>% filter(stringr::str_detect(group,start_group)) %>% pull(cell_id),
                                   start_n = 1,
                                   end_n = 1,
                                   groups_id = groups_id,
                                   groups_network = groups_network
  )
  ## OPTIONAL: add grouping  information
  dataset <- add_grouping(
    dataset,
    data$group
  )
  # ## optional: add dimensionality reduction information
  # dataset <- add_dimred(
  #   dataset,
  #   dataset$dimred
  # )
  
  ## Loop through markers and methods
  for (i in 1:length(methods_selected)){
    method_name <- methods_selected[[i]]
    print(marker)
    print(method_name)
    print("----")
    tryCatch({
      ## Generate model
      model <- infer_trajectory(dataset, 
                                methods_selected[[i]],
                                give_priors = c("start_id","start_n","end_n","groups_id","groups_network"),
                                seed=0)
      # Root model
      tryCatch({
        if(start_group %in% model$milestone_ids)
          model <- model %>% add_root(root_milestone_id = start_group)
        else if("milestone_begin" %in% model$milestone_ids)
          model <- model %>% add_root(root_milestone_id = "milestone_begin")
      },
      error = function(err) print(err)
      )
      # Simplify model
      tryCatch({
        if(length(model$milestone_ids) > 2)
          model <- simplify_trajectory(model)
      },
      error = function(err)  print(err)
      )
      print(model$trajectory_type)
      # Calculate pseudotime
      pseudotime_data <- calculate_pseudotime(model)
      # Create linear projection of model
      model_linear <- add_linear_trajectory(dataset,pseudotime_data)
      model_linear <- model_linear %>% add_root(root_milestone_id = "milestone_begin")
      # If trajectory topology is not linear, focus on the branch from desired start to end state
      if(model$trajectory_type != "linear"){
        idx = model$progressions$to==end_group # if multiple paths, exclude those not leading to end group
        cells_to_remove = model_linear$cell_ids[!idx]
        idx_percentages_remove = as.numeric(flatten(as.list(sapply(cells_to_remove,FUN=function(X) which(model_linear$milestone_percentages$cell_id %in% X)))))
        model_linear$cell_info <- model_linear$cell_info[idx,]
        model_linear$cell_ids <- model_linear$cell_ids[idx]
        model_linear$counts <- model_linear$counts[idx,]
        model_linear$expression <- model_linear$expression[idx,]
        model_linear$grouping <- model_linear$grouping[idx]
        model_linear$milestone_percentages <- model_linear$milestone_percentages[-idx_percentages_remove,]
        model_linear$milestone_percentages$percentage <- rescale(model_linear$milestone_percentages$percentage)
        model_linear$progressions <- model_linear$progressions[idx,]
        model_linear$progressions$percentage <- rescale(model_linear$progressions$percentage)
        model_linear$dimread <- model_linear$dimread[idx,]
        model_linear$pseudotime <- rescale(model_linear$pseudotime[idx])
      }
      
      ## Plot trajectories
      # Plot original model (if extraneous branches were removed, cells assigned to them are listed as "NA")
      #   Color by pseudotime
      h = plot_dimred(model,"pseudotime", pseudotime = model_linear$pseudotime,plot_trajectory=FALSE) 
      saveTrajectoryPlot(h, "Default - model", file.path(save_dir,method_name,"Pseudotime.png"))
      #   Color by groups
      h = plot_dimred(model,grouping = model_linear$grouping,plot_trajectory=FALSE)
      saveTrajectoryPlot(h, "Default - model", file.path(save_dir,method_name,"Pseudotime_grouping.png"))
      # If trajectory topology is not linear, plot version with all samples
      if(model$trajectory_type != "linear"){
        h = plot_dimred(model,"pseudotime", pseudotime = pseudotime_data,plot_trajectory=FALSE) 
        saveTrajectoryPlot(h, "Default - model", file.path(save_dir,method_name,"Pseudotime_all.png"))
        h = plot_dimred(model,grouping = data$group,plot_trajectory=FALSE) 
        saveTrajectoryPlot(h, "Default - model", file.path(save_dir,method_name,"Grouping_all.png"))
      }
      # Plot simplified linear projection
      h = plot_onedim(model_linear, grouping = model_linear$grouping)
      saveTrajectoryPlot(h, "Default - linear model - one dim", file.path(save_dir,method_name,"Onedim_Linear.png"))
      
      # Plot a heat map of the expression
      h = plot_heatmap(model_linear,features_oi=20,expression_source = model_linear$expression,label_milestones = FALSE)
      saveTrajectoryPlot(h, "Default - trajectory", file.path(save_dir,method_name,"Trajectory_Heatmap_20.png"))
      
      # Plot with selected features
      features_oi <- c("locs_density","log_voronoi_density_mean","log_interior_dbscan_cluster_density_mean")
      h_selected = plot_heatmap(model_linear,features_oi=features_oi,expression_source = model_linear$expression,label_milestones = FALSE)
      saveTrajectoryPlot(h_selected, "Default - trajectory", file.path(save_dir,method_name,"Trajectory_Heatmap_selected.png"))
      
      ## Save data
      save_name = paste(c(method_name,marker),collapse="_")
      pseudotime_table<-data.frame(cell=model_linear$grouping,pseudotime=model_linear$pseudotime)
      write.csv(pseudotime_table,
                file = file.path(save_dir,method_name,paste(c(save_name,".csv"),collapse="")));
      if(model$trajectory_type != "linear"){
        pseudotime_table<-data.frame(cell=data$group,pseudotime=pseudotime_data)
        write.csv(pseudotime_table,
                  file = file.path(save_dir,method_name,paste(c(save_name,"_all.csv"),collapse="")));
      }
      save(dataset,list = c("model", "model_linear"),file = file.path(save_dir,method_name,paste(c(save_name,"_trajectory.RData"),collapse="")))
      
      ## Clean up
      graphics.off() # Remove 
      rm(model_linear, model)
    }, 
    error = function(err) print(err)
    )
  }
}

print("Complete")
