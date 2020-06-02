boot_community <- function(nets,
                           lag = 0,
                           simluations = 50,
                           weight = F,
                           log.weight = F,
                           seed = 123)
{
  require('igraph')
  require('purr')
  require('Matrix')
  if(is.list(nets) == F)
  {
    stop("Nets needs to be a list of igraph networks. ")
  } else {
    if (is.null(names(nets)))
    {
      stop("Network list needs names for indexing.")
    } else {
      if (lag == 0)
      {
        cat("\r No lag term entered.")
      }

      ## Extracting network years
      net_years <- map(nets, extract_yearly_nets) ## pull out the network years

      ## Creating year range
      years <- min(sapply(net_years, find_range)):max(sapply(net_years, find_range))

      boot_sims_by_range <- list() ## return_data list
      pct_complete <- round(length((lag + 1):length(years)) * seq(0.1, 0.9, 0.1))
        for (i in (lag + 1):length(years))
        {
          if (i %in% pct_complete)
          {
            cat("\r Pct complete: ", round(i/length((lag + 1):length(years)), 2))
          }
          if (lag != 0)
          {
            range <- years[i - lag]:(years[i - 1])
          } else {
            range = years[i]
          }

          check_range <- sapply(net_years, function(x) {sum(x %in% range)})
          net_names <- names(nets)

          ## Which networks have data in the five year window?
          reduced_net_names <- net_names[check_range > 0]

          ## This look indexes the names of the network to keep
          ## only networks with data in the range
          reduced_net <- list()
          for (l in 1:length(reduced_net_names)) {
            reduced_net[[l]] <- nets[[reduced_net_names[l]]]
          }

          ## The networks are named for indexing
          names(reduced_net) <- reduced_net_names

          ## find the locations in the list of the networks
          ## that fall in the sliding 5 year window
          nets_to_keep <- list()
          for (h in 1:length(reduced_net)) {
            nets_to_keep[[h]] <- which(
              sapply(reduced_net[[h]], graph_attr, "year") %in% range
            )
          }

          ## Keep the networks in the 5 year window based on position
          for (k in 1:length(nets_to_keep)) {
            reduced_net[[k]] <- reduced_net[[k]][nets_to_keep[[k]]]
          }

          ## Name the networks by year for indexing
          for (year_name in 1:length(reduced_net)) {
            names(reduced_net[[year_name]]) <- sapply(
              reduced_net[[year_name]], graph_attr, "year"
            )
          }

          edge_attr_names_list <- list()
          for (edge_attr_loop in 1:length(reduced_net))
          {
            if (sum(sapply(reduced_net[[edge_attr_loop]], edge_density) == 1) == length(sapply(reduced_net[[edge_attr_loop]], edge_density)))
            {
              edge_attr_names_list[[edge_attr_loop]] <- unique(sapply(reduced_net[[edge_attr_loop]], edge_attr_names))
            }
          }
          potential_attrs <- unlist(edge_attr_names_list)


          ## Create adjacency matricies
          adj_mat <- list()
          for (adj_loop in 1:length(reduced_net)) {
            ## Check if the network is weighted
            edge_ats <- names(edge_attr(reduced_net[[adj_loop]][[1]]))
            edge_ats <- edge_ats[edge_ats != "na"]
            edge_ats <- edge_ats[edge_ats != "year"]
            edge_ats <- edge_ats[edge_ats %in% potential_attrs]
            edge_ats <- edge_ats[1]
            if (sum(edge_ats %in% potential_attrs) > 0) {
              ## If edge is weighted then pull out weighted edges
              adj_mat[[adj_loop]] <- lapply(
                reduced_net[[adj_loop]], as_adjacency_matrix, attr = edge_ats
              )
            } else {
              ## If it is unweighted then get dyads
              adj_mat[[adj_loop]] <- lapply(reduced_net[[adj_loop]], as_adjacency_matrix)
            }
          }

          for(make_mat in 1:length(adj_mat))
          {
            for (index_data in 1:length(adj_mat[[make_mat]]))
            {
              adj_mat[[make_mat]][[index_data]] <- as.matrix(adj_mat[[make_mat]][[index_data]])
              adj_mat[[make_mat]][[index_data]][is.na(adj_mat[[make_mat]][[index_data]])] <- 0
              adj_mat[[make_mat]][[index_data]][is.nan(adj_mat[[make_mat]][[index_data]])] <- 0
            }
          }

          ## Log the overdispersed networks

          if (log.weight == T)
          {
            for (log_weights in 1:length(adj_mat)) {
              for (check_nets in 1:length(adj_mat[[log_weights]])) {
                if (length(unique(
                  as.vector(adj_mat[[log_weights]][[check_nets]]))) > 20) {
                  if (sum(adj_mat[[log_weights]][[check_nets]] < 0) > 0)
                  {
                    adj_mat[[log_weights]][[check_nets]] <- log(
                      adj_mat[[log_weights]][[check_nets]] + 1 + abs(min(adj_mat[[log_weights]][[check_nets]]))
                    )
                  } else
                    {
                  adj_mat[[log_weights]][[check_nets]] <- log(
                    adj_mat[[log_weights]][[check_nets]] + 1
                  )
                    }
                  }
                }
              }
          }
        if (lag > 0)
        {
          reduce_lag <- rev(order((range - min(range) + 1)))

          range_data <- data.frame(range = range, influence = reduce_lag)
          range_data$reduce_influence_older <- 1/range_data$influence
          for (weight_loop in 1:nrow(range_data))
          {
            for (weight_nets in 1:length(adj_mat)) {
              for(time_since_t in 1:length(adj_mat[[weight_nets]])) {
                if (names(adj_mat[[weight_nets]])[time_since_t] == range_data$range[weight_loop]) {
                  adj_mat[[weight_nets]][[time_since_t]] <- adj_mat[[weight_nets]][[time_since_t]] * range_data$reduce[weight_loop]
                }
              }
            }
          }
        }

          ## find missing states in adjacency matricies
          all_nodes <- c()
          for (node_names in 1:length(adj_mat)) {
            for (check_index in 1:length(adj_mat[[node_names]])) {
              all_nodes <- unique(c(all_nodes, colnames(adj_mat[[node_names]][[check_index]])))
            }
          }



          ## Identify is a adj matrix is missing a state then adding and
          ## sorting all other states.
          for (add_nodes in 1:length(adj_mat)) {
            for (index_time in 1:length(adj_mat[[add_nodes]])) {
              missing_nodes <- all_nodes[!(
                all_nodes %in% colnames(adj_mat[[add_nodes]][[index_time]])
              )]
              if (length(missing_nodes) > 0) {
                new_cols <- data.frame(matrix(
                  nrow = nrow(adj_mat[[add_nodes]][[index_time]]),
                  ncol = length(missing_nodes), 0
                ))
                colnames(new_cols) <- missing_nodes
                adj_mat[[add_nodes]][[index_time]] <- cbind(
                  adj_mat[[add_nodes]][[index_time]], new_cols
                )
                new_rows <- t(data.frame(matrix(
                  nrow = ncol(adj_mat[[add_nodes]][[index_time]]),
                  ncol = length(missing_nodes), 0
                )))
                colnames(new_rows) <- colnames(adj_mat[[add_nodes]][[index_time]])
                rownames(new_rows) <- missing_nodes
                adj_mat[[add_nodes]][[index_time]] <- rbind(
                  adj_mat[[add_nodes]][[index_time]], new_rows
                )
                adj_mat[[add_nodes]][[index_time]] <- as.matrix(
                  adj_mat[[add_nodes]][[index_time]]
                )
              }
              adj_mat[[add_nodes]][[index_time]] <-  adj_mat[[add_nodes]][[index_time]][order(rownames( adj_mat[[add_nodes]][[index_time]])),]
              adj_mat[[add_nodes]][[index_time]] <-  adj_mat[[add_nodes]][[index_time]][,order(colnames( adj_mat[[add_nodes]][[index_time]]))]
            }
          }

          ## Create the sliding window
          sliding_adj_mat <- list()
          for (window in 1:length(adj_mat)) {
            sliding_adj_mat[[window]] <- Reduce("+", adj_mat[[window]])
          }



          ## Scale those networks between 0 and 1
          trans_mats <- map(sliding_adj_mat, zero_one_mat)

          if (weight == T)
          {
            density_weight <- c()
            for (density_check in 1:length(sliding_adj_mat)) {
              density_weight[density_check] <- 1 - sum(sliding_adj_mat[[density_check]]==0)/(nrow(sliding_adj_mat[[density_check]])^2)
            }
            density_weight[density_weight == 0] <- 1e-6
            for (density_set in 1:length(density_weight))
            {
              trans_mats[[density_set]] <- trans_mats[[density_set]]*density_weight[density_set]
            }
          }

          boot_sims <- simluations
          ## Bootstrap community detection to find different communities

          ## Create empty boot data.frame for community detection
          boot_trans_mats <-
            data.frame(
              name = all_nodes,
              time = paste0(years[i - lag], "-", years[i - 1])
            )
          boot_trans_mats <-
            plyr::arrange(
              boot_trans_mats,
              name
            )
          set.seed(seed)
          ## start bootstrap
          for (boot in 1:boot_sims) {
            boot_data <- sample(length(trans_mats), length(trans_mats), replace = TRUE)
            final_mat <- apply(simplify2array(trans_mats[boot_data]), 1:2, mean)
            net <- graph_from_adjacency_matrix(
              final_mat,
              mode = "directed",
              weighted = TRUE,
              diag = FALSE,
              add.rownames = TRUE
            )

            boot_community <- cluster_walktrap(net)
            boot_output <- data.frame(name = names(membership(boot_community)),
                                      time = paste0(years[i - lag], "-", years[i - 1]),
                                      community = as.vector(membership(boot_community)))

            boot_output$community <- as.character(boot_output$community)
            boot_output$community <- as.factor(boot_output$community)
            boot_output$community <- as.numeric(boot_output$community)
            colnames(boot_output)[3] <- paste0("community_sim_", boot)

            boot_trans_mats <- merge(
              boot_trans_mats,
              boot_output,
              by = c("name", "time"),
              all.x = T
            )
          }

          boot_sims_by_range[[i]] <- boot_trans_mats
        }
      return.item <- do.call(rbind, boot_sims_by_range)
      return(return.item)
      }
    }
}

