## ---------------------------
##
## Generate classroom networks
##
## Authors: Katia Bulekova, Wenrui Li
##          Boston University
##
## Date Created: 2020-07-31
##
## Email: ktrn@bu.edu
##
## ----------------------------

library(Matrix)
suppressMessages( library(igraph) )
suppressMessages( library(dplyr) )


#' Create igraph object from a sparse matrix.
#' 
#' @param netMatrix sparse matrix defining relationship
#' @return igraph graph object.
#'
buildGraph <- function ( netMatrix ) {
  
  NumOfStudents <- nrow(netMatrix)
  NumOfGroups <- ncol(netMatrix)
  Net <- vector("list", NumOfGroups)
  names(Net) <- colnames(netMatrix)
  Net[1:NumOfGroups] <- list(rsparsematrix(NumOfStudents ,NumOfStudents, 0))
  Size <- colSums(netMatrix)
  
  for ( i in seq_len(NumOfGroups) ) {
    if( Size[i] > 1 ){
      nonzero_entry <- which( netMatrix[,i] > 0 )
      Net[[i]][nonzero_entry,nonzero_entry] <- 1 - diag(length(nonzero_entry))
    }
  }
  
  graph <- Reduce("+", Net)
  graph <- graph_from_adjacency_matrix(graph,mode = "undirected",
                                       weighted = TRUE,
                                       diag = FALSE)
  
  vertex_attr(graph, "full_info_id") <- row.names(netMatrix)
  graph
}

##=================================================================

## output path
output_path = "../Data/networks/"

## the maximum number of students per a household
## baseline - number of students per bathroom
## 10 as intervention
dorm_df <- read.csv("../Data/input/housing_info.csv", 
                      stringsAsFactors = FALSE)


## ---------------------
## Baseline household
## ---------------------
## household can be defined by the group of students using the same bathroom
household <- dorm_df %>%
  group_by(buildingID, floor) %>%
  mutate(householdID = 
           paste(buildingID,
                 floor, 
                 (room - 1) %% 2 + 1, 
                 sep="_") ) %>%   # 2 households per floor
  mutate(smallHouseholdID = 
           paste(buildingID, 
                 floor, 
                 (room - 1) %% 5 + 1, 
                 sep="_") ) %>%   # 5 households per floor
  ungroup() %>% arrange(householdID, smallHouseholdID)
           
## construct sparse matrix for the baseline household
netMatrix <- household %>%
  select(studentID, householdID) %>%
  table() %>%
  as.data.frame.matrix() %>%
  as.matrix() %>%
  as("sparseMatrix")

graph <- buildGraph(netMatrix)
write_graph(graph,
            paste0( output_path, "HouseholdNet_baseline.graphml"),"graphml" )



## ------------------------------------
## Intervention (use smaller household)
## ------------------------------------
netMatrix <- household %>%
  select(studentID, smallHouseholdID) %>%
  table() %>%
  as.data.frame.matrix() %>%
  as.matrix() %>%
  as("sparseMatrix")


graph <- buildGraph(netMatrix)
write_graph(graph,
            paste0( output_path, "HouseholdNet_10.graphml"),"graphml" )


## -------------------------------------
## Construct room network
## -------------------------------------
netMatrix <- household %>%
  mutate( roomID = paste(buildingID, room, sep="_")) %>%
  select(studentID, roomID) %>%
  table() %>%
  as.data.frame.matrix() %>%
  as.matrix() %>%
  as("sparseMatrix")


graph <- buildGraph(netMatrix)
write_graph(graph,
            paste0( output_path, "RoomNet.graphml"),"graphml" )

## -------------------------------------
## Construct floor network
## -------------------------------------
netMatrix <- household %>%
  mutate( floorID = paste(buildingID, floor, sep="_")) %>%
  select(studentID, floorID) %>%
  table() %>%
  as.data.frame.matrix() %>%
  as.matrix() %>%
  as("sparseMatrix")


graph <- buildGraph(netMatrix)
write_graph(graph,
            paste0( output_path, "FloorNet.graphml"),"graphml" )

## --------------------------------------------------
##  Construct "friends" network inside the building
## --------------------------------------------------
friends <- function(rowID, roomV, StudentV){
  
  
  out <- lapply(rowID, function(x){
    # remove roomates and yourself
    rommates <-  which(roomV == roomV[x])
    net = StudentV[-rommates]
    
    if(length(net) > 2) {
      group <- sample(net, 3) 
    } else if (length(net) > 1)  {
      group <- c( sample(net, 2), NA)
      
    } else if (length(net) == 1)  {
      group <- c( net[1], NA, NA)
    } else {
      group <- c(NA, NA, NA)
    }
    group
    
  } )
  out
}

## Assign a bathroom ID to a student based on a number of floorbaths
set.seed(42)

friends.df <- household %>%
  group_by( buildingID ) %>%
  mutate( friends = friends(row_number(), room, studentID) ) %>%
  ungroup() %>% tidyr:: separate(friends, 
                                 into = c(NA, "f1", "f2", "f3"),
                                 extra = "drop") %>%
  mutate(f1 = as.numeric(f1), f2 = as.numeric(f2), f3 = as.numeric(f3)) %>%
  select(studentID, f1, f2, f3) %>% mutate(vid = row_number()) 


friends.df$fr1 = sapply(friends.df$f1, 
                        function(x)friends.df$vid[which(x == friends.df$studentID)])
friends.df$fr2 = sapply(friends.df$f2, 
                        function(x)if(!is.na(x))friends.df$vid[which(x == friends.df$studentID)] else NA)
friends.df$fr3 = sapply(friends.df$f3, 
                        function(x)if(!is.na(x))friends.df$vid[which(x == friends.df$studentID)] else NA)

g <- make_empty_graph(directed = FALSE) + 
  vertices(friends.df$vid) 

vertex_attr(g, "full_info_id") <- friends.df$studentID


edges_list1 <- as.vector( t(na.omit(cbind(friends.df$vid, friends.df$fr1))))
edges_list2 <- as.vector( t(na.omit(cbind(friends.df$vid, friends.df$fr2))))
edges_list3 <- as.vector( t(na.omit(cbind(friends.df$vid, friends.df$fr3))))

#add edges
g <- g + edges( c(edges_list1, edges_list2, edges_list3) )
write_graph(g,
            paste0(output_path, "Building_net.graphml"),"graphml" )
