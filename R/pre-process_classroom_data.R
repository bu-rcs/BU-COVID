library(Matrix)
suppressMessages( library(igraph) )

set.seed(42)

## Set the paths to input and output directories
input_path = "../Data/input/"
output_path = "../Data/networks/"

## Read input Datasets
faculty_courses <- read.csv( file = paste0(input_path, "faculty_courses.csv"), 
                             stringsAsFactors = F)

student_courses <- 
  read.csv( file = paste0(input_path, "student_courses.csv"), 
           stringsAsFactors = F)

course_info <- read.csv( file = paste0(input_path, "course_info.csv"), 
                         stringsAsFactors = F)

## combine faculty and student course lists
faculty_students_course <- data.frame(ID = c(student_courses$studentID, 
                                             faculty_courses$facultyID ),
                                      courseID = c(student_courses$courseID, 
                                                   faculty_courses$courseID ), 
                                      stringsAsFactors=F)


#' --------------------------------------------------------------------
#' Get contacts per class (k-reg network) for specified classes
#' 
#' @param k_reg 
#' @param AffiMat 
#' @return ClassNet - list object
#'
getClassNet <- function( k_reg, AffiMat) {
  NumOfPop <- nrow(AffiMat)
  NumofClass <- ncol(AffiMat)
  ClassNet <- vector("list", NumofClass)
  names(ClassNet) <- colnames(AffiMat)
  ClassNet[1:NumofClass] <- list(rsparsematrix(NumOfPop,NumOfPop,0))
  ClassSize <- colSums(AffiMat)
  
  for (i in 1:NumofClass) {
    
    if(ClassSize[i]<=k_reg & ClassSize[i]>1) {
      nonzero_entry <- which(AffiMat[,i]>0 )
      ClassNet[[i]][nonzero_entry,nonzero_entry] <- 1 - diag(length(nonzero_entry))
      
    }else if(ClassSize[i]>k_reg) {
      
      nonzero_entry <- which(AffiMat[,i]>0 )
      ClassNet[[i]][nonzero_entry,nonzero_entry] <- 
        get.adjacency(sample_k_regular(ClassSize[i],k_reg))
    }
  }
  return(ClassNet)
}



#' --------------------------------------------------------------------
#' Get daily network for a week
#' 
#' @param ClassSch - list with class schedule
#' @param ClassNet - list with class info 
#' @return DailyNet - list object with Daily network information
#'
getDailyNet <- function(ClassSch,ClassNet) {
  
  DailyNet <- vector("list", length(ClassSch))
  names(DailyNet) <- names(ClassSch)
  
  for (i in 1:length(ClassSch)) {
    
    ClassNet_Day <-  ClassNet[names(ClassNet) %in%  ClassSch[[i]]]
    DailyNet[[i]] <- Reduce("+",ClassNet_Day)
    
  }
  return(DailyNet)
  
}


#' --------------------------------------------------------------------
#' Get daily full network for a week (everyone )
#' 
#' @param ClassSch - list with class schedule
#' @param AffiMat - Sparce Matrix
#' @return DailyNetFull - list object
#'
getDailyNetFull <- function(ClassSch, AffiMat) {
  
  DailyNetFull <- vector("list", length(ClassSch))
  names(DailyNetFull) <- names(ClassSch)
  
  for (i in seq_along(ClassSch) ) {
    
    AffiMat_Day <-  AffiMat[,colnames(AffiMat) %in%  ClassSch[[i]]]
    DailyNetFull[[i]] <- tcrossprod(AffiMat_Day,AffiMat_Day)
    diag(DailyNetFull[[i]]) <- 0
    
  }
  return(DailyNetFull)
}

#' --------------------------------------------------------------------
#' Get daily full network for building Platoon netwoork
#' 
#' @param k_reg 
#' @param AffiMat 
#' @param ClassSch - list with class schedule
#' @param instructorList - dataframe with instructor information
#' @return DailyNetSparse - list object
#'
getDailyNetSparse <- function(k_reg, AffiMat, ClassSch, instructorList)
{
  DailyNetSparse <- vector("list", length(ClassSch))
  names(DailyNetSparse) <- names(ClassSch)
  NameOfClass <- colnames(AffiMat)
  NumOfPop <- nrow(AffiMat)
  DailyNetSparse[1:length(ClassSch)] <- list(rsparsematrix(NumOfPop,NumOfPop,0))
  ClassSize <- colSums(AffiMat)
  for (i in seq_len( ncol(AffiMat )) ){
    
    instructorList_ind <- which(instructorList[,1] == NameOfClass[i])
    classDay <- sapply(lapply(ClassSch, 
                              function(ch) grep(NameOfClass[i], ch)), 
                       function(x) length(x) > 0) 
    
    NumOfStudent <- ClassSize[i] - length(instructorList_ind)
    nonzero_entry <- which(AffiMat[,i] > 0 )
    
    ## case 1
    if(sum(classDay) ==1) {
      if(ClassSize[i] <= k_reg & ClassSize[i] > 1) {
        
        DailyNetSparse[[ seq_len(7)[classDay]]][nonzero_entry,nonzero_entry] <-
          DailyNetSparse[[ seq_len(7)[classDay] ]][nonzero_entry,nonzero_entry]+
          1 - diag(length(nonzero_entry))
        
      }else if(ClassSize[i]>k_reg){
        DailyNetSparse[[ seq_len(7)[classDay] ]][nonzero_entry,nonzero_entry] <-  
          DailyNetSparse[[ seq_len(7)[classDay] ]][nonzero_entry,nonzero_entry]+ 
          get.adjacency(sample_k_regular(ClassSize[i], k_reg))
      }
    
    
    ## case 2
    }else if(NumOfStudent == 1 & sum(classDay) > 1 ) {
      
      day_ind <- seq_len(7)[classDay][sample(sum(classDay),1)]
      if(ClassSize[i] <= k_reg & ClassSize[i] > 1) {
        
        DailyNetSparse[[ day_ind ]][nonzero_entry,nonzero_entry] <-  
          DailyNetSparse[[ day_ind ]][nonzero_entry,nonzero_entry] + 1 - 
          diag(length(nonzero_entry))
      }else if(ClassSize[i]>k_reg) {
        DailyNetSparse[[ day_ind ]][nonzero_entry,nonzero_entry] <-  
          DailyNetSparse[[ day_ind ]][nonzero_entry,nonzero_entry]+ 
          get.adjacency(sample_k_regular(ClassSize[i], k_reg))
      }
    
    
    ## case 3
    } else if(NumOfStudent>1 & sum(classDay) >1) {
      
      student_day <- split(1:NumOfStudent, cut(1:NumOfStudent, sum(classDay)))
      random_order <- sample(NumOfStudent,NumOfStudent)
      
      if(length(instructorList_ind) == 0) {
        
        nonzero_student <- nonzero_entry
        nonzero_instructor <- NULL
        
      }else if(length(instructorList_ind) >= 1) {
        
        instructor_id <- instructorList[instructorList_ind, 2]
        nonzero_instructor <- 
          nonzero_entry[as.integer(names(nonzero_entry)) %in% instructor_id]
        nonzero_student <- setdiff(nonzero_entry, nonzero_instructor)
        
      }
      
      nonzero_student <- nonzero_student[random_order]
      for(j in 1:sum(classDay)) {
        num_check <- length(instructorList_ind)+lengths(student_day)[j]
        nonzero_day <- c(nonzero_instructor,nonzero_student[student_day[[j]]])
        
        if(num_check <= k_reg & num_check > 1) {
          DailyNetSparse[[ seq_len(7)[classDay][j]]][nonzero_day,nonzero_day] <- 
            DailyNetSparse[[seq_len(7)[classDay][j] ]][nonzero_day,nonzero_day]+ 
            1 - diag(length(nonzero_day))
          
        } else if(num_check > k_reg ) {
          
          DailyNetSparse[[seq_len(7)[classDay][j] ]][nonzero_day,nonzero_day] <- 
            DailyNetSparse[[ seq_len(7)[classDay][j]]][nonzero_day,nonzero_day]+ 
            get.adjacency(sample_k_regular(num_check,k_reg))
        }
      }
    }
    
  }
  return(DailyNetSparse)
}

#' --------------------------------------------------------------------
#' Generate a graphml file from the network object
#' 
#' @param DailyNet - network object 
#' @param full_info_id - vector of IDs 
#' @param output_path - character string
#' @param NetName - Network name used for the output file name
#' @return None
#'
getgrapgml <- function(DailyNet, full_info_id, output_path, NetName)
{
  for(i in 1:length(DailyNet)){
    temp <- graph_from_adjacency_matrix(DailyNet[[i]],mode = "undirected",weighted = TRUE,diag = FALSE) 
    vertex_attr(temp, "full_info_id" ) <- as.character(full_info_id)
    write_graph(temp,paste0(output_path, NetName, Day[i],".graphml"),"graphml")
  }
}


## =======================================================================
## Main code

## Create Affiliation Sparse Matrix
AffiMat <- as.data.frame.matrix( table(faculty_students_course) )
AffiMat <- as(as.matrix(AffiMat),"sparseMatrix") 

## Get Daily Schedule
Day <- c('M','T','W','R','F','S','U')
ClassSch <- list()
for (i in 1:length(Day)) {
  # select courses in specified Day
  day <- Day[i]
  ClassSch[[i]] <- unique(course_info$courseID[course_info[,1+i]])
} 

k_reg <- 10
ClassNet <- getClassNet(k_reg, AffiMat)
DailyNet <- getDailyNet(ClassSch, ClassNet) 
DailyNetFull <- getDailyNetFull(ClassSch,AffiMat)
DailyNetSparse <- getDailyNetSparse(k_reg, AffiMat, ClassSch, faculty_courses)

full_info_id <-  as.numeric(rownames( AffiMat ))

## Output networks
getgrapgml(DailyNet, full_info_id, output_path, NetName = "ClassNetwork")
getgrapgml(DailyNetSparse, full_info_id, output_path, NetName = "ClassNetPlatoon")





