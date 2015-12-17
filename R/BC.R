#' Biogeographic connectedness
#'
#' Biogeographic connectedness (sensu Sidor et al.).
#'
#' @param taxon_locality_matrix Presence-absence (1-0) matrix of taxa (rows) in localities (columns).
#' @param tree Optional time-scaled phylogenetic tree of taxa in presence-absence matrix.
#' @param count.nodes Option to count nodes rather than use branch-lengths.
#' @param permute.tree Option to permute random trees by shuffling tips.
#' @param bootstrap Option to bootstrap occurrences.
#' @param jackknife Option to jacknife occurrences.
#' @param resample.replicates Number of replicates to use for bootstrapping or jacknifing.
#' @param tree.replicates Number of replicates to use if permuting trees.
#' @param k Value to use to shorten tree and avoid increasing phylogenetic diversity artefact.
#'
#' @return
#'
#' \item{BC_observed}{Empirical results for the: number of localities (\code{N_localities}), number of taxa (\code{N_taxa}), number of links (\code{N_links}), number of endemics (\code{N_endemics}), mean occurrences per locality (\code{Mean_occurrences_per_locality}), and biogeographic connectedness (\code{Biogeographic_connectedness}).}
#' \item{BC_resampled}{Output for the resampled (bootstrap or jacknife) data if requested.}
#' \item{BC_permuted}{Output for the permuted data if requested.}
#'
#' @details
#'
#' Nothing yet.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com} and Richard J. Butler \email{butler.richard.j@@gmail.com}
#'
#' @export
#'
#' @examples
#'
#' # Nothing yet.
BC <- function(taxon_locality_matrix, tree = NULL, count.nodes = FALSE, permute.tree = TRUE, bootstrap = FALSE, jackknife = FALSE, resample.replicates = 1000, tree.replicates = 1000, k = 1000) {
    
    # Function to get BC for a single set of data:
    BC_single <- function(taxon_locality_matrix) {
        
        # Calculate L (number of localities):
        L <- ncol(taxon_locality_matrix)
        
        # Calculate O (number of links, i.e. sum of matrix):
        O <- sum(taxon_locality_matrix)
        
        # Calculate N (number of taxa):
        N <- nrow(taxon_locality_matrix)
        
        # Calculate BC:
        BC <- (O - N) / ((L * N) - N)
        
        # Return BC:
        return(BC)
        
    }
    
    # Function to generate phylogenetically-modified taxon-locality matrix:
    PhyloTLMatrix <- function(taxon_locality_matrix, phylo_dist_matrix) {
        
        # Make copy of binary presence-absence matrix:
        pa_matrix <- taxon_locality_matrix
        
        # If not all taxa are present at all localitites (i.e., if BC is not 1):
        if(!(ncol(taxon_locality_matrix) * nrow(taxon_locality_matrix)) == sum(taxon_locality_matrix)) {
            
            # For each locality:
            for(i in 1:ncol(pa_matrix)) {
                
                # Get list of taxa not at ith locality:
                taxa_not_at_locality <- rownames(pa_matrix)[grep(TRUE, pa_matrix[, i] == 0)]
                
                # Only continue if locality has missing taxa:
                if(!is.null(taxa_not_at_locality)) {
                    
                    # For each taxon not at the ith locality:
                    for(j in taxa_not_at_locality) {
                        
                        # Find taxa at localities:
                        taxa_at_locality <- rownames(pa_matrix)[grep(TRUE, pa_matrix[, i] == 1)]
                        
                        # Get phylogenetic similarity of most closely related taxon from another locality (will be 1 if same taxon is found elsewhere):
                        taxon_locality_matrix[j, i] <- max(phylo_dist_matrix[j, taxa_at_locality])
                        
                    }
                    
                }
                
            }
            
        }
        
        # Return modified taxon-locality matrix:
        return(taxon_locality_matrix)
        
    }
    
    # Function to convert phylogeny to 0 to 1 similarity matrix:
    RescaledPhylogeneticSimilarity <- function(tree, k = k) {
        
        # Make root age maximu path length:
        tree$root.time <- max(diag(vcv(tree)))
        
        # Get node ages (expressed as time from tips of tree):
        node.ages <- GetNodeAges(tree)
        
        # Find nodes older than k:
        too.old.nodes <- as.numeric(names(which(node.ages > k)))
        
        # If there are nodes that are too old:
        if(length(too.old.nodes) > k) {
            
            # Identify edges with at least one (>1) too old node:
            edges.to.change <- apply(cbind(node.ages[tree$edge[, 1]], node.ages[tree$edge[, 2]]) > k, 1, sum)
            
            # Isolate Zero-Length Branches (ZLBs):
            ZLBs <- which(edges.to.change == 2)
            
            # Make ZLBs ZLBs:
            if(length(ZLBs > 0)) tree$edge.length[ZLBs] <- 0
            
            # Isolate branches to shorten:
            shorten <- which(edges.to.change == 1)
            
            # Shorten branches:
            tree$edge.length[shorten] <- k - node.ages[tree$edge[shorten, 2]]

        }
        
        # Get phylogenetic distance matrix:
        phylo_dist <- cophenetic.phylo(tree)
        
        # Rescale from 0 to 1:
        phylo_dist <- phylo_dist / max(phylo_dist)
        
        # Invert to make dissimilarity matrix a similarity matrix:
        phylo_dist <- 1 - phylo_dist
        
        # Return output:
        return(phylo_dist)
        
    }
    
    # If counting nodes:
    if(count.nodes && !is.null(tree)) {
        
        # Start by making all branches length 1:
        tree$edge.length <- rep(1, nrow(tree$edge))
        
        # Now make just the terminals 0.5:
        tree$edge.length[match(1:Ntip(tree), tree$edge[, 2])] <- 0.5
        
    }
    
    # Get unique lengths:
    path.lengths <- unique(diag(vcv(tree)))
    
    # If there is more than one different path length (possible non-ultrametric tree):
    if(length(path.lengths) > 1) {
        
        # For each first value:
        for(i in 1:(length(path.lengths) - 1)) {
            
            # For each second value:
            for(j in (i + 1):length(path.lengths)) {
                
                # Check to see if within floating point error and stop if not:
                if(!all.equal(path.lengths[i], path.lengths[j])) stop("Tree must be ultrametric.")
                
            }
            
        }
        
    }

# THIS BIT COULD PROBABLY BE IMPROVED! (CREATES SLIGHT, NON SIGNIFICANT DOWNWARD TREND.) ALSO, MAYBE THIS IS REMOVING REAL TEMPORAL TREND, I.E., OVER-CORRECTING.
    
    # If not counting nodes (and not a star tree):
    if(!count.nodes && !is.null(tree) && var(vcv(tree)[lower.tri(vcv(tree))]) / mean(diag(vcv(tree))) > 0) {
        
        # Rescale tree using delta transformation to avoid time-bias toward higher BC later:
        tree <- rescale(tree, "delta",  var(vcv(tree)[lower.tri(vcv(tree))]) / mean(diag(vcv(tree))))
        
    }
    
    
    
    
    
    
    # Ensure taxon-locality matrix only contains taxa shared by tree and matrix (if there is a tree):
    if(!is.null(tree)) taxon_locality_matrix <- taxon_locality_matrix[intersect(tree$tip.label, rownames(taxon_locality_matrix)), ]
    
# IS THERE A REASON TO RETAIN ZEROES? IT SEEMS LIKE SAMPLING MIGHT BE AN ISSUE HERE.
    
    # Remove zero value rows and columns:
    taxon_locality_matrix <- taxon_locality_matrix[which(apply(taxon_locality_matrix, 1, sum) > 0), which(apply(taxon_locality_matrix, 2, sum) > 0)]
    
    # If not permuting random trees set output to NULL (will be overwritten later if permuting):
    phyloBC.distribution <- NULL
    
    # If using a tree:
    if(!is.null(tree)) {
        
        # Ensure tree only includes taxa in matrix:
        if(length(setdiff(tree$tip.label, rownames(taxon_locality_matrix))) > 0) tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(taxon_locality_matrix)))
        
        # If permuting random trees:
        if(permute.tree) {
            
            # Create list variable to store random trees:
            rand.trees <- list()
            
            # For each replicate:
            for(i in 1:tree.replicates) {
                
                # Store observed tree as ith random tree:
                rand.trees[[i]] <- tree
                
                # Shuffle taxon names on ith tree:
                rand.trees[[i]]$tip.label <- sample(rand.trees[[i]]$tip.label)
                
            }
            
            # Make variable a multiPhylo object:
            class(rand.trees) <- "multiPhylo"
            
            # Get phylogenetic distance matrices for each random tree:
            rand.distances <- lapply(rand.trees, RescaledPhylogeneticSimilarity)
            
            # Create list to store phylogenetically modified presence-absence matrices:
            tl.matrices <- list()
            
            # For the ith tree store the phylogenetically-modified presence-absence matrix:
            for(i in 1:tree.replicates) tl.matrices[[i]] <- PhyloTLMatrix(taxon_locality_matrix, rand.distances[[i]])
            
            # Get distribution of phylogenetic BCs for randomly-shuffled trees:
            phyloBC.distribution <- unlist(lapply(tl.matrices, BC_single))
            
        }
        
        # Get phylogenetic similarity:
        phylo_dist <- RescaledPhylogeneticSimilarity(tree)
        
        # Get phylogenetically modified taxon-locality matrix:
        taxon_locality_matrix <- PhyloTLMatrix(taxon_locality_matrix, phylo_dist)
        
    }
    
    # If requested both bootstrapping and jackknifing:
    if(bootstrap && jackknife) {
        
        # Print warning:
        cat("WARNING: Can't perform both boostrapping and jacknifing (can lead to empty presence-absence matrix) so just doing latter.\n")
        
        # Set bootstrap to FALSE:
        bootstrap <- FALSE
        
    }
    
    # If performing either bootstrapping or jackknifing:
    if(bootstrap || jackknife) {
        
        # Create vectors to store results of each replicate:
        BC_outs <- Ns <- Os <- Ls <- Mean_occurrences <- N_endemics <- BC_values <- rep(NA, resample.replicates)
        
        # For each replicate:
        for(i in 1:resample.replicates) {
            
            # Create modifiable taxon locality matrix from inputted matrix:
            taxon_locality_matrix2 <- taxon_locality_matrix
            
            # Bootstrap:
            if(bootstrap) {
                
                # Bootstrap matrix:
                taxon_locality_matrix2 <- taxon_locality_matrix2[, sample(ncol(taxon_locality_matrix2), replace = TRUE)]
                
                # Remove empty rows/columns:
                taxon_locality_matrix2 <- taxon_locality_matrix2[which(apply(taxon_locality_matrix2, 1, sum) > 0), which(apply(taxon_locality_matrix2, 2, sum) > 0)]
                
            }
            
            # Jackknife:
            if(jackknife) {
                
                # Jacknife matrix:
                taxon_locality_matrix2 <- taxon_locality_matrix2[sample(nrow(taxon_locality_matrix2), replace = TRUE), ]
                
                # Remove empty rows/columns:
                taxon_locality_matrix2 <- taxon_locality_matrix2[which(apply(taxon_locality_matrix2, 1, sum) > 0), which(apply(taxon_locality_matrix2, 2, sum) > 0)]
                
            }
            
            # If there is still data and at least two taxa and two localitites:
            if(is.matrix(taxon_locality_matrix2) && length(taxon_locality_matrix2) > 0) {
                
                # Calculate and store the BC:
                BC_values[i] <- BC_single(taxon_locality_matrix2)
                
                # Record number of endemics:
                N_endemics[i] <- sum(apply(taxon_locality_matrix2, 1, sum) == 1)
                
                # Record mean number of occurrences per locality:
                Mean_occurrences[i] <- mean(apply(taxon_locality_matrix2, 2, sum))
                
                # Calculate L (number of localitites):
                Ls[i] <- ncol(taxon_locality_matrix2)
                
                # Calculate O (number of links, i.e. sum of matrix):
                Os[i] <- sum(taxon_locality_matrix2)
                
                # Calculate N (number of taxa):
                Ns[i] <- nrow(taxon_locality_matrix2)
                
            }
            
        }
        
        # Compile output:
        BC_outs <- round(cbind(c(1:resample.replicates), Ls, Ns, Os, N_endemics, Mean_occurrences, BC_values), 2)
        
        # Add names:
        colnames(BC_outs) <- c("Replicate_number", "N_localities", "N_taxa", "N_links", "N_endemics", "Mean_occurrences_per_locality", "Biogeographic_connectedness")
        
        # Case if not boostrapping or jacknifing:
    } else {
        
        # Set output bootstrap/jackknife variable to NULL:
        BC_outs <- NULL
        
    }
    
    # Calculate L (number of localitites):
    L <- ncol(taxon_locality_matrix)
    
    # Calculate O (number of links, i.e. sum of matrix):
    O <- sum(taxon_locality_matrix)
    
    # Calculate N (number of taxa):
    N <- nrow(taxon_locality_matrix)
    
    # Record number of endemics:
    N_endemic <- sum(apply(taxon_locality_matrix, 1, sum) == 1)
    
    # Record mean occurrences:
    Mean_occurrence <- mean(apply(taxon_locality_matrix, 2, sum))
    
    # Calculate BC:
    BC <- BC_single(taxon_locality_matrix)
    
    # Compile output:
    BC_out <- round(c(L, N, O, N_endemic, Mean_occurrence, BC), 2)
    
    # Add names:
    names(BC_out) <- c("N_localities", "N_taxa", "N_links", "N_endemics", "Mean_occurrences_per_locality", "Biogeographic_connectedness")
    
    # Record BC for observed and boostrapped and/or jackknifed data:
    out <- list(BC_out, BC_outs, phyloBC.distribution)
    
    # Add names to output:
    names(out) <- c("BC_observed", "BC_resampled", "BC_permuted")
    
    # Return output:
    return(out)
    
}
