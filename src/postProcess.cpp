#include <RcppArmadillo.h>
#include <set>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::IntegerMatrix create_co_occurrence_matrix_cpp(Rcpp::IntegerMatrix clustering_labels_matrix) {

  // Convert Rcpp::IntegerMatrix to Armadillo integer matrix (imat) for efficient operations.
  // This creates a fast, constant-time reference if possible, or a copy if necessary.
  const arma::imat labels_mat = Rcpp::as<arma::imat>(clustering_labels_matrix);

  const arma::uword N = labels_mat.n_rows; // Number of observations
  const arma::uword K = labels_mat.n_cols; // Number of clustering runs

  // Initialize the N x N co-occurrence matrix with zeros (using integer type for counts).
  // Note: Rcpp::IntegerMatrix is used for the return type, but we use arma::imat
  // for high-speed computation.
  arma::imat co_occurrence_matrix(N, N, arma::fill::zeros);

  // Iterate through each clustering run (column k)
  for (arma::uword k = 0; k < K; ++k) {

    // 1. Extract the current column (cluster labels for the run)
    //    arma::icolvec is the Armadillo integer column vector type
    arma::icolvec current_labels = labels_mat.col(k);

    // 2. Find unique cluster IDs in this run using std::set
    std::set<int> unique_cluster_ids;
    for (arma::uword i = 0; i < N; ++i) {
      unique_cluster_ids.insert(current_labels(i));
    }

    // 3. For each unique cluster, find members and update the co-occurrence matrix
    for (int cluster_id : unique_cluster_ids) {

      // Find all observation indices that belong to the current cluster_id.
      // arma::find returns a uvec (unsigned integer vector) of zero-based indices.
      arma::uvec obs_indices = arma::find(current_labels == cluster_id);

      // Optimization: Skip clusters with zero or one member if they are not
      // the focus, but we MUST include the single member case to increment the diagonal (self-co-occurrence).
      if (obs_indices.n_elem > 0) {

        // Key Efficiency Step: Increment the submatrix defined by obs_indices
        // in both rows and columns by 1.
        // This is a single, highly vectorized Armadillo operation.
        // The .submat(indices, indices) syntax is correct for non-contiguous indices.
        co_occurrence_matrix.submat(obs_indices, obs_indices) += 1;
      }
    }
  }

  // Convert the final Armadillo matrix back to an Rcpp::IntegerMatrix to return to R.
  return Rcpp::wrap(co_occurrence_matrix);
}


// [[Rcpp::export]]
Int32 find_ls_optimal_partition(Rcpp::IntegerMatrix co_occurrence_matrix,
                                              Rcpp::IntegerMatrix clustering_labels_matrix) {

  // Convert inputs to Armadillo types for fast computation
  const arma::imat CoocMat = Rcpp::as<arma::imat>(co_occurrence_matrix);
  const arma::imat Z = Rcpp::as<arma::imat>(clustering_labels_matrix);

  const arma::uword N = CoocMat.n_rows;  // Number of observations
  const arma::uword M = Z.n_cols;        // Number of candidate partitions

  if (N != Z.n_rows) {
    Rcpp::stop("The co-occurrence matrix rows (N) must match candidate partitions rows.");
  }

  // Determine the total number of runs (K) used to calculate the Co-occurrence Matrix.
  // This is equal to the maximum value in the matrix (which occurs on the diagonal).
  // The diag_vec() function is faster than searching the entire matrix.
  const double K = (double)CoocMat.diag().max();
  if (K < 1.0) {
    Rcpp::stop("Co-occurrence counts are zero. Cannot calculate PSM.");
  }

  // 1. Normalize the CoocMat to get the Posterior Similarity Matrix (PSM), C.
  //    C[i, j] = P(i, j together).
  arma::mat PSM = arma::conv_to<arma::mat>::from(CoocMat) / K;

  double min_loss = std::numeric_limits<double>::max();
  arma::uword optimal_partition_index = 0; // 0-based index

  // 2. Iterate through each candidate partition (column m)
  for (arma::uword m = 0; m < M; ++m) {

    const arma::icolvec current_partition = Z.col(m);
    double current_loss = 0.0;

    // 3. Calculate the LS loss for the current partition
    //    We only iterate over the unique pairs (i < j), corresponding to the
    //    upper triangle of the similarity matrix.
    for (arma::uword i = 0; i < N; ++i) {
      for (arma::uword j = i + 1; j < N; ++j) {

        // C_ij is the Posterior Similarity (from PSM)
        const double C_ij = PSM(i, j);

        // I_ij is the Indicator for the current candidate partition (0 or 1)
        // 1.0 if subjects i and j are in the same cluster, 0.0 otherwise.
        const double I_ij = (current_partition(i) == current_partition(j)) ? 1.0 : 0.0;

        // Loss is the squared difference: (I_ij - C_ij)^2
        current_loss += std::pow(I_ij - C_ij, 2.0);
      }
    }

    // 4. Check for minimum loss
    if (current_loss < min_loss) {
      min_loss = current_loss;
      optimal_partition_index = m;
    }
  }
  // Return the 1-based index (for R compatibility) and the minimum loss
  return  optimal_partition_index + 1.0;
}



