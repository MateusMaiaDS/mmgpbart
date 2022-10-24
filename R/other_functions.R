# Function to create the the function K that will be used
# in a Gaussian process (Andrew's Version)
kernel_function <- function(squared_distance_matrix_phi, nu) {

        # Calculating the square matrix
        kernel_matrix <- (exp(-squared_distance_matrix_phi)) / nu

        # Case nu = 0
        if(nu == 0 || nu > 1e13){
                kernel_matrix <- matrix(0, nrow = dim(squared_distance_matrix_phi)[1],
                                        ncol = dim(squared_distance_matrix_phi)[2])
        }
        # Getting the kernel matrix
        return(kernel_matrix)
}
