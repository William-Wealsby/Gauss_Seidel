#include<stdio.h>
#include<stdlib.h>

// runs 1 iteration of the gauss seidel method
void gs_iteration(double* upper,double* lower,int max_i,double* bvect,double* return_vect){
            
    double temp;
    //x_n+1 = lower*x_n+1 + upper*x + d
    for (int i=0;i<max_i;i++){
        
        temp =0;
        
        // loop to perform matrix multiplication        
        for (int j=0;j<max_i;j++){
            temp += (-1)*upper[j+max_i*i]*return_vect[j] + (-1)*lower[j+max_i*i]*return_vect[j];  
        }
        return_vect[i] = temp + bvect[i];

    } 

}

/* "return_matrix" is filled with zeros the lower indices
 * and the elements of "matrix" that are in the upper indices 
 * from "matrix" are copied over to "return_matrix"
 */
void upper(double* matrix,int max_i,double* return_matrix){
    for (int i=0;i<max_i;i++){
        for (int j=0;j<max_i;j++){
            if (i<j){return_matrix[j+i*max_i]=matrix[j+max_i*i];}
            else {return_matrix[j+max_i*i]=0;}
        }
    }
}

/* return matrix "lower" with same dimensions as "matrix"
 * "return_matrix" has only the lower indices from "matrix"
 * with zeros in the upper indices 
*/
void lower(double* matrix,int max_i,double* return_matrix){
    for (int i=0;i<max_i;i++){
        for (int j=0;j<max_i;j++){
            if (i>j){return_matrix[j+max_i*i]=matrix[j+max_i*i];}
            else {return_matrix[j+max_i*i]=0;}
        }
    }

}

// normalise the matrix and bvect so the diagonal elements of the matrix are all 1
void normalise_system(double* matrix,int max_i,double* bvect){
    for (int i=0;i<max_i;i++){
        double divisor = matrix[i+max_i*i];
        bvect[i] = bvect[i]/divisor; 
        for (int j=0;j<max_i;j++){matrix[j+max_i*i]/=divisor;}
    }
}
 
// main logic of gauss seidel interation method, allocates memory then calls gs_interation
void gauss_seidel(double* matrix,int max_i,double* bvect,double* return_vect,int itrs){
    
    double* upper_mat;
    double* lower_mat;
    upper_mat = (double*)malloc(max_i*max_i*sizeof(double));
    lower_mat = (double*)malloc(max_i*max_i*sizeof(double));

    for (int n=0;n<max_i;n++){return_vect[n]=bvect[n]/matrix[n*max_i+n];}

    normalise_system(matrix,max_i,bvect);
    upper(matrix,max_i,upper_mat);
    lower(matrix,max_i,lower_mat);

    printf("upper: %f \n",upper_mat[1]);
    printf("lower: %f \n",lower_mat[2]);
    printf("bvect[0]: %f \n", bvect[0]);
    printf("bvect[1]: %f \n", bvect[1]);

    // loop for n iterations of the gs methods
    for (int n=0;n<itrs;n++){
        gs_iteration(upper_mat,lower_mat,max_i,bvect,return_vect);
    }

    free(upper_mat);
    free(lower_mat);

}

int main(int argc, char** argv){
    
    //
    int size = 2;
    double* matrix;
    double* bvect;
    double* return_vect;
    int itrs = 5;
    
    // allocate memory, for vectors and matrices
    matrix = (double*)malloc(size*size*sizeof(double));
    bvect = (double*)malloc(size*sizeof(double));
    return_vect = (double*)malloc(size*sizeof(double));
    
    // set values for bvect and matrix
    matrix[0] = 2;
    matrix[1] = 1;
    matrix[2] = 1;
    matrix[3] = 2;

    bvect[0] = 4;
    bvect[1] = 5;

    // run guass_seidel algorithm for calculating the vector
    gauss_seidel(matrix,size,bvect,return_vect,itrs);
  
    // output return vector
    printf("final return\n");
    for (int i=0;i<size;i++){
        printf("%f ",return_vect[i]);
    }
    printf("\n");



    // free memory used by vectors and matrices 
    free(matrix);
    free(bvect);
    free(return_vect);

    return 0;
}
