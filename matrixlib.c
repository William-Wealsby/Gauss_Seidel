#include<stdio.h>

double Det2(double* matrix){ 
    return  matrix[0]*matrix[3]-matrix[1]*matrix[2]; 
}

int sign(int n){
    return 1-2*(n%2);
}

void submatrix(double* matrix,double* sub_mat,int n,int i,int j){ //matrix dimensions n^2, sub mat dimensions (n-1)^2
    int m=0;             
    for (int k=0;k<n*n;k++){        //remove ith row and jth column
        if (k%n!=i&&(k<n*j||k>=n*(j+1))){
            sub_mat[m]=matrix[k];
            m++;
        }
    }
}

void outerProduct(double* vect1,double* vect2,double* product,int n,int m){//vect1 length n, vect2 length m
    for (int i=0;i<n;i++){
        for (int j=0;j<m;j++){
            product[i*n+j] = vect1[i]*vect2[j];
        }
    } 
}

double innerProduct(double* vect1, double* vect2, int n){
    double product=0;
    for (int i=0;i<n;i++){
        product += vect1[i]*vect2[i];
    }
    return product;
}

void Inv2(double* matrix,double* inv_mat){
    double det = Det2(matrix);
    inv_mat[0] = matrix[3]/det;
    inv_mat[1] = -matrix[1]/det;
    inv_mat[2] = -matrix[2]/det;
    inv_mat[3] = matrix[0]/det;
    //for (int i=0;i<4;i++){printf(" %f", inv_mat[i]);}
}

void transpose(double* matrix,double* t_mat,int m,int n){
    for (int i=0;i<n;i++){
        for (int j=0;j<m;j++){
            t_mat[i*m+j] = matrix[j*n+i];
        }
    }
}

