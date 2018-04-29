// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <stdio.h>
#include <iostream>
#include <RcppArmadillo.h>
#include <math.h>  // for -INFINITY, NAN, isnan()

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
//#include <Rcpp.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R


/////////////////////////
// [[Rcpp::export]]
arma::vec vectorize(arma::mat mat, int axis){
    int row = mat.n_rows;
    int col = mat.n_cols;
    int i,j;
    arma::vec result(row*col);
    for(i=0;i<row;i++){
        for(j=0;j<col;j++){
            if(axis==0)  result(i+j*row) = mat(i,j);  // by column
            else result(j+i*col) = mat(i,j); //by row
        }
    }
    return result;
}


// [[Rcpp::export]]
arma::mat vec2mat(arma::vec vec, int row, int col, int axis){
    int i,j;
    arma::mat result(row, col);
    for(i=0; i<row; i++){
        for(j=0; j<col; j++){
            if(axis==0)  result(i,j) = vec(i+j*row) ; //by column
            else  result(i,j) = vec(j+i*col);     // by row
        }
    }
    return result;
}

///////////////////////////////////////
// [[Rcpp::export]]
arma::vec colsum(arma::mat matrix){
    long row = matrix.n_rows;
    int col = matrix.n_cols;
    long i;
    int j;
    arma::vec columnsum(col);
    
    for(j=0; j<col; j++) columnsum(j) = 0;
    for(i=0; i<row; i++){
        for(j=0; j<col; j++){
            columnsum(j) += matrix(i,j);
        }
    }
    return columnsum;
}

///////////////////////////////////////
// [[Rcpp::export]]
arma::vec rowsum(arma::mat matrix){
    long row = matrix.n_rows;
    int col = matrix.n_cols;
    long i;
    int j;
    arma::vec rowsum(row);
    
    for(i=0; i<row; i++) rowsum(i) = 0;
    for(j=0; j<col; j++){
        for(i=0; i<row; i++){
            rowsum(i) += matrix(i,j);
        }
    }
    return rowsum;
}

//////////////////////////////////////////////////////////////////////////////////
//subset a big matrix by the beginning and ending of row/col index
// [[Rcpp::export]]
arma::mat subsetmatrix(arma::mat rawmat, arma::vec rowindex, arma::vec colindex){
    int row = rowindex(1) - rowindex(0) + 1;
    int col = colindex(1) - colindex(0) + 1;
    arma::mat result(row,col);
    int i,j;
    for(i=0;i<row;i++){
        for(j=0;j<col;j++){
            result(i,j) = rawmat(rowindex(0)+i,colindex(0)+j);
        }
    }
    return result;
}

/////////////////////////////
//compute matrix power
// [[Rcpp::export]]
arma::mat matrixpower(arma::mat oldmat, int power){
    arma::mat newmat = oldmat;
    int i;
    
    for(i=1;i<power;i++) newmat *= oldmat;
    
    return newmat;
}

///////////////////////////
//compute matrix exponential
// [[Rcpp::export]]
arma::mat matrixexp(arma::mat oldmat, double t){
    arma::mat temp = t * oldmat;
    arma::mat newmat = arma::expmat(temp);
    return newmat;
}

/////////////
// [[Rcpp::export]]
double matrixsum(arma::mat mat1, arma::mat mat2){
    double tempsum=0;
    arma::mat temp = mat1 % mat2;
    int row = mat1.n_rows;
    int col = mat1.n_cols;
    int i,j;
    for(i=0;i<row;i++){
        for(j=0;j<col;j++){
            tempsum += temp(i,j);
        }
    }
    return(tempsum);
}

////////////
//integral of matrix exponential
// x and y starts from 0
// [[Rcpp::export]]
arma::mat matrixintegral(arma::mat Q, double interval, int x, int y){
    int M = Q.n_cols;
    arma::mat temp(2*M, 2*M);
    temp.zeros();
    int i,j;
    
    for(i=0;i<M;i++)
        for(j=0;j<M;j++)
            temp(i,j) = Q(i,j);
    
    for(i=M;i<2*M;i++)
        for(j=M;j<2*M;j++)
            temp(i,j) = Q(i-M,j-M);
    
    temp(x, M+y) = 1;
    
    
    arma::mat temp1 = matrixexp(temp, interval);
    //Rcpp::Rcout<<temp1<<std::endl;
    arma::vec rowindex(2);
    arma::vec colindex(2);
    rowindex(0)=0; rowindex(1) = M-1; colindex(0) = M; colindex(1) = 2*M-1;
    arma::mat result = subsetmatrix(temp1,rowindex,colindex);
    return result;
}


////////////////////////////////////////////////////////////
//' pmf for zero-inflated poisson
//' @param p proportion of structural zero's
//' @param theta the poisson mean
//' @param y the observed value
//' @param loga Logical. Whether to return the log probability or not.
//' @return the probability mass of the zero-inflated poisson distribution
//' @export
// [[Rcpp::export]]
double dzip(double p, double theta, int y, bool loga){
    double result;
    if(loga==FALSE){
        if(y==0) result = p + (1-p) * exp(-theta);
        else result = (1-p) * R::dpois(y, theta, loga);
    }
    else{
        if(y==0) result = log(p + (1-p) * exp(-theta));
        else result = log(1-p) + R::dpois(y, theta, loga);
    }
    return result;
}

////////////////////////////////////////////////////////////
//' generate zero-inflated poisson random variables
//' @param n length of the random series
//' @param p proportion of structural zero's
//' @param theta the poisson mean
//' @return a series of zero-inflated poisson random variables
//' @export
// [[Rcpp::export]]
arma::vec rzip(int n, double p, double theta){
    arma::vec result(n);
    int i;
    double u;
    
    for(i=0;i<n;i++){
        u =  Rcpp::runif(1,0,1)(0);
        if(u <= p) result(i) = 0;
        else result(i) = Rcpp::rpois(1, theta)(0);
    }
    return result;
}

//aft model with exponential base: pmf. xval starting from 1.
// [[Rcpp::export]]
double pmf_expbase(double lden, double eb, int xval){
    double result;
    result = exp( -eb*(xval-1)/lden ) - exp(-eb*xval/lden);
    return result;
}

//aft model with exponential base: cdf
// [[Rcpp::export]]
double cdf_expbase(double lden, double eb, int xval){
    double result;
    result = 1 - exp( -eb*xval/lden );
    return result;
}

//aft model with exponential base: random
// [[Rcpp::export]]
int random_expbase(double lden, double eb, int maxt){
    double u = Rcpp::runif(1,0,1)(0);
    int i;
    for(i=1; i<maxt; i++){
        if(cdf_expbase(lden,eb, i)>u) break;
    }
    return i;
}

//shift the poisson to the right by shift
// [[Rcpp::export]]
arma::vec rshiftpois(int n, int theta, int shift){
    arma::vec result;
    result = Rcpp::rpois(n, theta) + shift; /*everything plus one to avoid zero for dwell time dist.*/
    return(result);
}

//shift the poisson to the right by shift
// [[Rcpp::export]]
double dshiftpois(int x, int theta, int shift, bool loga){
    double result;
    int temp = x - shift;
    result = R::dpois(temp, theta, loga); /*everything plus one to avoid zero for dwell time dist.*/
    return result;
}



// [[Rcpp::export]]
double dlogp(int x, double p, bool loga){
    double result;
    if(loga==FALSE) result = -exp(x*log(p))/(x*log(1-p));
    else result = log(-exp(x*log(p))/(x*log(1-p)));
    return result;
}


// [[Rcpp::export]]
arma::vec rlogp(int n, double p){
    arma::vec result(n);
    double u,cumsum;
    int i,j;
    
    for(i=0;i<n;i++){
        u = Rcpp::runif(1,0,1)(0);
        cumsum = 0;
        j = 1;
        while(cumsum<u){
            cumsum += dlogp(j,p, FALSE);
            j++;
        }
        result(i) = j-1;
    }
    
    return result;
}

/*
arma::vec randnorm(int n, double mu, double sd){
    arma::vec result(n);
    double pdf;
    result = Rcpp::rnorm(n,mu,sd);
    pdf = R::dnorm(0,mu,sd,FALSE);
    //Rcpp::Rcout<<pdf<<std::endl;
    return(result);
}
*/


// [[Rcpp::export]]
arma::vec multinomrand(int n, int k, arma::vec prob, arma::vec label){
    arma::vec result(n);
    arma::vec cumprob(k);
    int i,j;
    double u;
    
    cumprob(0) = prob(0);
    for(i=1; i<k;i++){
        cumprob(i) = cumprob(i-1) + prob(i);
    }
    
    
    for(j=0;j<n;j++){
        u = Rcpp::runif(1,0,1)(0);   //to make type match
        for(i=0; i<k;i++){
            if(u<cumprob(i)) {
                result(j) = label(i);
                break;
            }
        }
    }
    return(result);
}

/////////////
// [[Rcpp::export]]
arma::mat hsmm_hmm (arma::mat omega, arma::mat dm, arma::vec mv){
    //each row in dm is a dwell time pmf
    //mv is vector of the length until the first zero in each dwell time distribution
    int m = omega.n_rows;
    int dmrow = dm.n_rows;
    int dmcol = dm.n_cols;
    int dim = arma::sum(mv); // dimension of the final result
    
    int i,j,p,q,mi,rowsum,colsum;
    //double tempsum;
    arma::mat temp(dmrow,dmcol);
    arma::mat ci(dmrow,dmcol);
    arma::mat cim(dmrow,dmcol);
    arma::mat gamma(dim,dim);
    
    //Rcpp::Rcout << dim << std::endl;
    
    for(i=0;i<m;i++){
        mi = mv[i];
        
        for(j=0;j<mi;j++){
            if(j==0) temp(i,j) = 0;
            else temp(i,j) = temp(i,j-1) + dm(i,j-1);
        }
        
        for(j=0;j<mi;j++){
            if(std::abs(1-temp(i,j))>0.000000001) ci(i,j) = dm(i,j)/(1-temp(i,j));
            else ci(i,j) = 1;
            if(1-ci(i,j)>0) cim(i,j)=1-ci(i,j);
            else cim(i,j) = 0;
        }
    }
    
    rowsum = 0;
    
    
    for(i=0; i<m; i++){
        colsum = 0;
        for(j=0; j<m; j++){
            if(i==j){
                if(mv[i]==1) gamma(rowsum,colsum) = cim(i,0);
                else{
                    for(p=0; p<mv[i]; p++){
                        for(q=0; q<mv[j]; q++){
                            if((q-p)==1) gamma(rowsum+p,colsum+q)=cim(i,p);
                            else if((p==mv[i]-1) & (q==mv[j]-1)) gamma(rowsum+p,colsum+q)=cim(i,p);
                            else gamma(rowsum+p,colsum+q)=0;
                        }
                    }
                }
            }
            else{
                for(p=0; p<mv[i]; p++){
                    for(q=0; q<mv[j]; q++){
                        if(q==0) gamma(rowsum+p, colsum+q)=omega(i,j)*ci(i,p);
                        else gamma(rowsum+p, colsum+q)=0;
                    }
                }
                
            }
            colsum += mv[j];
        }
        rowsum += mv[i];
    }
    
    
    return(gamma);
    
}



/////////////

// [[Rcpp::export]]
arma::mat hmm_gen (int dim, int M, int ntimes, arma::vec pi, arma::mat a, arma::vec theta,
                   arma::vec zeroprop){
    
    int i,m,n, prev, curr;
    arma::vec label(M);
    double u;
    
    arma::mat result(dim, 2*ntimes); /*first column for x, second column for state*/
    
    for(m=0;m<M;m++) label(m) = m+1;
    
    
    //starting generation
    for(i=0; i<ntimes; i++){
             result(0,2*i+1) = multinomrand(1, M, pi, label)(0);
        
             curr = result(0,2*i+1) - 1;
             u = Rcpp::runif(1,0,1)(0);
             if(u<=zeroprop(curr)) result(0,2*i) = 0;
             else result(0,2*i) = Rcpp::rpois(1, theta(curr))(0);
     
    }
    
    //iteration
    for(n=1; n<dim; n++){
        
        for(i=0; i<ntimes; i++){
            prev = result(n-1, 2*i+1) - 1;
            result(n,2*i+1) = multinomrand(1, M, a.row(prev).t(), label)(0);  //row to column vetor
            
            curr = result(n,2*i+1) - 1;
            u = Rcpp::runif(1,0,1)(0);
            if(u<=zeroprop(curr)) result(n,2*i) = 0;
            else result(n,2*i) = Rcpp::rpois(1, theta(curr))(0);
            
        }
    }
    
    return(result);
}

////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat getnodeprob_nocov(arma::vec ystar, arma::mat emit){
    long ns = ystar.n_rows; //original numeric series starting from 1
    int d = emit.n_rows;
    
    arma::mat nodeprob(ns, d);
    
    long i;
    int j;
    
    for(i=0; i<ns; i++){
        for(j=0; j<d; j++){
            
                //emit(0,0) is zeroprop for state 1, emit(0,1) is lambda for state 1
                nodeprob(i,j) = dzip(emit(j,0), emit(j,1), ystar(i), FALSE);
          
        }
    }
    return nodeprob;
}

//this function is not correct, not used
////////////////////////////////
/*
arma::mat getnodeprob_cov(arma::vec y, arma::mat x, arma::vec zerovec, arma::mat emitmat, int m){
    long ns = y.n_rows; //original numeric series starting from 1
    arma::mat nodeprob(ns, m);
    arma::vec eta;
    double zeroprop,lambda;
    long i;
    int j;
    
    for(i=0;i<ns;i++){
         //Rcpp::Rcout<<i<<std::endl;
        for(j=0;j<m;j++){
            if(j==0){
                eta = x.row(i) * zerovec ;
                zeroprop = exp(eta(0)) / (1+exp(eta(0)));
                eta = x.row(i) * emitmat.row(i).t();
                lambda = exp(eta(0));
                nodeprob(i,j) = dzip(zeroprop, lambda, y(i), FALSE);
            }
            else{
                eta = x.row(i) * emitmat.row(i).t();
                lambda = exp(eta(0));
                nodeprob(i,j) = dzip(zeroprop, lambda, y(i), FALSE);
            }
          }
    }
    return nodeprob;
}
*/

///////////////////////////////////////////
// [[Rcpp::export]]
arma::mat hmm_cov_gen (arma::vec parm, int M, long dim, int ncolcovpi, arma::mat covpi,
                       int ncolcovtrans, arma::mat covtrans, int ncolcovp1, arma::mat covp1,
                       int ncolcovpois, arma::mat covpois, arma::vec zeroindex){
    
    int i,j,m,n,nextindex,prev,curr;
    double tempsum, u;
    arma::vec label(M);
    
    arma::mat result(dim, 2); /*first column for x, second column for state*/
    
    
    arma::vec pi(M);
    arma::mat a(M,M);
    arma::vec theta(M);
    arma::vec zeroprop(M);

    //retrieve the full parameters
    //prior
    nextindex=0;
    pi(0) = 1;
    tempsum = 1;
    for(m=1; m<M; m++){
        pi(m) = parm(nextindex);
        
        for(j=0;j<ncolcovpi;j++)
            pi(m) += parm(nextindex+j+1)*covpi(0,j) ;
        
        pi(m) = exp(pi(m));
        
        tempsum += pi(m);
        nextindex = nextindex + ncolcovpi + 1;
    }
    
    for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
    
    
    //transition
    for(i=0; i<M; i++){
        
        for(j=0; j<M; j++){
            if(j==0) {
                a(i,j) = 1;
                tempsum = 1;
            }
            else {
                a(i,j) = parm(nextindex);
                
                for(m=0;m<ncolcovtrans;m++)
                    a(i,j) += parm(nextindex+m+1)*covtrans(0,m);
                
                a(i,j) = exp(a(i,j));
                tempsum += a(i,j);
                nextindex = nextindex + ncolcovtrans + 1;
            }
            //Rcpp::Rcout<<a(i,j)<<std::endl;
        }
        
        for(n=0;n<M; n++) a(i,n) = a(i,n) / tempsum;
        
    }
    
    
    for(i=0;i<M;i++){
        if(zeroindex(i)==0) zeroprop(i)=0;
        else{
            zeroprop(i) = parm(nextindex);
            
            for(m=0;m<ncolcovp1;m++)
                zeroprop(i) += parm(nextindex+m+1) * covp1(0,m);
            
            zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
            nextindex = nextindex + ncolcovp1 + 1;
        }
    }
    
    //Rcpp::Rcout << "p1=" << p1 << std::endl;
    
    for(m=0; m<M; m++){
        theta(m) = parm(nextindex);
        for(j=0; j<ncolcovpois;j++){
            theta(m) += parm(nextindex+j+1) * covpois(0,j);
        }
        theta(m) = exp(theta(m));
        nextindex = nextindex + ncolcovpois + 1;
    }
    
    
    //start generation
    
    for(m=0;m<M;m++) label(m) = m+1;
  
    result(0,1) = multinomrand(1, M, pi, label)(0);
    curr = result(0,1) - 1;
    u = Rcpp::runif(1,0,1)(0);
    if(u<=zeroprop(curr)) result(0,0) = 0;
    else result(0,0) = Rcpp::rpois(1, theta(curr))(0);
    
    
    
    //iteration steps
    for(n=1; n<dim; n++){
        
        //still need to retrieve the natural parameters
        nextindex=0;
        pi(0) = 1;
        tempsum = 1;
        for(m=1; m<M; m++){
            pi(m) = parm(nextindex);
            
            for(j=0;j<ncolcovpi;j++)
                pi(m) += parm(nextindex+j+1)*covpi(n,j) ;
            
            pi(m) = exp(pi(m));
            
            tempsum += pi(m);
            nextindex = nextindex + ncolcovpi + 1;
        }
        
        for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
        
        
        //transition
        for(i=0; i<M; i++){
            
            for(j=0; j<M; j++){
                if(j==0) {
                    a(i,j) = 1;
                    tempsum = 1;
                }
                else {
                    a(i,j) = parm(nextindex);
                    
                    for(m=0;m<ncolcovtrans;m++)
                        a(i,j) += parm(nextindex+m+1)*covtrans(n,m);
                    
                    a(i,j) = exp(a(i,j));
                    tempsum += a(i,j);
                    nextindex = nextindex + ncolcovtrans + 1;
                }
                //Rcpp::Rcout<<a(i,j)<<std::endl;
            }
            
            for(m=0;m<M; m++) a(i,m) = a(i,m) / tempsum;
            
        }
        
        
        for(i=0;i<M;i++){
            if(zeroindex(i)==0) zeroprop(i)=0;
            else{
                zeroprop(i) = parm(nextindex);
                
                for(m=0;m<ncolcovp1;m++)
                    zeroprop(i) += parm(nextindex+m+1) * covp1(n,m);
                
                zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
                nextindex = nextindex + ncolcovp1 + 1;
            }
        }
        //Rcpp::Rcout << "p1=" << p1 << std::endl;
        
        for(m=0; m<M; m++){
            theta(m) = parm(nextindex);
            for(j=0; j<ncolcovpois;j++){
                theta(m) += parm(nextindex+j+1) * covpois(n,j);
            }
            theta(m) = exp(theta(m));
            nextindex = nextindex + ncolcovpois + 1;
        }

        
            prev = result(n-1, 1) - 1;
            result(n,1) = multinomrand(1, M, a.row(prev).t(), label)(0);  //row to column vetor
            
            curr = result(n,1) - 1;
            u = Rcpp::runif(1,0,1)(0);
            if(u<=zeroprop(curr)) result(n,0) = 0;
            else result(n,0) = Rcpp::rpois(1, theta(curr))(0);
            
        
    }
  
    return result;
    
}


/////////////

//[[Rcpp::export]]
arma::mat hsmm_gen(int dim, int M, arma::vec pi, arma::vec theta, arma::vec zeroprop,
                   arma::mat omega, arma::vec p, std::string dt_dist){
    
    
    int j,count;
    
    int m,n,prev,curr;
    arma::vec label(M);
    double u;
    arma::mat result(dim, 2); /*first column for x, second column for state*/
    
    for(m=0;m<M;m++) label(m) = m+1;
    
    
    //starting generation
    curr = multinomrand(1, M, pi, label)(0);
    
    //check dt_dist
    if(dt_dist=="log"){
      count = rlogp(1,p(curr-1))(0);
    }else{
        count = rshiftpois(1, p(curr-1), 1)(0);
    }

         for(j=0;j<count;j++) {
             result(j,1) = curr;
             u = Rcpp::runif(1,0,1)(0);
             if(u<=zeroprop(curr-1)) result(j,0)=0;
             else result(j,0) = Rcpp::rpois(1, theta(curr-1))(0);
         }

    n = count;
    
    //iteration
    while(n<dim){
        
        prev = result(n-1,1) - 1;
        curr = multinomrand(1, M, omega.row(prev).t(), label)(0);
        
        //check dt_dist
        if(dt_dist=="log"){
          count = rlogp(1, p(curr-1))(0);
        }else{
          count = rshiftpois(1, p(curr-1), 1)(0);
        }
            
             for(j=0;j<count and n+j<dim; j++) {
                 result(n+j,1) = curr;
                 u = Rcpp::runif(1,0,1)(0);
                 if(u<=zeroprop(curr-1)) result(n+j,0)=0;
                 else result(n+j,0) = Rcpp::rpois(1, theta(curr-1))(0);
             }

        n += count;
    }
    
    return(result);
    
}


//////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat hsmm_cov_gen(arma::vec parm, int M, long dim, std::string dt_dist, arma::vec zeroindex,
                       int ncolcovp, arma::mat covp, int ncolcovpi, arma::mat covpi,
                       int ncolcovomega, arma::mat covomega, int ncolcovp1, arma::mat covp1,
                       int ncolcovpois, arma::mat covpois){
    
        /*parm = logitp from logarithm dist---M, logit(pi)---M-1, logit(p1), log(theta)---M,
         omega---M*(M-2): (1,3)...(1,n)(2,3)...(2,n)(3,2)(3,4)...(3,n).......*/
        
        //intercept column is not included in the covariates, but should be in the parm
        
        //parameters are given alternatively in the order of beta0 and beta1
    
        int i,j,k,m,nextindex,prev,curr;
        long count,n;
        double tempsum,u;
        arma::vec label(M);
    
        arma::mat result(dim, 2); /*first column for x, second column for state*/
    
    
        arma::vec p(M);  //dt_parm
        arma::vec pi(M);
        arma::vec theta(M);
        arma::vec zeroprop(M);

        
        arma::mat omega(M,M);
        //arma::mat a(totalmv,totalmv);
    
        //////
    //retrieve some of the parameters
    nextindex = 0;
    //dwell time
    
    //check dt_dist
    if(dt_dist=="log"){
        for(i=0;i<M;i++){
                p(i) = parm(nextindex);
                for(m=0;m<ncolcovp;m++)
                    p(i) += parm(nextindex+m+1)*covp(0,m);
                p(i) = exp(p(i)) / (1+exp(p(i)));
                
            nextindex = nextindex + ncolcovp + 1;
        }
    }else{
        for(i=0;i<M;i++){
        
                p(i) = parm(nextindex);
                for(m=0;m<ncolcovp;m++)
                    p(i) += parm(nextindex+m+1)*covp(0,m);
                p(i) = exp(p(i)) ;
            
                
            }
            nextindex = nextindex + ncolcovp + 1;
        }
    
    
    
    //recover pi
    pi(0) = 1;
    tempsum = 1;
    for(m=1; m<M; m++){
        pi(m) = parm(nextindex);
        
        for(j=0;j<ncolcovpi;j++)
            pi(m) += parm(nextindex+j+1)*covpi(0,j) ;
        
        pi(m) = exp(pi(m));
        
        tempsum += pi(m);
        nextindex = nextindex + ncolcovpi + 1;
    }
    
    for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
    
 
    
    //////////////////
    //start generation
    for(m=0;m<M;m++) label(m) = m+1;
    
    
    //starting generation
    curr = multinomrand(1, M, pi, label)(0);
    
    //check dt_dist
    if(dt_dist=="log"){
        count = rlogp(1,p(curr-1))(0);
    }else{
        count = rshiftpois(1, p(curr-1), 1)(0);
    }
    
    
    
    /////////////////
    for(k=0;k<count;k++) {
        
        //get some of the parameters in each iteration
        nextindex = 0;
        //dwell time
       
        //check dt_dist
        if(dt_dist=="log"){
            for(i=0;i<M;i++){
                
                    p(i) = parm(nextindex);
                    for(m=0;m<ncolcovp;m++)
                        p(i) += parm(nextindex+m+1)*covp(k,m);
                    p(i) = exp(p(i)) / (1+exp(p(i)));
                    
                nextindex = nextindex + ncolcovp + 1;
            }
        }else{
            for(i=0;i<M;i++){
                
                
                    p(i) = parm(nextindex);
                    for(m=0;m<ncolcovp;m++)
                        p(i) += parm(nextindex+m+1)*covp(k,m);
                    p(i) = exp(p(i)) ;
                   
                nextindex = nextindex + ncolcovp + 1;
            }
            
        }
        
        //recover newtheta,p1, newpi
        pi(0) = 1;
        tempsum = 1;
        for(m=1; m<M; m++){
            pi(m) = parm(nextindex);
            
            for(j=0;j<ncolcovpi;j++)
                pi(m) += parm(nextindex+j+1)*covpi(k,j) ;
            
            pi(m) = exp(pi(m));
            
            tempsum += pi(m);
            nextindex = nextindex + ncolcovpi + 1;
        }
        
        for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
        
        //zeroprops
        
        for(i=0;i<M;i++){
            if(zeroindex(i)==0) zeroprop(i)=0;
            else{
                zeroprop(i) = parm(nextindex);
                
                for(m=0;m<ncolcovp1;m++)
                    zeroprop(i) += parm(nextindex+m+1) * covp1(k,m);
                
                zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
                nextindex = nextindex + ncolcovp1 + 1;
            }
        }
        
        
        //theta
        for(m=0; m<M; m++){
            theta(m) = parm(nextindex);
            for(j=0; j<ncolcovpois;j++){
                theta(m) += parm(nextindex+j+1) * covpois(k,j);
            }
            theta(m) = exp(theta(m));
            nextindex = nextindex + ncolcovpois + 1;
        }
        
        
        //recover omega:
        //if M=2 then assign off-diagonal with 1; else
        //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
        
        omega.zeros();
        if(M==2) {omega(0,1)=1; omega(1,0)=1;}
        else{
            for(i=0;i<M;i++){
                tempsum = 1;
                for(j=2;j<M;j++){
                    omega(i,j) = parm(nextindex); //new
                    for(m=0;m<ncolcovomega;m++)
                        omega(i,j) += parm(nextindex+m+1)*covomega(k,m);
                    //omega(i,j) = exp(omega(i,j))/(1+exp(omega(i,j)));
                    omega(i,j) = exp(omega(i,j)); //new
                    tempsum += omega(i,j); //new
                    nextindex = nextindex + ncolcovomega + 1;
                }
                
                for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //new
                if(i==0) omega(i,1) = 1 / tempsum;
                else omega(i,0) = 1 / tempsum;
                //tempsum = 1 - arma::sum(omega.row(i)); //prevent negative values
                //if(i==0) omega(i,1) = MAX(0,tempsum);
                //else omega(i,0) = MAX(0,tempsum);
            }
            
            for(i=2;i<M;i++){
                for(j=2;j<=i;j++){
                    omega(i,j-1) = omega(i,j);
                    if(i==j) omega(i,j)=0;
                }
            }
            
            
        }//end section for zeroinflated poisson distribution
        
        
    
        result(k,1) = curr;
        u = Rcpp::runif(1,0,1)(0);
        if(u<=zeroprop(curr-1)) result(k,0)=0;
        else result(k,0) = Rcpp::rpois(1, theta(curr-1))(0);
    }
    
    n = count;
  
  

    //////////////////////////////////////////////////////
    //iteration
    while(n<dim){
        
        //retrieve the parameters
        nextindex = 0;
        //dwell time
      
        //check dt_dist
        if(dt_dist=="log"){
            for(i=0;i<M;i++){
                
                    p(i) = parm(nextindex);
                    for(m=0;m<ncolcovp;m++)
                        p(i) += parm(nextindex+m+1)*covp(n,m);
                    p(i) = exp(p(i)) / (1+exp(p(i)));
                   
                nextindex = nextindex + ncolcovp + 1;
            }
        }else{
            for(i=0;i<M;i++){
                
                    p(i) = parm(nextindex);
                    for(m=0;m<ncolcovp;m++)
                        p(i) += parm(nextindex+m+1)*covp(n,m);
                    p(i) = exp(p(i)) ;
                   
                nextindex = nextindex + ncolcovp + 1;
            }
            
        }
        
        //recover newtheta,p1, newpi
        pi(0) = 1;
        tempsum = 1;
        for(m=1; m<M; m++){
            pi(m) = parm(nextindex);
            
            for(j=0;j<ncolcovpi;j++)
                pi(m) += parm(nextindex+j+1)*covpi(n,j) ;
            
            pi(m) = exp(pi(m));
            
            tempsum += pi(m);
            nextindex = nextindex + ncolcovpi + 1;
        }
        
        for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
        
        //zeroprops
        
        for(i=0;i<M;i++){
            if(zeroindex(i)==0) zeroprop(i)=0;
            else{
                zeroprop(i) = parm(nextindex);
                
                for(m=0;m<ncolcovp1;m++)
                    zeroprop(i) += parm(nextindex+m+1) * covp1(n,m);
                
                zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
                nextindex = nextindex + ncolcovp1 + 1;
            }
        }
        
        
        
        //theta
        for(m=0; m<M; m++){
            theta(m) = parm(nextindex);
            for(j=0; j<ncolcovpois;j++){
                theta(m) += parm(nextindex+j+1) * covpois(n,j);
            }
            theta(m) = exp(theta(m));
            nextindex = nextindex + ncolcovpois + 1;
        }
        
        //then recover omega for this iteration:
        //if M=2 then assign off-diagonal with 1; else
        //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
        
        omega.zeros();
        if(M==2) {omega(0,1)=1; omega(1,0)=1;}
        else{
            for(i=0;i<M;i++){
                tempsum = 1;
                for(j=2;j<M;j++){
                    omega(i,j) = parm(nextindex); //new
                    for(m=0;m<ncolcovomega;m++)
                        omega(i,j) += parm(nextindex+m+1)*covomega(n,m);
                    //omega(i,j) = exp(omega(i,j))/(1+exp(omega(i,j)));
                    omega(i,j) = exp(omega(i,j)); //new
                    tempsum += omega(i,j); //new
                    nextindex = nextindex + ncolcovomega + 1;
                }
                
                for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //new
                if(i==0) omega(i,1) = 1 / tempsum;
                else omega(i,0) = 1 / tempsum;
                //tempsum = 1 - arma::sum(omega.row(i)); //prevent negative values
                //if(i==0) omega(i,1) = MAX(0,tempsum);
                //else omega(i,0) = MAX(0,tempsum);
            }
            
            for(i=2;i<M;i++){
                for(j=2;j<=i;j++){
                    omega(i,j-1) = omega(i,j);
                    if(i==j) omega(i,j)=0;
                }
            }
            
            
            
        }//end section for zeroinflated poisson distribution
        
        prev = result(n-1,1) - 1;
        //sample from the previous omega
        curr = multinomrand(1, M, omega.row(prev).t(), label)(0);
        
        //check dt_dist
        if(dt_dist=="log"){
            count = rlogp(1, p(curr-1))(0);
        }else{
            count = rshiftpois(1, p(curr-1), 1)(0);
        }
        
        
       
        

        
        
        ////////////////
        for(k=0;k<count and n+k<dim; k++) {
            
            
            //get some of the parameters in each iteration
            nextindex = 0;
            //dwell time
            
            //check dt_dist
            if(dt_dist=="log"){
                for(i=0;i<M;i++){
                  
                        p(i) = parm(nextindex);
                        for(m=0;m<ncolcovp;m++)
                            p(i) += parm(nextindex+m+1)*covp(n+k,m);
                        p(i) = exp(p(i)) / (1+exp(p(i)));
                       
                    nextindex = nextindex + ncolcovp + 1;
                }
            }else{
                for(i=0;i<M;i++){
                    
                    p(i) = parm(nextindex);
                        for(m=0;m<ncolcovp;m++)
                            p(i) += parm(nextindex+m+1)*covp(n+k,m);
                        p(i) = exp(p(i)) ;
                        
                    nextindex = nextindex + ncolcovp + 1;
                }
                
            }
            
            //recover newtheta,p1, newpi
            pi(0) = 1;
            tempsum = 1;
            for(m=1; m<M; m++){
                pi(m) = parm(nextindex);
                
                for(j=0;j<ncolcovpi;j++)
                    pi(m) += parm(nextindex+j+1)*covpi(n+k,j) ;
                
                pi(m) = exp(pi(m));
                
                tempsum += pi(m);
                nextindex = nextindex + ncolcovpi + 1;
            }
            
            for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
            
            //zeroprops
            
            for(i=0;i<M;i++){
                if(zeroindex(i)==0) zeroprop(i)=0;
                else{
                    zeroprop(i) = parm(nextindex);
                    
                    for(m=0;m<ncolcovp1;m++)
                        zeroprop(i) += parm(nextindex+m+1) * covp1(n+k,m);
                    
                    zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
                    nextindex = nextindex + ncolcovp1 + 1;
                }
            }
            
            
            //theta
            for(m=0; m<M; m++){
                theta(m) = parm(nextindex);
                for(j=0; j<ncolcovpois;j++){
                    theta(m) += parm(nextindex+j+1) * covpois(n+k,j);
                }
                theta(m) = exp(theta(m));
                nextindex = nextindex + ncolcovpois + 1;
            }
            
            
            //recover omega:
            //if M=2 then assign off-diagonal with 1; else
            //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
            
            omega.zeros();
            if(M==2) {omega(0,1)=1; omega(1,0)=1;}
            else{
                for(i=0;i<M;i++){
                    tempsum = 1;
                    for(j=2;j<M;j++){
                        omega(i,j) = parm(nextindex); //new
                        for(m=0;m<ncolcovomega;m++)
                            omega(i,j) += parm(nextindex+m+1)*covomega(n+k,m);
                        //omega(i,j) = exp(omega(i,j))/(1+exp(omega(i,j)));
                        omega(i,j) = exp(omega(i,j)); //new
                        tempsum += omega(i,j); //new
                        nextindex = nextindex + ncolcovomega + 1;
                    }
                    
                    for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //new
                    if(i==0) omega(i,1) = 1 / tempsum;
                    else omega(i,0) = 1 / tempsum;
                    //tempsum = 1 - arma::sum(omega.row(i)); //prevent negative values
                    //if(i==0) omega(i,1) = MAX(0,tempsum);
                    //else omega(i,0) = MAX(0,tempsum);
                }
                
            for(i=2;i<M;i++){
                for(j=2;j<=i;j++){
                    omega(i,j-1) = omega(i,j);
                    if(i==j) omega(i,j)=0;
                    }
                }
     
            }//end section for zeroinflated poisson distribution
            
            
            result(n+k,1) = curr;
            u = Rcpp::runif(1,0,1)(0);
            if(u<=zeroprop(curr-1)) result(n+k,0)=0;
            else result(n+k,0) = Rcpp::rpois(1, theta(curr-1))(0);
        }
        
        n += count;
    }

     
    return(result);
     
}


//////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::vec hmm_nllk(arma::vec parm, int M, arma::vec y, arma::vec zeroindex){
    
    arma::vec loglik;
    arma::vec negloglik;
    long dim = y.n_rows;
    int i,j,m,n;
    double tempsum;
    int nextindex;
    
    arma::vec pi(M);
    arma::mat a(M,M);
    arma::vec theta(M);
    arma::vec zeroprop(M);
    
    
    arma::vec forward(M);
    arma::vec meanvec(M);
    arma::vec forwardsum;  //forwardsum is a vec, forwardsum(0) is double
    arma::mat tempmat(1,M);
    
    //retrieve the full parameters
    
    
    pi(0) = 1;
    tempsum = 1;
    for(m=1; m<M; m++) {
        pi(m) = exp(parm(m-1));
        tempsum += pi(m);
    }
    for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
    
    for(i=0; i<M; i++){
        
        for(j=0; j<M; j++){
            if(j==0) {
                a(i,j) = 1;
                tempsum = 1;
            }
            else {
                a(i,j) = exp(parm((i+1)*(M-1)+j-1));
                tempsum += a(i,j);
            }
        }
        
        for(n=0;n<M; n++) a(i,n) = a(i,n) / tempsum;
        
    }
    
    nextindex = M*M-1;
    for(i=0;i<M;i++){
        if(zeroindex(i)==0) zeroprop(i) = 0;
        else{
           zeroprop(i) = exp(parm(nextindex))/(1+exp(parm(nextindex)));
            nextindex += 1;
        }
        
     //Rcpp::Rcout << "p1=" << p1 << std::endl;
    }
    
     for(m=0; m<M; m++) theta(m) = exp(parm(nextindex+m));
    
    
    
    //initialize the forward variable

     for(m=0; m<M; m++) meanvec(m) = dzip(zeroprop(m), theta(m), y(0), FALSE);
    
    forward = pi % meanvec; //element wise multiplication
    forwardsum = pi.t() * meanvec;
    loglik = log(forwardsum);
    //Rcpp::Rcout << "loglik=" << loglik << std::endl;
    
    for(m=0; m<M; m++) forward(m) = forward(m) / forwardsum(0);
    
    //recursion
    for(n=1; n<dim; n++){
  
        for(m=0; m<M; m++) meanvec(m) = dzip(zeroprop(m), theta(m), y(n), FALSE);
    
        tempmat = forward.t() * a; //row matrix
        forward = tempmat.t() % meanvec;
        forwardsum = tempmat * meanvec;
        loglik += log(forwardsum);
        for(m=0; m<M; m++) forward(m) = forward(m) / forwardsum(0);
    }
    
    negloglik = -loglik;
    return(negloglik);
    
}

// [[Rcpp::export]]
double hmm_common_nocov_negloglik(arma::vec allparm, int M, arma::vec ally, arma::vec ntimes,
                                  arma::vec zeroindex){
    
    //wrapper function of the zip_negloglik
    
    
    double negloglik;
    long i,cumtimes;
    long timeindex = ntimes.n_rows;
    
    negloglik = 0;
    cumtimes = 0;
    
    
    for(i=0; i<timeindex; i++){
        
        negloglik += hmm_nllk(allparm, M, ally.subvec(cumtimes, cumtimes + ntimes[i] - 1),
                               zeroindex)(0);
        //use vec.subvec(a,b) to extract elements from a to b in col vec
        cumtimes += ntimes[i];
        //Rcpp::Rcout<<i<<std::endl;
    }
    //Rcpp::Rcout<<negloglik<<std::endl;
    return negloglik;
    
}

////////////////////

// [[Rcpp::export]]
arma::vec hmm_cov_negloglik(arma::vec parm, int M, arma::vec y, int ncolcovpi, arma::mat covpi,
                            int ncolcovtrans, arma::mat covtrans, int ncolcovp1, arma::mat covp1,
                            int ncolcovpois, arma::mat covpois, arma::vec zeroindex){
    
    //intercept column is not included in the covariates, but should be in the parm
    
    //parameters are given alternatively in the order of beta0 and beta1
    arma::vec loglik;
    arma::vec negloglik;
    long dim = y.n_rows;
    int i,j,k,m,n,nextindex;
    double tempsum;
    
    
    arma::vec pi(M);
    arma::mat a(M,M);
    arma::vec theta(M);
    arma::vec zeroprop(M);
    
    
    arma::vec forward(M);
    arma::vec meanvec(M);
    arma::vec forwardsum;  //forwardsum is a vec, forwardsum(0) is double
    arma::mat tempmat(1,M);
    
    //retrieve the full parameters
    //prior
    nextindex=0;
    pi(0) = 1;
    tempsum = 1;
    for(m=1; m<M; m++){
        pi(m) = parm(nextindex);
        
        for(j=0;j<ncolcovpi;j++)
            pi(m) += parm(nextindex+j+1)*covpi(0,j) ;
        
        pi(m) = exp(pi(m));
        
        tempsum += pi(m);
        nextindex = nextindex + ncolcovpi + 1;
    }
    
    for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
    
    
    //transition
    for(i=0; i<M; i++){
        
        for(j=0; j<M; j++){
            if(j==0) {
                a(i,j) = 1;
                tempsum = 1;
            }
            else {
                a(i,j) = parm(nextindex);
                
                for(m=0;m<ncolcovtrans;m++)
                    a(i,j) += parm(nextindex+m+1)*covtrans(0,m);
                
                a(i,j) = exp(a(i,j));
                tempsum += a(i,j);
                nextindex = nextindex + ncolcovtrans + 1;
            }
            //Rcpp::Rcout<<a(i,j)<<std::endl;
        }
        
        for(n=0;n<M; n++) a(i,n) = a(i,n) / tempsum;
        
    }
    
    
    for(i=0;i<M;i++){
        if(zeroindex(i)==0) zeroprop(i)=0;
        else{
            zeroprop(i) = parm(nextindex);
            
            for(m=0;m<ncolcovp1;m++)
                zeroprop(i) += parm(nextindex+m+1) * covp1(0,m);
            
            zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
            nextindex = nextindex + ncolcovp1 + 1;
        }
    }
    
     //Rcpp::Rcout << "p1=" << p1 << std::endl;
    
     for(m=0; m<M; m++){
         theta(m) = parm(nextindex);
         for(j=0; j<ncolcovpois;j++){
             theta(m) += parm(nextindex+j+1) * covpois(0,j);
         }
         theta(m) = exp(theta(m));
         nextindex = nextindex + ncolcovpois + 1;
     }
     //end of zero-inflated poisson
    
    //Rcpp::Rcout << nextindex << std::endl;
    
    
    //initialize the forward variable

     for(m=0; m<M; m++) meanvec(m) = dzip(zeroprop(m), theta(m), y(0), FALSE);
    
    
    forward = pi % meanvec; //element wise multiplication
    forwardsum = pi.t() * meanvec;
    loglik = log(forwardsum);
    //Rcpp::Rcout << "loglik=" << loglik << std::endl;
    
    for(m=0; m<M; m++) forward(m) = forward(m) / forwardsum(0);
    
    //recursion
    for(k=1; k<dim; k++){
        
        if(ncolcovpi + ncolcovtrans>0){
            
        nextindex=0;
        pi(0) = 1;
        tempsum = 1;
        for(m=1; m<M; m++){
            pi(m) = parm(nextindex);
            
            for(j=0;j<ncolcovpi;j++)
                pi(m) += parm(nextindex+j+1)*covpi(k,j) ;
            
            pi(m) = exp(pi(m));
            
            tempsum += pi(m);
            nextindex = nextindex + ncolcovpi + 1;
        }
        
        for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
        
        
        //transition
        for(i=0; i<M; i++){
            
            for(j=0; j<M; j++){
                if(j==0) {
                    a(i,j) = 1;
                    tempsum = 1;
                }
                else {
                    a(i,j) = parm(nextindex);
                    
                    for(m=0;m<ncolcovtrans;m++)
                        a(i,j) += parm(nextindex+m+1)*covtrans(k,m);
                    
                    a(i,j) = exp(a(i,j));
                    tempsum += a(i,j);
                    nextindex = nextindex + ncolcovtrans + 1;
                }
                //Rcpp::Rcout<<a(i,j)<<std::endl;
            }
            
            for(n=0;n<M; n++) a(i,n) = a(i,n) / tempsum;
            
        }
        
        
        for(i=0;i<M;i++){
            if(zeroindex(i)==0) zeroprop(i)=0;
            else{
                zeroprop(i) = parm(nextindex);
                
                for(m=0;m<ncolcovp1;m++)
                    zeroprop(i) += parm(nextindex+m+1) * covp1(k,m);
                
                zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
                nextindex = nextindex + ncolcovp1 + 1;
            }
        }
        //Rcpp::Rcout << "p1=" << p1 << std::endl;
        
         for(m=0; m<M; m++){
             theta(m) = parm(nextindex);
             for(j=0; j<ncolcovpois;j++){
                 theta(m) += parm(nextindex+j+1) * covpois(k,j);
             }
             theta(m) = exp(theta(m));
             nextindex = nextindex + ncolcovpois + 1;
         }
        
        }
        else{
            nextindex = M-1 + M*(M-1);
            for(i=0;i<M;i++){
                if(zeroindex(i)==0) zeroprop(i)=0;
                else{
                    zeroprop(i) = parm(nextindex);
                    
                    for(m=0;m<ncolcovp1;m++)
                        zeroprop(i) += parm(nextindex+m+1) * covp1(k,m);
                    
                    zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
                    nextindex = nextindex + ncolcovp1 + 1;
                }
            }
            //Rcpp::Rcout << "p1=" << p1 << std::endl;
            
            for(m=0; m<M; m++){
                theta(m) = parm(nextindex);
                for(j=0; j<ncolcovpois;j++){
                    theta(m) += parm(nextindex+j+1) * covpois(k,j);
                }
                theta(m) = exp(theta(m));
                nextindex = nextindex + ncolcovpois + 1;
            }

        }
        
         for(m=0; m<M; m++) meanvec(m) = dzip(zeroprop(m), theta(m), y(k), FALSE);
        
        
        tempmat = forward.t() * a; //row matrix
        forward = tempmat.t() % meanvec;
        forwardsum = tempmat * meanvec;
        loglik += log(forwardsum);
        for(m=0; m<M; m++) forward(m) = forward(m) / forwardsum(0);
    }
    
    negloglik = -loglik;
    //Rcpp::Rcout << "loglik=" << loglik << std::endl;
    return(negloglik);
}




//////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double hmm_common_negloglik(arma::vec allparm, int M, arma::vec ally, arma::vec ntimes,
                            int ncolcovpi, arma::mat allcovpi,
                            int ncolcovtrans, arma::mat allcovtrans, int ncolcovp1, arma::mat allcovp1,
                            int ncolcovpois, arma::mat allcovpois, arma::vec zeroindex){
    
    //wrapper function of the zip_negloglik
    
    
    double negloglik;
    long i,cumtimes;
    long timeindex = ntimes.n_rows;
    
    negloglik = 0;
    cumtimes = 0;
    
    
    for(i=0; i<timeindex; i++){
        
        negloglik += hmm_cov_negloglik(allparm, M, ally.subvec(cumtimes, cumtimes + ntimes[i] - 1),
                                       ncolcovpi,allcovpi.rows(cumtimes, cumtimes + ntimes[i] - 1),
                                       ncolcovtrans,allcovtrans.rows(cumtimes, cumtimes + ntimes[i] - 1),
                                       ncolcovp1,allcovp1.rows(cumtimes, cumtimes + ntimes[i] - 1),
                                       ncolcovpois,allcovpois.rows(cumtimes,cumtimes+ntimes[i]-1),
                                       zeroindex)(0);
        //use vec.subvec(a,b) to extract elements from a to b in col vec
        cumtimes += ntimes[i];
        //Rcpp::Rcout<<i<<std::endl;
    }
    //Rcpp::Rcout<<negloglik<<std::endl;
    return negloglik;
    
}


/*
double hmm_mixed_negloglik(arma::vec allparm, int M, arma::vec ally, arma::vec ntimes,
                           int ncolcovpi, arma::mat allcovpi,
                           int ncolcovtrans, arma::mat allcovtrans, int ncolcovp1, arma::mat allcovp1,
                           int ncolcovpois, arma::mat allcovpois, int nrand, double B,
                           arma::vec zeroindex){
    
    //wrapper function of the zip_negloglik
    //last items in allparm are sd's for random effects with mean 0
    
    double b,negloglik,tempsum,factor,first;
    long i,j,cumtimes;
    long timeindex = ntimes.n_rows;
    arma::vec sd(nrand);
    //the dimension of parm has to be modified (AT LEAST M*M+M when no covariates)
    int parmdim = M * M + M + ncolcovpi * (M-1) + ncolcovtrans * M * (M-1) + ncolcovp1 + ncolcovpois * M;
    arma::vec parm(parmdim);
    
    negloglik = 0;
    cumtimes = 0;
    
    //Rcpp::Rcout<<parmdim<<std::endl;
    
    //random effects: currently testing for 1 for the last theta
    for(j=0;j<nrand;j++) sd(j) = allparm(parmdim+j);
    
    for(i=0; i<timeindex; i++){
        
        //importance sampling for each subject's log likelihood
        //L(beta)=f(y) = 1/B * sum(f(y|ui,beta)*f(ui,beta)/g(ui))
        //let g(ui)=f(ui,beta) then LogLi = mean(f(y|ui,beta)); then sum up for all subjects
        //use the trick of sum of logs to avoid underflow
        //logLi = -logB + log(f(y|u_1,beta)) + log[1+...+exp[log(f(y|u_B,beta))-log(f(y|u_1,beta))]]
        //NLLK = logB + negloglik1 - log(1+...+exp(negloglik1-negloglikB))
        
        //fixed effects
        for(j=0;j<parmdim;j++) parm(j) = allparm(j);
        parm(parmdim - ncolcovpois - 1) = parm(parmdim - ncolcovpois - 1) + Rcpp::rnorm(1,0,sd(0))(0);
        first = hmm_cov_negloglik(parm, M, ally.subvec(cumtimes, cumtimes + ntimes[i] - 1),
                                  ncolcovpi,allcovpi.rows(cumtimes, cumtimes + ntimes[i] - 1),
                                  ncolcovtrans,allcovtrans.rows(cumtimes, cumtimes + ntimes[i] - 1),
                                  ncolcovp1,allcovp1.rows(cumtimes, cumtimes + ntimes[i] - 1),
                                  ncolcovpois,allcovpois.rows(cumtimes,cumtimes+ntimes[i]-1),
                                  zeroindex)(0);
        tempsum = 1;
        
        for(b=0; b<B; b++){
            for(j=0;j<parmdim;j++) parm(j) = allparm(j);
            parm(parmdim - ncolcovpois - 1) = parm(parmdim - ncolcovpois - 1) + Rcpp::rnorm(1,0,sd(0))(0);
            factor = hmm_cov_negloglik(parm, M, ally.subvec(cumtimes, cumtimes + ntimes[i] - 1),
                                       ncolcovpi,allcovpi.rows(cumtimes, cumtimes + ntimes[i] - 1),
                                       ncolcovtrans,allcovtrans.rows(cumtimes, cumtimes + ntimes[i] - 1),
                                       ncolcovp1,allcovp1.rows(cumtimes, cumtimes + ntimes[i] - 1),
                                       ncolcovpois,allcovpois.rows(cumtimes,cumtimes+ntimes[i]-1),
                                       zeroindex)(0);
            tempsum = tempsum + exp(first-factor);
            
        }
        
        
        negloglik += log(B) + first - log(tempsum);
        
        cumtimes += ntimes[i];
        //Rcpp::Rcout<<i<<std::endl;
    }
    
    return negloglik;
    
}
*/


//////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::vec hsmm_nllk(arma::vec parm, int M, arma::vec trunc, arma::vec y,
                    std::string dt_dist, arma::vec zeroindex){
    
    /*parm = logitp from logarithm dist---M, logit(pi)---M-1, logit(p1), log(theta)---M,
     omega---M*(M-2): (1,3)...(1,n)(2,3)...(2,n)(3,2)(3,4)...(3,n).......*/
    
    arma::vec negloglik;
    arma::vec loglik;
    double tempsum;
    int nextindex;
    arma::vec zeroprop(M);
    long dim = y.n_rows;
    int i,j,m,n;
    int dmcol = arma::max(trunc); //get the maximum
    int totalmv = arma::sum(trunc);
    arma::vec pi(M);
    arma::vec logtheta(M);
    arma::vec newtheta(totalmv);
    arma::vec newpi(totalmv);
    arma::vec newzeroprop(totalmv);
    //Rcpp::Rcout<<dmcol<<std::endl;
    arma::mat dm(M,dmcol);
    
    arma::mat omega(M,M);
    arma::mat a(totalmv,totalmv);
    
    arma::vec forward(totalmv);
    arma::vec meanvec(totalmv);
    arma::vec forwardsum;  //forwardsum is a vec, forwardsum(0) is double
    arma::mat tempmat(1,totalmv);
    
    //recover dm: dwell time
    dm.zeros();
    //check dt_dist
    if(dt_dist=="log"){
     for(i=0;i<M;i++){
         for(j=0;j<trunc(i);j++){
             dm(i,j)=dlogp(j+1, exp(parm(i))/(1+exp(parm(i))), FALSE);
         }
      }
    }else{
        for(i=0;i<M;i++){
            for(j=0;j<trunc(i);j++){
                dm(i,j)=dshiftpois(j+1, exp(parm(i)), 1, FALSE);
            }
        }
    }
    
    //recover newtheta,p1, newpi
    
    pi(0) = 1;
    tempsum = 1;
    for(m=1; m<M; m++) {
        pi(m) = exp(parm(M+m-1));
        tempsum += pi(m);
    }
    for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
    
    nextindex = M + M-1;
    for(i=0;i<M;i++){
        if(zeroindex(i)==0) zeroprop(i) = 0;
        else{
            zeroprop(i) = exp(parm(nextindex))/(1+exp(parm(nextindex)));
            nextindex += 1;
        }
        
        //Rcpp::Rcout << "p1=" << p1 << std::endl;
    }
    
    
     for(m=0; m<M; m++) logtheta(m) = parm(nextindex+m);
    
    nextindex = nextindex + M;
    
     //get newlogtheta,newpi
     tempsum = 0;
     for(i=0;i<M;i++){
         for(j=0;j<trunc[i];j++){
             newtheta(tempsum+j) = exp(logtheta(i));
             newpi(tempsum+j) = pi(i)/trunc(i);
             newzeroprop(tempsum+j) = zeroprop(i);
            //Rcpp::Rcout << "theta=" << newtheta(tempsum+j) << "pi=" << newpi(tempsum+j) << std::endl;
         }
         tempsum += trunc(i);
     }
    
    //recover omega:   from M*(M-2) [3M ~ M*M+M]
    //if M=2 then assign off-diagonal with 1; else
    //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
    
     omega.zeros();
     if(M==2) {omega(0,1)=1; omega(1,0)=1;}
     else{
         for(i=0;i<M;i++){
             tempsum = 1;
             for(j=2;j<M;j++){
                 omega(i,j) = exp(parm(nextindex+(M-2)*i+j-2)); //updated
                 tempsum += omega(i,j); //updated
                 //omega(i,j) = exp(parm(3*M+(M-2)*i+j-2))/(1+exp(parm(3*M+(M-2)*i+j-2))); //old
             }
             
             for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //updated
             
             if(i==0) omega(i,1) = 1 / tempsum;
             else omega(i,0) = 1 / tempsum;
             //tempsum = 1 - arma::sum(omega.row(i)); //old
             //if(i==0) omega(i,1) = MAX(0,tempsum); //old
             //else omega(i,0) = MAX(0,tempsum); //old
         }
        
         for(i=2;i<M;i++){
             for(j=2;j<=i;j++){
                 omega(i,j-1) = omega(i,j);
                 if(i==j) omega(i,j)=0;
             }
         }
        
      }
    
    
    // recover transition matrix
    
    a = hsmm_hmm (omega, dm, trunc);
    
    //forward algorithm
    
    //initialize the forward variable

     for(m=0; m<totalmv; m++) meanvec(m) = dzip(newzeroprop(m), newtheta(m), y(0), FALSE);
    
    forward = newpi % meanvec; //element wise multiplication
    forwardsum = newpi.t() * meanvec;
    loglik = log(forwardsum);
    //Rcpp::Rcout << "loglik=" << loglik << std::endl;
    
    for(m=0; m<totalmv; m++) forward(m) = forward(m) / forwardsum(0);
    
    //recursion
    for(n=1; n<dim; n++){
    
         for(m=0; m<totalmv; m++) meanvec(m) = dzip(newzeroprop(m), newtheta(m), y(n), FALSE);
        
        
        tempmat = forward.t() * a; //row matrix
        forward = tempmat.t() % meanvec;
        forwardsum = tempmat * meanvec;
        loglik += log(forwardsum);
        for(m=0; m<totalmv; m++) forward(m) = forward(m) / forwardsum(0);
    }
    
    negloglik = -loglik;
    return(negloglik);
    
    
}


//////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double hsmm_common_nocov_nllk(arma::vec allparm, int M, arma::vec ally, arma::vec trunc,arma::vec ntimes,
                              std::string dt_dist, arma::vec zeroprop){
    
    //wrapper function of the hsmm_negloglik
    
    double negloglik;
    long i,cumtimes;
    long timeindex = ntimes.n_rows;
    
    negloglik = 0;
    cumtimes = 0;
    
    
    for(i=0; i<timeindex; i++){
        
        negloglik += hsmm_nllk(allparm, M, trunc, ally.subvec(cumtimes, cumtimes + ntimes[i] - 1), 
                                dt_dist, zeroprop)(0);
        
        //use vec.subvec(a,b) to extract elements from a to b in col vec
        cumtimes += ntimes[i];
        //Rcpp::Rcout<<i<<std::endl;
    }
    //Rcpp::Rcout<<negloglik<<std::endl;
    return negloglik;
    
}


////////////////////
// [[Rcpp::export]]
arma::mat hsmm_cov_nllk(arma::vec parm, int M, arma::vec y,arma::vec trunc,
                        std::string dt_dist, arma::vec zeroindex,
                        int ncolcovp, arma::mat covp, int ncolcovpi, arma::mat covpi,
                        int ncolcovomega, arma::mat covomega, int ncolcovp1, arma::mat covp1,
                        int ncolcovpois, arma::mat covpois){
    /*parm = logitp from logarithm dist---M, logit(pi)---M-1, logit(p1), log(theta)---M,
     omega---M*(M-2): (1,3)...(1,n)(2,3)...(2,n)(3,2)(3,4)...(3,n).......*/
    
    //intercept column is not included in the covariates, but should be in the parm
    
    //parameters are given alternatively in the order of beta0 and beta1
    arma::vec negloglik;
    arma::vec loglik;
    double tempsum;
    arma::vec zeroprop(M);
    long dim = y.n_rows;
    int i,j,k,m,nextindex;
    int dmcol = arma::max(trunc); //get the maximum
    int totalmv = arma::sum(trunc);
    arma::vec pi(M);
    arma::vec theta(M);
    arma::vec newtheta(totalmv);
    arma::vec newpi(totalmv);
    arma::vec newzeroprop(totalmv);
    //Rcpp::Rcout<<dmcol<<std::endl;
    arma::mat dm(M,dmcol);
    
    arma::mat omega(M,M);
    arma::mat a(totalmv,totalmv);
    
    arma::vec forward(totalmv);
    arma::vec meanvec(totalmv);
    arma::vec forwardsum;  //forwardsum is a vec, forwardsum(0) is double
    arma::mat tempmat(1,totalmv);
    
    //retrieve the full parameters
    nextindex = 0;
    //dwell time
    dm.zeros();
    //check dt_dist
    if(dt_dist=="log"){
        for(i=0;i<M;i++){
            
            for(j=0;j<trunc(i);j++){
                tempsum = parm(nextindex);
                for(m=0;m<ncolcovp;m++)
                    tempsum += parm(nextindex+m+1)*covp(0,m);
                tempsum = exp(tempsum) / (1+exp(tempsum));
                dm(i,j)=dlogp(j+1, tempsum, FALSE);
                
            }
            nextindex = nextindex + ncolcovp + 1;
        }
    }else{
        for(i=0;i<M;i++){
            
            for(j=0;j<trunc(i);j++){
                tempsum = parm(nextindex);
                for(m=0;m<ncolcovp;m++)
                    tempsum += parm(nextindex+m+1)*covp(0,m);
                tempsum = exp(tempsum) ;
                dm(i,j)=dshiftpois(j+1, tempsum, 1, FALSE);
                
            }
            nextindex = nextindex + ncolcovp + 1;
        }
        
    }
    
    
    //recover newtheta,p1, newpi
    pi(0) = 1;
    tempsum = 1;
    for(m=1; m<M; m++){
        pi(m) = parm(nextindex);
        
        for(j=0;j<ncolcovpi;j++)
            pi(m) += parm(nextindex+j+1)*covpi(0,j) ;
        
        pi(m) = exp(pi(m));
        
        tempsum += pi(m);
        nextindex = nextindex + ncolcovpi + 1;
    }
    
    for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
    
    //zero proportions
    for(i=0;i<M;i++){
        if(zeroindex(i)==0) zeroprop(i)=0;
        else{
            zeroprop(i) = parm(nextindex);
            
            for(m=0;m<ncolcovp1;m++)
                zeroprop(i) += parm(nextindex+m+1) * covp1(0,m);
            
            zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
            nextindex = nextindex + ncolcovp1 + 1;
        }
    }
    
    
    
    //theta
    for(m=0; m<M; m++){
        theta(m) = parm(nextindex);
        for(j=0; j<ncolcovpois;j++){
            theta(m) += parm(nextindex+j+1) * covpois(0,j);
        }
        theta(m) = exp(theta(m));
        nextindex = nextindex + ncolcovpois + 1;
    }
    
    //Rcpp::Rcout << nextindex << std::endl;
    //get newlogtheta,newpi
    tempsum = 0;
    for(i=0;i<M;i++){
        for(j=0;j<trunc[i];j++){
            newtheta(tempsum+j) = theta(i);
            newpi(tempsum+j) = pi(i)/trunc(i);
            newzeroprop(tempsum+j) = zeroprop(i);
            //Rcpp::Rcout << "theta=" << newtheta(tempsum+j) << "pi=" << newpi(tempsum+j) << std::endl;
        }
        tempsum += trunc(i);
    }
    
    //recover omega:
    //if M=2 then assign off-diagonal with 1; else
    //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
    
    omega.zeros();
    if(M==2) {omega(0,1)=1; omega(1,0)=1;}
    else{
        for(i=0;i<M;i++){
            tempsum = 1;
            for(j=2;j<M;j++){
                omega(i,j) = parm(nextindex); //new
                for(m=0;m<ncolcovomega;m++)
                    omega(i,j) += parm(nextindex+m+1)*covomega(0,m);
                //omega(i,j) = exp(omega(i,j))/(1+exp(omega(i,j)));
                omega(i,j) = exp(omega(i,j)); //new
                tempsum += omega(i,j); //new
                nextindex = nextindex + ncolcovomega + 1;
            }
            
            for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //new
            if(i==0) omega(i,1) = 1 / tempsum;
            else omega(i,0) = 1 / tempsum;
            //tempsum = 1 - arma::sum(omega.row(i)); //prevent negative values
            //if(i==0) omega(i,1) = MAX(0,tempsum);
            //else omega(i,0) = MAX(0,tempsum);
        }
        
        for(i=2;i<M;i++){
            for(j=2;j<=i;j++){
                omega(i,j-1) = omega(i,j);
                if(i==j) omega(i,j)=0;
            }
        }
        
        
    }  //end of poisson and zerofl part
    
    
    
    a = hsmm_hmm (omega, dm, trunc);
    
    
    //initialize the forward variable
    for(m=0; m<totalmv; m++) meanvec(m) = dzip(newzeroprop(m), newtheta(m), y(0), FALSE);
    
    forward = newpi % meanvec; //element wise multiplication
    forwardsum = newpi.t() * meanvec;
    loglik = log(forwardsum);
    //Rcpp::Rcout << "loglik=" << loglik << std::endl;
    
    for(m=0; m<totalmv; m++) forward(m) = forward(m) / forwardsum(0);
    
    //recursion
    for(k=1; k<dim; k++){
        
        if(ncolcovp+ncolcovpi+ncolcovomega>0){
            
            nextindex = 0;
            //dwell time
            dm.zeros();
            //check dt_dist
            if(dt_dist=="log"){
                for(i=0;i<M;i++){
                    
                    for(j=0;j<trunc(i);j++){
                        tempsum = parm(nextindex);
                        for(m=0;m<ncolcovp;m++)
                            tempsum += parm(nextindex+m+1)*covp(k,m);
                        tempsum = exp(tempsum) / (1+exp(tempsum));
                        dm(i,j)=dlogp(j+1, tempsum, FALSE);
                        
                    }
                    nextindex = nextindex + ncolcovp + 1;
                }
            }else{
                for(i=0;i<M;i++){
                    
                    for(j=0;j<trunc(i);j++){
                        tempsum = parm(nextindex);
                        for(m=0;m<ncolcovp;m++)
                            tempsum += parm(nextindex+m+1)*covp(k,m);
                        tempsum = exp(tempsum) ;
                        dm(i,j)=dshiftpois(j+1, tempsum, 1, FALSE);
                        
                    }
                    nextindex = nextindex + ncolcovp + 1;
                }
                
            }
            
            //recover newtheta,p1, newpi
            pi(0) = 1;
            tempsum = 1;
            for(m=1; m<M; m++){
                pi(m) = parm(nextindex);
                
                for(j=0;j<ncolcovpi;j++)
                    pi(m) += parm(nextindex+j+1)*covpi(k,j) ;
                
                pi(m) = exp(pi(m));
                
                tempsum += pi(m);
                nextindex = nextindex + ncolcovpi + 1;
            }
            
            for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
            
            //zeroprops
            
            for(i=0;i<M;i++){
                if(zeroindex(i)==0) zeroprop(i)=0;
                else{
                    zeroprop(i) = parm(nextindex);
                    
                    for(m=0;m<ncolcovp1;m++)
                        zeroprop(i) += parm(nextindex+m+1) * covp1(k,m);
                    
                    zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
                    nextindex = nextindex + ncolcovp1 + 1;
                }
            }
            
            
            
            //theta
            for(m=0; m<M; m++){
                theta(m) = parm(nextindex);
                for(j=0; j<ncolcovpois;j++){
                    theta(m) += parm(nextindex+j+1) * covpois(k,j);
                }
                theta(m) = exp(theta(m));
                nextindex = nextindex + ncolcovpois + 1;
            }
            
            //Rcpp::Rcout << nextindex << std::endl;
            //get newlogtheta,newpi
            tempsum = 0;
            for(i=0;i<M;i++){
                for(j=0;j<trunc[i];j++){
                    newtheta(tempsum+j) = theta(i);
                    newpi(tempsum+j) = pi(i)/trunc(i);
                    newzeroprop(tempsum+j) = zeroprop(i);
                    //Rcpp::Rcout << "theta=" << newtheta(tempsum+j) << "pi=" << newpi(tempsum+j) << std::endl;
                }
                tempsum += trunc(i);
            }
            
            
            //recover omega:
            //if M=2 then assign off-diagonal with 1; else
            //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
            
            omega.zeros();
            if(M==2) {omega(0,1)=1; omega(1,0)=1;}
            else{
                for(i=0;i<M;i++){
                    tempsum = 1;
                    for(j=2;j<M;j++){
                        omega(i,j) = parm(nextindex); //new
                        for(m=0;m<ncolcovomega;m++)
                            omega(i,j) += parm(nextindex+m+1)*covomega(k,m);
                        //omega(i,j) = exp(omega(i,j))/(1+exp(omega(i,j)));
                        omega(i,j) = exp(omega(i,j)); //new
                        tempsum += omega(i,j); //new
                        nextindex = nextindex + ncolcovomega + 1;
                    }
                    
                    for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //new
                    if(i==0) omega(i,1) = 1 / tempsum;
                    else omega(i,0) = 1 / tempsum;
                    //tempsum = 1 - arma::sum(omega.row(i)); //prevent negative values
                    //if(i==0) omega(i,1) = MAX(0,tempsum);
                    //else omega(i,0) = MAX(0,tempsum);
                }
                
                for(i=2;i<M;i++){
                    for(j=2;j<=i;j++){
                        omega(i,j-1) = omega(i,j);
                        if(i==j) omega(i,j)=0;
                    }
                }
                
            }
            
            // recover transition matrix
            
            a = hsmm_hmm (omega, dm, trunc);
            
        }else{
            //only update emit, zeroprop covariates
            
            nextindex = M + M - 1;
            
            for(i=0;i<M;i++){
                if(zeroindex(i)==0) zeroprop(i)=0;
                else{
                    zeroprop(i) = parm(nextindex);
                    
                    for(m=0;m<ncolcovp1;m++)
                        zeroprop(i) += parm(nextindex+m+1) * covp1(k,m);
                    
                    zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
                    nextindex = nextindex + ncolcovp1 + 1;
                }
            }
            
            
            
            //theta
            for(m=0; m<M; m++){
                theta(m) = parm(nextindex);
                for(j=0; j<ncolcovpois;j++){
                    theta(m) += parm(nextindex+j+1) * covpois(k,j);
                }
                theta(m) = exp(theta(m));
                nextindex = nextindex + ncolcovpois + 1;
            }
            
            //Rcpp::Rcout << nextindex << std::endl;
            //get newlogtheta,newpi
            tempsum = 0;
            for(i=0;i<M;i++){
                for(j=0;j<trunc[i];j++){
                    newtheta(tempsum+j) = theta(i);
                    //newpi(tempsum+j) = pi(i)/trunc(i);
                    newzeroprop(tempsum+j) = zeroprop(i);
                    //Rcpp::Rcout << "theta=" << newtheta(tempsum+j) << "pi=" << newpi(tempsum+j) << std::endl;
                }
                tempsum += trunc(i);
            }
            
        }
        
        
        for(m=0; m<totalmv; m++) meanvec(m) = dzip(newzeroprop(m), newtheta(m), y(k), FALSE);
        
        tempmat = forward.t() * a; //row matrix
        forward = tempmat.t() % meanvec;
        forwardsum = tempmat * meanvec;
        loglik += log(forwardsum);
        for(m=0; m<totalmv; m++) forward(m) = forward(m) / forwardsum(0);
    }
    
    negloglik = -loglik;
    
    
    //Rcpp::Rcout << "loglik=" << loglik << std::endl;
    return(negloglik);
}


//////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double hsmm_common_negloglik(arma::vec allparm, int M, arma::vec ally, arma::vec trunc,arma::vec ntimes,
                             std::string dt_dist, arma::vec zeroindex,
                             int ncolcovp, arma::mat allcovp, int ncolcovpi, arma::mat allcovpi,
                             int ncolcovomega, arma::mat allcovomega, int ncolcovp1, arma::mat allcovp1,
                             int ncolcovpois, arma::mat allcovpois){
    
    //wrapper function of the zip_negloglik
    
    
    double negloglik;
    long i,cumtimes;
    long timeindex = ntimes.n_rows;
    
    negloglik = 0;
    cumtimes = 0;
    
    
    for(i=0; i<timeindex; i++){
        
        negloglik += hsmm_cov_nllk(allparm, M, ally.subvec(cumtimes, cumtimes + ntimes[i] - 1), trunc,
                                    dt_dist, zeroindex,
                                    ncolcovp, allcovp.rows(cumtimes, cumtimes + ntimes[i] - 1),
                                    ncolcovpi,allcovpi.rows(cumtimes, cumtimes + ntimes[i] - 1),
                                    ncolcovomega,allcovomega.rows(cumtimes, cumtimes + ntimes[i] - 1),
                                    ncolcovp1,allcovp1.rows(cumtimes, cumtimes + ntimes[i] - 1),
                                    ncolcovpois,allcovpois.rows(cumtimes, cumtimes + ntimes[i] - 1))(0);
        
        
        //use vec.subvec(a,b) to extract elements from a to b in col vec
        cumtimes += ntimes[i];
        //Rcpp::Rcout<<i<<std::endl;
    }
    //Rcpp::Rcout<<negloglik<<std::endl;
    return negloglik;
    
}


//////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::vec hmm_viterbi(arma::vec pi, arma::mat a, arma::vec theta, int M, arma::vec y,
                      arma::vec zeroprop){
    
    long dim = y.n_rows;
    int i,j,m,n;
    double colmax;
    
    
    arma::vec meanvec(M);
    arma::vec forward(M);
    arma::vec forwardsum;
    arma::vec state(dim);
    arma::vec tempmax(M);
    arma::mat xi(dim, M);
    arma::mat tempmat(M,M);
    
    
    //initialize the forward variable
    
     for(m=0; m<M; m++) meanvec(m) = dzip(zeroprop(m), theta(m), y(0), FALSE);
    
        
    forward = pi % meanvec; //element wise multiplication
    forwardsum = pi.t() * meanvec;
    for(m=0; m<M; m++) xi(0,m) = forward(m) / forwardsum(0);
    
    //recursion
    for(n=1; n<dim; n++){

        for(m=0; m<M; m++) meanvec(m) = dzip(zeroprop(m), theta(m), y(n), FALSE);
        
        //difficult part
        for(i=0; i<M; i++){
            for(j=0;j<M;j++){
                tempmat(i,j) = xi(n-1, i) * a(i,j);
            }
        }
        
        for(j=0; j<M; j++){
            for(i=0;i<M;i++){
                if(i==0) colmax = tempmat(i,j);
                else colmax = MAX(colmax, tempmat(i,j));
            }
            tempmax(j) = colmax;
        }
        
        
        forward = tempmax % meanvec;
        forwardsum = tempmax.t() * meanvec;
        for(m=0; m<M; m++) xi(n,m) = forward(m) / forwardsum(0);
    }
    
    //termination
    for(m=0; m<M; m++){
        if(m==0) {
            colmax = xi(dim-1, m);
            state(dim-1) = m;
        }
        else{
            if(xi(dim-1,m) >= colmax) {
                state(dim-1) = m;
                colmax = xi(dim-1,m);
            }
        }
    }
    
    //backtracking
    for(n=dim-2; n>=0; n--){
        
        for(m=0;m<M;m++){
            if(m==0) {
                colmax = xi(n,m) * a(m, state(n+1));
                state(n) = m;
            }
            else{
                if(xi(n,m) * a(m, state(n+1))>=colmax) {
                    state(n) = m;
                    colmax = xi(n,m) * a(m, state(n+1));
                }
            }
        }
    }
    
    for(n=0; n<dim; n++) state(n) = state(n) + 1;
    return(state);
}



//////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::vec hsmm_viterbi(arma::vec y, int M, arma::vec pi, arma::vec theta,
                           arma::mat omega, arma::vec p, arma::vec trunc,
                       std::string dt_dist, arma::vec zeroprop){
    
    long dim = y.n_rows;
    double tempsum;
    arma::vec state(dim);
    int i,j;
    int dmcol = arma::max(trunc); //get the maximum
    int totalmv = arma::sum(trunc);
    arma::vec newtheta(totalmv);
    arma::vec newpi(totalmv);
    arma::vec newzeroprop(totalmv);
    //Rcpp::Rcout<<dmcol<<std::endl;
    arma::mat dm(M,dmcol);
    
    arma::mat a(totalmv,totalmv);
    
    
    
    //recover dm: dwell time
    dm.zeros();
    //check dt_dist
    if(dt_dist=="log"){
        for(i=0;i<M;i++){
            for(j=0;j<trunc(i);j++){
                dm(i,j)=dlogp(j+1, p(i), FALSE);
            }
        }
    }else{
        for(i=0;i<M;i++){
            for(j=0;j<trunc(i);j++){
                dm(i,j)=dshiftpois(j+1, p(i), 1, FALSE);
            }
        }
    }
    // recover transition matrix
    
    a = hsmm_hmm (omega, dm, trunc);
    
    //recover newtheta,newpi
    
    tempsum = 0;
    for(i=0;i<M;i++){
        for(j=0;j<trunc[i];j++){
            
              newtheta(tempsum+j) = theta(i);
            newzeroprop(tempsum+j) = zeroprop(i);
            newpi(tempsum+j) = pi(i)/trunc(i);
        }
        tempsum += trunc(i);
    }
    
    state = hmm_viterbi(newpi, a, newtheta, totalmv, y, newzeroprop);
    
    for(i=0;i<dim;i++){
        tempsum=0;
        for(j=0;j<M;j++){
            if(state(i)<=tempsum+trunc(j)){
                state(i)=j+1;
                break;
            }
            tempsum += trunc(j);
        }
    }
    
    return state;
}







/////////////////////////

// [[Rcpp::export]]
arma::vec hmm_cov_viterbi(arma::vec parm, int M, arma::vec y, int ncolcovpi, arma::mat covpi,
                          int ncolcovtrans, arma::mat covtrans, int ncolcovp1, arma::mat covp1,
                          int ncolcovpois, arma::mat covpois, arma::vec zeroindex){
    
    long dim = y.n_rows;
    int i,j,m,n,nextindex;
    double tempsum, colmax;
    
    
    arma::vec pi(M);
    arma::mat a(M,M);
    arma::vec theta(M);
    arma::vec zeroprop(M);
    
    arma::vec meanvec(M);
    arma::vec forward(M);
    arma::vec forwardsum;
    arma::vec state(dim);
    arma::vec tempmax(M);
    arma::mat xi(dim, M);
    arma::mat tempmat(M,M);
    
    //retrieve the full parameters
    //prior
    nextindex=0;
    pi(0) = 1;
    tempsum = 1;
    for(m=1; m<M; m++){
        pi(m) = parm(nextindex);
        
        for(j=0;j<ncolcovpi;j++)
            pi(m) += parm(nextindex+j+1)*covpi(0,j) ;
        
        pi(m) = exp(pi(m));
        
        tempsum += pi(m);
        nextindex = nextindex + ncolcovpi + 1;
    }
    
    for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
    
    
    //transition
    for(i=0; i<M; i++){
        
        for(j=0; j<M; j++){
            if(j==0) {
                a(i,j) = 1;
                tempsum = 1;
            }
            else {
                a(i,j) = parm(nextindex);
                
                for(m=0;m<ncolcovtrans;m++)
                    a(i,j) += parm(nextindex+m+1)*covtrans(0,m);
                
                a(i,j) = exp(a(i,j));
                tempsum += a(i,j);
                nextindex = nextindex + ncolcovtrans + 1;
            }
            //Rcpp::Rcout<<a(i,j)<<std::endl;
        }
        
        for(n=0;n<M; n++) a(i,n) = a(i,n) / tempsum;
        
    }
    
   
    for(i=0;i<M;i++){
        if(zeroindex(i)==0) zeroprop(i)=0;
        else{
            zeroprop(i) = parm(nextindex);
            
            for(m=0;m<ncolcovp1;m++)
                zeroprop(i) += parm(nextindex+m+1) * covp1(0,m);
            
            zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
            nextindex = nextindex + ncolcovp1 + 1;
        }
    }

    //Rcpp::Rcout << "p1=" << p1 << std::endl;
    
     for(m=0; m<M; m++){
         theta(m) = parm(nextindex);
         for(j=0; j<ncolcovpois;j++){
             theta(m) += parm(nextindex+j+1) * covpois(0,j);
         }
         theta(m) = exp(theta(m));
         nextindex = nextindex + ncolcovpois + 1;
      }
    
    //Rcpp::Rcout << nextindex << std::endl;
    
    //initialize the forward variable
    //check emit_dist

     for(m=0; m<M; m++) meanvec(m) = dzip(zeroprop(m), theta(m), y(0), FALSE);
    
    
    forward = pi % meanvec; //element wise multiplication
    forwardsum = pi.t() * meanvec;
    for(m=0; m<M; m++) xi(0,m) = forward(m) / forwardsum(0);
    
    //recursion
    for(n=1; n<dim; n++){
        
        //retrieve the full parameters
        //prior
        nextindex=0;
        pi(0) = 1;
        tempsum = 1;
        for(m=1; m<M; m++){
            pi(m) = parm(nextindex);
            
            for(j=0;j<ncolcovpi;j++)
                pi(m) += parm(nextindex+j+1)*covpi(n,j) ;
            
            pi(m) = exp(pi(m));
            
            tempsum += pi(m);
            nextindex = nextindex + ncolcovpi + 1;
        }
        
        for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
        
        
        //transition
        for(i=0; i<M; i++){
            
            for(j=0; j<M; j++){
                if(j==0) {
                    a(i,j) = 1;
                    tempsum = 1;
                }
                else {
                    a(i,j) = parm(nextindex);
                    
                    for(m=0;m<ncolcovtrans;m++)
                        a(i,j) += parm(nextindex+m+1)*covtrans(n,m);
                    
                    a(i,j) = exp(a(i,j));
                    tempsum += a(i,j);
                    nextindex = nextindex + ncolcovtrans + 1;
                }
                //Rcpp::Rcout<<a(i,j)<<std::endl;
            }
            
            for(m=0;m<M; m++) a(i,m) = a(i,m) / tempsum;
            
        }
        
        
        for(i=0;i<M;i++){
            if(zeroindex(i)==0) zeroprop(i)=0;
            else{
                zeroprop(i) = parm(nextindex);
                
                for(m=0;m<ncolcovp1;m++)
                    zeroprop(i) += parm(nextindex+m+1) * covp1(n,m);
                
                zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
                nextindex = nextindex + ncolcovp1 + 1;
            }
        }
       
         //Rcpp::Rcout << "p1=" << p1 << std::endl;
        
         for(m=0; m<M; m++){
             theta(m) = parm(nextindex);
             for(j=0; j<ncolcovpois;j++){
                 theta(m) += parm(nextindex+j+1) * covpois(n,j);
             }
             theta(m) = exp(theta(m));
             nextindex = nextindex + ncolcovpois + 1;
         }
        
        
        for(m=0; m<M; m++) meanvec(m) = dzip(zeroprop(m), theta(m), y(n), FALSE);
     
        
        //difficult part
        for(i=0; i<M; i++){
            for(j=0;j<M;j++){
                tempmat(i,j) = xi(n-1, i) * a(i,j);
            }
        }
        
        for(j=0; j<M; j++){
            for(i=0;i<M;i++){
                if(i==0) colmax = tempmat(i,j);
                else colmax = MAX(colmax, tempmat(i,j));
            }
            tempmax(j) = colmax;
        }
        
        
        forward = tempmax % meanvec;
        forwardsum = tempmax.t() * meanvec;
        for(m=0; m<M; m++) xi(n,m) = forward(m) / forwardsum(0);
    }
    
    //termination
    for(m=0; m<M; m++){
        if(m==0) {
            colmax = xi(dim-1, m);
            state(dim-1) = m;
        }
        else{
            if(xi(dim-1,m) >= colmax) {
                state(dim-1) = m;
                colmax = xi(dim-1,m);
            }
        }
    }
    
    //backtracking
    for(n=dim-2; n>=0; n--){
        
        
        nextindex=0;
        pi(0) = 1;
        tempsum = 1;
        for(m=1; m<M; m++){
            pi(m) = parm(nextindex);
            
            for(j=0;j<ncolcovpi;j++)
                pi(m) += parm(nextindex+j+1)*covpi(n,j) ;
            
            pi(m) = exp(pi(m));
            
            tempsum += pi(m);
            nextindex = nextindex + ncolcovpi + 1;
        }
        
        for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
        
        //transition
        for(i=0; i<M; i++){
            
            for(j=0; j<M; j++){
                if(j==0) {
                    a(i,j) = 1;
                    tempsum = 1;
                }
                else {
                    a(i,j) = parm(nextindex);
                    
                    for(m=0;m<ncolcovtrans;m++)
                        a(i,j) += parm(nextindex+m+1)*covtrans(n,m);
                    
                    a(i,j) = exp(a(i,j));
                    tempsum += a(i,j);
                    nextindex = nextindex + ncolcovtrans + 1;
                }
                //Rcpp::Rcout<<a(i,j)<<std::endl;
            }
            
            for(m=0;m<M; m++) a(i,m) = a(i,m) / tempsum;
            
        }
        
        
        for(m=0;m<M;m++){
            if(m==0) {
                colmax = xi(n,m) * a(m, state(n+1));
                state(n) = m;
            }
            else{
                if(xi(n,m) * a(m, state(n+1))>=colmax) {
                    state(n) = m;
                    colmax = xi(n,m) * a(m, state(n+1));
                }
            }
        }
    }
    
    for(n=0; n<dim; n++) state(n) = state(n) + 1;
    return(state);
}


/////////////////////////

// [[Rcpp::export]]
arma::vec hsmm_cov_viterbi(arma::vec parm, int M, arma::vec y,arma::vec trunc, arma::vec zeroindex,
                           int ncolcovp, arma::mat covp, int ncolcovpi, arma::mat covpi,
                           int ncolcovomega, arma::mat covomega, int ncolcovp1, arma::mat covp1,
                           int ncolcovpois, arma::mat covpois, std::string dt_dist){
    
    double tempsum,colmax;
    
    long dim = y.n_rows;
    int i,j,m,n,nextindex;
    int dmcol = arma::max(trunc); //get the maximum
    int totalmv = arma::sum(trunc);
    arma::vec pi(M);
    arma::vec theta(M);
    arma::vec newtheta(totalmv);
    arma::vec newpi(totalmv);
    arma::vec zeroprop(M);
    arma::vec newzeroprop(totalmv);
    //Rcpp::Rcout<<dmcol<<std::endl;
    arma::mat dm(M,dmcol);
    
    arma::mat omega(M,M);
    arma::mat a(totalmv,totalmv);
    
    
    arma::vec meanvec(totalmv);
    arma::vec forward(totalmv);
    arma::vec forwardsum;
    arma::vec state(dim);
    arma::vec tempmax(totalmv);
    arma::mat xi(dim, totalmv);
    arma::mat tempmat(totalmv,totalmv);
    
    //retrieve the full parameters
    nextindex = 0;
    //dwell time
    dm.zeros();
    //check dt_dist
    if(dt_dist=="log"){
        for(i=0;i<M;i++){
            
            for(j=0;j<trunc(i);j++){
                tempsum = parm(nextindex);
                for(m=0;m<ncolcovp;m++)
                    tempsum += parm(nextindex+m+1)*covp(0,m);
                tempsum = exp(tempsum) / (1+exp(tempsum));
                dm(i,j)=dlogp(j+1, tempsum, FALSE);
                
            }
            nextindex = nextindex + ncolcovp + 1;
        }
    }else{
        for(i=0;i<M;i++){
            
            for(j=0;j<trunc(i);j++){
                tempsum = parm(nextindex);
                for(m=0;m<ncolcovp;m++)
                    tempsum += parm(nextindex+m+1)*covp(0,m);
                tempsum = exp(tempsum) ;
                dm(i,j)=dshiftpois(j+1, tempsum, 1, FALSE);
                
            }
            nextindex = nextindex + ncolcovp + 1;
        }
        
    }
    
    //recover newtheta,p1, newpi
    pi(0) = 1;
    tempsum = 1;
    for(m=1; m<M; m++){
        pi(m) = parm(nextindex);
        
        for(j=0;j<ncolcovpi;j++)
            pi(m) += parm(nextindex+j+1)*covpi(0,j) ;
        
        pi(m) = exp(pi(m));
        
        tempsum += pi(m);
        nextindex = nextindex + ncolcovpi + 1;
    }
    
    for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
    
    //p1
    
    //zero proportions
    for(i=0;i<M;i++){
        if(zeroindex(i)==0) zeroprop(i)=0;
        else{
            zeroprop(i) = parm(nextindex);
            
            for(m=0;m<ncolcovp1;m++)
                zeroprop(i) += parm(nextindex+m+1) * covp1(0,m);
            
            zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
            nextindex = nextindex + ncolcovp1 + 1;
        }
    }
        
        //theta
        for(m=0; m<M; m++){
            theta(m) = parm(nextindex);
            for(j=0; j<ncolcovpois;j++){
                theta(m) += parm(nextindex+j+1) * covpois(0,j);
            }
            theta(m) = exp(theta(m));
            nextindex = nextindex + ncolcovpois + 1;
        }
        
        //Rcpp::Rcout << nextindex << std::endl;
        //get newlogtheta,newpi
        tempsum = 0;
        for(i=0;i<M;i++){
            for(j=0;j<trunc[i];j++){
                newtheta(tempsum+j) = theta(i);
                newpi(tempsum+j) = pi(i)/trunc(i);
                newzeroprop(tempsum+j) = zeroprop(i);
                //Rcpp::Rcout << "theta=" << newtheta(tempsum+j) << "pi=" << newpi(tempsum+j) << std::endl;
            }
            tempsum += trunc(i);
        }
        
        //recover omega:
        //if M=2 then assign off-diagonal with 1; else
        //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
        
        omega.zeros();
        if(M==2) {omega(0,1)=1; omega(1,0)=1;}
        else{
            for(i=0;i<M;i++){
                tempsum = 1;
                for(j=2;j<M;j++){
                    omega(i,j) = parm(nextindex); //new
                    for(m=0;m<ncolcovomega;m++)
                        omega(i,j) += parm(nextindex+m+1)*covomega(0,m);
                    //omega(i,j) = exp(omega(i,j))/(1+exp(omega(i,j)));
                    omega(i,j) = exp(omega(i,j)); //new
                    tempsum += omega(i,j); //new
                    nextindex = nextindex + ncolcovomega + 1;
                }
                
                for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //new
                if(i==0) omega(i,1) = 1 / tempsum;
                else omega(i,0) = 1 / tempsum;
                //tempsum = 1 - arma::sum(omega.row(i)); //prevent negative values
                //if(i==0) omega(i,1) = MAX(0,tempsum);
                //else omega(i,0) = MAX(0,tempsum);
            }
            
            for(i=2;i<M;i++){
                for(j=2;j<=i;j++){
                    omega(i,j-1) = omega(i,j);
                    if(i==j) omega(i,j)=0;
                }
            }
            
        
    }  //end of poisson and zerofl part
    
    
    
    a = hsmm_hmm (omega, dm, trunc);
    
    //Rcpp::Rcout << nextindex << std::endl;
    
    
    //initialize the forward variable
    
    for(m=0; m<totalmv; m++) meanvec(m) = dzip(newzeroprop(m),newtheta(m),y(0),FALSE);
    
    
    forward = newpi % meanvec; //element wise multiplication
    forwardsum = newpi.t() * meanvec;
    
    //Rcpp::Rcout << "loglik=" << loglik << std::endl;
    
    for(m=0; m<totalmv; m++) xi(0,m) = forward(m) / forwardsum(0);
    
    
    
    //recursion
    for(n=1; n<dim; n++){
        
        nextindex = 0;
        //dwell time
        dm.zeros();
        if(dt_dist=="log"){
            for(i=0;i<M;i++){
                
                for(j=0;j<trunc(i);j++){
                    tempsum = parm(nextindex);
                    for(m=0;m<ncolcovp;m++)
                        tempsum += parm(nextindex+m+1)*covp(n,m);
                    tempsum = exp(tempsum) / (1+exp(tempsum));
                    dm(i,j)=dlogp(j+1, tempsum, FALSE);
                    
                }
                nextindex = nextindex + ncolcovp + 1;
            }
        }else{
            for(i=0;i<M;i++){
                
                for(j=0;j<trunc(i);j++){
                    tempsum = parm(nextindex);
                    for(m=0;m<ncolcovp;m++)
                        tempsum += parm(nextindex+m+1)*covp(n,m);
                    tempsum = exp(tempsum) ;
                    dm(i,j)=dshiftpois(j+1, tempsum, 1, FALSE);
                    
                }
                nextindex = nextindex + ncolcovp + 1;
            }
            
        }

        //recover newtheta,p1, newpi
        pi(0) = 1;
        tempsum = 1;
        for(m=1; m<M; m++){
            pi(m) = parm(nextindex);
            
            for(j=0;j<ncolcovpi;j++)
                pi(m) += parm(nextindex+j+1)*covpi(n,j) ;
            
            pi(m) = exp(pi(m));
            
            tempsum += pi(m);
            nextindex = nextindex + ncolcovpi + 1;
        }
        
        for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
        
        //p1
        for(i=0;i<M;i++){
            if(zeroindex(i)==0) zeroprop(i)=0;
            else{
                zeroprop(i) = parm(nextindex);
                
                for(m=0;m<ncolcovp1;m++)
                    zeroprop(i) += parm(nextindex+m+1) * covp1(n,m);
                
                zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
                nextindex = nextindex + ncolcovp1 + 1;
            }
        }
        
            //theta
            for(m=0; m<M; m++){
                theta(m) = parm(nextindex);
                for(j=0; j<ncolcovpois;j++){
                    theta(m) += parm(nextindex+j+1) * covpois(n,j);
                }
                theta(m) = exp(theta(m));
                nextindex = nextindex + ncolcovpois + 1;
            }
            
            //Rcpp::Rcout << nextindex << std::endl;
            //get newlogtheta,newpi
            tempsum = 0;
            for(i=0;i<M;i++){
                for(j=0;j<trunc[i];j++){
                    newtheta(tempsum+j) = theta(i);
                    newpi(tempsum+j) = pi(i)/trunc(i);
                    newzeroprop(tempsum+j) = zeroprop(i);
                    //Rcpp::Rcout << "theta=" << newtheta(tempsum+j) << "pi=" << newpi(tempsum+j) << std::endl;
                }
                tempsum += trunc(i);
            }
            
            
            //recover omega:
            //if M=2 then assign off-diagonal with 1; else
            //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
            
            omega.zeros();
            if(M==2) {omega(0,1)=1; omega(1,0)=1;}
            else{
                for(i=0;i<M;i++){
                    tempsum = 1;
                    for(j=2;j<M;j++){
                        omega(i,j) = parm(nextindex); //new
                        for(m=0;m<ncolcovomega;m++)
                            omega(i,j) += parm(nextindex+m+1)*covomega(n,m);
                        //omega(i,j) = exp(omega(i,j))/(1+exp(omega(i,j)));
                        omega(i,j) = exp(omega(i,j)); //new
                        tempsum += omega(i,j); //new
                        nextindex = nextindex + ncolcovomega + 1;
                    }
                    
                    for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //new
                    if(i==0) omega(i,1) = 1 / tempsum;
                    else omega(i,0) = 1 / tempsum;
                    //tempsum = 1 - arma::sum(omega.row(i)); //prevent negative values
                    //if(i==0) omega(i,1) = MAX(0,tempsum);
                    //else omega(i,0) = MAX(0,tempsum);
                }
                
                
                for(i=2;i<M;i++){
                    for(j=2;j<=i;j++){
                        omega(i,j-1) = omega(i,j);
                        if(i==j) omega(i,j)=0;
                    }
                }
                
            }
        //end section for zeroinflated poisson distribution
    
    
        //Rcpp::Rcout<<"nextindex="<<nextindex<<std::endl;
        // recover transition matrix
        
        a = hsmm_hmm (omega, dm, trunc);
        
        
         for(m=0; m<totalmv; m++) meanvec(m) = dzip(newzeroprop(m),newtheta(m),y(n),FALSE);
        
        
        
        //difficult part
        for(i=0; i<totalmv; i++){
            for(j=0;j<totalmv;j++){
                tempmat(i,j) = xi(n-1, i) * a(i,j);
            }
        }
        
        for(j=0; j<totalmv; j++){
            for(i=0;i<totalmv;i++){
                if(i==0) colmax = tempmat(i,j);
                else colmax = MAX(colmax, tempmat(i,j));
            }
            tempmax(j) = colmax;
        }
        
        
        forward = tempmax % meanvec;
        forwardsum = tempmax.t() * meanvec;
        for(m=0; m<totalmv; m++) xi(n,m) = forward(m) / forwardsum(0);
    }
    
    //termination
    for(m=0; m<totalmv; m++){
        if(m==0) {
            colmax = xi(dim-1, m);
            state(dim-1) = m;
        }
        else{
            if(xi(dim-1,m) >= colmax) {
                state(dim-1) = m;
                colmax = xi(dim-1,m);
            }
        }
    }
    
    //backtracking
    for(n=dim-2; n>=0; n--){
        
        
        //retrieve the full parameters
        nextindex = 0;
        //dwell time
        dm.zeros();
        if(dt_dist=="log"){
            for(i=0;i<M;i++){
                
                for(j=0;j<trunc(i);j++){
                    tempsum = parm(nextindex);
                    for(m=0;m<ncolcovp;m++)
                        tempsum += parm(nextindex+m+1)*covp(n,m);
                    tempsum = exp(tempsum) / (1+exp(tempsum));
                    dm(i,j)=dlogp(j+1, tempsum, FALSE);
                    
                }
                nextindex = nextindex + ncolcovp + 1;
            }
        }else{
            for(i=0;i<M;i++){
                
                for(j=0;j<trunc(i);j++){
                    tempsum = parm(nextindex);
                    for(m=0;m<ncolcovp;m++)
                        tempsum += parm(nextindex+m+1)*covp(n,m);
                    tempsum = exp(tempsum) ;
                    dm(i,j)=dshiftpois(j+1, tempsum, 1, FALSE);
                    
                }
                nextindex = nextindex + ncolcovp + 1;
            }
            
        }
        
        
        //recover newtheta,p1, newpi
        pi(0) = 1;
        tempsum = 1;
        for(m=1; m<M; m++){
            pi(m) = parm(nextindex);
            
            for(j=0;j<ncolcovpi;j++)
                pi(m) += parm(nextindex+j+1)*covpi(n,j) ;
            
            pi(m) = exp(pi(m));
            
            tempsum += pi(m);
            nextindex = nextindex + ncolcovpi + 1;
        }
        
        for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
        
        //p1
        for(i=0;i<M;i++){
            if(zeroindex(i)==0) zeroprop(i)=0;
            else{
                zeroprop(i) = parm(nextindex);
                
                for(m=0;m<ncolcovp1;m++)
                    zeroprop(i) += parm(nextindex+m+1) * covp1(n,m);
                
                zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
                nextindex = nextindex + ncolcovp1 + 1;
            }
        }
        
        //theta
        for(m=0; m<M; m++){
            theta(m) = parm(nextindex);
            for(j=0; j<ncolcovpois;j++){
                theta(m) += parm(nextindex+j+1) * covpois(n,j);
            }
            theta(m) = exp(theta(m));
            nextindex = nextindex + ncolcovpois + 1;
        }
        
        //Rcpp::Rcout << nextindex << std::endl;
        //get newlogtheta,newpi
        tempsum = 0;
        for(i=0;i<M;i++){
            for(j=0;j<trunc[i];j++){
                newtheta(tempsum+j) = theta(i);
                newpi(tempsum+j) = pi(i)/trunc(i);
                newzeroprop(tempsum+j) = zeroprop(i);
                //Rcpp::Rcout << "theta=" << newtheta(tempsum+j) << "pi=" << newpi(tempsum+j) << std::endl;
            }
            tempsum += trunc(i);
        }
        
        
        //recover omega:
        //if M=2 then assign off-diagonal with 1; else
        //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
        
        omega.zeros();
        if(M==2) {omega(0,1)=1; omega(1,0)=1;}
        else{
            for(i=0;i<M;i++){
                tempsum = 1;
                for(j=2;j<M;j++){
                    omega(i,j) = parm(nextindex); //new
                    for(m=0;m<ncolcovomega;m++)
                        omega(i,j) += parm(nextindex+m+1)*covomega(n,m);
                    //omega(i,j) = exp(omega(i,j))/(1+exp(omega(i,j)));
                    omega(i,j) = exp(omega(i,j)); //new
                    tempsum += omega(i,j); //new
                    nextindex = nextindex + ncolcovomega + 1;
                }
                
                for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //new
                if(i==0) omega(i,1) = 1 / tempsum;
                else omega(i,0) = 1 / tempsum;
                //tempsum = 1 - arma::sum(omega.row(i)); //prevent negative values
                //if(i==0) omega(i,1) = MAX(0,tempsum);
                //else omega(i,0) = MAX(0,tempsum);
            }
            
            for(i=2;i<M;i++){
                for(j=2;j<=i;j++){
                    omega(i,j-1) = omega(i,j);
                    if(i==j) omega(i,j)=0;
                }
            }
            
        }

        
        //Rcpp::Rcout<<"nextindex="<<nextindex<<std::endl;
        // recover transition matrix
        
        a = hsmm_hmm (omega, dm, trunc);
        
        
        for(m=0;m<totalmv;m++){
            if(m==0) {
                colmax = xi(n,m) * a(m, state(n+1));
                state(n) = m;
            }
            else{
                if(xi(n,m) * a(m, state(n+1))>=colmax) {
                    state(n) = m;
                    colmax = xi(n,m) * a(m, state(n+1));
                }
            }
        }
    }
    
    for(n=0; n<dim; n++) state(n) = state(n) + 1;
    
    for(i=0;i<dim;i++){
        tempsum=0;
        for(j=0;j<M;j++){
            if(state(i)<=tempsum+trunc(j)){
                state(i)=j+1;
                break;
            }
            tempsum += trunc(j);
        }
    }
    
    
    return(state);
}

///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double ziploglik_nocov(arma::vec delta, arma::mat gamma, double theta, arma::vec lambda,
                       arma::vec y, arma::vec ntimes){
    double loglik = 0;
    //long dim = y.n_rows;
    int M = lambda.n_rows;
    int n = ntimes.n_rows;
    int j,m,t;
    double tempsum;
    
    arma::vec forward(M);
    arma::vec meanvec(M);
    arma::vec forwardsum;  //then forwardsum(0) is double
    arma::mat tempmat(1,M);
    
    int count = 0;
    
    for(j=0; j<n; j++){
        
        tempsum = 0;
        //initialize the forward variable
        meanvec(0) = dzip(theta, lambda(0), y(count), FALSE);
        for(m=1; m<M; m++) meanvec(m) = R::dpois(y(count), lambda(m), FALSE);
        
        forward = delta % meanvec; //element wise multiplication
        forwardsum = delta.t() * meanvec;
        tempsum = log(forwardsum(0));
        //Rcpp::Rcout << "loglik=" << loglik << std::endl;
        
        for(m=0; m<M; m++) forward(m) = forward(m) / forwardsum(0);
        
        //recursion
        for(t=1; t<ntimes(j); t++){
            
            meanvec(0) = dzip(theta, lambda(0), y(count+t), FALSE);
            for(m=1; m<M; m++) meanvec(m) = R::dpois(y(count+t), lambda(m), FALSE);
            
            tempmat = forward.t() * gamma; //row matrix
            forward = tempmat.t() % meanvec;
            forwardsum = tempmat * meanvec;
            tempsum += log(forwardsum(0));
            for(m=0; m<M; m++) forward(m) = forward(m) / forwardsum(0);
        }
        
        loglik += tempsum;
        count += ntimes(j);
        
    }
    
    return loglik;
    
}

///////////////////////////////////
// [[Rcpp::export]]
Rcpp::List retrieve_nocov(arma::vec parm, int M){
    int nextindex = 0;
    arma::vec delta(M);
    arma::mat gamma(M,M);
    double theta;
    arma::vec lambda(M);
    int i, j, m;
    double tempsum;
    
    //retrieve the full parameters
    delta(0) = 1;
    tempsum = 1;
    for(m=1; m<M; m++) {
        delta(m) = exp(parm(m-1));
        tempsum += delta(m);
    }
    for(m=0; m<M; m++) delta(m) = delta(m) / tempsum;
    
    for(i=0; i<M; i++){
        for(j=0; j<M; j++){
            if(j==0) {
                gamma(i,j) = 1;
                tempsum = 1;
            }
            else {
                gamma(i,j) = exp(parm((i+1)*(M-1)+j-1));
                tempsum += gamma(i,j);
            }
        }
        
        for(m=0;m<M; m++) gamma(i,m) = gamma(i,m) / tempsum;
        
    }
    
    nextindex = M*M-1;
    theta = exp(parm(nextindex))/(1+exp(parm(nextindex)));
    
    for(m=0; m<M; m++) lambda(m) = exp(parm(nextindex+m+1));
    
    return Rcpp::List::create(Rcpp::Named("delta")=delta,
                              Rcpp::Named("gamma")=gamma,
                              Rcpp::Named("theta")=theta,
                              Rcpp::Named("lambda")=lambda);
    
}




//////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double zipnegloglik_nocov(arma::vec parm, int M, arma::vec y, arma::vec ntimes){
    
    Rcpp::List mod = retrieve_nocov(parm,M);
    arma::vec delta = mod("delta");
    arma::mat gamma = mod("gamma");
    double theta = mod("theta");
    arma::vec lambda = mod("lambda");
    
    double negloglik = - ziploglik_nocov(delta, gamma, theta, lambda, y, ntimes);
    return negloglik;
}


////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List forwardbackward(arma::vec Pi, arma::mat P, arma::mat nodeprob,
                           arma::vec y, arma::vec ntimes){
    
    long dim = y.n_rows;
    int M = nodeprob.n_cols;
    int n = ntimes.n_rows;
    int j,m,t,i,k;
    double tempsum;
    arma::vec tempval;
    arma::mat tempmat(1,M);
    arma::mat tempmat2(M,M);
    
    arma::mat alpha(dim, M);
    arma::vec scale(dim);
    arma::mat beta(dim, M);
    arma::mat Gamma(dim, M);
    arma::mat xi(dim-1, M*M);
    
    arma::vec colsumgamma(M);
    arma::vec tempsumxi(M*M);
    arma::mat colsumxi(M,M);
    
    int count = 0;
    for(j=0; j<n; j++){
        
        //forward probs
        alpha.row(count) = Pi.t() % nodeprob.row(count);
        tempval = nodeprob.row(count) * Pi;
        scale(count) = tempval(0);//double
        for(m=0; m<M; m++) alpha(count, m) = alpha(count,m) / scale(count);
        
        for(t=1; t<ntimes(j); t++){
            tempmat = alpha.row(count+t-1) * P; //row matrix
            alpha.row(count+t) = tempmat % nodeprob.row(count+t);
            tempval = tempmat * nodeprob.row(count+t).t();
            scale(count+t) = tempval(0);//double
            for(m=0; m<M; m++) alpha(count+t, m) = alpha(count+t,m) / scale(count+t);
        }
        
        //backward probs and state conditional probs
        for(m=0; m<M; m++) beta(count + ntimes(j) - 1,m) = 1 / (M * scale(count + ntimes(j) - 1));
        Gamma.row(count+ntimes(j)-1) = alpha.row(count+ntimes(j)-1) % beta.row(count+ntimes(j)-1);
        tempval = alpha.row(count+ntimes(j)-1) * beta.row(count+ntimes(j)-1).t();
        for(m=0; m<M; m++) Gamma(count+ntimes(j)-1,m)=Gamma(count+ntimes(j)-1,m)/tempval(0); //double
        
        for(t=ntimes(j)-2; t>=0; t--){
            tempmat = beta.row(count+t+1) % nodeprob.row(count+t+1);
            beta.row(count+t) = tempmat * P.t();
            for(m=0; m<M; m++) beta(count+t,m) = beta(count+t,m) / scale(count + t);
            
            Gamma.row(count+t) = alpha.row(count+t) % beta.row(count+t);
            tempval = alpha.row(count+t) * beta.row(count+t).t();
            for(m=0; m<M; m++) Gamma(count+t,m)=Gamma(count+t,m)/tempval(0); //double
        }
        
        //transition probs
        for(t=0; t<ntimes(j)-1; t++){
            tempmat = beta.row(count+t+1) % nodeprob.row(count+t+1);
            tempmat2 = P % (alpha.row(count+t).t() * tempmat);
            
            tempsum = 0;
            for(i=0; i<M; i++){
                for(k=0; k<M; k++){
                    xi(count+t, i + k*M) = tempmat2(i,k);
                    tempsum += xi(count+t, i + k*M);
                }
            }
            
            for(m=0; m<M*M; m++) xi(count+t, m) = xi(count+t, m) / tempsum;
        }
        
        count += ntimes(j);
    }
    
    
    //get the column sums
    colsumgamma = colsum(Gamma);
    tempsumxi = colsum(xi);
    
    for(i=0; i<M; i++)
        for(k=0; k<M; k++)
            colsumxi(i,k) = tempsumxi(i+k*M);
    
    
    return Rcpp::List::create(Rcpp::Named("colsumgamma")=colsumgamma,
                              Rcpp::Named("colsumxi")=colsumxi,
                              Rcpp::Named("Gamma")=Gamma,
                              Rcpp::Named("xi")=xi);
}

////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::vec grad_ziploglik_nocov(arma::vec delta, arma::mat gamma, double theta, arma::vec lambda,
                               arma::vec y, arma::vec ntimes){
    
    long dim = y.n_rows;
    int M = lambda.n_rows;
    int n = ntimes.n_rows;
    int i,j,m,t,k;
    arma::mat nodeprob(dim, M);
    
    
    int count = 0;
    for(j=0; j<n; j++){
        for(t=0; t<ntimes(j); t++){
            nodeprob(count+t,0) = dzip(theta, lambda(0), y(count+t), FALSE);
            for(m=1; m<M; m++) nodeprob(count+t,m) = R::dpois(y(count+t), lambda(m), FALSE);
        }
        count += ntimes(j);
    }
    
    Rcpp::List fb = forwardbackward(delta, gamma, nodeprob, y, ntimes);
    arma::mat Gamma = fb("Gamma");
    arma::vec colsumgamma = fb("colsumgamma");
    arma::mat colsumxi = fb("colsumxi");
    //arma::mat xi = fb("xi");
    
    //get the derivatives
    //prior
    arma::vec temppi(M);
    arma::vec dlogitpi(M);
    
    for(m=0; m<M; m++) temppi(m) = 0;
    count = 0;
    for(j=0; j<n; j++){
        for(m=0; m<M; m++) temppi(m) += Gamma(count,m);
        count += ntimes(j);
    }
    
    for(m=0; m<M; m++) {
        temppi(m) = temppi(m) / n;
        dlogitpi(m) = temppi(m) - delta(m);
    }
    
    
    //tpm
    arma::vec rowsumtrans = rowsum(colsumxi);
    
    arma::mat dlogtrans(M,M);
    for(i=0; i<M; i++)
        for(k=0; k<M; k++)
            dlogtrans(i,k) = colsumxi(i,k) - rowsumtrans(i) * gamma(i,k);
    
    //theta and lambda
    double dlogtheta = 0;
    arma::vec dloglambda(M);
    for(m=0; m<M; m++) dloglambda(m) = 0;
    
    for(t=0; t<dim; t++){
        
        if(y(t)>0){
            dlogtheta -= theta * Gamma(t,0);
            dloglambda(0) += Gamma(t,0) * (y(t) - lambda(0));
        }
        else{
            dlogtheta -= Gamma(t,0)*theta*(1-theta)*(exp(-lambda(0))-1) / (theta+(1-theta)*exp(-lambda(0)));
            dloglambda(0) += Gamma(t,0)*(theta-1)*lambda(0)*exp(-lambda(0)) / (theta+(1-theta)*exp(-lambda(0)));
        }
        
        for(m=1; m<M; m++){
            dloglambda(m) += Gamma(t,m)*(y(t) - lambda(m));
        }
    }
    
    
    arma::vec result(M-1+M*(M-1)+1+M);
    int nextindex = 0;
    for(j=0; j<M-1; j++) result(j) = dlogitpi(j+1);
    nextindex = M - 1;
    for(i=0; i<M; i++){
        for(j=1; j<M; j++){
            result(nextindex + i*(M-1) + j - 1) = dlogtrans(i,j);
        }
    }
    nextindex += M * (M-1);
    result(nextindex) = dlogtheta;
    nextindex += 1;
    for(i=0; i<M; i++) result(nextindex+i) = dloglambda(i);
    
    
    return result;
    
    
}


//////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::vec grad_zipnegloglik_nocov(arma::vec parm, int M, arma::vec y, arma::vec ntimes){
    
    Rcpp::List mod = retrieve_nocov(parm,M);
    arma::vec delta = mod("delta");
    arma::mat gamma = mod("gamma");
    double theta = mod("theta");
    arma::vec lambda = mod("lambda");
    
    ////change to negative !!!!!!
    arma::vec result = grad_ziploglik_nocov(delta, gamma, theta, lambda, y, ntimes);
    int m;
    for(m=0; m<M-1+M*(M-1)+1+M; m++) result(m) = -result(m);
    return result;
}




////////////////////
//////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double ziploglik_cov(arma::vec delta, arma::mat gamma, arma::vec thetaparm, arma::mat lambdaparm,
                     arma::vec y, arma::mat x, arma::vec ntimes){
    //the first column of x should be 1 for intercept
    double loglik = 0;
    //long dim = y.n_rows;
    int M = gamma.n_rows;
    int n = ntimes.n_rows;
    int j,m,t;
    arma::vec tempval;
    double tempsum;
    double theta;
    arma::vec lambda(M);
    
    arma::vec forward(M);
    arma::vec meanvec(M);
    arma::vec forwardsum;  //then forwardsum(0) is double
    arma::mat tempmat(1,M);
    
    int count = 0;
    
    for(j=0; j<n; j++){
        
        tempsum = 0;
        
        //initialize the forward variable
        tempval = x.row(count) * thetaparm;
        theta = exp(tempval(0))/(1+exp(tempval(0)));
        for(m=0;m<M;m++) {
            tempval = x.row(count) * lambdaparm.row(m).t();
            lambda(m) = exp(tempval(0));
        }
        
        meanvec(0) = dzip(theta, lambda(0), y(count), FALSE);
        for(m=1; m<M; m++) meanvec(m) = R::dpois(y(count), lambda(m), FALSE);
        
        forward = delta % meanvec; //element wise multiplication
        forwardsum = delta.t() * meanvec;
        tempsum = log(forwardsum(0));
        //Rcpp::Rcout << "loglik=" << loglik << std::endl;
        
        for(m=0; m<M; m++) forward(m) = forward(m) / forwardsum(0);
        
        //recursion
        for(t=1; t<ntimes(j); t++){
            
            tempval = x.row(count+t) * thetaparm;
            theta = exp(tempval(0))/(1+exp(tempval(0)));
            for(m=0;m<M;m++) {
                tempval = x.row(count+t) * lambdaparm.row(m).t();
                lambda(m) = exp(tempval(0));
            }
            
            meanvec(0) = dzip(theta, lambda(0), y(count+t), FALSE);
            for(m=1; m<M; m++) meanvec(m) = R::dpois(y(count+t), lambda(m), FALSE);
            
            tempmat = forward.t() * gamma; //row matrix
            forward = tempmat.t() % meanvec;
            forwardsum = tempmat * meanvec;
            tempsum += log(forwardsum(0));
            for(m=0; m<M; m++) forward(m) = forward(m) / forwardsum(0);
        }
        
        loglik += tempsum;
        count += ntimes(j);
        
    }
    
    return loglik;
    
}

////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////
// [[Rcpp::export]]
Rcpp::List retrieve_cov(arma::vec parm, int M, int ncolx){
    int nextindex = 0;
    arma::vec delta(M);
    arma::mat gamma(M,M);
    arma::vec thetaparm(ncolx);
    arma::mat lambdaparm(M,ncolx);
    int i, j, m;
    double tempsum;
    
    //retrieve the full parameters
    delta(0) = 1;
    tempsum = 1;
    for(m=1; m<M; m++) {
        delta(m) = exp(parm(m-1));
        tempsum += delta(m);
    }
    for(m=0; m<M; m++) delta(m) = delta(m) / tempsum;
    
    for(i=0; i<M; i++){
        for(j=0; j<M; j++){
            if(j==0) {
                gamma(i,j) = 1;
                tempsum = 1;
            }
            else {
                gamma(i,j) = exp(parm((i+1)*(M-1)+j-1));
                tempsum += gamma(i,j);
            }
        }
        
        for(m=0;m<M; m++) gamma(i,m) = gamma(i,m) / tempsum;
        
    }
    
    nextindex = M*M-1;
    for(i=0; i<ncolx; i++) thetaparm(i) = parm(nextindex+i);
    nextindex += ncolx;
    
    for(i=0; i<M; i++){
        for(j=0; j<ncolx; j++){
            lambdaparm(i,j) = parm(nextindex + i*ncolx + j);
        }
    }
    
    return Rcpp::List::create(Rcpp::Named("delta")=delta,
                              Rcpp::Named("gamma")=gamma,
                              Rcpp::Named("thetaparm")=thetaparm,
                              Rcpp::Named("lambdaparm")=lambdaparm);
    
}


//////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double zipnegloglik_cov(arma::vec parm, arma::vec y, arma::mat covariates, int M, arma::vec ntimes){
    int ncolx = covariates.n_cols;
    Rcpp::List mod = retrieve_cov(parm,M,ncolx);
    arma::vec delta = mod("delta");
    arma::mat gamma = mod("gamma");
    arma::vec thetaparm = mod("thetaparm");
    arma::mat lambdaparm = mod("lambdaparm");
    
    double negloglik = - ziploglik_cov(delta, gamma, thetaparm, lambdaparm, y, covariates, ntimes);
    return negloglik;
}


////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::vec grad_ziploglik_cov(arma::vec delta, arma::mat gamma, arma::vec thetaparm, arma::mat lambdaparm,
                             arma::vec y, arma::mat x, arma::vec ntimes){
    
    long dim = y.n_rows;
    int ncolx = x.n_cols;
    int M = gamma.n_rows;
    int n = ntimes.n_rows;
    int i,j,m,t,k;
    double theta;
    arma::vec tempval;
    arma::vec lambda(M);
    arma::mat nodeprob(dim, M);
    
    
    int count = 0;
    for(j=0; j<n; j++){
        for(t=0; t<ntimes(j); t++){
            tempval = x.row(count+t) * thetaparm;
            theta = exp(tempval(0))/(1+exp(tempval(0)));
            for(m=0;m<M;m++) {
                tempval = x.row(count+t) * lambdaparm.row(m).t();
                lambda(m) = exp(tempval(0));
            }
            
            nodeprob(count+t,0) = dzip(theta, lambda(0), y(count+t), FALSE);
            for(m=1; m<M; m++) nodeprob(count+t,m) = R::dpois(y(count+t), lambda(m), FALSE);
        }
        count += ntimes(j);
    }
    
    Rcpp::List fb = forwardbackward(delta, gamma, nodeprob, y, ntimes);
    arma::mat Gamma = fb("Gamma");
    arma::vec colsumgamma = fb("colsumgamma");
    arma::mat colsumxi = fb("colsumxi");
    //arma::mat xi = fb("xi");
    
    //get the derivatives
    //prior
    arma::vec temppi(M);
    arma::vec dlogitpi(M);
    
    for(m=0; m<M; m++) temppi(m) = 0;
    count = 0;
    for(j=0; j<n; j++){
        for(m=0; m<M; m++) temppi(m) += Gamma(count,m);
        count += ntimes(j);
    }
    
    for(m=0; m<M; m++) {
        temppi(m) = temppi(m) / n;
        dlogitpi(m) = temppi(m) - delta(m);
    }
    
    
    //tpm
    arma::vec rowsumtrans = rowsum(colsumxi);
    
    arma::mat dlogtrans(M,M);
    for(i=0; i<M; i++)
        for(k=0; k<M; k++)
            dlogtrans(i,k) = colsumxi(i,k) - rowsumtrans(i) * gamma(i,k);
    
    //theta and lambda
    arma::vec dlogtheta(ncolx);
    arma::mat dloglambda(M,ncolx);
    for(m=0; m<ncolx; m++) dlogtheta(m) = 0;
    for(i=0; i<M; i++)
        for(j=0; j<ncolx; j++)
            dloglambda(i,j) = 0;
    
    for(t=0; t<dim; t++){
        
        tempval = x.row(t) * thetaparm;
        theta = exp(tempval(0))/(1+exp(tempval(0)));
        for(m=0;m<M;m++) {
            tempval = x.row(t) * lambdaparm.row(m).t();
            lambda(m) = exp(tempval(0));
        }
        
        if(y(t)>0){
            for(i=0;i<ncolx;i++){
                dlogtheta(i) -= theta * Gamma(t,0) * x(t,i) ;
                dloglambda(0,i) += Gamma(t,0) * x(t,i) * (y(t) - lambda(0));
            }
        }
        else{
            for(i=0;i<ncolx;i++){
                dlogtheta(i) -= Gamma(t,0)*x(t,i)*theta*(1-theta)*(exp(-lambda(0))-1) / (theta+(1-theta)*exp(-lambda(0)));
                dloglambda(0,i) += Gamma(t,0)*x(t,i)*(theta-1)*lambda(0)*exp(-lambda(0)) / (theta+(1-theta)*exp(-lambda(0)));
            }
        }
        
        for(m=1; m<M; m++){
            for(i=0;i<ncolx;i++){
                dloglambda(m,i) += Gamma(t,m)*x(t,i)*(y(t) - lambda(m));
            }
        }
    }
    
    
    arma::vec result(M-1+M*(M-1)+ncolx*(1+M));
    int nextindex = 0;
    for(j=0; j<M-1; j++) result(j) = dlogitpi(j+1);
    nextindex = M - 1;
    for(i=0; i<M; i++){
        for(j=1; j<M; j++){
            result(nextindex + i*(M-1) + j - 1) = dlogtrans(i,j);
        }
    }
    nextindex += M * (M-1);
    for(i=0; i<ncolx; i++)
        result(nextindex+i) = dlogtheta(i);
    
    nextindex += ncolx;
    for(i=0; i<M; i++)
        for(j=0; j<ncolx; j++)
            result(nextindex+i*ncolx+j) = dloglambda(i,j);
    
    
    return result;
    
    
}


//////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::vec grad_zipnegloglik_cov(arma::vec parm, arma::vec y, arma::mat covariates, int M, arma::vec ntimes){
    int ncolx = covariates.n_cols;
    Rcpp::List mod = retrieve_cov(parm,M,ncolx);
    arma::vec delta = mod("delta");
    arma::mat gamma = mod("gamma");
    arma::vec thetaparm = mod("thetaparm");
    arma::mat lambdaparm = mod("lambdaparm");
    
    ////change to negative !!!!!!
    arma::vec result = grad_ziploglik_cov(delta, gamma, thetaparm, lambdaparm, y, covariates, ntimes);
    int m;
    for(m=0; m<M-1+M*(M-1)+ncolx*(1+M); m++) result(m) = -result(m);
    return result;
}


////////////////////////////////////////////////
///////////////convert nupdate to the matrix of nbetween
// [[Rcpp::export]]
arma::mat nupdate_to_nbetween(int nupdate, int nmin, long dim, arma::vec ntimes){
    int n = ntimes.n_rows;
    int unit;
    long cumsum;
    arma::mat nbetween(n, nupdate);
    int i,j;
    
    for(i=0; i<n; i++){
        cumsum = 0;
        unit = (ntimes(i) - nmin) / (nupdate - 1);
        for(j=0; j<nupdate; j++){
            if(j==0){
                nbetween(i,j) = nmin;
                cumsum = nmin;
            }
            else if(j < nupdate-1) {
                nbetween(i, j) = unit;
                cumsum += unit;
            }
            else{
                nbetween(i, j) = ntimes(i) - cumsum;
            }
        }
    }
    
    return nbetween;
    
}


/////////////////////////////////////////////////////////
//stochastic gradient descent for zip_nocov
////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat sgd_zip_nocov(arma::vec parm, arma::vec y, int M, arma::vec ntimes,
                        int nmin, int nupdate, double power, arma::vec rate){
    
    long dim = y.n_rows;
    int n = ntimes.n_rows;
    int i,j,t;
    double ifloat = 0;
    //double coef;
    arma::vec coef;
    int cumsum1, cumsum2;
    
    arma::mat nbetween = nupdate_to_nbetween(nupdate, nmin, dim, ntimes);
    //only 3 different lengths of y
    arma::vec length = colsum(nbetween);
    arma::vec yhead(length(0));
    arma::vec ybody(length(2));
    arma::vec ytail(length(nupdate-1));
    
    arma::vec pointer(n);
    for(i=0; i<n; i++) pointer(i) = 0;
    
    arma::vec neggradient(M-1+M*(M-1)+1+M);
    //y.subvec(cumtimes, cumtimes + ntimes[i] - 1)
    
    arma::mat result(nupdate, M-1+M*(M-1)+1+M);
    
    
    //stochastic gradient descent
    for(i=0; i<nupdate; i++){
        
        if(i==0){
            //assign y
            cumsum1 = 0;
            cumsum2 = 0;
            for(j=0; j<n; j++){
                for(t=0; t<nbetween(j,i); t++){
                    yhead(cumsum1+t) = y(cumsum2+t);
                }
                cumsum1 += nbetween(j,i);
                cumsum2 += ntimes(j);
            }
            pointer += nbetween.col(i);
            
            //get the forward probability
            //forward = forward_nocov(parm, M, yhead, nbetween.col(i));
            
            //get the gradient
            neggradient = grad_zipnegloglik_nocov(parm, M, yhead, nbetween.col(i));
            
        }
        else if(i<nupdate - 1){
            cumsum1 = 0;
            cumsum2 = 0;
            for(j=0; j<n; j++){
                for(t=0; t<nbetween(j,i); t++){
                    ybody(cumsum1+t) = y(cumsum2+pointer(j)+t);
                }
                cumsum1 += nbetween(j,i);
                cumsum2 += ntimes(j);
            }
            pointer += nbetween.col(i);
            
            //replace the delta in the parm for gradient calculation
            //for(m=0; m<M-1; m++) parm(m) = log(forward(m+1)) - log(forward(0));
            
            //get the forward probability
            //forward = forward_nocov(parm, M, ybody, nbetween.col(i));
            
            //get the gradient
            neggradient = grad_zipnegloglik_nocov(parm, M, ybody, nbetween.col(i));
            
        }
        else{
            cumsum1 = 0;
            cumsum2 = 0;
            for(j=0; j<n; j++){
                for(t=0; t<nbetween(j,i); t++){
                    ytail(cumsum1+t) = y(cumsum2+pointer(j)+t);
                }
                cumsum1 += nbetween(j,i);
                cumsum2 += ntimes(j);
            }
            pointer += nbetween.col(i);
            
            //replace the delta in the parm for gradient calculation
            //for(m=0; m<M-1; m++) parm(m) = log(forward(m+1)) - log(forward(0));
            
            //get the forward probability
            //forward = forward_nocov(parm, M, ytail, nbetween.col(i));
            
            neggradient = grad_zipnegloglik_nocov(parm, M, ytail, nbetween.col(i));
            
        }
        
        //coef = (i+1)^(-power)
        coef = exp(-power * log(ifloat+1)) * rate / cumsum1;
        ifloat ++;
        
        parm -= coef % neggradient;
        //Rcpp::Rcout<<i<<std::endl;
        //Rcpp::Rcout<<coef<<std::endl;
        //Rcpp::Rcout<<cumsum1<<std::endl;
        result.row(i) = parm.t();
        //result.row(i) = neggradient.t();
    }
    
    return result;
    
}


/////////////////////////////////////////////////////////
//stochastic gradient descent for zip_cov
////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat sgd_zip_cov(arma::vec parm, arma::vec y, arma::mat x, int M, arma::vec ntimes,
                      int nmin, int nupdate, double power, arma::vec rate){
    
    long dim = y.n_rows;
    int n = ntimes.n_rows;
    int ncolx = x.n_cols;
    //Rcpp::List mod = retrieve_cov(parm,M,ncolx);
    int i,j,t;
    //double coef;
    arma::vec coef;
    int cumsum1, cumsum2;
    double ifloat = 0;
    arma::mat nbetween = nupdate_to_nbetween(nupdate, nmin, dim, ntimes);
    //only 3 different lengths of y
    arma::vec length = colsum(nbetween);
    arma::vec yhead(length(0));
    arma::vec ybody(length(2));
    arma::vec ytail(length(nupdate-1));
    arma::mat xhead(length(0), ncolx);
    arma::mat xbody(length(2), ncolx);
    arma::mat xtail(length(nupdate-1), ncolx);
    
    arma::vec pointer(n);
    for(i=0; i<n; i++) pointer(i) = 0;
    
    arma::vec neggradient(M-1+M*(M-1)+ncolx*(1+M));
    //y.subvec(cumtimes, cumtimes + ntimes[i] - 1)
    
    arma::mat result(nupdate, M-1+M*(M-1)+ncolx*(1+M));
    
    
    //stochastic gradient descent
    for(i=0; i<nupdate; i++){
        
        if(i==0){
            //assign y
            cumsum1 = 0;
            cumsum2 = 0;
            for(j=0; j<n; j++){
                for(t=0; t<nbetween(j,i); t++){
                    yhead(cumsum1+t) = y(cumsum2+t);
                    xhead.row(cumsum1+t) = x.row(cumsum2+t);
                }
                cumsum1 += nbetween(j,i);
                cumsum2 += ntimes(j);
            }
            pointer += nbetween.col(i);
            
            //get the forward probability
            //forward = forward_nocov(parm, M, yhead, nbetween.col(i));
            
            //get the gradient
            neggradient = grad_zipnegloglik_cov(parm, yhead, xhead, M, nbetween.col(i));
            
        }
        else if(i<nupdate - 1){
            cumsum1 = 0;
            cumsum2 = 0;
            for(j=0; j<n; j++){
                for(t=0; t<nbetween(j,i); t++){
                    ybody(cumsum1+t) = y(cumsum2+pointer(j)+t);
                    xbody.row(cumsum1+t) = x.row(cumsum2+pointer(j)+t);
                }
                cumsum1 += nbetween(j,i);
                cumsum2 += ntimes(j);
            }
            pointer += nbetween.col(i);
            
            //replace the delta in the parm for gradient calculation
            //for(m=0; m<M-1; m++) parm(m) = log(forward(m+1)) - log(forward(0));
            
            //get the forward probability
            //forward = forward_nocov(parm, M, ybody, nbetween.col(i));
            
            //get the gradient
            neggradient = grad_zipnegloglik_cov(parm, ybody, xbody, M, nbetween.col(i));
            
        }
        else{
            cumsum1 = 0;
            cumsum2 = 0;
            for(j=0; j<n; j++){
                for(t=0; t<nbetween(j,i); t++){
                    ytail(cumsum1+t) = y(cumsum2+pointer(j)+t);
                    xtail.row(cumsum1+t) = x.row(cumsum2+pointer(j)+t);
                }
                cumsum1 += nbetween(j,i);
                cumsum2 += ntimes(j);
            }
            pointer += nbetween.col(i);
            
            //replace the delta in the parm for gradient calculation
            //for(m=0; m<M-1; m++) parm(m) = log(forward(m+1)) - log(forward(0));
            
            //get the forward probability
            //forward = forward_nocov(parm, M, ytail, nbetween.col(i));
            
            neggradient = grad_zipnegloglik_cov(parm, ytail, xtail, M, nbetween.col(i));
            
        }
        
        //coef = (i+1)^(-power)
        coef = exp(-power * log(ifloat+1)) * rate / cumsum1;
        ifloat++;
        parm -= coef % neggradient;
        
        //Rcpp::Rcout<<i<<std::endl;
        //Rcpp::Rcout<<coef<<std::endl;
        //Rcpp::Rcout<<cumsum1<<std::endl;
        result.row(i) = parm.t();
        //result.row(i) = neggradient.t();
        //if(i==nupdate-2) break;
    }
    
    return result;
    
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//now begins the specialized nonparametric hidden semi-markov model
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////

//[[Rcpp::export]]
arma::mat hsmm_gen_np(int dim, int M, arma::vec pi, arma::vec theta, arma::vec zeroprop,
                      arma::mat omega, arma::mat dt){
    
    
    int j,count;
    
    int m,n,prev,curr;
    arma::vec label(M); //for number of states
    
    double u;
    arma::mat result(dim, 2); /*first column for x, second column for state*/
    //int dtrow = dt.n_rows;
    int dtcol = dt.n_cols;
    arma::vec dtlabel(dtcol);
    
    for(m=0;m<M;m++) label(m) = m+1;
    for(m=0;m<dtcol;m++) dtlabel(m) = m + 1;
    ///multinomrand(int n, int k, arma::vec prob, arma::vec label)
    
    //starting generation
    curr = multinomrand(1, M, pi, label)(0);
    
    count = multinomrand(1, dtcol, dt.row(curr-1).t(), dtlabel)(0);
    
    for(j=0;j<count;j++) {
        result(j,1) = curr;
        u = Rcpp::runif(1,0,1)(0);
        if(u<=zeroprop(curr-1)) result(j,0)=0;
        else result(j,0) = Rcpp::rpois(1, theta(curr-1))(0);
    }
    
    n = count;
    
    //iteration
    while(n<dim){
        
        prev = result(n-1,1) - 1;
        curr = multinomrand(1, M, omega.row(prev).t(), label)(0);
        
        count = multinomrand(1, dtcol, dt.row(curr-1).t(), dtlabel)(0);
        
        for(j=0;j<count and n+j<dim; j++) {
            result(n+j,1) = curr;
            u = Rcpp::runif(1,0,1)(0);
            if(u<=zeroprop(curr-1)) result(n+j,0)=0;
            else result(n+j,0) = Rcpp::rpois(1, theta(curr-1))(0);
        }
        
        n += count;
    }
    
    return(result);
    
}


//////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat hsmm_cov_gen_np(arma::vec parm, arma::mat dt, int M, long dim, arma::vec zeroindex,
                          int ncolcovpi, arma::mat covpi,
                          int ncolcovomega, arma::mat covomega, int ncolcovp1, arma::mat covp1,
                          int ncolcovpois, arma::mat covpois){
    
    /*parm =  logit(pi)---M-1, logit(p1), log(theta)---M,
     omega---M*(M-2): (1,3)...(1,n)(2,3)...(2,n)(3,2)(3,4)...(3,n).......*/
    //dt  is the matrix of state duration probabilities
    
    //intercept column is not included in the covariates, but should be in the parm
    
    //parameters are given alternatively in the order of beta0 and beta1
    
    int i,j,k,m,nextindex,prev,curr;
    long count,n;
    double tempsum,u;
    int dtcol = dt.n_cols;
    arma::vec label(M);
    arma::vec dtlabel(dtcol);
    arma::mat result(dim, 2); /*first column for x, second column for state*/
    
    
    arma::vec pi(M);
    arma::vec theta(M);
    arma::vec zeroprop(M);
    arma::mat omega(M,M);
    //arma::mat a(totalmv,totalmv);
    
    //////
    //retrieve some of the parameters
    nextindex = 0;
    //recover pi
    pi(0) = 1;
    tempsum = 1;
    for(m=1; m<M; m++){
        pi(m) = parm(nextindex);
        
        for(j=0;j<ncolcovpi;j++)
        pi(m) += parm(nextindex+j+1)*covpi(0,j) ;
        
        pi(m) = exp(pi(m));
        
        tempsum += pi(m);
        nextindex = nextindex + ncolcovpi + 1;
    }
    
    for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
    
    
    
    //////////////////
    //start generation
    for(m=0;m<M;m++) label(m) = m+1;
    for(m=0;m<dtcol;m++) dtlabel(m) = m + 1;
    
    
    //starting generation
    curr = multinomrand(1, M, pi, label)(0);
    
    count = multinomrand(1, dtcol, dt.row(curr-1).t(), dtlabel)(0);
    
    
    /////////////////
    for(k=0;k<count;k++) {
        
        //get some of the parameters in each iteration
        nextindex = 0;
        
        //recover newtheta,p1, newpi
        pi(0) = 1;
        tempsum = 1;
        for(m=1; m<M; m++){
            pi(m) = parm(nextindex);
            
            for(j=0;j<ncolcovpi;j++)
            pi(m) += parm(nextindex+j+1)*covpi(k,j) ;
            
            pi(m) = exp(pi(m));
            
            tempsum += pi(m);
            nextindex = nextindex + ncolcovpi + 1;
        }
        
        for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
        
        //zeroprops
        
        for(i=0;i<M;i++){
            if(zeroindex(i)==0) zeroprop(i)=0;
            else{
                zeroprop(i) = parm(nextindex);
                
                for(m=0;m<ncolcovp1;m++)
                zeroprop(i) += parm(nextindex+m+1) * covp1(k,m);
                
                zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
                nextindex = nextindex + ncolcovp1 + 1;
            }
        }
        
        
        //theta
        for(m=0; m<M; m++){
            theta(m) = parm(nextindex);
            for(j=0; j<ncolcovpois;j++){
                theta(m) += parm(nextindex+j+1) * covpois(k,j);
            }
            theta(m) = exp(theta(m));
            nextindex = nextindex + ncolcovpois + 1;
        }
        
        
        //recover omega:
        //if M=2 then assign off-diagonal with 1; else
        //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
        
        omega.zeros();
        if(M==2) {omega(0,1)=1; omega(1,0)=1;}
        else{
            for(i=0;i<M;i++){
                tempsum = 1;
                for(j=2;j<M;j++){
                    omega(i,j) = parm(nextindex); //new
                    for(m=0;m<ncolcovomega;m++)
                    omega(i,j) += parm(nextindex+m+1)*covomega(k,m);
                    //omega(i,j) = exp(omega(i,j))/(1+exp(omega(i,j)));
                    omega(i,j) = exp(omega(i,j)); //new
                    tempsum += omega(i,j); //new
                    nextindex = nextindex + ncolcovomega + 1;
                }
                
                for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //new
                if(i==0) omega(i,1) = 1 / tempsum;
                else omega(i,0) = 1 / tempsum;
                //tempsum = 1 - arma::sum(omega.row(i)); //prevent negative values
                //if(i==0) omega(i,1) = MAX(0,tempsum);
                //else omega(i,0) = MAX(0,tempsum);
            }
            
            for(i=2;i<M;i++){
                for(j=2;j<=i;j++){
                    omega(i,j-1) = omega(i,j);
                    if(i==j) omega(i,j)=0;
                }
            }
            
            
        }//end section for zeroinflated poisson distribution
        
        
        
        result(k,1) = curr;
        u = Rcpp::runif(1,0,1)(0);
        if(u<=zeroprop(curr-1)) result(k,0)=0;
        else result(k,0) = Rcpp::rpois(1, theta(curr-1))(0);
    }
    
    n = count;
    
    
    
    //////////////////////////////////////////////////////
    //iteration
    while(n<dim){
        
        //retrieve the parameters
        nextindex = 0;
        
        //recover newtheta,p1, newpi
        pi(0) = 1;
        tempsum = 1;
        for(m=1; m<M; m++){
            pi(m) = parm(nextindex);
            
            for(j=0;j<ncolcovpi;j++)
            pi(m) += parm(nextindex+j+1)*covpi(n,j) ;
            
            pi(m) = exp(pi(m));
            
            tempsum += pi(m);
            nextindex = nextindex + ncolcovpi + 1;
        }
        
        for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
        
        //zeroprops
        
        for(i=0;i<M;i++){
            if(zeroindex(i)==0) zeroprop(i)=0;
            else{
                zeroprop(i) = parm(nextindex);
                
                for(m=0;m<ncolcovp1;m++)
                zeroprop(i) += parm(nextindex+m+1) * covp1(n,m);
                
                zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
                nextindex = nextindex + ncolcovp1 + 1;
            }
        }
        
        
        
        //theta
        for(m=0; m<M; m++){
            theta(m) = parm(nextindex);
            for(j=0; j<ncolcovpois;j++){
                theta(m) += parm(nextindex+j+1) * covpois(n,j);
            }
            theta(m) = exp(theta(m));
            nextindex = nextindex + ncolcovpois + 1;
        }
        
        //then recover omega for this iteration:
        //if M=2 then assign off-diagonal with 1; else
        //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
        
        omega.zeros();
        if(M==2) {omega(0,1)=1; omega(1,0)=1;}
        else{
            for(i=0;i<M;i++){
                tempsum = 1;
                for(j=2;j<M;j++){
                    omega(i,j) = parm(nextindex); //new
                    for(m=0;m<ncolcovomega;m++)
                    omega(i,j) += parm(nextindex+m+1)*covomega(n,m);
                    //omega(i,j) = exp(omega(i,j))/(1+exp(omega(i,j)));
                    omega(i,j) = exp(omega(i,j)); //new
                    tempsum += omega(i,j); //new
                    nextindex = nextindex + ncolcovomega + 1;
                }
                
                for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //new
                if(i==0) omega(i,1) = 1 / tempsum;
                else omega(i,0) = 1 / tempsum;
                //tempsum = 1 - arma::sum(omega.row(i)); //prevent negative values
                //if(i==0) omega(i,1) = MAX(0,tempsum);
                //else omega(i,0) = MAX(0,tempsum);
            }
            
            for(i=2;i<M;i++){
                for(j=2;j<=i;j++){
                    omega(i,j-1) = omega(i,j);
                    if(i==j) omega(i,j)=0;
                }
            }
            
            
            
        }//end section for zeroinflated poisson distribution
        
        prev = result(n-1,1) - 1;
        //sample from the previous omega
        curr = multinomrand(1, M, omega.row(prev).t(), label)(0);
        
        count = multinomrand(1, dtcol, dt.row(curr-1).t(), dtlabel)(0);
        
        
        ////////////////
        for(k=0;k<count and n+k<dim; k++) {
            
            
            //get some of the parameters in each iteration
            nextindex = 0;
            
            //recover newtheta,p1, newpi
            pi(0) = 1;
            tempsum = 1;
            for(m=1; m<M; m++){
                pi(m) = parm(nextindex);
                
                for(j=0;j<ncolcovpi;j++)
                pi(m) += parm(nextindex+j+1)*covpi(n+k,j) ;
                
                pi(m) = exp(pi(m));
                
                tempsum += pi(m);
                nextindex = nextindex + ncolcovpi + 1;
            }
            
            for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
            
            //zeroprops
            
            for(i=0;i<M;i++){
                if(zeroindex(i)==0) zeroprop(i)=0;
                else{
                    zeroprop(i) = parm(nextindex);
                    
                    for(m=0;m<ncolcovp1;m++)
                    zeroprop(i) += parm(nextindex+m+1) * covp1(n+k,m);
                    
                    zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
                    nextindex = nextindex + ncolcovp1 + 1;
                }
            }
            
            
            //theta
            for(m=0; m<M; m++){
                theta(m) = parm(nextindex);
                for(j=0; j<ncolcovpois;j++){
                    theta(m) += parm(nextindex+j+1) * covpois(n+k,j);
                }
                theta(m) = exp(theta(m));
                nextindex = nextindex + ncolcovpois + 1;
            }
            
            
            //recover omega:
            //if M=2 then assign off-diagonal with 1; else
            //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
            
            omega.zeros();
            if(M==2) {omega(0,1)=1; omega(1,0)=1;}
            else{
                for(i=0;i<M;i++){
                    tempsum = 1;
                    for(j=2;j<M;j++){
                        omega(i,j) = parm(nextindex); //new
                        for(m=0;m<ncolcovomega;m++)
                        omega(i,j) += parm(nextindex+m+1)*covomega(n+k,m);
                        //omega(i,j) = exp(omega(i,j))/(1+exp(omega(i,j)));
                        omega(i,j) = exp(omega(i,j)); //new
                        tempsum += omega(i,j); //new
                        nextindex = nextindex + ncolcovomega + 1;
                    }
                    
                    for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //new
                    if(i==0) omega(i,1) = 1 / tempsum;
                    else omega(i,0) = 1 / tempsum;
                    //tempsum = 1 - arma::sum(omega.row(i)); //prevent negative values
                    //if(i==0) omega(i,1) = MAX(0,tempsum);
                    //else omega(i,0) = MAX(0,tempsum);
                }
                
                for(i=2;i<M;i++){
                    for(j=2;j<=i;j++){
                        omega(i,j-1) = omega(i,j);
                        if(i==j) omega(i,j)=0;
                    }
                }
                
            }//end section for zeroinflated poisson distribution
            
            
            result(n+k,1) = curr;
            u = Rcpp::runif(1,0,1)(0);
            if(u<=zeroprop(curr-1)) result(n+k,0)=0;
            else result(n+k,0) = Rcpp::rpois(1, theta(curr-1))(0);
        }
        
        n += count;
    }
    
    
    return(result);
    
}

/////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::vec hsmm_viterbi_np(arma::vec y, int M, arma::vec trunc, arma::vec pi, arma::vec theta,
                          arma::mat omega, arma::mat dm, arma::vec zeroprop){
    
    long dim = y.n_rows;
    double tempsum;
    arma::vec state(dim);
    int i,j;
    //int dmcol = arma::max(trunc); //get the maximum
    int totalmv = arma::sum(trunc);
    arma::vec newtheta(totalmv);
    arma::vec newpi(totalmv);
    arma::vec newzeroprop(totalmv);
    //Rcpp::Rcout<<dmcol<<std::endl;
    
    
    arma::mat a(totalmv,totalmv);
    
    a = hsmm_hmm (omega, dm, trunc);
    
    //recover newtheta,newpi
    
    tempsum = 0;
    for(i=0;i<M;i++){
        for(j=0;j<trunc[i];j++){
            
            newtheta(tempsum+j) = theta(i);
            newzeroprop(tempsum+j) = zeroprop(i);
            newpi(tempsum+j) = pi(i)/trunc(i);
        }
        tempsum += trunc(i);
    }
    
    state = hmm_viterbi(newpi, a, newtheta, totalmv, y, newzeroprop);
    
    for(i=0;i<dim;i++){
        tempsum=0;
        for(j=0;j<M;j++){
            if(state(i)<=tempsum+trunc(j)){
                state(i)=j+1;
                break;
            }
            tempsum += trunc(j);
        }
    }
    
    return state;
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////
/////////////////////////

// [[Rcpp::export]]
arma::vec hsmm_cov_viterbi_np(arma::vec parm, arma::mat dm, int M, arma::vec trunc, arma::vec y,arma::vec zeroindex,
                              int ncolcovpi, arma::mat covpi,
                              int ncolcovomega, arma::mat covomega, int ncolcovp1, arma::mat covp1,
                              int ncolcovpois, arma::mat covpois){
    
    double tempsum,colmax;
    
    long dim = y.n_rows;
    int i,j,m,n,nextindex;
    //int dmcol = arma::max(trunc); //get the maximum
    int totalmv = arma::sum(trunc);
    arma::vec pi(M);
    arma::vec theta(M);
    arma::vec newtheta(totalmv);
    arma::vec newpi(totalmv);
    arma::vec zeroprop(M);
    arma::vec newzeroprop(totalmv);
    //Rcpp::Rcout<<dmcol<<std::endl;
    
    
    arma::mat omega(M,M);
    arma::mat a(totalmv,totalmv);
    
    
    arma::vec meanvec(totalmv);
    arma::vec forward(totalmv);
    arma::vec forwardsum;
    arma::vec state(dim);
    arma::vec tempmax(totalmv);
    arma::mat xi(dim, totalmv);
    arma::mat tempmat(totalmv,totalmv);
    
    //retrieve the full parameters
    nextindex = 0;
    //dwell time
    
    //recover newtheta,p1, newpi
    pi(0) = 1;
    tempsum = 1;
    for(m=1; m<M; m++){
        pi(m) = parm(nextindex);
        
        for(j=0;j<ncolcovpi;j++)
        pi(m) += parm(nextindex+j+1)*covpi(0,j) ;
        
        pi(m) = exp(pi(m));
        
        tempsum += pi(m);
        nextindex = nextindex + ncolcovpi + 1;
    }
    
    for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
    
    //p1
    
    //zero proportions
    for(i=0;i<M;i++){
        if(zeroindex(i)==0) zeroprop(i)=0;
        else{
            zeroprop(i) = parm(nextindex);
            
            for(m=0;m<ncolcovp1;m++)
            zeroprop(i) += parm(nextindex+m+1) * covp1(0,m);
            
            zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
            nextindex = nextindex + ncolcovp1 + 1;
        }
    }
    
    //theta
    for(m=0; m<M; m++){
        theta(m) = parm(nextindex);
        for(j=0; j<ncolcovpois;j++){
            theta(m) += parm(nextindex+j+1) * covpois(0,j);
        }
        theta(m) = exp(theta(m));
        nextindex = nextindex + ncolcovpois + 1;
    }
    
    //Rcpp::Rcout << nextindex << std::endl;
    //get newlogtheta,newpi
    tempsum = 0;
    for(i=0;i<M;i++){
        for(j=0;j<trunc[i];j++){
            newtheta(tempsum+j) = theta(i);
            newpi(tempsum+j) = pi(i)/trunc(i);
            newzeroprop(tempsum+j) = zeroprop(i);
            //Rcpp::Rcout << "theta=" << newtheta(tempsum+j) << "pi=" << newpi(tempsum+j) << std::endl;
        }
        tempsum += trunc(i);
    }
    
    //recover omega:
    //if M=2 then assign off-diagonal with 1; else
    //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
    
    omega.zeros();
    if(M==2) {omega(0,1)=1; omega(1,0)=1;}
    else{
        for(i=0;i<M;i++){
            tempsum = 1;
            for(j=2;j<M;j++){
                omega(i,j) = parm(nextindex); //new
                for(m=0;m<ncolcovomega;m++)
                omega(i,j) += parm(nextindex+m+1)*covomega(0,m);
                //omega(i,j) = exp(omega(i,j))/(1+exp(omega(i,j)));
                omega(i,j) = exp(omega(i,j)); //new
                tempsum += omega(i,j); //new
                nextindex = nextindex + ncolcovomega + 1;
            }
            
            for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //new
            if(i==0) omega(i,1) = 1 / tempsum;
            else omega(i,0) = 1 / tempsum;
            //tempsum = 1 - arma::sum(omega.row(i)); //prevent negative values
            //if(i==0) omega(i,1) = MAX(0,tempsum);
            //else omega(i,0) = MAX(0,tempsum);
        }
        
        for(i=2;i<M;i++){
            for(j=2;j<=i;j++){
                omega(i,j-1) = omega(i,j);
                if(i==j) omega(i,j)=0;
            }
        }
        
        
    }  //end of poisson and zerofl part
    
    
    
    a = hsmm_hmm (omega, dm, trunc);
    
    //Rcpp::Rcout << nextindex << std::endl;
    
    
    //initialize the forward variable
    
    for(m=0; m<totalmv; m++) meanvec(m) = dzip(newzeroprop(m),newtheta(m),y(0),FALSE);
    
    
    forward = newpi % meanvec; //element wise multiplication
    forwardsum = newpi.t() * meanvec;
    
    //Rcpp::Rcout << "loglik=" << loglik << std::endl;
    
    for(m=0; m<totalmv; m++) xi(0,m) = forward(m) / forwardsum(0);
    
    
    
    //recursion
    for(n=1; n<dim; n++){
        
        nextindex = 0;
        
        //recover newtheta,p1, newpi
        pi(0) = 1;
        tempsum = 1;
        for(m=1; m<M; m++){
            pi(m) = parm(nextindex);
            
            for(j=0;j<ncolcovpi;j++)
            pi(m) += parm(nextindex+j+1)*covpi(n,j) ;
            
            pi(m) = exp(pi(m));
            
            tempsum += pi(m);
            nextindex = nextindex + ncolcovpi + 1;
        }
        
        for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
        
        //p1
        for(i=0;i<M;i++){
            if(zeroindex(i)==0) zeroprop(i)=0;
            else{
                zeroprop(i) = parm(nextindex);
                
                for(m=0;m<ncolcovp1;m++)
                zeroprop(i) += parm(nextindex+m+1) * covp1(n,m);
                
                zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
                nextindex = nextindex + ncolcovp1 + 1;
            }
        }
        
        //theta
        for(m=0; m<M; m++){
            theta(m) = parm(nextindex);
            for(j=0; j<ncolcovpois;j++){
                theta(m) += parm(nextindex+j+1) * covpois(n,j);
            }
            theta(m) = exp(theta(m));
            nextindex = nextindex + ncolcovpois + 1;
        }
        
        //Rcpp::Rcout << nextindex << std::endl;
        //get newlogtheta,newpi
        tempsum = 0;
        for(i=0;i<M;i++){
            for(j=0;j<trunc[i];j++){
                newtheta(tempsum+j) = theta(i);
                newpi(tempsum+j) = pi(i)/trunc(i);
                newzeroprop(tempsum+j) = zeroprop(i);
                //Rcpp::Rcout << "theta=" << newtheta(tempsum+j) << "pi=" << newpi(tempsum+j) << std::endl;
            }
            tempsum += trunc(i);
        }
        
        
        //recover omega:
        //if M=2 then assign off-diagonal with 1; else
        //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
        
        omega.zeros();
        if(M==2) {omega(0,1)=1; omega(1,0)=1;}
        else{
            for(i=0;i<M;i++){
                tempsum = 1;
                for(j=2;j<M;j++){
                    omega(i,j) = parm(nextindex); //new
                    for(m=0;m<ncolcovomega;m++)
                    omega(i,j) += parm(nextindex+m+1)*covomega(n,m);
                    //omega(i,j) = exp(omega(i,j))/(1+exp(omega(i,j)));
                    omega(i,j) = exp(omega(i,j)); //new
                    tempsum += omega(i,j); //new
                    nextindex = nextindex + ncolcovomega + 1;
                }
                
                for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //new
                if(i==0) omega(i,1) = 1 / tempsum;
                else omega(i,0) = 1 / tempsum;
                //tempsum = 1 - arma::sum(omega.row(i)); //prevent negative values
                //if(i==0) omega(i,1) = MAX(0,tempsum);
                //else omega(i,0) = MAX(0,tempsum);
            }
            
            
            for(i=2;i<M;i++){
                for(j=2;j<=i;j++){
                    omega(i,j-1) = omega(i,j);
                    if(i==j) omega(i,j)=0;
                }
            }
            
        }
        //end section for zeroinflated poisson distribution
        
        
        //Rcpp::Rcout<<"nextindex="<<nextindex<<std::endl;
        // recover transition matrix
        
        a = hsmm_hmm (omega, dm, trunc);
        
        
        for(m=0; m<totalmv; m++) meanvec(m) = dzip(newzeroprop(m),newtheta(m),y(n),FALSE);
        
        
        
        //difficult part
        for(i=0; i<totalmv; i++){
            for(j=0;j<totalmv;j++){
                tempmat(i,j) = xi(n-1, i) * a(i,j);
            }
        }
        
        for(j=0; j<totalmv; j++){
            for(i=0;i<totalmv;i++){
                if(i==0) colmax = tempmat(i,j);
                else colmax = MAX(colmax, tempmat(i,j));
            }
            tempmax(j) = colmax;
        }
        
        
        forward = tempmax % meanvec;
        forwardsum = tempmax.t() * meanvec;
        for(m=0; m<totalmv; m++) xi(n,m) = forward(m) / forwardsum(0);
    }
    
    //termination
    for(m=0; m<totalmv; m++){
        if(m==0) {
            colmax = xi(dim-1, m);
            state(dim-1) = m;
        }
        else{
            if(xi(dim-1,m) >= colmax) {
                state(dim-1) = m;
                colmax = xi(dim-1,m);
            }
        }
    }
    
    //backtracking
    for(n=dim-2; n>=0; n--){
        
        
        //retrieve the full parameters
        nextindex = 0;
        
        //recover newtheta,p1, newpi
        pi(0) = 1;
        tempsum = 1;
        for(m=1; m<M; m++){
            pi(m) = parm(nextindex);
            
            for(j=0;j<ncolcovpi;j++)
            pi(m) += parm(nextindex+j+1)*covpi(n,j) ;
            
            pi(m) = exp(pi(m));
            
            tempsum += pi(m);
            nextindex = nextindex + ncolcovpi + 1;
        }
        
        for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
        
        //p1
        for(i=0;i<M;i++){
            if(zeroindex(i)==0) zeroprop(i)=0;
            else{
                zeroprop(i) = parm(nextindex);
                
                for(m=0;m<ncolcovp1;m++)
                zeroprop(i) += parm(nextindex+m+1) * covp1(n,m);
                
                zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
                nextindex = nextindex + ncolcovp1 + 1;
            }
        }
        
        //theta
        for(m=0; m<M; m++){
            theta(m) = parm(nextindex);
            for(j=0; j<ncolcovpois;j++){
                theta(m) += parm(nextindex+j+1) * covpois(n,j);
            }
            theta(m) = exp(theta(m));
            nextindex = nextindex + ncolcovpois + 1;
        }
        
        //Rcpp::Rcout << nextindex << std::endl;
        //get newlogtheta,newpi
        tempsum = 0;
        for(i=0;i<M;i++){
            for(j=0;j<trunc[i];j++){
                newtheta(tempsum+j) = theta(i);
                newpi(tempsum+j) = pi(i)/trunc(i);
                newzeroprop(tempsum+j) = zeroprop(i);
                //Rcpp::Rcout << "theta=" << newtheta(tempsum+j) << "pi=" << newpi(tempsum+j) << std::endl;
            }
            tempsum += trunc(i);
        }
        
        
        //recover omega:
        //if M=2 then assign off-diagonal with 1; else
        //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
        
        omega.zeros();
        if(M==2) {omega(0,1)=1; omega(1,0)=1;}
        else{
            for(i=0;i<M;i++){
                tempsum = 1;
                for(j=2;j<M;j++){
                    omega(i,j) = parm(nextindex); //new
                    for(m=0;m<ncolcovomega;m++)
                    omega(i,j) += parm(nextindex+m+1)*covomega(n,m);
                    //omega(i,j) = exp(omega(i,j))/(1+exp(omega(i,j)));
                    omega(i,j) = exp(omega(i,j)); //new
                    tempsum += omega(i,j); //new
                    nextindex = nextindex + ncolcovomega + 1;
                }
                
                for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //new
                if(i==0) omega(i,1) = 1 / tempsum;
                else omega(i,0) = 1 / tempsum;
                //tempsum = 1 - arma::sum(omega.row(i)); //prevent negative values
                //if(i==0) omega(i,1) = MAX(0,tempsum);
                //else omega(i,0) = MAX(0,tempsum);
            }
            
            for(i=2;i<M;i++){
                for(j=2;j<=i;j++){
                    omega(i,j-1) = omega(i,j);
                    if(i==j) omega(i,j)=0;
                }
            }
            
        }
        
        
        //Rcpp::Rcout<<"nextindex="<<nextindex<<std::endl;
        // recover transition matrix
        
        a = hsmm_hmm (omega, dm, trunc);
        
        
        for(m=0;m<totalmv;m++){
            if(m==0) {
                colmax = xi(n,m) * a(m, state(n+1));
                state(n) = m;
            }
            else{
                if(xi(n,m) * a(m, state(n+1))>=colmax) {
                    state(n) = m;
                    colmax = xi(n,m) * a(m, state(n+1));
                }
            }
        }
    }
    
    for(n=0; n<dim; n++) state(n) = state(n) + 1;
    
    for(i=0;i<dim;i++){
        tempsum=0;
        for(j=0;j<M;j++){
            if(state(i)<=tempsum+trunc(j)){
                state(i)=j+1;
                break;
            }
            tempsum += trunc(j);
        }
    }
    
    
    return(state);
}

///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//now begins the (Stochastic) gradient descent for the specicalized hidden semi Markov model
//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double hsmmloglik_nocov(arma::mat dt, arma::vec delta, arma::mat gamma, double theta, arma::vec lambda,
                        arma::vec y, arma::vec ntimes, arma::vec trunc){
    double loglik = 0;
    //long dim = y.n_rows;
    int M = lambda.n_rows;
    int n = ntimes.n_rows;
    int j,m,t;
    double tempsum;
    //int dmcol = arma::max(trunc); //get the maximum
    int totalmv = arma::sum(trunc);
    
    arma::vec newlambda(totalmv);
    arma::vec newpi(totalmv);
    arma::mat a(totalmv,totalmv);
    arma::vec newtheta(totalmv);
    
    arma::vec forward(totalmv);
    arma::vec meanvec(totalmv);
    arma::vec forwardsum;  //then forwardsum(0) is double
    arma::mat tempmat(1,totalmv);
    
    
    tempsum = 0;
    for(m=0; m<M; m++){
        for(j=0; j<trunc(m); j++){
            newlambda(tempsum+j) = lambda(m);
            newpi(tempsum+j) = delta(m) / trunc(m);
            if(m==0){
                newtheta(tempsum+j) = theta;
            }
            else{
                newtheta(tempsum+j) = 0;
            }
        }
        tempsum += trunc(m);
    }
    
    a = hsmm_hmm(gamma, dt, trunc);
    
    int count = 0;
    
    for(j=0; j<n; j++){
        
        tempsum = 0;
        //initialize the forward variable
        for(m=0; m<totalmv; m++) meanvec(m) = dzip(newtheta(m),newlambda(m), y(count), FALSE);
        
        forward = newpi % meanvec; //element wise multiplication
        forwardsum = newpi.t() * meanvec;
        tempsum = log(forwardsum(0));
        //Rcpp::Rcout << "loglik=" << loglik << std::endl;
        
        for(m=0; m<totalmv; m++) forward(m) = forward(m) / forwardsum(0);
        
        //recursion
        for(t=1; t<ntimes(j); t++){
            
            for(m=0; m<totalmv; m++) meanvec(m) = dzip(newtheta(m),newlambda(m), y(count+t), FALSE);
            
            tempmat = forward.t() * a; //row matrix
            forward = tempmat.t() % meanvec;
            forwardsum = tempmat * meanvec;
            tempsum += log(forwardsum(0));
            for(m=0; m<totalmv; m++) forward(m) = forward(m) / forwardsum(0);
            
        }
        
        loglik += tempsum;
        count += ntimes(j);
        
    }
    
    return loglik;
    
}


///////////////////////////////////////////////////////////////////
///////////////////////////////////
// [[Rcpp::export]]
Rcpp::List hsmmretrieve_nocov(arma::vec parm, int M, arma::vec trunc){
    int nextindex = 0;
    arma::vec delta(M);
    arma::mat omega(M,M);
    double theta;
    arma::vec lambda(M);
    int i, j, m;
    double tempsum;
    
    int dtcol = arma::max(trunc); //get the maximum
    arma::mat dt(M, dtcol);
    //int totalmv = arma::sum(trunc);
    //retrieve the full parameters
    
    dt.zeros();
    
    for(i=0;i<M;i++){
        for(j=0;j<trunc(i);j++){
            if(j==0){
                tempsum = 1;
                dt(i,j) = 1;
            }
            else{
                dt(i,j) = exp(parm(nextindex + j - 1));
                tempsum += dt(i,j);
            }
            //automatically 0 when j>trunc(i);
        }
        
        for(m=0;m<trunc(i); m++) dt(i,m) = dt(i,m) / tempsum;
        nextindex += trunc(i)-1;
        //Rcpp::Rcout<<nextindex<<std::endl;
    }
    
    delta(0) = 1;
    tempsum = 1;
    for(m=1; m<M; m++) {
        delta(m) = exp(parm(nextindex+m-1));
        tempsum += delta(m);
    }
    for(m=0; m<M; m++) delta(m) = delta(m) / tempsum;
    
    nextindex += M - 1;
    theta = exp(parm(nextindex)) / (1+exp(parm(nextindex)));
    nextindex ++;
    
    for(m=0; m<M; m++) lambda(m) = exp(parm(nextindex+m));
    
    nextindex += M;
    
    
    //recover omega:   from M*(M-2) [3M ~ M*M+M]
    //if M=2 then assign off-diagonal with 1; else
    //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
    
    omega.zeros();
    if(M==2) {omega(0,1)=1; omega(1,0)=1;}
    else{
        for(i=0;i<M;i++){
            tempsum = 1;
            for(j=2;j<M;j++){
                omega(i,j) = exp(parm(nextindex+(M-2)*i+j-2)); //updated
                tempsum += omega(i,j); //updated
                //omega(i,j) = exp(parm(3*M+(M-2)*i+j-2))/(1+exp(parm(3*M+(M-2)*i+j-2))); //old
            }
            
            for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //updated
            
            if(i==0) omega(i,1) = 1 / tempsum;
            else omega(i,0) = 1 / tempsum;
            //tempsum = 1 - arma::sum(omega.row(i)); //old
            //if(i==0) omega(i,1) = MAX(0,tempsum); //old
            //else omega(i,0) = MAX(0,tempsum); //old
        }
        
        for(i=2;i<M;i++){
            for(j=2;j<=i;j++){
                omega(i,j-1) = omega(i,j);
                if(i==j) omega(i,j)=0;
            }
        }
        
    }
    
    
    return Rcpp::List::create(Rcpp::Named("dt")=dt,
                              Rcpp::Named("delta")=delta,
                              Rcpp::Named("gamma")=omega,
                              Rcpp::Named("theta")=theta,
                              Rcpp::Named("lambda")=lambda);
    
}

////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double hsmmnegloglik_nocov(arma::vec parm, int M, arma::vec y, arma::vec ntimes, arma::vec trunc){
    
    Rcpp::List mod = hsmmretrieve_nocov(parm,M,trunc);
    arma::mat dt = mod("dt");
    arma::vec delta = mod("delta");
    arma::mat gamma = mod("gamma");
    double theta = mod("theta");
    arma::vec lambda = mod("lambda");
    
    double negloglik = - hsmmloglik_nocov(dt, delta, gamma, theta, lambda,
                                          y, ntimes, trunc);
    
    return negloglik;
}

/////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::vec gradhsmmloglik_nocov(arma::mat dt, arma::vec delta, arma::mat gamma, double theta, arma::vec lambda,
                               arma::vec y, arma::vec ntimes, arma::vec trunc){
    //double loglik = 0;
    long dim = y.n_rows;
    int M = lambda.n_rows;
    int n = ntimes.n_rows;
    int i,j,m,t;
    double tempsum;
    int dmcol = arma::max(trunc); //get the maximum
    int totalmv = arma::sum(trunc);
    
    arma::vec newlambda(totalmv);
    arma::vec newpi(totalmv);
    arma::mat a(totalmv,totalmv);
    arma::vec newtheta(totalmv);
    
    arma::mat nodeprob(dim, totalmv);
    arma::vec cumld(M);
    
    tempsum = 0;
    for(m=0; m<M; m++){
        for(j=0; j<trunc(m); j++){
            newlambda(tempsum+j) = lambda(m);
            newpi(tempsum+j) = delta(m) / trunc(m);
            if(m==0){
                newtheta(tempsum+j) = theta;
            }
            else{
                newtheta(tempsum+j) = 0;
            }
        }
        cumld(m) = tempsum;
        tempsum += trunc(m);
        
    }
    
    a = hsmm_hmm(gamma, dt, trunc);
    
    for(t=0; t<dim; t++)
    for(m=0; m<totalmv; m++)
    nodeprob(t,m) = dzip(newtheta(m),newlambda(m), y(t), FALSE);
    
    Rcpp::List mod = forwardbackward(newpi, a, nodeprob, y, ntimes);
    
    //arma::vec newgammasum = mod("colsumgamma");
    arma::mat newgamma = mod("Gamma");
    arma::mat newcolsumxi = mod("colsumxi");
    
    //get back old Gamma
    arma::mat oldgamma(dim,M);
    oldgamma.zeros();
    
    
    for(t=0; t<dim; t++){
        tempsum = 0;
        for(m=0; m<M; m++){
            for(j=0; j<trunc(m); j++){
                oldgamma(t,m) += newgamma(t,tempsum+j);
            }
            tempsum += trunc(m);
        }
    }
    
    //get back old tpm
    arma::mat oldcolsumxi(M,M);
    oldcolsumxi.zeros();
    arma::mat nondiag(totalmv, M-1);
    nondiag.zeros();
    arma::vec tempindex(M-1); //record other states
    arma::vec label(M);
    for(m=0;m<M;m++) label(m) = m;
    
    tempsum = 0;
    for(m=0; m<M; m++){
        
        if(m==0){
            tempindex = label.subvec(1,M-1);
        }
        else if(m==M-1){
            tempindex = label.subvec(0,M-2);
        }
        else{
            tempindex.subvec(0,m-1) = label.subvec(0,m-1);
            tempindex.subvec(m,M-2) = label.subvec(m+1,M-1);
        }
        
        //Rcpp::Rcout<<tempindex<<std::endl;
        //Rcpp::Rcout<<cumld<<std::endl;
        for(j=0; j<trunc(m); j++){
            for(i=0; i<M-1; i++){
                nondiag(tempsum+j,i) = newcolsumxi(tempsum+j, cumld( tempindex(i) ));
                oldcolsumxi(m,tempindex(i)) += nondiag(tempsum+j, i);
            }
        }
        tempsum += trunc(m);
    }
    
    
    //get back the state durations
    arma::vec tempeta = rowsum(nondiag);
    arma::mat eta(M,dmcol);
    arma::vec cumsums(M);
    
    eta.zeros();
    
    tempsum = 0;
    for(m=0;m<M;m++){
        cumsums(m) = 0;
        for(j=0;j<trunc(m);j++){
            eta(m,j) = tempeta(tempsum+j);
            cumsums(m) += eta(m,j);
        }
        tempsum += trunc(m);
    }
    //Rcpp::Rcout<<cumsums<<std::endl;
    
    ///////////////////////////////////////////////////////////////////////
    //now compute the derivatives
    //for state duration
    arma::mat dlogdur(M,dmcol);
    for(m=0;m<M;m++){
        for(j=0;j<trunc(m);j++){
            dlogdur(m,j) = eta(m,j) - cumsums(m)*dt(m,j);
        }
    }
    
    //for glogit prior
    arma::vec temppi(M);
    arma::vec dlogitpi(M);
    
    for(m=0; m<M; m++) temppi(m) = 0;
    int count = 0;
    for(j=0; j<n; j++){
        for(m=0; m<M; m++) temppi(m) += oldgamma(count,m);
        count += ntimes(j);
    }
    
    for(m=0; m<M; m++) {
        temppi(m) = temppi(m) / n;
        dlogitpi(m) = temppi(m) - delta(m);
    }
    
    //for glogit tpm
    
    arma::mat dlogtrans(M,M);
    cumsums = rowsum(oldcolsumxi);
    for(m=0;m<M;m++){
        for(j=0;j<M;j++){
            dlogtrans(m,j) = oldcolsumxi(m,j) - cumsums(m)*gamma(m,j);
        }
    }
    
    //theta and lambda
    double dlogtheta = 0;
    arma::vec dloglambda(M);
    for(m=0; m<M; m++) dloglambda(m) = 0;
    
    for(t=0; t<dim; t++){
        
        if(y(t)>0){
            dlogtheta -= theta * oldgamma(t,0);
            dloglambda(0) += oldgamma(t,0) * (y(t) - lambda(0));
        }
        else{
            dlogtheta -= oldgamma(t,0)*theta*(1-theta)*(exp(-lambda(0))-1) / (theta+(1-theta)*exp(-lambda(0)));
            dloglambda(0) += oldgamma(t,0)*(theta-1)*lambda(0)*exp(-lambda(0)) / (theta+(1-theta)*exp(-lambda(0)));
        }
        
        for(m=1; m<M; m++){
            dloglambda(m) += oldgamma(t,m)*(y(t) - lambda(m));
        }
    }
    
    arma::vec result(totalmv-M+M-1+1+M+M*(M-2));
    
    int nextindex = 0;
    for(m=0;m<M;m++){
        for(j=1;j<trunc(m);j++){
            result(nextindex+j-1) = dlogdur(m,j);
        }
        nextindex += trunc(m)-1;
    }
    
    
    for(j=0; j<M-1; j++) result(nextindex+j) = dlogitpi(j+1);
    nextindex += M - 1;
    
    result(nextindex) = dlogtheta;
    nextindex += 1;
    for(i=0; i<M; i++) result(nextindex+i) = dloglambda(i);
    nextindex += M;
    
    //if M=2, gamma = matrix(0,1,1,0); no derivatives for tpm is needed
    if(M>2){
        for(i=0; i<M; i++){
            if(i==0){
                tempindex = label.subvec(1,M-1);
            }
            else if(i==M-1){
                tempindex = label.subvec(0,M-2);
            }
            else{
                tempindex.subvec(0,i-1) = label.subvec(0,i-1);
                tempindex.subvec(i,M-2) = label.subvec(i+1,M-1);
            }
            
            for(j=1; j<M-1; j++){
                result(nextindex + i*(M-2) + j - 1) = dlogtrans(i,tempindex(j));
            }
        }
    }
    
    
    return result;
    
    
    /*
     return Rcpp::List::create(Rcpp::Named("dlogdur")=dlogdur,
     Rcpp::Named("dlogitpi")=dlogitpi,
     Rcpp::Named("dlogtheta")=dlogtheta,
     Rcpp::Named("dloglambda")=dloglambda,
     Rcpp::Named("dlogtrans")=dlogtrans);
     */
}

/////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::vec gradhsmmnegloglik_nocov(arma::vec parm, int M, arma::vec y, arma::vec ntimes, arma::vec trunc){
    
    Rcpp::List mod = hsmmretrieve_nocov(parm,M,trunc);
    arma::mat dt = mod("dt");
    arma::vec delta = mod("delta");
    arma::mat gamma = mod("gamma");
    double theta = mod("theta");
    arma::vec lambda = mod("lambda");
    
    arma::vec result = - gradhsmmloglik_nocov(dt, delta, gamma, theta, lambda,
                                              y, ntimes, trunc);
    
    return result;
}

/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double hsmmloglik_cov(arma::mat dt, arma::vec delta, arma::mat gamma, arma::vec thetaparm, arma::mat lambdaparm,
                      arma::vec y, arma::mat x, arma::vec ntimes, arma::vec trunc){
    double loglik = 0;
    long dim = y.n_rows;
    int M = dt.n_rows;
    int n = ntimes.n_rows;
    int j,m,t;
    double tempsum;
    //int dmcol = arma::max(trunc); //get the maximum
    int totalmv = arma::sum(trunc);
    
    arma::vec lambda(M);
    double theta;
    arma::vec newlambda(totalmv);
    arma::vec newpi(totalmv);
    arma::mat a(totalmv,totalmv);
    arma::vec newtheta(totalmv);
    
    arma::vec forward(totalmv);
    arma::vec meanvec(totalmv);
    arma::vec forwardsum;  //then forwardsum(0) is double
    arma::mat tempmat(1,totalmv);
    arma::vec tempval;
    arma::mat nodeprob(dim,totalmv);
    
    for(t=0; t<dim; t++){
        tempval = x.row(t) * thetaparm;
        theta = exp(tempval(0)) / (1+exp(tempval(0)));
        for(m=0; m<M; m++){
            tempval = x.row(t) * lambdaparm.row(m).t();
            lambda(m) = exp(tempval(0));
        }
        
        tempsum = 0;
        for(m=0; m<M; m++){
            for(j=0; j<trunc(m); j++){
                newlambda(tempsum+j) = lambda(m);
                newpi(tempsum+j) = delta(m) / trunc(m);
                if(m==0){
                    newtheta(tempsum+j) = theta;
                }
                else{
                    newtheta(tempsum+j) = 0;
                }
            }
            tempsum += trunc(m);
        }
        
        for(m=0; m<totalmv; m++) nodeprob(t,m) = dzip(newtheta(m),newlambda(m), y(t), FALSE);
        
    }
    
    
    a = hsmm_hmm(gamma, dt, trunc);
    
    int count = 0;
    
    for(j=0; j<n; j++){
        
        tempsum = 0;
        //initialize the forward variable
        forward = newpi % nodeprob.row(count).t(); //element wise multiplication
        forwardsum = newpi.t() * nodeprob.row(count).t();
        tempsum = log(forwardsum(0));
        //Rcpp::Rcout << "loglik=" << loglik << std::endl;
        
        for(m=0; m<totalmv; m++) forward(m) = forward(m) / forwardsum(0);
        
        //recursion
        for(t=1; t<ntimes(j); t++){
            
            tempmat = forward.t() * a; //row matrix
            forward = tempmat.t() % nodeprob.row(count+t).t();
            forwardsum = tempmat * nodeprob.row(count+t).t();
            tempsum += log(forwardsum(0));
            for(m=0; m<totalmv; m++) forward(m) = forward(m) / forwardsum(0);
            
        }
        
        loglik += tempsum;
        count += ntimes(j);
        
    }
    
    return loglik;
    
}

///////////////////////////////////////////////////////////////////
///////////////////////////////////
// [[Rcpp::export]]
Rcpp::List hsmmretrieve_cov(arma::vec parm, int M, int ncolx, arma::vec trunc){
    int nextindex = 0;
    arma::vec delta(M);
    arma::mat omega(M,M);
    arma::vec thetaparm(ncolx);
    arma::mat lambdaparm(M,ncolx);
    int i, j, m;
    double tempsum;
    
    int dtcol = arma::max(trunc); //get the maximum
    arma::mat dt(M, dtcol);
    
    //retrieve the full parameters
    
    dt.zeros();
    
    for(i=0;i<M;i++){
        for(j=0;j<trunc(i);j++){
            if(j==0){
                tempsum = 1;
                dt(i,j) = 1;
            }
            else{
                dt(i,j) = exp(parm(nextindex + j - 1));
                tempsum += dt(i,j);
            }
            //automatically 0 when j>trunc(i);
        }
        
        for(m=0;m<trunc(i); m++) dt(i,m) = dt(i,m) / tempsum;
        nextindex += trunc(i)-1;
        //Rcpp::Rcout<<nextindex<<std::endl;
    }
    
    
    delta(0) = 1;
    tempsum = 1;
    for(m=1; m<M; m++) {
        delta(m) = exp(parm(nextindex+m-1));
        tempsum += delta(m);
    }
    for(m=0; m<M; m++) delta(m) = delta(m) / tempsum;
    
    
    
    nextindex += M - 1;
    for(i=0; i<ncolx; i++) thetaparm(i) = parm(nextindex+i);
    nextindex += ncolx;
    
    for(i=0; i<M; i++){
        for(j=0; j<ncolx; j++){
            lambdaparm(i,j) = parm(nextindex + i*ncolx + j);
        }
    }
    nextindex += M*(ncolx);
    
    
    //recover omega:   from M*(M-2) [3M ~ M*M+M]
    //if M=2 then assign off-diagonal with 1; else
    //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
    
    omega.zeros();
    if(M==2) {omega(0,1)=1; omega(1,0)=1;}
    else{
        for(i=0;i<M;i++){
            tempsum = 1;
            for(j=2;j<M;j++){
                omega(i,j) = exp(parm(nextindex+(M-2)*i+j-2)); //updated
                tempsum += omega(i,j); //updated
                //omega(i,j) = exp(parm(3*M+(M-2)*i+j-2))/(1+exp(parm(3*M+(M-2)*i+j-2))); //old
            }
            
            for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //updated
            
            if(i==0) omega(i,1) = 1 / tempsum;
            else omega(i,0) = 1 / tempsum;
            //tempsum = 1 - arma::sum(omega.row(i)); //old
            //if(i==0) omega(i,1) = MAX(0,tempsum); //old
            //else omega(i,0) = MAX(0,tempsum); //old
        }
        
        for(i=2;i<M;i++){
            for(j=2;j<=i;j++){
                omega(i,j-1) = omega(i,j);
                if(i==j) omega(i,j)=0;
            }
        }
        
    }
    
    
    return Rcpp::List::create(Rcpp::Named("dt")=dt,
                              Rcpp::Named("delta")=delta,
                              Rcpp::Named("gamma")=omega,
                              Rcpp::Named("thetaparm")=thetaparm,
                              Rcpp::Named("lambdaparm")=lambdaparm);
    
}

//////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double hsmmnegloglik_cov(arma::vec parm, arma::vec y, arma::mat covariates, int M, arma::vec ntimes, arma::vec trunc){
    int ncolx = covariates.n_cols;
    Rcpp::List mod = hsmmretrieve_cov(parm,M,ncolx,trunc);
    arma::mat dt = mod("dt");
    arma::vec delta = mod("delta");
    arma::mat gamma = mod("gamma");
    arma::vec thetaparm = mod("thetaparm");
    arma::mat lambdaparm = mod("lambdaparm");
    
    double negloglik = - hsmmloglik_cov(dt, delta, gamma, thetaparm, lambdaparm, y, covariates, ntimes, trunc);
    return negloglik;
}


/////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::vec gradhsmmloglik_cov(arma::mat dt, arma::vec delta, arma::mat gamma, arma::vec thetaparm, arma::mat lambdaparm,
                             arma::vec y, arma::mat x, arma::vec ntimes, arma::vec trunc){
    //double loglik = 0;
    long dim = y.n_rows;
    int ncolx = x.n_cols;
    int M = dt.n_rows;
    int n = ntimes.n_rows;
    int i,j,m,t;
    double tempsum;
    int dmcol = arma::max(trunc); //get the maximum
    int totalmv = arma::sum(trunc);
    
    double theta;
    arma::vec lambda(M);
    arma::vec newlambda(totalmv);
    arma::vec newpi(totalmv);
    arma::mat a(totalmv,totalmv);
    arma::vec newtheta(totalmv);
    
    arma::mat nodeprob(dim, totalmv);
    arma::vec cumld(M);
    arma::vec tempval;
    
    tempsum = 0;
    for(m=0; m<M; m++){
        cumld(m) = tempsum;
        tempsum += trunc(m);
    }
    
    for(t=0; t<dim; t++){
        tempval = x.row(t) * thetaparm;
        theta = exp(tempval(0)) / (1+exp(tempval(0)));
        for(m=0; m<M; m++){
            tempval = x.row(t) * lambdaparm.row(m).t();
            lambda(m) = exp(tempval(0));
        }
        
        tempsum = 0;
        for(m=0; m<M; m++){
            for(j=0; j<trunc(m); j++){
                newlambda(tempsum+j) = lambda(m);
                newpi(tempsum+j) = delta(m) / trunc(m);
                if(m==0){
                    newtheta(tempsum+j) = theta;
                }
                else{
                    newtheta(tempsum+j) = 0;
                }
            }
            tempsum += trunc(m);
        }
        
        for(m=0; m<totalmv; m++) nodeprob(t,m) = dzip(newtheta(m),newlambda(m), y(t), FALSE);
    }
    
    a = hsmm_hmm(gamma, dt, trunc);
    
    
    Rcpp::List mod = forwardbackward(newpi, a, nodeprob, y, ntimes);
    
    //arma::vec newgammasum = mod("colsumgamma");
    arma::mat newgamma = mod("Gamma");
    arma::mat newcolsumxi = mod("colsumxi");
    
    //get back old Gamma
    arma::mat oldgamma(dim,M);
    oldgamma.zeros();
    
    
    for(t=0; t<dim; t++){
        tempsum = 0;
        for(m=0; m<M; m++){
            for(j=0; j<trunc(m); j++){
                oldgamma(t,m) += newgamma(t,tempsum+j);
            }
            tempsum += trunc(m);
        }
    }
    
    //get back old tpm
    arma::mat oldcolsumxi(M,M);
    oldcolsumxi.zeros();
    arma::mat nondiag(totalmv, M-1);
    nondiag.zeros();
    arma::vec tempindex(M-1); //record other states
    arma::vec label(M);
    for(m=0;m<M;m++) label(m) = m;
    
    tempsum = 0;
    for(m=0; m<M; m++){
        
        if(m==0){
            tempindex = label.subvec(1,M-1);
        }
        else if(m==M-1){
            tempindex = label.subvec(0,M-2);
        }
        else{
            tempindex.subvec(0,m-1) = label.subvec(0,m-1);
            tempindex.subvec(m,M-2) = label.subvec(m+1,M-1);
        }
        
        //Rcpp::Rcout<<tempindex<<std::endl;
        //Rcpp::Rcout<<cumld<<std::endl;
        for(j=0; j<trunc(m); j++){
            for(i=0; i<M-1; i++){
                nondiag(tempsum+j,i) = newcolsumxi(tempsum+j, cumld( tempindex(i) ));
                oldcolsumxi(m,tempindex(i)) += nondiag(tempsum+j, i);
            }
        }
        tempsum += trunc(m);
    }
    
    
    //get back the state durations
    arma::vec tempeta = rowsum(nondiag);
    arma::mat eta(M,dmcol);
    arma::vec cumsums(M);
    
    eta.zeros();
    
    tempsum = 0;
    for(m=0;m<M;m++){
        cumsums(m) = 0;
        for(j=0;j<trunc(m);j++){
            eta(m,j) = tempeta(tempsum+j);
            cumsums(m) += eta(m,j);
        }
        tempsum += trunc(m);
    }
    //Rcpp::Rcout<<cumsums<<std::endl;
    
    ///////////////////////////////////////////////////////////////////////
    //now compute the derivatives
    //for state duration
    arma::mat dlogdur(M,dmcol);
    for(m=0;m<M;m++){
        for(j=0;j<trunc(m);j++){
            dlogdur(m,j) = eta(m,j) - cumsums(m)*dt(m,j);
        }
    }
    
    //for glogit prior
    arma::vec temppi(M);
    arma::vec dlogitpi(M);
    
    for(m=0; m<M; m++) temppi(m) = 0;
    int count = 0;
    for(j=0; j<n; j++){
        for(m=0; m<M; m++) temppi(m) += oldgamma(count,m);
        count += ntimes(j);
    }
    
    for(m=0; m<M; m++) {
        temppi(m) = temppi(m) / n;
        dlogitpi(m) = temppi(m) - delta(m);
    }
    
    //for glogit tpm
    
    arma::mat dlogtrans(M,M);
    cumsums = rowsum(oldcolsumxi);
    for(m=0;m<M;m++){
        for(j=0;j<M;j++){
            dlogtrans(m,j) = oldcolsumxi(m,j) - cumsums(m)*gamma(m,j);
        }
    }
    
    
    //theta and lambda
    arma::vec dlogtheta(ncolx);
    arma::mat dloglambda(M,ncolx);
    for(m=0; m<ncolx; m++) dlogtheta(m) = 0;
    for(i=0; i<M; i++)
    for(j=0; j<ncolx; j++)
    dloglambda(i,j) = 0;
    
    for(t=0; t<dim; t++){
        
        tempval = x.row(t) * thetaparm;
        theta = exp(tempval(0))/(1+exp(tempval(0)));
        for(m=0;m<M;m++) {
            tempval = x.row(t) * lambdaparm.row(m).t();
            lambda(m) = exp(tempval(0));
        }
        
        if(y(t)>0){
            for(i=0;i<ncolx;i++){
                dlogtheta(i) -= theta * oldgamma(t,0) * x(t,i) ;
                dloglambda(0,i) += oldgamma(t,0) * x(t,i) * (y(t) - lambda(0));
            }
        }
        else{
            for(i=0;i<ncolx;i++){
                dlogtheta(i) -= oldgamma(t,0)*x(t,i)*theta*(1-theta)*(exp(-lambda(0))-1) / (theta+(1-theta)*exp(-lambda(0)));
                dloglambda(0,i) += oldgamma(t,0)*x(t,i)*(theta-1)*lambda(0)*exp(-lambda(0)) / (theta+(1-theta)*exp(-lambda(0)));
            }
        }
        
        for(m=1; m<M; m++){
            for(i=0;i<ncolx;i++){
                dloglambda(m,i) += oldgamma(t,m)*x(t,i)*(y(t) - lambda(m));
            }
        }
    }
    
    
    /////////////////////
    arma::vec result(totalmv-M+M-1+ncolx*(1+M)+M*(M-2));
    
    int nextindex = 0;
    for(m=0;m<M;m++){
        for(j=1;j<trunc(m);j++){
            result(nextindex+j-1) = dlogdur(m,j);
        }
        nextindex += trunc(m)-1;
    }
    
    
    for(j=0; j<M-1; j++) result(nextindex+j) = dlogitpi(j+1);
    nextindex += M - 1;
    
    for(i=0; i<ncolx; i++) result(nextindex+i) = dlogtheta(i);
    nextindex += ncolx;
    
    for(i=0; i<M; i++)
    for(j=0; j<ncolx; j++)
    result(nextindex+i*ncolx+j) = dloglambda(i,j);
    nextindex += ncolx*M;
    
    //if M=2, gamma = matrix(0,1,1,0); no derivatives for tpm is needed
    if(M>2){
        for(i=0; i<M; i++){
            if(i==0){
                tempindex = label.subvec(1,M-1);
            }
            else if(i==M-1){
                tempindex = label.subvec(0,M-2);
            }
            else{
                tempindex.subvec(0,i-1) = label.subvec(0,i-1);
                tempindex.subvec(i,M-2) = label.subvec(i+1,M-1);
            }
            
            for(j=1; j<M-1; j++){
                result(nextindex + i*(M-2) + j - 1) = dlogtrans(i,tempindex(j));
            }
        }
    }
    
    //Rcpp::Rcout<<oldcolsumxi<<std::endl;
    return result;
    
    /*
     
     return Rcpp::List::create(Rcpp::Named("dlogdur")=dlogdur,
     Rcpp::Named("dlogitpi")=dlogitpi,
     Rcpp::Named("dlogtheta")=dlogtheta,
     Rcpp::Named("dloglambda")=dloglambda,
     Rcpp::Named("dlogtrans")=dlogtrans);
     */
}

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::vec gradhsmmnegloglik_cov(arma::vec parm, arma::vec y, arma::mat covariates, int M, arma::vec ntimes, arma::vec trunc){
    int ncolx = covariates.n_cols;
    int totalmv = arma::sum(trunc);
    Rcpp::List mod = hsmmretrieve_cov(parm,M,ncolx,trunc);
    arma::mat dt = mod("dt");
    arma::vec delta = mod("delta");
    arma::mat gamma = mod("gamma");
    arma::vec thetaparm = mod("thetaparm");
    arma::mat lambdaparm = mod("lambdaparm");
    
    
    arma::vec result = gradhsmmloglik_cov(dt, delta, gamma, thetaparm, lambdaparm, y, covariates, ntimes, trunc);
    int m;
    for(m=0; m<totalmv-M+M-1+ncolx*(1+M)+M*(M-2); m++) result(m) = -result(m);
    
    
    return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//stochastic gradient descent for hsmm_nocov
////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat sgd_hsmm_nocov(arma::vec parm, arma::vec y, int M, arma::vec ntimes,
                         int nmin, int nupdate, double power, arma::vec rate, arma::vec trunc){
    
    long dim = y.n_rows;
    int n = ntimes.n_rows;
    int i,j,t;
    //double coef;
    arma::vec coef;
    int cumsum1, cumsum2;
    int totalmv = arma::sum(trunc);
    double ifloat = 0;
    arma::mat nbetween = nupdate_to_nbetween(nupdate, nmin, dim, ntimes);
    //only 3 different lengths of y
    arma::vec length = colsum(nbetween);
    arma::vec yhead(length(0));
    arma::vec ybody(length(2));
    arma::vec ytail(length(nupdate-1));
    
    arma::vec pointer(n);
    for(i=0; i<n; i++) pointer(i) = 0;
    
    arma::vec neggradient(totalmv-M+M-1+M*(M-2)+1+M);
    //y.subvec(cumtimes, cumtimes + ntimes[i] - 1)
    
    arma::mat result(nupdate, totalmv-M+M-1+M*(M-2)+1+M);
    
    
    //stochastic gradient descent
    for(i=0; i<nupdate; i++){
        
        if(i==0){
            //assign y
            cumsum1 = 0;
            cumsum2 = 0;
            for(j=0; j<n; j++){
                for(t=0; t<nbetween(j,i); t++){
                    yhead(cumsum1+t) = y(cumsum2+t);
                }
                cumsum1 += nbetween(j,i);
                cumsum2 += ntimes(j);
            }
            pointer += nbetween.col(i);
            
            //get the forward probability
            //forward = forward_nocov(parm, M, yhead, nbetween.col(i));
            
            //get the gradient
            neggradient = gradhsmmnegloglik_nocov(parm, M, yhead, nbetween.col(i), trunc);
            
        }
        else if(i<nupdate - 1){
            cumsum1 = 0;
            cumsum2 = 0;
            for(j=0; j<n; j++){
                for(t=0; t<nbetween(j,i); t++){
                    ybody(cumsum1+t) = y(cumsum2+pointer(j)+t);
                }
                cumsum1 += nbetween(j,i);
                cumsum2 += ntimes(j);
            }
            pointer += nbetween.col(i);
            
            //replace the delta in the parm for gradient calculation
            //for(m=0; m<M-1; m++) parm(m) = log(forward(m+1)) - log(forward(0));
            
            //get the forward probability
            //forward = forward_nocov(parm, M, ybody, nbetween.col(i));
            
            //get the gradient
            neggradient = gradhsmmnegloglik_nocov(parm, M, ybody, nbetween.col(i), trunc);
            
        }
        else{
            cumsum1 = 0;
            cumsum2 = 0;
            for(j=0; j<n; j++){
                for(t=0; t<nbetween(j,i); t++){
                    ytail(cumsum1+t) = y(cumsum2+pointer(j)+t);
                }
                cumsum1 += nbetween(j,i);
                cumsum2 += ntimes(j);
            }
            pointer += nbetween.col(i);
            
            //replace the delta in the parm for gradient calculation
            //for(m=0; m<M-1; m++) parm(m) = log(forward(m+1)) - log(forward(0));
            
            //get the forward probability
            //forward = forward_nocov(parm, M, ytail, nbetween.col(i));
            
            neggradient = gradhsmmnegloglik_nocov(parm, M, ytail, nbetween.col(i), trunc);
            
        }
        
        //coef = (i+1)^(-power)
        coef = exp(-power * log(ifloat+1)) * rate / cumsum1;
        ifloat++;
        parm -= coef % neggradient;
        //Rcpp::Rcout<<i<<std::endl;
        //Rcpp::Rcout<<coef<<std::endl;
        //Rcpp::Rcout<<cumsum1<<std::endl;
        result.row(i) = parm.t();
        //result.row(i) = neggradient.t();
    }
    
    return result;
    
}


////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//stochastic gradient descent for zip_cov
////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat sgd_hsmm_cov(arma::vec parm, arma::vec y, arma::mat x, int M, arma::vec ntimes,
                       int nmin, int nupdate, double power, arma::vec rate, arma::vec trunc){
    
    long dim = y.n_rows;
    int n = ntimes.n_rows;
    int ncolx = x.n_cols;
    //Rcpp::List mod = retrieve_cov(parm,M,ncolx);
    int i,j,t;
    //double coef;
    arma::vec coef;
    int cumsum1, cumsum2;
    int totalmv = arma::sum(trunc);
    double ifloat = 0;
    
    arma::mat nbetween = nupdate_to_nbetween(nupdate, nmin, dim, ntimes);
    //only 3 different lengths of y
    arma::vec length = colsum(nbetween);
    arma::vec yhead(length(0));
    arma::vec ybody(length(2));
    arma::vec ytail(length(nupdate-1));
    arma::mat xhead(length(0), ncolx);
    arma::mat xbody(length(2), ncolx);
    arma::mat xtail(length(nupdate-1), ncolx);
    
    arma::vec pointer(n);
    for(i=0; i<n; i++) pointer(i) = 0;
    
    arma::vec neggradient(totalmv-M+M-1+M*(M-2)+ncolx*(1+M));
    //y.subvec(cumtimes, cumtimes + ntimes[i] - 1)
    
    arma::mat result(nupdate, totalmv-M+M-1+M*(M-2)+ncolx*(1+M));
    
    
    //stochastic gradient descent
    for(i=0; i<nupdate; i++){
        
        if(i==0){
            //assign y
            cumsum1 = 0;
            cumsum2 = 0;
            for(j=0; j<n; j++){
                for(t=0; t<nbetween(j,i); t++){
                    yhead(cumsum1+t) = y(cumsum2+t);
                    xhead.row(cumsum1+t) = x.row(cumsum2+t);
                }
                cumsum1 += nbetween(j,i);
                cumsum2 += ntimes(j);
            }
            pointer += nbetween.col(i);
            
            //get the forward probability
            //forward = forward_nocov(parm, M, yhead, nbetween.col(i));
            
            //get the gradient
            neggradient = gradhsmmnegloglik_cov(parm, yhead, xhead, M, nbetween.col(i),trunc);
            
        }
        else if(i<nupdate - 1){
            cumsum1 = 0;
            cumsum2 = 0;
            for(j=0; j<n; j++){
                for(t=0; t<nbetween(j,i); t++){
                    ybody(cumsum1+t) = y(cumsum2+pointer(j)+t);
                    xbody.row(cumsum1+t) = x.row(cumsum2+pointer(j)+t);
                }
                cumsum1 += nbetween(j,i);
                cumsum2 += ntimes(j);
            }
            pointer += nbetween.col(i);
            
            //replace the delta in the parm for gradient calculation
            //for(m=0; m<M-1; m++) parm(m) = log(forward(m+1)) - log(forward(0));
            
            //get the forward probability
            //forward = forward_nocov(parm, M, ybody, nbetween.col(i));
            
            //get the gradient
            neggradient = gradhsmmnegloglik_cov(parm, ybody, xbody, M, nbetween.col(i),trunc);
            
        }
        else{
            cumsum1 = 0;
            cumsum2 = 0;
            for(j=0; j<n; j++){
                for(t=0; t<nbetween(j,i); t++){
                    ytail(cumsum1+t) = y(cumsum2+pointer(j)+t);
                    xtail.row(cumsum1+t) = x.row(cumsum2+pointer(j)+t);
                }
                cumsum1 += nbetween(j,i);
                cumsum2 += ntimes(j);
            }
            pointer += nbetween.col(i);
            
            //replace the delta in the parm for gradient calculation
            //for(m=0; m<M-1; m++) parm(m) = log(forward(m+1)) - log(forward(0));
            
            //get the forward probability
            //forward = forward_nocov(parm, M, ytail, nbetween.col(i));
            
            neggradient = gradhsmmnegloglik_cov(parm, ytail, xtail, M, nbetween.col(i),trunc);
            
        }
        
        //coef = (i+1)^(-power)
        coef = exp(-power * log(ifloat+1)) * rate / cumsum1;
        ifloat++;
        parm -= coef % neggradient;
        
        //Rcpp::Rcout<<i<<std::endl;
        //Rcpp::Rcout<<coef<<std::endl;
        //Rcpp::Rcout<<cumsum1<<std::endl;
        result.row(i) = parm.t();
        //result.row(i) = neggradient.t();
        //if(i==nupdate-2) break;
    }
    
    return result;
    
}

//[[Rcpp::export]]
arma::cube getallexpm(arma::mat tpm, arma::vec udiff){
    int dim = udiff.n_rows;
    int M = tpm.n_rows;
    arma::cube result(M,M,dim);
    int i;
    for(i=0; i<dim; i++) result.slice(i) = matrixexp(tpm, udiff(i));
    return result;
}

//[[Rcpp::export]]
int locate(double interval, arma::vec udiff){
    int i, result=0;
    int dim = udiff.n_rows;
    for(i=0; i<dim; i++){
        if(interval==udiff(i)) {
            result=i;
            break;
        }
    }
    return result;
}


/////////////
// [[Rcpp::export]]
arma::mat get_c(arma::mat xi, arma::vec diffvec, arma::vec udiff){
    int ncol = xi.n_cols;
    int nrow = udiff.n_rows;
    int ns = xi.n_rows;
    arma::mat C(nrow, ncol);
    C.zeros();
    
    int i,j,k;
    for(i=0; i<ns; i++){
        for(j=0; j<ncol; j++){
            for(k=0; k<nrow; k++){
                if(diffvec(i) == udiff(k)){
                    C(k, j) += xi(i,j);
                    break;
                }
            }
        }
    }
    return C;
}


//////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat get_taun(arma::mat Q, arma::vec udiff, arma::mat C){
    int nstate = Q.n_cols;
    //int ns = diffvec.n_rows;
    int ndelta = udiff.n_rows;
    int i,j,k;
    
    arma::mat countmat(nstate, nstate);
    arma::mat tempint(nstate, nstate);
    arma::mat pkl(nstate, nstate);
    arma::mat di(nstate, nstate);
    arma::vec tempvec(nstate*nstate);
    arma::mat cmat(nstate,nstate);
    arma::mat tempmat(nstate, nstate);
    arma::mat result(nstate,nstate);
    countmat.zeros();
    result.zeros();
    
    //arma::mat identity = arma::eye<arma::mat>(nstate,nstate);
    //Rcpp::Rcout<<identity<<std::endl;
    //double tempsum;
    //int a,b;
    
    for(i=0; i<nstate; i++){
        for(j=0; j<nstate; j++){
            for(k=0; k<ndelta; k++){
                tempint = matrixintegral(Q, udiff(k), i, j)  ;
                pkl = matrixexp(Q, udiff(k));
                
                //tempvec = rowsum(pkl);
                //for(l=0;l<nstate;l++){
                //    for(m=0;m<nstate;m++){
                //        pkl(l,m) = pkl(l,m) / tempvec(l);
                //    }
                //}
                
                di = tempint / pkl;
                //Rcpp::Rcout<<di<<std::endl;
                tempvec = C.row(k).t();
                cmat = vec2mat(tempvec, nstate, nstate, 0);
                
                if(j==i){
                    //tau_i
                    countmat(i,j) += matrixsum(cmat,di);
                }
                else{
                    //n_ij
                    countmat(i,j) += Q(i,j) * matrixsum(cmat,di);
                }
            }
        }
    }
    
    /*
     for(i=0; i<nstate; i++){
     tempsum = 0;
     for(j=0; j<nstate; j++){
     if(i!=j){
     result(i,j) = countmat(i,j) / countmat(i,i);
     tempsum += result(i,j);
     }
     }
     for(j=0; j<nstate;j++){
     if(i==j) result(i,j) = -tempsum;
     }
     }
     */
    return countmat;
}

//////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List fb_cont(arma::vec Pi, arma::mat P, arma::mat nodeprob,
                   long dim, arma::vec ntimes, arma::vec timeindex,
                   arma::vec udiff, arma::cube expms){
    
    
    int M = nodeprob.n_cols;
    int n = ntimes.n_rows;
    int j,m,t,i,k,location;
    double interval;
    arma::mat transit(M,M);
    double tempsum;
    arma::vec tempval;
    arma::mat tempmat(1,M);
    arma::mat tempmat2(M,M);
    
    arma::mat alpha(dim, M);
    arma::vec scale(dim);
    arma::mat beta(dim, M);
    arma::mat Gamma(dim, M);
    /*when multiple series, rows in between will be zeros
     for safety reasons, avoid multiple series.*/
    arma::mat xi(dim-1, M*M);
    
    arma::vec colsumgamma(M);
    arma::vec tempsumxi(M*M);
    arma::mat colsumxi(M,M);
    double loglik = 0;
    
    int count = 0;
    for(j=0; j<n; j++){
        
        //forward probs
        alpha.row(count) = Pi.t() % nodeprob.row(count);
        tempval = nodeprob.row(count) * Pi;
        scale(count) = tempval(0);//double
        for(m=0; m<M; m++) alpha(count, m) = alpha(count,m) / scale(count);
        
        for(t=1; t<ntimes(j); t++){
            
            interval = timeindex(count+t) - timeindex(count+t-1);
            location = locate(interval, udiff);
            transit = expms.slice(location);
            
            tempmat = alpha.row(count+t-1) * transit; //row matrix
            alpha.row(count+t) = tempmat % nodeprob.row(count+t);
            tempval = tempmat * nodeprob.row(count+t).t();
            scale(count+t) = tempval(0);//double
            for(m=0; m<M; m++) alpha(count+t, m) = alpha(count+t,m) / scale(count+t);
        }
        
        //backward probs and state conditional probs
        for(m=0; m<M; m++) beta(count + ntimes(j) - 1,m) = 1 / (M * scale(count + ntimes(j) - 1));
        Gamma.row(count+ntimes(j)-1) = alpha.row(count+ntimes(j)-1) % beta.row(count+ntimes(j)-1);
        tempval = alpha.row(count+ntimes(j)-1) * beta.row(count+ntimes(j)-1).t();
        for(m=0; m<M; m++) Gamma(count+ntimes(j)-1,m)=Gamma(count+ntimes(j)-1,m)/tempval(0); //double
        
        for(t=ntimes(j)-2; t>=0; t--){
            interval = timeindex(count+t+1) - timeindex(count+t);
            location = locate(interval, udiff);
            transit = expms.slice(location);
            
            tempmat = beta.row(count+t+1) % nodeprob.row(count+t+1);
            beta.row(count+t) = tempmat * transit.t();
            for(m=0; m<M; m++) beta(count+t,m) = beta(count+t,m) / scale(count + t);
            
            Gamma.row(count+t) = alpha.row(count+t) % beta.row(count+t);
            tempval = alpha.row(count+t) * beta.row(count+t).t();
            for(m=0; m<M; m++) Gamma(count+t,m)=Gamma(count+t,m)/tempval(0); //double
        }
        
        //transition probs
        for(t=0; t<ntimes(j)-1; t++){
            interval = timeindex(count+t+1) - timeindex(count+t);
            location = locate(interval, udiff);
            transit = expms.slice(location);
            
            tempmat = beta.row(count+t+1) % nodeprob.row(count+t+1);
            tempmat2 = transit % (alpha.row(count+t).t() * tempmat);
            
            tempsum = 0;
            for(i=0; i<M; i++){
                for(k=0; k<M; k++){
                    xi(count+t, i + k*M) = tempmat2(i,k);
                    tempsum += xi(count+t, i + k*M);
                }
            }
            //Rcpp::Rcout<<count+t<<std::endl;
            for(m=0; m<M*M; m++) xi(count+t, m) = xi(count+t, m) / tempsum;
        }
        
        count += ntimes(j);
    }
    
    
    //get the column sums
    colsumgamma = colsum(Gamma);
    tempsumxi = colsum(xi);
    
    for(i=0; i<M; i++)
        for(k=0; k<M; k++)
            colsumxi(i,k) = tempsumxi(i+k*M);
    
    loglik = arma::sum(log(scale));

    return Rcpp::List::create(Rcpp::Named("colsumgamma")=colsumgamma,
                              Rcpp::Named("colsumxi")=colsumxi,
                              Rcpp::Named("Gamma")=Gamma,
                              Rcpp::Named("xi")=xi,
                              Rcpp::Named("loglik")=loglik);
}



//[[Rcpp::export]]
arma::mat hmm_gen_cont (int dim, int M, arma::vec pi, arma::mat a, arma::vec theta,
                        arma::vec zeroprop, arma::vec timeindex){
    
    int m,n, prev, curr;
    arma::vec label(M);
    double u;
    double interval;
    arma::mat transit(M,M);
    arma::mat result(dim, 2); /*first column for x, second column for state*/
    
    for(m=0;m<M;m++) label(m) = m+1;
    
    
    //starting generation
    
        result(0,1) = multinomrand(1, M, pi, label)(0);
        
        curr = result(0,1) - 1;
        u = Rcpp::runif(1,0,1)(0);
        if(u<=zeroprop(curr)) result(0,0) = 0;
        else result(0,0) = Rcpp::rpois(1, theta(curr))(0);
        
    
    
    //iteration
    for(n=1; n<dim; n++){
        interval = timeindex(n) - timeindex(n-1);
        transit = matrixexp(a, interval);
            prev = result(n-1, 1) - 1;
            result(n,1) = multinomrand(1, M, transit.row(prev).t(), label)(0);  //row to column vetor
            
            curr = result(n,1) - 1;
            u = Rcpp::runif(1,0,1)(0);
            if(u<=zeroprop(curr)) result(n,0) = 0;
            else result(n,0) = Rcpp::rpois(1, theta(curr))(0);
 
    }
    
    return(result);
}

////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double ziploglik_nocov_cont(arma::vec delta, arma::mat gamma, double theta, arma::vec lambda,
                            arma::vec y, arma::vec ntimes, arma::vec timeindex, arma::vec udiff,
                            arma::cube expms){
    double loglik = 0;
    //long dim = y.n_rows;
    int M = lambda.n_rows;
    int n = ntimes.n_rows;
    int j,m,t;
    double tempsum;
    
    arma::vec forward(M);
    arma::vec meanvec(M);
    arma::vec forwardsum;  //then forwardsum(0) is double
    arma::mat tempmat(1,M);
    arma::mat transit(M,M);
    int count = 0, location, interval;
    
    for(j=0; j<n; j++){
        
        tempsum = 0;
        //initialize the forward variable
        meanvec(0) = dzip(theta, lambda(0), y(count), FALSE);
        for(m=1; m<M; m++) meanvec(m) = R::dpois(y(count), lambda(m), FALSE);
        
        forward = delta % meanvec; //element wise multiplication
        forwardsum = delta.t() * meanvec;
        tempsum = log(forwardsum(0));
        //Rcpp::Rcout << "loglik=" << loglik << std::endl;
        
        for(m=0; m<M; m++) forward(m) = forward(m) / forwardsum(0);
        
        //recursion
        for(t=1; t<ntimes(j); t++){
            
            meanvec(0) = dzip(theta, lambda(0), y(count+t), FALSE);
            for(m=1; m<M; m++) meanvec(m) = R::dpois(y(count+t), lambda(m), FALSE);
            interval = timeindex(count+t) - timeindex(count+t-1);
            location = locate(interval, udiff);
            transit = expms.slice(location);
            
            tempmat = forward.t() * transit; //row matrix
            forward = tempmat.t() % meanvec;
            forwardsum = tempmat * meanvec;
            tempsum += log(forwardsum(0));
            for(m=0; m<M; m++) forward(m) = forward(m) / forwardsum(0);
        }
        
        loglik += tempsum;
        count += ntimes(j);
        
    }
    
    return loglik;
    
}

///////////////////////////////////
//' retrieve the natural parameters from working parameters for a continuous-time
//' zero-inflated Poisson hidden Markov model where zero-inflation only happens in state 1
//' @param parm working parameters
//' @param M number of hidden states
//' @return a list of natural parameters
//' @export
// [[Rcpp::export]]
Rcpp::List retrieve_nocov_cont(arma::vec parm, int M){
    int nextindex = 0;
    arma::vec delta(M);
    arma::mat gamma(M,M);
    double theta;
    arma::vec lambda(M);
    int i, j, m;
    double tempsum;
    
    //retrieve the full parameters
    delta(0) = 1;
    tempsum = 1;
    for(m=1; m<M; m++) {
        delta(m) = exp(parm(m-1));
        tempsum += delta(m);
    }
    for(m=0; m<M; m++) delta(m) = delta(m) / tempsum;
    
    nextindex = M-1;
    for(i=0; i<M; i++){
        tempsum = 0;
        for(j=0; j<M; j++){
            if(i!=j){
                gamma(i,j) = exp(parm(nextindex)) / (1+exp(parm(nextindex)));
                tempsum += gamma(i,j);
                nextindex ++;
            }
        }
        gamma(i,i) = -tempsum;
    }
    
    nextindex = M*M-1;
    theta = exp(parm(nextindex))/(1+exp(parm(nextindex)));
    
    for(m=0; m<M; m++) lambda(m) = exp(parm(nextindex+m+1));
    
    return Rcpp::List::create(Rcpp::Named("delta")=delta,
                              Rcpp::Named("gamma")=gamma,
                              Rcpp::Named("theta")=theta,
                              Rcpp::Named("lambda")=lambda);
    
}

//////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////
//' negative log likelihood function for zero-inflated Poisson hidden Markov model without covariates,
//' where zero-inflation only happens in state 1
//' @param parm working parameters
//' @param M number of hidden states
//' @param y observed series
//' @param ntimes length of the observed series
//' @param timeindex vector of observed time points
//' @param udiff unique time intervals
//' @return negative log likelihood
//' @export
// [[Rcpp::export]]
double zipnegloglik_nocov_cont(arma::vec parm, int M, arma::vec y, arma::vec ntimes,
                               arma::vec timeindex, arma::vec udiff){
    
    Rcpp::List mod = retrieve_nocov_cont(parm,M);
    arma::vec delta = mod("delta");
    arma::mat gamma = mod("gamma");
    double theta = mod("theta");
    arma::vec lambda = mod("lambda");
    
    arma::cube expms = getallexpm(gamma, udiff);
                                  
    double negloglik = - ziploglik_nocov_cont(delta, gamma, theta, lambda, y, ntimes,
                                              timeindex,udiff,expms);
    return negloglik;
}


////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::vec grad_ziploglik_nocov_cont(arma::vec delta, arma::mat gamma, double theta, arma::vec lambda,
                                    arma::vec y, arma::vec ntimes, arma::vec timeindex,
                                    arma::vec udiff, arma::cube expms){
    
    long dim = y.n_rows;
    int M = lambda.n_rows;
    int n = ntimes.n_rows;
    int i,j,m,t,k;
    arma::mat nodeprob(dim, M);

    int count = 0;
    for(j=0; j<n; j++){
        for(t=0; t<ntimes(j); t++){
            nodeprob(count+t,0) = dzip(theta, lambda(0), y(count+t), FALSE);
            for(m=1; m<M; m++) nodeprob(count+t,m) = R::dpois(y(count+t), lambda(m), FALSE);
        }
        count += ntimes(j);
    }
    
    Rcpp::List fb = fb_cont(delta, gamma, nodeprob, dim, ntimes, timeindex, udiff, expms);
    arma::mat Gamma = fb("Gamma");
    arma::vec colsumgamma = fb("colsumgamma");
    arma::mat xi = fb("xi");
    //arma::mat xi = fb("xi");
    
    //get the derivatives
    //prior
    arma::vec temppi(M);
    arma::vec dlogitpi(M);
    
    for(m=0; m<M; m++) temppi(m) = 0;
    count = 0;
    for(j=0; j<n; j++){
        for(m=0; m<M; m++) temppi(m) += Gamma(count,m);
        count += ntimes(j);
    }
    
    for(m=0; m<M; m++) {
        temppi(m) = temppi(m) / n;
        dlogitpi(m) = temppi(m) - delta(m);
    }
    
    
    //tpm
    arma::vec vdiff(dim-1);
    for(t=0; t<dim-1; t++) vdiff(t) = timeindex(t+1) - timeindex(t);
    arma::mat C = get_c(xi, vdiff, udiff );
    arma::mat countmat = get_taun(gamma, udiff, C);

    arma::mat dlogtrans(M,M);
    for(i=0; i<M; i++)
        for(k=0; k<M; k++)
            dlogtrans(i,k) = (1-gamma(i,k))*(countmat(i,k) - countmat(i,i)*gamma(i,k));
            //dlogtrans(i,k) = countmat(i,k) - countmat(i,i) * gamma(i,k);
    
    
    //theta and lambda
    double dlogtheta = 0;
    arma::vec dloglambda(M);
    for(m=0; m<M; m++) dloglambda(m) = 0;
    
    for(t=0; t<dim; t++){
        
        if(y(t)>0){
            dlogtheta -= theta * Gamma(t,0);
            dloglambda(0) += Gamma(t,0) * (y(t) - lambda(0));
        }
        else{
            dlogtheta -= Gamma(t,0)*theta*(1-theta)*(exp(-lambda(0))-1) / (theta+(1-theta)*exp(-lambda(0)));
            dloglambda(0) += Gamma(t,0)*(theta-1)*lambda(0)*exp(-lambda(0)) / (theta+(1-theta)*exp(-lambda(0)));
        }
        
        for(m=1; m<M; m++){
            dloglambda(m) += Gamma(t,m)*(y(t) - lambda(m));
        }
    }
    
    
    arma::vec result(M-1+M*(M-1)+1+M);
    int nextindex = 0;
    for(j=0; j<M-1; j++) result(j) = dlogitpi(j+1);
    nextindex = M - 1;
    for(i=0; i<M; i++){
        for(j=0; j<M; j++){
            if(j!=i){
                result(nextindex) = dlogtrans(i,j);
                nextindex++;
            }
        }
    }
    
    result(nextindex) = dlogtheta;
    nextindex += 1;
    for(i=0; i<M; i++) result(nextindex+i) = dloglambda(i);
    
    return result;
    
}


//////////////////////////////////////////////////////////////////////////////////
//' gradient for negative log likelihood function from zero-inflated Poisson hidden Markov model
//' without covariates, where zero-inflation only happens in state 1
//' @param parm working parameters
//' @param M number of hidden states
//' @param y observed series
//' @param ntimes length of the observed series
//' @param timeindex vector of observed time points
//' @param udiff unique time intervals
//' @return gradient for negative log likelihood
//' @export
// [[Rcpp::export]]
arma::vec grad_zipnegloglik_nocov_cont(arma::vec parm, int M, arma::vec y, arma::vec ntimes,
                                       arma::vec timeindex,arma::vec udiff){
    
    Rcpp::List mod = retrieve_nocov_cont(parm,M);
    arma::vec delta = mod("delta");
    arma::mat gamma = mod("gamma");
    double theta = mod("theta");
    arma::vec lambda = mod("lambda");
    
    arma::cube expms = getallexpm(gamma, udiff);
    
    ////change to negative !!!!!!
    arma::vec result = grad_ziploglik_nocov_cont(delta, gamma, theta, lambda, y, ntimes,
                                            timeindex,udiff, expms);
    int m;
    for(m=0; m<M-1+M*(M-1)+1+M; m++) result(m) = -result(m);
    return result;
}

////////////////////////////////////////
// [[Rcpp::export]]
arma::vec hmm_viterbi_cont(arma::vec pi, arma::mat a, arma::vec theta, int M, arma::vec y,
                           arma::vec zeroprop, arma::vec timeindex, arma::vec udiff, arma::cube expms){
    
    long dim = y.n_rows;
    int i,j,m,n, location;
    double colmax, interval;
    
    arma::vec meanvec(M);
    arma::vec forward(M);
    arma::vec forwardsum;
    arma::vec state(dim);
    arma::vec tempmax(M);
    arma::mat xi(dim, M);
    arma::mat tempmat(M,M);
    arma::mat transit(M,M);
    
    //initialize the forward variable
    
    for(m=0; m<M; m++) meanvec(m) = dzip(zeroprop(m), theta(m), y(0), FALSE);
    
    forward = pi % meanvec; //element wise multiplication
    forwardsum = pi.t() * meanvec;
    for(m=0; m<M; m++) xi(0,m) = forward(m) / forwardsum(0);
    
    //recursion
    for(n=1; n<dim; n++){
        
        for(m=0; m<M; m++) meanvec(m) = dzip(zeroprop(m), theta(m), y(n), FALSE);
        
        interval = timeindex(n) - timeindex(n-1);
        location = locate(interval, udiff);
        transit = expms.slice(location);
        
        //difficult part
        for(i=0; i<M; i++){
            for(j=0;j<M;j++){
                tempmat(i,j) = xi(n-1, i) * transit(i,j);
            }
        }
        
        for(j=0; j<M; j++){
            for(i=0;i<M;i++){
                if(i==0) colmax = tempmat(i,j);
                else colmax = MAX(colmax, tempmat(i,j));
            }
            tempmax(j) = colmax;
        }
        
        
        forward = tempmax % meanvec;
        forwardsum = tempmax.t() * meanvec;
        for(m=0; m<M; m++) xi(n,m) = forward(m) / forwardsum(0);
    }
    
    //termination
    for(m=0; m<M; m++){
        if(m==0) {
            colmax = xi(dim-1, m);
            state(dim-1) = m;
        }
        else{
            if(xi(dim-1,m) >= colmax) {
                state(dim-1) = m;
                colmax = xi(dim-1,m);
            }
        }
    }
    
    //backtracking
    for(n=dim-2; n>=0; n--){
        
        for(m=0;m<M;m++){
            interval = timeindex(n+1) - timeindex(n);
            location = locate(interval, udiff);
            transit = expms.slice(location);

            if(m==0) {
                colmax = xi(n,m) * transit(m, state(n+1));
                state(n) = m;
            }
            else{
                if(xi(n,m) * transit(m, state(n+1))>=colmax) {
                    state(n) = m;
                    colmax = xi(n,m) * transit(m, state(n+1));
                }
            }
        }
    }
    
    for(n=0; n<dim; n++) state(n) = state(n) + 1;
    return(state);
}


///////////////////////////////////////////
//[[Rcpp::export]]
arma::mat hmm_cov_gen_cont (arma::vec parm, int M, long dim, int ncolcovpi, arma::mat covpi,
                       int ncolcovtrans, arma::mat covtrans, int ncolcovp1, arma::mat covp1,
                       int ncolcovpois, arma::mat covpois, arma::vec zeroindex, arma::vec timeindex){
    
    int i,j,m,n,nextindex,prev,curr;
    double tempsum, u, interval;
    arma::vec label(M);
    
    arma::mat result(dim, 2); /*first column for x, second column for state*/
    
    arma::mat transit(M,M);
    arma::vec pi(M);
    arma::mat a(M,M);
    arma::vec theta(M);
    arma::vec zeroprop(M);
    
    //retrieve the full parameters
    //prior
    nextindex=0;
    pi(0) = 1;
    tempsum = 1;
    for(m=1; m<M; m++){
        pi(m) = parm(nextindex);
        
        for(j=0;j<ncolcovpi;j++)
            pi(m) += parm(nextindex+j+1)*covpi(0,j) ;
        
        pi(m) = exp(pi(m));
        
        tempsum += pi(m);
        nextindex = nextindex + ncolcovpi + 1;
    }
    
    for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
    
    
    //transition
    for(i=0; i<M; i++){
        tempsum = 0;
        for(j=0; j<M; j++){
            if(i!=j){
                a(i,j) = exp(parm(nextindex))/(1+exp(parm(nextindex)));
                tempsum += a(i,j);
                nextindex++;
            }
        }
        a(i,i) = -tempsum;
        
    }
    
    
    for(i=0;i<M;i++){
        if(zeroindex(i)==0) zeroprop(i)=0;
        else{
            zeroprop(i) = parm(nextindex);
            
            for(m=0;m<ncolcovp1;m++)
                zeroprop(i) += parm(nextindex+m+1) * covp1(0,m);
            
            zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
            nextindex = nextindex + ncolcovp1 + 1;
        }
    }
    
    //Rcpp::Rcout << "p1=" << p1 << std::endl;
    
    for(m=0; m<M; m++){
        theta(m) = parm(nextindex);
        for(j=0; j<ncolcovpois;j++){
            theta(m) += parm(nextindex+j+1) * covpois(0,j);
        }
        theta(m) = exp(theta(m));
        nextindex = nextindex + ncolcovpois + 1;
    }
    
    
    //start generation
    
    for(m=0;m<M;m++) label(m) = m+1;
    
    result(0,1) = multinomrand(1, M, pi, label)(0);
    curr = result(0,1) - 1;
    u = Rcpp::runif(1,0,1)(0);
    if(u<=zeroprop(curr)) result(0,0) = 0;
    else result(0,0) = Rcpp::rpois(1, theta(curr))(0);
    
    
    
    //iteration steps
    for(n=1; n<dim; n++){
        
        //still need to retrieve the natural parameters
        nextindex=0;
        pi(0) = 1;
        tempsum = 1;
        for(m=1; m<M; m++){
            pi(m) = parm(nextindex);
            
            for(j=0;j<ncolcovpi;j++)
                pi(m) += parm(nextindex+j+1)*covpi(n,j) ;
            
            pi(m) = exp(pi(m));
            
            tempsum += pi(m);
            nextindex = nextindex + ncolcovpi + 1;
        }
        
        for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
        
        
        //transition
        for(i=0; i<M; i++){
            tempsum = 0;
            for(j=0; j<M; j++){
                if(i!=j){
                    a(i,j) = exp(parm(nextindex))/(1+exp(parm(nextindex)));
                    tempsum += a(i,j);
                    nextindex++;
                }
            }
            a(i,i) = -tempsum;
            
        }
        
        for(i=0;i<M;i++){
            if(zeroindex(i)==0) zeroprop(i)=0;
            else{
                zeroprop(i) = parm(nextindex);
                
                for(m=0;m<ncolcovp1;m++)
                    zeroprop(i) += parm(nextindex+m+1) * covp1(n,m);
                
                zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
                nextindex = nextindex + ncolcovp1 + 1;
            }
        }
        //Rcpp::Rcout << "p1=" << p1 << std::endl;
        
        for(m=0; m<M; m++){
            theta(m) = parm(nextindex);
            for(j=0; j<ncolcovpois;j++){
                theta(m) += parm(nextindex+j+1) * covpois(n,j);
            }
            theta(m) = exp(theta(m));
            nextindex = nextindex + ncolcovpois + 1;
        }
        
        interval = timeindex(n) - timeindex(n-1);
        transit = matrixexp(a, interval);
        prev = result(n-1, 1) - 1;
        result(n,1) = multinomrand(1, M, transit.row(prev).t(), label)(0);  //row to column vetor
        
        curr = result(n,1) - 1;
        u = Rcpp::runif(1,0,1)(0);
        if(u<=zeroprop(curr)) result(n,0) = 0;
        else result(n,0) = Rcpp::rpois(1, theta(curr))(0);
  
    }
    
    return result;
    
}



////////////////////
//////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double ziploglik_cov_cont(arma::vec delta, arma::mat gamma, arma::vec thetaparm, arma::mat lambdaparm,
                          arma::vec y, arma::mat x, arma::vec ntimes,
                          arma::vec timeindex, arma::vec udiff, arma::cube expms){
    //the first column of x should be 1 for intercept
    double loglik = 0, interval;
    //long dim = y.n_rows;
    int M = gamma.n_rows;
    int n = ntimes.n_rows;
    int j,m,t,location;
    arma::vec tempval;
    double tempsum;
    double theta;
    arma::vec lambda(M);
    arma::mat transit(M,M);
    
    arma::vec forward(M);
    arma::vec meanvec(M);
    arma::vec forwardsum;  //then forwardsum(0) is double
    arma::mat tempmat(1,M);
    
    int count = 0;
    
    for(j=0; j<n; j++){
        
        tempsum = 0;
        
        //initialize the forward variable
        tempval = x.row(count) * thetaparm;
        theta = exp(tempval(0))/(1+exp(tempval(0)));
        for(m=0;m<M;m++) {
            tempval = x.row(count) * lambdaparm.row(m).t();
            lambda(m) = exp(tempval(0));
        }
        
        meanvec(0) = dzip(theta, lambda(0), y(count), FALSE);
        for(m=1; m<M; m++) meanvec(m) = R::dpois(y(count), lambda(m), FALSE);
        
        forward = delta % meanvec; //element wise multiplication
        forwardsum = delta.t() * meanvec;
        tempsum = log(forwardsum(0));
        //Rcpp::Rcout << "loglik=" << loglik << std::endl;
        
        for(m=0; m<M; m++) forward(m) = forward(m) / forwardsum(0);
        
        //recursion
        for(t=1; t<ntimes(j); t++){
            
            tempval = x.row(count+t) * thetaparm;
            theta = exp(tempval(0))/(1+exp(tempval(0)));
            for(m=0;m<M;m++) {
                tempval = x.row(count+t) * lambdaparm.row(m).t();
                lambda(m) = exp(tempval(0));
            }
            
            meanvec(0) = dzip(theta, lambda(0), y(count+t), FALSE);
            for(m=1; m<M; m++) meanvec(m) = R::dpois(y(count+t), lambda(m), FALSE);
            
            interval = timeindex(count+t) - timeindex(count+t-1);
            location = locate(interval, udiff);
            transit = expms.slice(location);
            
            tempmat = forward.t() * transit; //row matrix
            forward = tempmat.t() % meanvec;
            forwardsum = tempmat * meanvec;
            tempsum += log(forwardsum(0));
            for(m=0; m<M; m++) forward(m) = forward(m) / forwardsum(0);
        }
        
        loglik += tempsum;
        count += ntimes(j);
        
    }
    
    return loglik;
    
}

///////////////////////////////////
//' retrieve the natural parameters from the working parameters in zero-inflated Poisson
//' hidden Markov model with covariates, where zero-inflation only happens in state 1
//' @param parm working parameters
//' @param M number of hidden states
//' @param ncolx number of covariates including the intercept
//' @return a list of natural parameters
//' @export
// [[Rcpp::export]]
Rcpp::List retrieve_cov_cont(arma::vec parm, int M, int ncolx){
    int nextindex = 0;
    arma::vec delta(M);
    arma::mat gamma(M,M);
    arma::vec thetaparm(ncolx);
    arma::mat lambdaparm(M,ncolx);
    int i, j, m;
    double tempsum;
    
    //retrieve the full parameters
    delta(0) = 1;
    tempsum = 1;
    for(m=1; m<M; m++) {
        delta(m) = exp(parm(m-1));
        tempsum += delta(m);
    }
    for(m=0; m<M; m++) delta(m) = delta(m) / tempsum;
    
    nextindex = M-1;
    for(i=0; i<M; i++){
        tempsum = 0;
        for(j=0; j<M; j++){
            if(i!=j){
                gamma(i,j) = exp(parm(nextindex)) / (1+exp(parm(nextindex)));
                tempsum += gamma(i,j);
                nextindex ++;
            }
        }
        gamma(i,i) = -tempsum;
    }
    
    for(i=0; i<ncolx; i++) thetaparm(i) = parm(nextindex+i);
    nextindex += ncolx;
    
    for(i=0; i<M; i++){
        for(j=0; j<ncolx; j++){
            lambdaparm(i,j) = parm(nextindex + i*ncolx + j);
        }
    }
    
    return Rcpp::List::create(Rcpp::Named("delta")=delta,
                              Rcpp::Named("gamma")=gamma,
                              Rcpp::Named("thetaparm")=thetaparm,
                              Rcpp::Named("lambdaparm")=lambdaparm);
    
}


//////////////////////////////////////////////////////////////////////////////////
//' negative log likelihood function for zero-inflated Poisson hidden Markov model with covariates,
//' where zero-inflation only happens in state 1
//' @param parm working parameters
//' @param y observed series
//' @param covariates design matrix of covariates including the intercept
//' @param M number of hidden states
//' @param ntimes length of the observed series
//' @param timeindex vector of observed time points
//' @param udiff unique time intervals
//' @return negative log likelihood
//' @export
// [[Rcpp::export]]
double zipnegloglik_cov_cont(arma::vec parm, arma::vec y, arma::mat covariates, int M, arma::vec ntimes,
                             arma::vec timeindex, arma::vec udiff){
    int ncolx = covariates.n_cols;
    Rcpp::List mod = retrieve_cov_cont(parm,M,ncolx);
    arma::vec delta = mod("delta");
    arma::mat gamma = mod("gamma");
    arma::vec thetaparm = mod("thetaparm");
    arma::mat lambdaparm = mod("lambdaparm");
    arma::cube expms = getallexpm(gamma, udiff);
    
    double negloglik = - ziploglik_cov_cont(delta, gamma, thetaparm, lambdaparm, y, covariates, ntimes,
                                            timeindex, udiff, expms);
    return negloglik;
}


////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::vec grad_ziploglik_cov_cont(arma::vec delta, arma::mat gamma, arma::vec thetaparm, arma::mat lambdaparm,
                                  arma::vec y, arma::mat x, arma::vec ntimes, arma::vec timeindex,
                                  arma::vec udiff, arma::cube expms){
    
    long dim = y.n_rows;
    int ncolx = x.n_cols;
    int M = gamma.n_rows;
    int n = ntimes.n_rows;
    int i,j,m,t,k;
    double theta;
    arma::vec tempval;
    arma::vec lambda(M);
    arma::mat nodeprob(dim, M);
    
    
    int count = 0;
    for(j=0; j<n; j++){
        for(t=0; t<ntimes(j); t++){
            tempval = x.row(count+t) * thetaparm;
            theta = exp(tempval(0))/(1+exp(tempval(0)));
            for(m=0;m<M;m++) {
                tempval = x.row(count+t) * lambdaparm.row(m).t();
                lambda(m) = exp(tempval(0));
            }
            
            nodeprob(count+t,0) = dzip(theta, lambda(0), y(count+t), FALSE);
            for(m=1; m<M; m++) nodeprob(count+t,m) = R::dpois(y(count+t), lambda(m), FALSE);
        }
        count += ntimes(j);
    }
    
    
    Rcpp::List fb = fb_cont(delta, gamma, nodeprob, dim, ntimes, timeindex, udiff, expms);
    arma::mat Gamma = fb("Gamma");
    arma::vec colsumgamma = fb("colsumgamma");
    arma::mat colsumxi = fb("colsumxi");
    arma::mat xi = fb("xi");
    
    //get the derivatives
    //prior
    arma::vec temppi(M);
    arma::vec dlogitpi(M);
    
    for(m=0; m<M; m++) temppi(m) = 0;
    count = 0;
    for(j=0; j<n; j++){
        for(m=0; m<M; m++) temppi(m) += Gamma(count,m);
        count += ntimes(j);
    }
    
    for(m=0; m<M; m++) {
        temppi(m) = temppi(m) / n;
        dlogitpi(m) = temppi(m) - delta(m);
    }
    
    
    //tpm
    arma::vec vdiff(dim-1);
    for(t=0; t<dim-1; t++) vdiff(t) = timeindex(t+1) - timeindex(t);
    arma::mat C = get_c(xi, vdiff, udiff );
    arma::mat countmat = get_taun(gamma, udiff, C);
    
    arma::mat dlogtrans(M,M);
    for(i=0; i<M; i++)
        for(k=0; k<M; k++)
            dlogtrans(i,k) = (1-gamma(i,k))*(countmat(i,k) - countmat(i,i)*gamma(i,k));
            //dlogtrans(i,k) = countmat(i,k) - countmat(i,i) * gamma(i,k);
    
    //theta and lambda
    arma::vec dlogtheta(ncolx);
    arma::mat dloglambda(M,ncolx);
    for(m=0; m<ncolx; m++) dlogtheta(m) = 0;
    for(i=0; i<M; i++)
        for(j=0; j<ncolx; j++)
            dloglambda(i,j) = 0;
    
    for(t=0; t<dim; t++){
        
        tempval = x.row(t) * thetaparm;
        theta = exp(tempval(0))/(1+exp(tempval(0)));
        for(m=0;m<M;m++) {
            tempval = x.row(t) * lambdaparm.row(m).t();
            lambda(m) = exp(tempval(0));
        }
        
        if(y(t)>0){
            for(i=0;i<ncolx;i++){
                dlogtheta(i) -= theta * Gamma(t,0) * x(t,i) ;
                dloglambda(0,i) += Gamma(t,0) * x(t,i) * (y(t) - lambda(0));
            }
        }
        else{
            for(i=0;i<ncolx;i++){
                dlogtheta(i) -= Gamma(t,0)*x(t,i)*theta*(1-theta)*(exp(-lambda(0))-1) / (theta+(1-theta)*exp(-lambda(0)));
                dloglambda(0,i) += Gamma(t,0)*x(t,i)*(theta-1)*lambda(0)*exp(-lambda(0)) / (theta+(1-theta)*exp(-lambda(0)));
            }
        }
        
        for(m=1; m<M; m++){
            for(i=0;i<ncolx;i++){
                dloglambda(m,i) += Gamma(t,m)*x(t,i)*(y(t) - lambda(m));
            }
        }
    }
    
    ///////
    arma::vec result(M-1+M*(M-1)+ncolx*(1+M));
    int nextindex = 0;
    for(j=0; j<M-1; j++) result(j) = dlogitpi(j+1);
    nextindex = M - 1;
    nextindex = M - 1;
    for(i=0; i<M; i++){
        for(j=0; j<M; j++){
            if(j!=i){
                result(nextindex) = dlogtrans(i,j);
                nextindex++;
            }
        }
    }

    for(i=0; i<ncolx; i++)
        result(nextindex+i) = dlogtheta(i);
    
    nextindex += ncolx;
    for(i=0; i<M; i++)
        for(j=0; j<ncolx; j++)
            result(nextindex+i*ncolx+j) = dloglambda(i,j);
    
    
    return result;
    
    
}


//////////////////////////////////////////////////////////////////////////////////
//' gradient for negative log likelihood function in zero-inflated Poisson hidden Markov model
//' with covariates, where zero-inflation only happens in state 1
//' @param parm working parameters
//' @param y observed series
//' @param covariates design matrix of covariates including the intercept
//' @param M number of hidden states
//' @param ntimes length of the observed series
//' @param timeindex vector of observed time points
//' @param udiff unique time intervals
//' @return gradient for negative log likelihood
//' @export
// [[Rcpp::export]]
arma::vec grad_zipnegloglik_cov_cont(arma::vec parm, arma::vec y, arma::mat covariates, int M,
                                     arma::vec ntimes, arma::vec timeindex, arma::vec udiff){
    int ncolx = covariates.n_cols;
    Rcpp::List mod = retrieve_cov_cont(parm,M,ncolx);
    arma::vec delta = mod("delta");
    arma::mat gamma = mod("gamma");
    arma::vec thetaparm = mod("thetaparm");
    arma::mat lambdaparm = mod("lambdaparm");
    arma::cube expms = getallexpm(gamma, udiff);
    
    ////change to negative !!!!!!
    arma::vec result = grad_ziploglik_cov_cont(delta, gamma, thetaparm, lambdaparm, y, covariates, ntimes,
                                               timeindex, udiff, expms);
    int m;
    for(m=0; m<M-1+M*(M-1)+ncolx*(1+M); m++) result(m) = -result(m);
    return result;
}


///////////////////////////////////////////////////////////
/////////////////////////

// [[Rcpp::export]]
arma::vec hmm_cov_viterbi_cont(arma::vec parm, int M, arma::vec y, int ncolcovpi, arma::mat covpi,
                          int ncolcovtrans, arma::mat covtrans, int ncolcovp1, arma::mat covp1,
                          int ncolcovpois, arma::mat covpois, arma::vec zeroindex,
                          arma::vec timeindex, arma::vec udiff, arma::cube expms){
    
    long dim = y.n_rows;
    int i,j,m,n,nextindex, location;
    double tempsum, colmax, interval;
    
    arma::vec pi(M);
    arma::mat a(M,M);
    arma::vec theta(M);
    arma::vec zeroprop(M);
    
    arma::vec meanvec(M);
    arma::vec forward(M);
    arma::vec forwardsum;
    arma::vec state(dim);
    arma::vec tempmax(M);
    arma::mat xi(dim, M);
    arma::mat tempmat(M,M);
    arma::mat transit(M,M);
    //retrieve the full parameters
    //prior
    nextindex=0;
    pi(0) = 1;
    tempsum = 1;
    for(m=1; m<M; m++){
        pi(m) = parm(nextindex);
        
        for(j=0;j<ncolcovpi;j++)
            pi(m) += parm(nextindex+j+1)*covpi(0,j) ;
        
        pi(m) = exp(pi(m));
        
        tempsum += pi(m);
        nextindex = nextindex + ncolcovpi + 1;
    }
    
    for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;

    //transition
 
    for(i=0; i<M; i++){
        tempsum = 0;
        for(j=0; j<M; j++){
            if(i!=j){
                a(i,j) = exp(parm(nextindex))/(1+exp(parm(nextindex)));
                tempsum += a(i,j);
                nextindex++;
            }
        }
        a(i,i) = -tempsum;
        
    }
    
    for(i=0;i<M;i++){
        if(zeroindex(i)==0) zeroprop(i)=0;
        else{
            zeroprop(i) = parm(nextindex);
            
            for(m=0;m<ncolcovp1;m++)
                zeroprop(i) += parm(nextindex+m+1) * covp1(0,m);
            
            zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
            nextindex = nextindex + ncolcovp1 + 1;
        }
    }
    
    //Rcpp::Rcout << "p1=" << p1 << std::endl;
    
    for(m=0; m<M; m++){
        theta(m) = parm(nextindex);
        for(j=0; j<ncolcovpois;j++){
            theta(m) += parm(nextindex+j+1) * covpois(0,j);
        }
        theta(m) = exp(theta(m));
        nextindex = nextindex + ncolcovpois + 1;
    }
    
    //Rcpp::Rcout << nextindex << std::endl;
    
    //initialize the forward variable
    //check emit_dist
    
    for(m=0; m<M; m++) meanvec(m) = dzip(zeroprop(m), theta(m), y(0), FALSE);
    
    
    forward = pi % meanvec; //element wise multiplication
    forwardsum = pi.t() * meanvec;
    for(m=0; m<M; m++) xi(0,m) = forward(m) / forwardsum(0);
    
    //recursion
    for(n=1; n<dim; n++){
        
        nextindex = M * M - 1;
        for(i=0;i<M;i++){
            if(zeroindex(i)==0) zeroprop(i)=0;
            else{
                zeroprop(i) = parm(nextindex);
                
                for(m=0;m<ncolcovp1;m++)
                    zeroprop(i) += parm(nextindex+m+1) * covp1(n,m);
                
                zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
                nextindex = nextindex + ncolcovp1 + 1;
            }
        }
        
        //Rcpp::Rcout << "p1=" << p1 << std::endl;
        
        for(m=0; m<M; m++){
            theta(m) = parm(nextindex);
            for(j=0; j<ncolcovpois;j++){
                theta(m) += parm(nextindex+j+1) * covpois(n,j);
            }
            theta(m) = exp(theta(m));
            nextindex = nextindex + ncolcovpois + 1;
        }
        
        
        for(m=0; m<M; m++) meanvec(m) = dzip(zeroprop(m), theta(m), y(n), FALSE);
        
        
        //difficult part
        interval = timeindex(n) - timeindex(n-1);
        location = locate(interval, udiff);
        transit = expms.slice(location);
        
        for(i=0; i<M; i++){
            for(j=0;j<M;j++){
                tempmat(i,j) = xi(n-1, i) * transit(i,j);
            }
        }
        
        for(j=0; j<M; j++){
            for(i=0;i<M;i++){
                if(i==0) colmax = tempmat(i,j);
                else colmax = MAX(colmax, tempmat(i,j));
            }
            tempmax(j) = colmax;
        }
        
        
        forward = tempmax % meanvec;
        forwardsum = tempmax.t() * meanvec;
        for(m=0; m<M; m++) xi(n,m) = forward(m) / forwardsum(0);
    }
    
    //termination
    for(m=0; m<M; m++){
        if(m==0) {
            colmax = xi(dim-1, m);
            state(dim-1) = m;
        }
        else{
            if(xi(dim-1,m) >= colmax) {
                state(dim-1) = m;
                colmax = xi(dim-1,m);
            }
        }
    }
    
    //backtracking
    for(n=dim-2; n>=0; n--){
        
        
        interval = timeindex(n+1) - timeindex(n);
        location = locate(interval, udiff);
        transit = expms.slice(location);

        
        for(m=0;m<M;m++){
            if(m==0) {
                colmax = xi(n,m) * transit(m, state(n+1));
                state(n) = m;
            }
            else{
                if(xi(n,m) * transit(m, state(n+1))>=colmax) {
                    state(n) = m;
                    colmax = xi(n,m) * transit(m, state(n+1));
                }
            }
        }
    }
    
    for(n=0; n<dim; n++) state(n) = state(n) + 1;
    return(state);
}


////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List smooth_nocov_cont(arma::vec delta, arma::mat gamma, double theta, arma::vec lambda,
                                    arma::vec y, arma::vec ntimes, arma::vec timeindex,
                                    arma::vec udiff, arma::cube expms){
    
    long dim = y.n_rows;
    int M = lambda.n_rows;
    int n = ntimes.n_rows;
    int j,m,t;
    arma::mat nodeprob(dim, M);
    
    int count = 0;
    for(j=0; j<n; j++){
        for(t=0; t<ntimes(j); t++){
            nodeprob(count+t,0) = dzip(theta, lambda(0), y(count+t), FALSE);
            for(m=1; m<M; m++) nodeprob(count+t,m) = R::dpois(y(count+t), lambda(m), FALSE);
        }
        count += ntimes(j);
    }
    
    Rcpp::List fb = fb_cont(delta, gamma, nodeprob, dim, ntimes, timeindex, udiff, expms);
    //arma::mat Gamma = fb("Gamma");
    //arma::vec colsumgamma = fb("colsumgamma");
    //arma::mat xi = fb("xi");
    //arma::mat xi = fb("xi");
    return fb;
}

////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List smooth_cov_cont(arma::vec delta, arma::mat gamma, arma::vec thetaparm, arma::mat lambdaparm,
                                  arma::vec y, arma::mat x, arma::vec ntimes, arma::vec timeindex,
                                  arma::vec udiff, arma::cube expms){
    
    long dim = y.n_rows;
    //int ncolx = x.n_cols;
    int M = gamma.n_rows;
    int n = ntimes.n_rows;
    int j,m,t;
    double theta;
    arma::vec tempval;
    arma::vec lambda(M);
    arma::mat nodeprob(dim, M);
    
    
    int count = 0;
    for(j=0; j<n; j++){
        for(t=0; t<ntimes(j); t++){
            tempval = x.row(count+t) * thetaparm;
            theta = exp(tempval(0))/(1+exp(tempval(0)));
            for(m=0;m<M;m++) {
                tempval = x.row(count+t) * lambdaparm.row(m).t();
                lambda(m) = exp(tempval(0));
            }
            
            nodeprob(count+t,0) = dzip(theta, lambda(0), y(count+t), FALSE);
            for(m=1; m<M; m++) nodeprob(count+t,m) = R::dpois(y(count+t), lambda(m), FALSE);
        }
        count += ntimes(j);
    }
    
    
    Rcpp::List fb = fb_cont(delta, gamma, nodeprob, dim, ntimes, timeindex, udiff, expms);
    return fb;
}

///////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
// [[Rcpp::export]]

arma::mat hsmm_cov_gen_aft(arma::vec delta, arma::vec lden, arma::mat dtparm, arma::vec tpmparm,
                           arma::vec zeroparm, arma::mat emitparm, arma::vec trunc,
                           int M, long dim, int ncolcovomega,
                           arma::mat covp,arma::mat covomega,
                           arma::mat covp1, arma::mat covpois){
    
    /*parm = logitp from logarithm dist---M, logit(pi)---M-1, logit(p1), log(theta)---M,
     omega---M*(M-2): (1,3)...(1,n)(2,3)...(2,n)(3,2)(3,4)...(3,n).......*/
    
    //intercept column is not included in the covariates, but should be in the parm
    
    //parameters are given alternatively in the order of beta0 and beta1
    
    int i,j,k,m,nextindex,prev,curr;
    long count,n;
    double tempsum;
    arma::vec label(M);
    
    arma::mat result(dim, 2); /*first column for x, second column for state*/
    arma::vec theta(M);
    double zeroprop;
    
    arma::mat omega(M,M);
    //arma::mat a(totalmv,totalmv);
    
    //////
    //retrieve some of the parameters
    nextindex = 0;
    //dwell time
    
    //check dt_dist
    arma::vec eb(M);
    arma::vec temp(M);
    arma::vec tempval(1);
    temp = dtparm * covp.row(0).t(); //the first column is not intercept
    //Rcpp::Rcout<<temp<<std::endl;
    
    for(i=0;i<M;i++)
        eb(i) = exp(temp(i));
    
    
    //////////////////
    //start generation
    for(m=0;m<M;m++) label(m) = m+1;
    
    
    //starting generation
    curr = multinomrand(1, M, delta, label)(0);
    
    //check dt_dist
    count = random_expbase(lden(curr-1), eb(curr-1), trunc(curr-1));
    
    
    /////////////////
    for(k=0;k<count;k++) {
        nextindex = 0;
        //zeroprops
        tempval = zeroparm.t() * covp1.row(k).t(); //the first column is intercept
        zeroprop = exp(tempval(0)) / (1+exp(tempval(0)));
        //theta
        temp = emitparm * covpois.row(k).t(); // the first column is intercept
        for(i=0;i<M;i++)
            theta(i) = exp(temp(i));
        
        
        //recover omega:
        //if M=2 then assign off-diagonal with 1; else
        //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
        
        omega.zeros();
       //tpmparm is different because it is a long vector in the input
        if(M==2) {omega(0,1)=1; omega(1,0)=1;}
        else{
            for(i=0;i<M;i++){
                tempsum = 1;
                
                for(j=2;j<M;j++){
                    omega(i,j) = tpmparm(nextindex); //new
                    for(m=0;m<ncolcovomega;m++)
                        omega(i,j) += tpmparm(nextindex+m+1)*covomega(k,m);
                    //omega(i,j) = exp(omega(i,j))/(1+exp(omega(i,j)));
                    omega(i,j) = exp(omega(i,j)); //new
                    tempsum += omega(i,j); //new
                    nextindex = nextindex + ncolcovomega + 1;
                }
                
                for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //new
                if(i==0) omega(i,1) = 1 / tempsum;
                else omega(i,0) = 1 / tempsum;
            }
            
            for(i=2;i<M;i++){
                for(j=2;j<=i;j++){
                    omega(i,j-1) = omega(i,j);
                    if(i==j) omega(i,j)=0;
                }
            }
            
            
        }
        
        
        result(k,1) = curr;
        if(curr==1) result(k,0) = rzip(1, zeroprop, theta(0))(0);
        else result(k,0) = Rcpp::rpois(1, theta(curr-1))(0);
        //Rcpp::Rcout<<omega<<std::endl;
    }
    
    n = count;

    //////////////////////////////////////////////////////
    //iteration
    while(n<dim){
        
        //retrieve the parameters
        nextindex = 0;
        temp = dtparm * covp.row(n).t(); //the first column is not intercept
        //Rcpp::Rcout<<temp<<std::endl;
        
        for(i=0;i<M;i++)
            eb(i) = exp(temp(i));
        
        
        omega.zeros();
        //tpmparm is different because it is a long vector in the input
        if(M==2) {omega(0,1)=1; omega(1,0)=1;}
        else{
            for(i=0;i<M;i++){
                tempsum = 1;
                
                for(j=2;j<M;j++){
                    omega(i,j) = tpmparm(nextindex); //new
                    for(m=0;m<ncolcovomega;m++)
                        omega(i,j) += tpmparm(nextindex+m+1)*covomega(n,m);
                    //omega(i,j) = exp(omega(i,j))/(1+exp(omega(i,j)));
                    omega(i,j) = exp(omega(i,j)); //new
                    tempsum += omega(i,j); //new
                    nextindex = nextindex + ncolcovomega + 1;
                }
                
                for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //new
                if(i==0) omega(i,1) = 1 / tempsum;
                else omega(i,0) = 1 / tempsum;
            }
            
            for(i=2;i<M;i++){
                for(j=2;j<=i;j++){
                    omega(i,j-1) = omega(i,j);
                    if(i==j) omega(i,j)=0;
                }
            }
            
        }

        
        //////////////////
        //start generation
        for(m=0;m<M;m++) label(m) = m+1;
        
        
        //starting generation
        prev = result(n-1,1) - 1;
        //sample from the previous omega
        curr = multinomrand(1, M, omega.row(prev).t(), label)(0);
        
        //check dt_dist
        count = random_expbase(lden(curr-1), eb(curr-1), trunc(curr-1));
     
        ////////////////
        for(k=0;k<count and n+k<dim; k++) {
            
            //get some of the parameters in each iteration
            nextindex = 0;
            //zeroprops
            tempval = zeroparm.t() * covp1.row(n+k).t(); //the first column is intercept
            zeroprop = exp(tempval(0)) / (1+exp(tempval(0)));
            //theta
            temp = emitparm * covpois.row(n+k).t(); // the first column is intercept
            for(i=0;i<M;i++)
                theta(i) = exp(temp(i));
            
            
            //recover omega:
            //if M=2 then assign off-diagonal with 1; else
            //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
            
            omega.zeros();
            //tpmparm is different because it is a long vector in the input
            if(M==2) {omega(0,1)=1; omega(1,0)=1;}
            else{
                for(i=0;i<M;i++){
                    tempsum = 1;
                    
                    for(j=2;j<M;j++){
                        omega(i,j) = tpmparm(nextindex); //new
                        for(m=0;m<ncolcovomega;m++)
                            omega(i,j) += tpmparm(nextindex+m+1)*covomega(n+k,m);
                        //omega(i,j) = exp(omega(i,j))/(1+exp(omega(i,j)));
                        omega(i,j) = exp(omega(i,j)); //new
                        tempsum += omega(i,j); //new
                        nextindex = nextindex + ncolcovomega + 1;
                    }
                    
                    for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //new
                    if(i==0) omega(i,1) = 1 / tempsum;
                    else omega(i,0) = 1 / tempsum;
                }
                
                for(i=2;i<M;i++){
                    for(j=2;j<=i;j++){
                        omega(i,j-1) = omega(i,j);
                        if(i==j) omega(i,j)=0;
                    }
                }
            
            }
            
            /////////////////////
            result(n+k,1) = curr;
            if(curr==1) result(n+k,0) = rzip(1, zeroprop, theta(0))(0);
            else result(n+k,0) = Rcpp::rpois(1, theta(curr-1))(0);

        }
        
        n += count;
    }
    
    return(result);
    
}

////////////////////
// [[Rcpp::export]]
arma::mat hsmm_cov_loglik_aft(arma::vec y, arma::vec delta, arma::vec lden, arma::mat dtparm,
                              arma::vec tpmparm, arma::vec zeroparm, arma::mat emitparm,
                              arma::vec trunc, int M, int ncolcovomega,
                              arma::mat covp,arma::mat covomega,
                              arma::mat covp1, arma::mat covpois){
  
    arma::vec loglik;
    double tempsum;
    arma::vec zeroprop(M);
    long dim = y.n_rows;
    int i,j,k,m,nextindex;
    int dmcol = arma::max(trunc); //get the maximum
    int totalmv = arma::sum(trunc);
    arma::vec theta(M);
    arma::vec newtheta(totalmv);
    arma::vec newpi(totalmv);
    arma::vec newzeroprop(totalmv);
    //Rcpp::Rcout<<dmcol<<std::endl;
    arma::mat dm(M,dmcol);
    
    arma::mat omega(M,M);
    arma::mat a(totalmv,totalmv);
    
    arma::vec forward(totalmv);
    arma::vec meanvec(totalmv);
    arma::vec forwardsum;  //forwardsum is a vec, forwardsum(0) is double
    arma::mat tempmat(1,totalmv);
    arma::vec temp(M);
    arma::vec eb(M);
    arma::vec tempval(1);
    
    //retrieve the full parameters
    
    temp = dtparm * covp.row(0).t(); //the first column is not intercept
    for(i=0;i<M;i++)
        eb(i) = exp(temp(i));

    //dwell time
    dm.zeros();
    //check dt_dist
    for(i=0;i<M;i++){
            for(j=0;j<trunc(i);j++){
                dm(i,j) = pmf_expbase(lden(i),eb(i),j+1);
            }
        }

    
    //zero proportions
    zeroprop.zeros();
    tempval = zeroparm.t() * covp1.row(0).t(); //the first column is intercept
    zeroprop(0) = exp(tempval(0)) / (1+exp(tempval(0)));
    
    //theta
    temp = emitparm * covpois.row(0).t(); // the first column is intercept
    for(i=0;i<M;i++)
        theta(i) = exp(temp(i));
    
    //get newlogtheta,newpi
    tempsum = 0;
    for(i=0;i<M;i++){
        for(j=0;j<trunc[i];j++){
            newtheta(tempsum+j) = theta(i);
            newpi(tempsum+j) = delta(i)/trunc(i);
            newzeroprop(tempsum+j) = zeroprop(i);
            
        }
        tempsum += trunc(i);
    }
    
    //recover omega:
    //if M=2 then assign off-diagonal with 1; else
    //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
    
    nextindex = 0;
    omega.zeros();
    if(M==2) {omega(0,1)=1; omega(1,0)=1;}
    else{
        for(i=0;i<M;i++){
            tempsum = 1;
            for(j=2;j<M;j++){
                omega(i,j) = tpmparm(nextindex); //new
                for(m=0;m<ncolcovomega;m++)
                    omega(i,j) += tpmparm(nextindex+m+1)*covomega(0,m);
                //omega(i,j) = exp(omega(i,j))/(1+exp(omega(i,j)));
                omega(i,j) = exp(omega(i,j)); //new
                tempsum += omega(i,j); //new
                nextindex = nextindex + ncolcovomega + 1;
            }
            
            for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //new
            if(i==0) omega(i,1) = 1 / tempsum;
            else omega(i,0) = 1 / tempsum;
            //tempsum = 1 - arma::sum(omega.row(i)); //prevent negative values
            //if(i==0) omega(i,1) = MAX(0,tempsum);
            //else omega(i,0) = MAX(0,tempsum);
        }
        
        for(i=2;i<M;i++){
            for(j=2;j<=i;j++){
                omega(i,j-1) = omega(i,j);
                if(i==j) omega(i,j)=0;
            }
        }
        
        
    }  //end of poisson and zerofl part
    
    a = hsmm_hmm (omega, dm, trunc);
    
    //initialize the forward variable
    for(m=0; m<totalmv; m++) meanvec(m) = dzip(newzeroprop(m), newtheta(m), y(0), FALSE);
    
    forward = newpi % meanvec; //element wise multiplication
    forwardsum = newpi.t() * meanvec;
    loglik = log(forwardsum);
    //Rcpp::Rcout << "loglik=" << loglik << std::endl;
    
    for(m=0; m<totalmv; m++) forward(m) = forward(m) / forwardsum(0);
    
    //recursion
    for(k=1; k<dim; k++){
        
        temp = dtparm * covp.row(k).t(); //the first column is not intercept
        for(i=0;i<M;i++)
            eb(i) = exp(temp(i));
        
        //dwell time
        dm.zeros();
        //check dt_dist
        for(i=0;i<M;i++){
            for(j=0;j<trunc(i);j++){
                dm(i,j) = pmf_expbase(lden(i),eb(i),j+1);
            }
        }
        
        //zero proportions
        zeroprop.zeros();
        tempval = zeroparm.t() * covp1.row(k).t(); //the first column is intercept
        zeroprop(0) = exp(tempval(0)) / (1+exp(tempval(0)));
        
        //theta
        temp = emitparm * covpois.row(k).t(); // the first column is intercept
        for(i=0;i<M;i++)
            theta(i) = exp(temp(i));
        
        //Rcpp::Rcout << nextindex << std::endl;
        //get newlogtheta,newpi
        tempsum = 0;
        for(i=0;i<M;i++){
            for(j=0;j<trunc[i];j++){
                newtheta(tempsum+j) = theta(i);
                newpi(tempsum+j) = delta(i)/trunc(i);
                newzeroprop(tempsum+j) = zeroprop(i);
                //Rcpp::Rcout << "theta=" << newtheta(tempsum+j) << "pi=" << newpi(tempsum+j) << std::endl;
            }
            tempsum += trunc(i);
        }
        
        //recover omega:
        //if M=2 then assign off-diagonal with 1; else
        //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
        
        nextindex = 0;
        omega.zeros();
        if(M==2) {omega(0,1)=1; omega(1,0)=1;}
        else{
            for(i=0;i<M;i++){
                tempsum = 1;
                for(j=2;j<M;j++){
                    omega(i,j) = tpmparm(nextindex); //new
                    for(m=0;m<ncolcovomega;m++)
                        omega(i,j) += tpmparm(nextindex+m+1)*covomega(k,m);
                    //omega(i,j) = exp(omega(i,j))/(1+exp(omega(i,j)));
                    omega(i,j) = exp(omega(i,j)); //new
                    tempsum += omega(i,j); //new
                    nextindex = nextindex + ncolcovomega + 1;
                }
                
                for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //new
                if(i==0) omega(i,1) = 1 / tempsum;
                else omega(i,0) = 1 / tempsum;
                //tempsum = 1 - arma::sum(omega.row(i)); //prevent negative values
                //if(i==0) omega(i,1) = MAX(0,tempsum);
                //else omega(i,0) = MAX(0,tempsum);
            }
            
            for(i=2;i<M;i++){
                for(j=2;j<=i;j++){
                    omega(i,j-1) = omega(i,j);
                    if(i==j) omega(i,j)=0;
                }
            }
            
            
        }  //end of poisson and zerofl part
        
        a = hsmm_hmm (omega, dm, trunc);
        
        
        for(m=0; m<totalmv; m++) meanvec(m) = dzip(newzeroprop(m), newtheta(m), y(k), FALSE);
        
        tempmat = forward.t() * a; //row matrix
        forward = tempmat.t() % meanvec;
        forwardsum = tempmat * meanvec;
        loglik += log(forwardsum);
        for(m=0; m<totalmv; m++) forward(m) = forward(m) / forwardsum(0);
    }
    
      return(loglik);
}

////////////////////
// [[Rcpp::export]]
Rcpp::List retrieve_hsmm_aft(arma::vec allparm,
                             arma::vec trunc, int M, int ncolcovomega,
                             int ncolcovp, int ncolcovp1, int ncolcovpois){
    
    int i,j,nextindex = 0;
    double tempsum;
    arma::vec delta(M);
    arma::vec lden(M);
    int ncoldt, ncoltpm, ncolzero, ncolemit;
    
    if(ncolcovp==0) ncoldt = 1;
    else ncoldt = ncolcovp;
    
    ncoltpm = ncolcovomega + 1;
    
    ncolzero = ncolcovp1+1;
    
    ncolemit = ncolcovpois+1;
    
    arma::mat dtparm(M,ncoldt);
    arma::vec tpmparm(ncoltpm * M * (M-2));
    arma::vec zeroparm(ncolzero);
    arma::mat emitparm(M, ncolemit);
    
    for(i=0;i<M;i++){
        lden(i) = exp(allparm(i));
    }
    nextindex = M;
 
    dtparm.zeros();
    if(ncolcovp>0){
        for(i=0;i<M;i++){
            for(j=0; j<ncolcovp; j++){
                dtparm(i,j) = allparm(nextindex + ncolcovp*i + j);
            }
        }
    }
    nextindex += ncolcovp * M;
    
    delta(0) = 1;
    tempsum = 1;
    for(i=1; i<M; i++) {
        delta(i) = exp(allparm(nextindex+i-1));
        tempsum += delta(i);
    }
    for(i=0; i<M; i++) delta(i) /= tempsum;
    nextindex += M-1;
  
    zeroparm.zeros();
    for(i=0; i<ncolzero; i++)  zeroparm(i) = allparm(nextindex+i) ;
    nextindex += ncolzero;
 
    emitparm.zeros();
    for(i=0; i<M; i++){
        for(j=0; j<ncolemit; j++){
            emitparm(i,j) = allparm(nextindex + i*ncolemit + j);
        }
    }
    nextindex += ncolemit * M;

    tpmparm.zeros();
    if(M>2){
        for(i=0;i<M*(M-2);i++){
            for(j=0; j<ncoltpm; j++){
                tpmparm(i*ncoltpm + j) = allparm(nextindex + i*ncoltpm + j);
            }
        }
    }
    
    return Rcpp::List::create(Rcpp::Named("delta")=delta,
                              Rcpp::Named("lden")=lden,
                              Rcpp::Named("dtparm")=dtparm,
                              Rcpp::Named("tpmparm")=tpmparm,
                              Rcpp::Named("zeroparm")=zeroparm,
                              Rcpp::Named("emitparm")=emitparm);
}


////////////////////
// [[Rcpp::export]]
double hsmm_cov_nllk_aft(arma::vec y, arma::vec allparm,
                              arma::vec trunc, int M, int ncolcovomega,
                              int ncolcovp, int ncolcovp1, int ncolcovpois,
                              arma::mat covp, arma::mat covomega,
                              arma::mat covp1, arma::mat covpois){
    //order in parm:  dtrate, dtparm, delta, zeroparm,  emitparm, tpmparm
 
    Rcpp::List mod = retrieve_hsmm_aft(allparm,trunc, M, ncolcovomega,
                                       ncolcovp, ncolcovp1, ncolcovpois);
    arma::vec delta = mod("delta");
    arma::vec lden = mod("lden");
    arma::mat dtparm = mod("dtparm");
    arma::vec tpmparm = mod("tpmparm");
    arma::vec zeroparm = mod("zeroparm");
    arma::mat emitparm = mod("emitparm");
    
    arma::mat loglik = hsmm_cov_loglik_aft(y, delta, lden, dtparm,
                        tpmparm, zeroparm, emitparm, trunc, M, ncolcovomega,
                                        covp,covomega, covp1, covpois);
    return -loglik(0);
}


/////////////////////////////////////////////////////////
/////////////////////////

// [[Rcpp::export]]
arma::vec hsmm_cov_viterbi_exp(arma::vec y,arma::vec trunc,int M,
                               arma::vec delta, arma::vec lden, arma::mat dtparm, arma::vec tpmparm,
                               arma::vec zeroparm, arma::mat emitparm,  int ncolcovomega,
                               arma::mat covp,arma::mat covomega,
                               arma::mat covp1, arma::mat covpois){

    double tempsum, colmax;
    arma::vec zeroprop(M);
    long dim = y.n_rows;
    int i,j,m,n,nextindex;
    int dmcol = arma::max(trunc); //get the maximum
    int totalmv = arma::sum(trunc);
    arma::vec theta(M);
    arma::vec newtheta(totalmv);
    arma::vec newpi(totalmv);
    arma::vec newzeroprop(totalmv);
    //Rcpp::Rcout<<dmcol<<std::endl;
    arma::mat dm(M,dmcol);
    
    arma::mat omega(M,M);
    arma::mat a(totalmv,totalmv);
    
    arma::vec forward(totalmv);
    arma::vec meanvec(totalmv);
    arma::vec forwardsum;  //forwardsum is a vec, forwardsum(0) is double
    arma::vec temp(M);
    arma::vec eb(M);
    arma::vec tempval(1);
    
    arma::vec state(dim);
    arma::vec tempmax(totalmv);
    arma::mat xi(dim, totalmv);
    arma::mat tempmat(totalmv,totalmv);
    
    //retrieve the full parameters
    
    temp = dtparm * covp.row(0).t(); //the first column is not intercept
    for(i=0;i<M;i++)
        eb(i) = exp(temp(i));
    
    //dwell time
    dm.zeros();
    //check dt_dist
    for(i=0;i<M;i++){
        for(j=0;j<trunc(i);j++){
            dm(i,j) = pmf_expbase(lden(i),eb(i),j+1);
        }
    }
    
    
    //zero proportions
    zeroprop.zeros();
    tempval = zeroparm.t() * covp1.row(0).t(); //the first column is intercept
    zeroprop(0) = exp(tempval(0)) / (1+exp(tempval(0)));
    
    //theta
    temp = emitparm * covpois.row(0).t(); // the first column is intercept
    for(i=0;i<M;i++)
        theta(i) = exp(temp(i));
    
    //get newlogtheta,newpi
    tempsum = 0;
    for(i=0;i<M;i++){
        for(j=0;j<trunc[i];j++){
            newtheta(tempsum+j) = theta(i);
            newpi(tempsum+j) = delta(i)/trunc(i);
            newzeroprop(tempsum+j) = zeroprop(i);
            
        }
        tempsum += trunc(i);
    }
    
    //recover omega:
    //if M=2 then assign off-diagonal with 1; else
    //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
    
    nextindex = 0;
    omega.zeros();
    if(M==2) {omega(0,1)=1; omega(1,0)=1;}
    else{
        for(i=0;i<M;i++){
            tempsum = 1;
            for(j=2;j<M;j++){
                omega(i,j) = tpmparm(nextindex); //new
                for(m=0;m<ncolcovomega;m++)
                    omega(i,j) += tpmparm(nextindex+m+1)*covomega(0,m);
                //omega(i,j) = exp(omega(i,j))/(1+exp(omega(i,j)));
                omega(i,j) = exp(omega(i,j)); //new
                tempsum += omega(i,j); //new
                nextindex = nextindex + ncolcovomega + 1;
            }
            
            for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //new
            if(i==0) omega(i,1) = 1 / tempsum;
            else omega(i,0) = 1 / tempsum;
            //tempsum = 1 - arma::sum(omega.row(i)); //prevent negative values
            //if(i==0) omega(i,1) = MAX(0,tempsum);
            //else omega(i,0) = MAX(0,tempsum);
        }
        
        for(i=2;i<M;i++){
            for(j=2;j<=i;j++){
                omega(i,j-1) = omega(i,j);
                if(i==j) omega(i,j)=0;
            }
        }
        
        
    }  //end of poisson and zerofl part
    
    a = hsmm_hmm (omega, dm, trunc);

    
    //Rcpp::Rcout << nextindex << std::endl;
    
    
    //initialize the forward variable
    
    for(m=0; m<totalmv; m++) meanvec(m) = dzip(newzeroprop(m),newtheta(m),y(0),FALSE);
    
    
    forward = newpi % meanvec; //element wise multiplication
    forwardsum = newpi.t() * meanvec;
    
    //Rcpp::Rcout << "loglik=" << loglik << std::endl;
    
    for(m=0; m<totalmv; m++) xi(0,m) = forward(m) / forwardsum(0);
    
    
    
    //recursion
    for(n=1; n<dim; n++){
        
        nextindex = 0;
        //dwell time
        temp = dtparm * covp.row(n).t(); //the first column is not intercept
        for(i=0;i<M;i++)
            eb(i) = exp(temp(i));
        
        //dwell time
        dm.zeros();
        //check dt_dist
        for(i=0;i<M;i++){
            for(j=0;j<trunc(i);j++){
                dm(i,j) = pmf_expbase(lden(i),eb(i),j+1);
            }
        }
        
        //zero proportions
        zeroprop.zeros();
        tempval = zeroparm.t() * covp1.row(n).t(); //the first column is intercept
        zeroprop(0) = exp(tempval(0)) / (1+exp(tempval(0)));
        
        //theta
        temp = emitparm * covpois.row(n).t(); // the first column is intercept
        for(i=0;i<M;i++)
            theta(i) = exp(temp(i));
        
        //Rcpp::Rcout << nextindex << std::endl;
        //get newlogtheta,newpi
        tempsum = 0;
        for(i=0;i<M;i++){
            for(j=0;j<trunc[i];j++){
                newtheta(tempsum+j) = theta(i);
                newpi(tempsum+j) = delta(i)/trunc(i);
                newzeroprop(tempsum+j) = zeroprop(i);
                //Rcpp::Rcout << "theta=" << newtheta(tempsum+j) << "pi=" << newpi(tempsum+j) << std::endl;
            }
            tempsum += trunc(i);
        }
        
        //recover omega:
        //if M=2 then assign off-diagonal with 1; else
        //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
        
        nextindex = 0;
        omega.zeros();
        if(M==2) {omega(0,1)=1; omega(1,0)=1;}
        else{
            for(i=0;i<M;i++){
                tempsum = 1;
                for(j=2;j<M;j++){
                    omega(i,j) = tpmparm(nextindex); //new
                    for(m=0;m<ncolcovomega;m++)
                        omega(i,j) += tpmparm(nextindex+m+1)*covomega(n,m);
                    //omega(i,j) = exp(omega(i,j))/(1+exp(omega(i,j)));
                    omega(i,j) = exp(omega(i,j)); //new
                    tempsum += omega(i,j); //new
                    nextindex = nextindex + ncolcovomega + 1;
                }
                
                for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //new
                if(i==0) omega(i,1) = 1 / tempsum;
                else omega(i,0) = 1 / tempsum;
                //tempsum = 1 - arma::sum(omega.row(i)); //prevent negative values
                //if(i==0) omega(i,1) = MAX(0,tempsum);
                //else omega(i,0) = MAX(0,tempsum);
            }
            
            for(i=2;i<M;i++){
                for(j=2;j<=i;j++){
                    omega(i,j-1) = omega(i,j);
                    if(i==j) omega(i,j)=0;
                }
            }
            
            
        }  //end of poisson and zerofl part
        
        a = hsmm_hmm (omega, dm, trunc);
        
        
        for(m=0; m<totalmv; m++) meanvec(m) = dzip(newzeroprop(m),newtheta(m),y(n),FALSE);
        
        //difficult part
        for(i=0; i<totalmv; i++){
            for(j=0;j<totalmv;j++){
                tempmat(i,j) = xi(n-1, i) * a(i,j);
            }
        }
        
        for(j=0; j<totalmv; j++){
            for(i=0;i<totalmv;i++){
                if(i==0) colmax = tempmat(i,j);
                else colmax = MAX(colmax, tempmat(i,j));
            }
            tempmax(j) = colmax;
        }
        
        
        forward = tempmax % meanvec;
        forwardsum = tempmax.t() * meanvec;
        for(m=0; m<totalmv; m++) xi(n,m) = forward(m) / forwardsum(0);
    }
    
    //termination
    for(m=0; m<totalmv; m++){
        if(m==0) {
            colmax = xi(dim-1, m);
            state(dim-1) = m;
        }
        else{
            if(xi(dim-1,m) >= colmax) {
                state(dim-1) = m;
                colmax = xi(dim-1,m);
            }
        }
    }
    
    //backtracking
    for(n=dim-2; n>=0; n--){
        
        
        //retrieve the full parameters
        nextindex = 0;
        //dwell time
        temp = dtparm * covp.row(n).t(); //the first column is not intercept
        for(i=0;i<M;i++)
            eb(i) = exp(temp(i));
        
        //dwell time
        dm.zeros();
        //check dt_dist
        for(i=0;i<M;i++){
            for(j=0;j<trunc(i);j++){
                dm(i,j) = pmf_expbase(lden(i),eb(i),j+1);
            }
        }
        
        //zero proportions
        zeroprop.zeros();
        tempval = zeroparm.t() * covp1.row(n).t(); //the first column is intercept
        zeroprop(0) = exp(tempval(0)) / (1+exp(tempval(0)));
        
        //theta
        temp = emitparm * covpois.row(n).t(); // the first column is intercept
        for(i=0;i<M;i++)
            theta(i) = exp(temp(i));
        
        //Rcpp::Rcout << nextindex << std::endl;
        //get newlogtheta,newpi
        tempsum = 0;
        for(i=0;i<M;i++){
            for(j=0;j<trunc[i];j++){
                newtheta(tempsum+j) = theta(i);
                newpi(tempsum+j) = delta(i)/trunc(i);
                newzeroprop(tempsum+j) = zeroprop(i);
                //Rcpp::Rcout << "theta=" << newtheta(tempsum+j) << "pi=" << newpi(tempsum+j) << std::endl;
            }
            tempsum += trunc(i);
        }
        
        //recover omega:
        //if M=2 then assign off-diagonal with 1; else
        //1.fill columns 3 ~ M (rectangles); 2. shift the lower triangle to the left; 3.reassigne diagonal=0;
        
        nextindex = 0;
        omega.zeros();
        if(M==2) {omega(0,1)=1; omega(1,0)=1;}
        else{
            for(i=0;i<M;i++){
                tempsum = 1;
                for(j=2;j<M;j++){
                    omega(i,j) = tpmparm(nextindex); //new
                    for(m=0;m<ncolcovomega;m++)
                        omega(i,j) += tpmparm(nextindex+m+1)*covomega(n,m);
                    //omega(i,j) = exp(omega(i,j))/(1+exp(omega(i,j)));
                    omega(i,j) = exp(omega(i,j)); //new
                    tempsum += omega(i,j); //new
                    nextindex = nextindex + ncolcovomega + 1;
                }
                
                for(j=2; j<M; j++) omega(i,j) = omega(i,j) / tempsum; //new
                if(i==0) omega(i,1) = 1 / tempsum;
                else omega(i,0) = 1 / tempsum;
                //tempsum = 1 - arma::sum(omega.row(i)); //prevent negative values
                //if(i==0) omega(i,1) = MAX(0,tempsum);
                //else omega(i,0) = MAX(0,tempsum);
            }
            
            for(i=2;i<M;i++){
                for(j=2;j<=i;j++){
                    omega(i,j-1) = omega(i,j);
                    if(i==j) omega(i,j)=0;
                }
            }
            
            
        }  //end of poisson and zerofl part
        
        a = hsmm_hmm (omega, dm, trunc);
        
        
        for(m=0;m<totalmv;m++){
            if(m==0) {
                colmax = xi(n,m) * a(m, state(n+1));
                state(n) = m;
            }
            else{
                if(xi(n,m) * a(m, state(n+1))>=colmax) {
                    state(n) = m;
                    colmax = xi(n,m) * a(m, state(n+1));
                }
            }
        }
    }
    
    for(n=0; n<dim; n++) state(n) = state(n) + 1;
    
    for(i=0;i<dim;i++){
        tempsum=0;
        for(j=0;j<M;j++){
            if(state(i)<=tempsum+trunc(j)){
                state(i)=j+1;
                break;
            }
            tempsum += trunc(j);
        }
    }
    
    
    return(state);
}


/////////////////////////////////////////////////////////////////////////////////////////
//////////add covariates in the transition rate of continuous time hidden Markov model
///////////////////////////////////////////
//[[Rcpp::export]]
arma::mat hmm_cov_gen_cont3 (arma::vec parm, int M, long dim, int ncolcovpi, arma::mat covpi,
                            int ncolcovtrans, arma::mat covtrans, int ncolcovp1, arma::mat covp1,
                            int ncolcovpois, arma::mat covpois, arma::vec zeroindex, arma::vec timeindex){
    
    int i,j,m,n,nextindex,prev,curr;
    double tempsum, u, interval;
    arma::vec label(M);
    
    arma::mat result(dim, 2); /*first column for x, second column for state*/
    
    arma::mat transit(M,M);
    arma::vec pi(M);
    arma::mat a(M,M);
    arma::vec theta(M);
    arma::vec zeroprop(M);
    
    //retrieve the full parameters
    //prior
    nextindex=0;
    pi(0) = 1;
    tempsum = 1;
    for(m=1; m<M; m++){
        pi(m) = parm(nextindex);
        
        for(j=0;j<ncolcovpi;j++)
            pi(m) += parm(nextindex+j+1)*covpi(0,j) ;
        
        pi(m) = exp(pi(m));
        
        tempsum += pi(m);
        nextindex = nextindex + ncolcovpi + 1;
    }
    
    for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
    
    
    //transition
    for(i=0; i<M; i++){
        tempsum = 0;
        for(j=0; j<M; j++){
            if(i!=j){
                a(i,j) = parm(nextindex);
                for(m=0; m<ncolcovtrans; m++)
                    a(i,j) += parm(nextindex + m + 1) * covtrans(0,m);
                //Rcpp::Rcout << "aij=" << a(i,j) << std::endl;
                a(i,j) = exp(a(i,j))/(1+exp(a(i,j)));
                tempsum += a(i,j);
                nextindex = nextindex + ncolcovtrans + 1;
                //a(i,j) = exp(parm(nextindex))/(1+exp(parm(nextindex)));
                //tempsum += a(i,j);
                //nextindex++;
            }
        }
        a(i,i) = -tempsum;
        
    }
    //Rcpp::Rcout << "rate" << a << std::endl;
    
    for(i=0;i<M;i++){
        if(zeroindex(i)==0) zeroprop(i)=0;
        else{
            zeroprop(i) = parm(nextindex);
            
            for(m=0;m<ncolcovp1;m++)
                zeroprop(i) += parm(nextindex+m+1) * covp1(0,m);
            
            zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
            nextindex = nextindex + ncolcovp1 + 1;
        }
    }
    
    //Rcpp::Rcout << "p1=" << p1 << std::endl;
    
    for(m=0; m<M; m++){
        theta(m) = parm(nextindex);
        for(j=0; j<ncolcovpois;j++){
            theta(m) += parm(nextindex+j+1) * covpois(0,j);
        }
        theta(m) = exp(theta(m));
        nextindex = nextindex + ncolcovpois + 1;
    }
    
    
    //start generation
    
    for(m=0;m<M;m++) label(m) = m+1;
    
    result(0,1) = multinomrand(1, M, pi, label)(0);
    curr = result(0,1) - 1;
    u = Rcpp::runif(1,0,1)(0);
    if(u<=zeroprop(curr)) result(0,0) = 0;
    else result(0,0) = Rcpp::rpois(1, theta(curr))(0);
    
    
    
    //iteration steps
    for(n=1; n<dim; n++){
        
        //still need to retrieve the natural parameters
        nextindex=0;
        pi(0) = 1;
        tempsum = 1;
        for(m=1; m<M; m++){
            pi(m) = parm(nextindex);
            
            for(j=0;j<ncolcovpi;j++)
                pi(m) += parm(nextindex+j+1)*covpi(n,j) ;
            
            pi(m) = exp(pi(m));
            
            tempsum += pi(m);
            nextindex = nextindex + ncolcovpi + 1;
        }
        
        for(m=0; m<M; m++) pi(m) = pi(m) / tempsum;
        
        
        //transition
        for(i=0; i<M; i++){
            tempsum = 0;
            for(j=0; j<M; j++){
                if(i!=j){
                    a(i,j) = parm(nextindex);
                    for(m=0; m<ncolcovtrans; m++)
                        a(i,j) += parm(nextindex + m + 1) * covtrans(n,m);
                    a(i,j) = exp(a(i,j))/(1+exp(a(i,j)));
                    tempsum += a(i,j);
                    nextindex = nextindex + ncolcovtrans + 1;
                    //a(i,j) = exp(parm(nextindex))/(1+exp(parm(nextindex)));
                    //tempsum += a(i,j);
                    //nextindex++;
                }
            }
            a(i,i) = -tempsum;
        }
        
        for(i=0;i<M;i++){
            if(zeroindex(i)==0) zeroprop(i)=0;
            else{
                zeroprop(i) = parm(nextindex);
                
                for(m=0;m<ncolcovp1;m++)
                    zeroprop(i) += parm(nextindex+m+1) * covp1(n,m);
                
                zeroprop(i) = exp(zeroprop(i))/(1+exp(zeroprop(i)));
                nextindex = nextindex + ncolcovp1 + 1;
            }
        }
        //Rcpp::Rcout << "p1=" << p1 << std::endl;
        
        for(m=0; m<M; m++){
            theta(m) = parm(nextindex);
            for(j=0; j<ncolcovpois;j++){
                theta(m) += parm(nextindex+j+1) * covpois(n,j);
            }
            theta(m) = exp(theta(m));
            nextindex = nextindex + ncolcovpois + 1;
        }
        
        interval = timeindex(n) - timeindex(n-1);
        transit = matrixexp(a, interval);
        //Rcpp::Rcout << "tpm=" << transit << std::endl;
        prev = result(n-1, 1) - 1;
        result(n,1) = multinomrand(1, M, transit.row(prev).t(), label)(0);  //row to column vetor
        
        curr = result(n,1) - 1;
        u = Rcpp::runif(1,0,1)(0);
        if(u<=zeroprop(curr)) result(n,0) = 0;
        else result(n,0) = Rcpp::rpois(1, theta(curr))(0);
        
    }
    
    return result;
    
}


////////////////////
//////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double ziploglik_cov_cont3(arma::vec delta, arma::mat gammaparm, arma::vec thetaparm, arma::mat lambdaparm,
                          arma::vec y, arma::mat x, arma::vec ntimes,
                          arma::vec timeindex){
    //the first column of x should be 1 for intercept
    double loglik = 0;
    //long dim = y.n_rows;
    int M = delta.n_rows;
    int n = ntimes.n_rows;
    int j,k,m,t;
    arma::vec tempval;
    double tempsum;
    double theta;
    arma::vec lambda(M);
    
    arma::mat rate(M,M);
    arma::mat transit(M,M);
    int gammarow;
    double tempgammasum;
    double interval;
    
    arma::vec forward(M);
    arma::vec meanvec(M);
    arma::vec forwardsum;  //then forwardsum(0) is double
    arma::mat tempmat(1,M);
    
    int count = 0;
    
    for(j=0; j<n; j++){
        
        tempsum = 0;
        
        //initialize the forward variable
        tempval = x.row(count) * thetaparm;
        theta = exp(tempval(0))/(1+exp(tempval(0)));
        for(m=0;m<M;m++) {
            tempval = x.row(count) * lambdaparm.row(m).t();
            lambda(m) = exp(tempval(0));
        }
        
        //the gammaparm is stacked
        gammarow = 0;
        for(k=0; k<M; k++){
            tempgammasum = 0;
            for(m=0; m<M; m++){
                if(k!=m){
                    tempval = x.row(count) * gammaparm.row(gammarow).t();
                    rate(k,m) = exp(tempval(0))/(1+exp(tempval(0)));
                    
                    tempgammasum += rate(k,m);
                    gammarow += 1;
                }
                
            }
            rate(k,k) = -tempgammasum;
        }
        //Rcpp::Rcout << rate << std::endl;
        
        meanvec(0) = dzip(theta, lambda(0), y(count), FALSE);
        for(m=1; m<M; m++) meanvec(m) = R::dpois(y(count), lambda(m), FALSE);
        
        forward = delta % meanvec; //element wise multiplication
        forwardsum = delta.t() * meanvec;
        tempsum = log(forwardsum(0));
        //Rcpp::Rcout << "loglik=" << loglik << std::endl;
        
        for(m=0; m<M; m++) forward(m) = forward(m) / forwardsum(0);
        
        //recursion
        for(t=1; t<ntimes(j); t++){
            
            tempval = x.row(count+t) * thetaparm;
            theta = exp(tempval(0))/(1+exp(tempval(0)));
            for(m=0;m<M;m++) {
                tempval = x.row(count+t) * lambdaparm.row(m).t();
                lambda(m) = exp(tempval(0));
            }
            
            //the gammaparm is stacked
            gammarow = 0;
            for(k=0; k<M; k++){
                tempgammasum = 0;
                for(m=0; m<M; m++){
                    if(k!=m){
                        tempval = x.row(count+t) * gammaparm.row(gammarow).t();
                        rate(k,m) = exp(tempval(0))/(1+exp(tempval(0)));
                        tempgammasum += rate(k,m);
                        gammarow += 1;
                    }
                    
                }
                rate(k,k) = -tempgammasum;
            }
            //Rcpp::Rcout << rate << std::endl;
            
            //forward variables
            meanvec(0) = dzip(theta, lambda(0), y(count+t), FALSE);
            for(m=1; m<M; m++) meanvec(m) = R::dpois(y(count+t), lambda(m), FALSE);

            interval = timeindex(count+t) - timeindex(count+t-1);
            transit = matrixexp(rate, interval);
            
            tempmat = forward.t() * transit; //row matrix
            forward = tempmat.t() % meanvec;
            forwardsum = tempmat * meanvec;
            tempsum += log(forwardsum(0));
            for(m=0; m<M; m++) forward(m) = forward(m) / forwardsum(0);
        }
        
        loglik += tempsum;
        count += ntimes(j);
        
    }
    
    return loglik;
    
}

//////////////////////////////////////////////////////////////////////
//[[Rcpp::export]]
arma::cube getallexpm3(int M, arma::mat gammaparm, arma::mat tpmcov, arma::vec timeindex){
    int dim = timeindex.n_rows - 1;
 
    arma::cube result(M,M,dim);
    int i,k,m,gammarow;
    double tempgammasum;
    arma::vec tempval;
    arma::mat rate(M,M);
    
    for(i=0; i<dim; i++) {
        //Rcpp::Rcout << "i=" << i << std::endl;
        gammarow = 0;
        for(k=0; k<M; k++){
            tempgammasum = 0;
            for(m=0; m<M; m++){
                if(k!=m){
                    //covariate starting from the second row
                    tempval = tpmcov.row(i+1) * gammaparm.row(gammarow).t();
                    rate(k,m) = exp(tempval(0))/(1+exp(tempval(0)));
                    tempgammasum += rate(k,m);
                    gammarow += 1;
                }
                
            }
            rate(k,k) = -tempgammasum;
        }
        
        result.slice(i) = matrixexp(rate, timeindex(i+1)-timeindex(i));

    }
    return result;
}


////////////////////////////////
// [[Rcpp::export]]
arma::mat getnodeprob_cov_cont(arma::vec y, arma::mat x, arma::vec thetaparm, arma::mat lambdaparm, int m){
    long ns = y.n_rows; //original numeric series starting from 1
    arma::mat nodeprob(ns, m);
    arma::vec eta;
    double zeroprop,lambda;
    long i;
    int j;
    
    for(i=0;i<ns;i++){
        //Rcpp::Rcout<<i<<std::endl;
        for(j=0;j<m;j++){
            if(j==0){
                eta = x.row(i) * thetaparm ;
                zeroprop = exp(eta(0)) / (1+exp(eta(0)));
                eta = x.row(i) * lambdaparm.row(j).t();
                lambda = exp(eta(0));
                nodeprob(i,j) = dzip(zeroprop, lambda, y(i), FALSE);
            }
            else{
                zeroprop = 0;
                eta = x.row(i) * lambdaparm.row(j).t();
                lambda = exp(eta(0));
                nodeprob(i,j) = dzip(zeroprop, lambda, y(i), FALSE);
            }
        }
    }
    return nodeprob;
}


///////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List fb_cont3(arma::vec Pi, arma::cube expms, arma::mat nodeprob,
                   long dim, arma::vec ntimes, arma::vec timeindex){
    
    
    int M = nodeprob.n_cols;
    int n = ntimes.n_rows;
    int j,m,t,i,k;

    arma::mat transit(M,M);
    double tempsum;
    arma::vec tempval;
    arma::mat tempmat(1,M);
    arma::mat tempmat2(M,M);
    
    arma::mat alpha(dim, M);
    arma::vec scale(dim);
    arma::mat beta(dim, M);
    arma::mat Gamma(dim, M);
    /*when multiple series, rows in between will be zeros
     for safety reasons, avoid multiple series.*/
    arma::mat xi(dim-1, M*M);
    
    arma::vec colsumgamma(M);
    arma::vec tempsumxi(M*M);
    arma::mat colsumxi(M,M);
    double loglik = 0;
    
    int count = 0;
    for(j=0; j<n; j++){
        
        //forward probs
        alpha.row(count) = Pi.t() % nodeprob.row(count);
        tempval = nodeprob.row(count) * Pi;
        scale(count) = tempval(0);//double
        for(m=0; m<M; m++) alpha(count, m) = alpha(count,m) / scale(count);
        
        for(t=1; t<ntimes(j); t++){

            transit = expms.slice(t-1); //expms only length-1 row
            
            tempmat = alpha.row(count+t-1) * transit; //row matrix
            alpha.row(count+t) = tempmat % nodeprob.row(count+t);
            tempval = tempmat * nodeprob.row(count+t).t();
            scale(count+t) = tempval(0);//double
            for(m=0; m<M; m++) alpha(count+t, m) = alpha(count+t,m) / scale(count+t);
        }
        
        //backward probs and state conditional probs
        for(m=0; m<M; m++) beta(count + ntimes(j) - 1,m) = 1 / (M * scale(count + ntimes(j) - 1));
        Gamma.row(count+ntimes(j)-1) = alpha.row(count+ntimes(j)-1) % beta.row(count+ntimes(j)-1);
        tempval = alpha.row(count+ntimes(j)-1) * beta.row(count+ntimes(j)-1).t();
        for(m=0; m<M; m++) Gamma(count+ntimes(j)-1,m)=Gamma(count+ntimes(j)-1,m)/tempval(0); //double
        
        for(t=ntimes(j)-2; t>=0; t--){
            
            transit = expms.slice(t);
            tempmat = beta.row(count+t+1) % nodeprob.row(count+t+1);
            beta.row(count+t) = tempmat * transit.t();
            for(m=0; m<M; m++) beta(count+t,m) = beta(count+t,m) / scale(count + t);
            
            Gamma.row(count+t) = alpha.row(count+t) % beta.row(count+t);
            tempval = alpha.row(count+t) * beta.row(count+t).t();
            for(m=0; m<M; m++) Gamma(count+t,m)=Gamma(count+t,m)/tempval(0); //double
        }
        
        //transition probs
        for(t=0; t<ntimes(j)-1; t++){
            
            transit = expms.slice(t);
            
            tempmat = beta.row(count+t+1) % nodeprob.row(count+t+1);
            tempmat2 = transit % (alpha.row(count+t).t() * tempmat);
            
            tempsum = 0;
            for(i=0; i<M; i++){
                for(k=0; k<M; k++){
                    xi(count+t, i + k*M) = tempmat2(i,k);
                    tempsum += xi(count+t, i + k*M);
                }
            }
            //Rcpp::Rcout<<count+t<<std::endl;
            for(m=0; m<M*M; m++) xi(count+t, m) = xi(count+t, m) / tempsum;
        }
        
        count += ntimes(j);
    }
    
    
    //get the column sums
    colsumgamma = colsum(Gamma);
    tempsumxi = colsum(xi);
    
    for(i=0; i<M; i++)
        for(k=0; k<M; k++)
            colsumxi(i,k) = tempsumxi(i+k*M);
    
    loglik = arma::sum(log(scale));
    
    return Rcpp::List::create(Rcpp::Named("colsumgamma")=colsumgamma,
                              Rcpp::Named("colsumxi")=colsumxi,
                              Rcpp::Named("Gamma")=Gamma,
                              Rcpp::Named("xi")=xi,
                              Rcpp::Named("loglik")=loglik);
}


///////////////////////////////////
//' retrieve the natural parameters from the working parameters in zero-inflated Poisson
//' hidden Markov model with covariates in state-dependent parameters and transition rates
//' @param parm working parameters
//' @param M number of hidden states
//' @param ncolx number of covariates including the intercept
//' @return a list of natural parameters
//' @export
// [[Rcpp::export]]
Rcpp::List retrieve_cov_cont3(arma::vec parm, int M, int ncolx){
    int nextindex = 0;
    arma::vec delta(M);
    arma::mat gammaparm(M*(M-1),ncolx);
    arma::vec thetaparm(ncolx);
    arma::mat lambdaparm(M,ncolx);
    int i, j, m;
    double tempsum;
    
    //retrieve the full parameters
    delta(0) = 1;
    tempsum = 1;
    for(m=1; m<M; m++) {
        delta(m) = exp(parm(m-1));
        tempsum += delta(m);
    }
    for(m=0; m<M; m++) delta(m) = delta(m) / tempsum;
    
    nextindex = M-1;
    
    for(i=0; i<M*(M-1); i++){
        for(j=0; j<ncolx; j++){
            gammaparm(i,j) = parm(nextindex);
            nextindex += 1;
        }
    }
    
    for(i=0; i<ncolx; i++) thetaparm(i) = parm(nextindex+i);
    nextindex += ncolx;
    
    for(i=0; i<M; i++){
        for(j=0; j<ncolx; j++){
            lambdaparm(i,j) = parm(nextindex + i*ncolx + j);
        }
    }
    
    return Rcpp::List::create(Rcpp::Named("delta")=delta,
                              Rcpp::Named("gammaparm")=gammaparm,
                              Rcpp::Named("thetaparm")=thetaparm,
                              Rcpp::Named("lambdaparm")=lambdaparm);
    
}


//////////////////////////////////////////////////////////////////////////////////
//' negative log likelihood function for zero-inflated Poisson hidden Markov model
//' with covariates in state-dependent parameters and transition rates
//' @param parm working parameters
//' @param y observed series
//' @param covariates design matrix of covariates including the intercept
//' @param M number of hidden states
//' @param ntimes length of the observed series
//' @param timeindex vector of observed time points
//' @return negative log likelihood
//' @export
// [[Rcpp::export]]
double zipnegloglik_cov_cont3(arma::vec parm, arma::vec y, arma::mat covariates, int M, arma::vec ntimes,
                             arma::vec timeindex){
    int ncolx = covariates.n_cols;
    Rcpp::List mod = retrieve_cov_cont3(parm,M,ncolx);
    arma::vec delta = mod("delta");
    arma::mat gammaparm = mod("gammaparm");
    arma::vec thetaparm = mod("thetaparm");
    arma::mat lambdaparm = mod("lambdaparm");

    double negloglik = - ziploglik_cov_cont3(delta, gammaparm, thetaparm, lambdaparm, y, covariates, ntimes,
                                            timeindex);
    return negloglik;
}

/////////////////////////////////////////////
////////////////////////////////////////////
//////////////////////////////////////////////
//' Convolution of two real vectors of the same length.
//' @param vec1 the first vector
//' @param vec2 the second vector
//' @return a vector of full convolution
//' @export
// [[Rcpp::export]]

arma::vec convolution(arma::vec vec1, arma::vec vec2) {
    return(arma::conv(vec1, vec2));
}



/////////
/*
// a good strategy to get all unique intervals and compute all transitions first!
// the below function is wrong !!! should solve equation since pij appears on both sides
// [[Rcpp::export]]

double gamma_conv(double shape1, double scale1, double shape2, double scale2,
                     int n, double upper){
    //mean = shape * scale;
    double step = upper / n;
    double cdf = 0;
    
    for(int i=0; i<n; i++){
        cdf += R::dgamma((n-i)*step,shape1,scale1,FALSE) * R::dgamma((i+1)*step,shape2,scale2,FALSE);
    }
    //R::pgamma((n-i)*step,shape1,scale1,TRUE,FALSE)
    //Rcpp::Rcout<<Rcpp::rgamma(1, shape2, scale2)<<std::endl;
    //cdf = R::dgamma(upper, shape1, scale1, FALSE);

    return cdf/n;
}
*/






////////////////////////////////////////////
///important utility function to compute transition probabilities for GAMMA SMP
///each row of parmmat is the shape and scale parameters for the off diagonal transitions

// [[Rcpp::export]]

arma::mat smprcpp(arma::mat parmmat, int M, double time, int ngrid, double lower){
    int i,j,k,m,s,t,nextindex;
    arma::cube pp(M,M,ngrid);
    pp.zeros();
    for(i=0; i<M; i++)
        pp(i,i,0) = 1;
    
    arma::vec grids(ngrid);
    double h = (time - lower) / (ngrid - 1);
    for(i=0; i<ngrid; i++)
        grids(i) = lower + i * h;
    
    arma::cube densities(M,M,ngrid);
    densities.zeros();
    arma::cube cdfs(M,M,ngrid);
    cdfs.zeros();
    
    
    for(m=0; m<ngrid; m++){
        nextindex = 0;
        for(i=0; i<M; i++){
            for(j=0; j<M; j++){
                if(i!=j){
                    densities(i,j,m) = R::dgamma(grids(m),parmmat(nextindex,0),parmmat(nextindex,1),FALSE);
                    cdfs(i,j,m) = R::pgamma(grids(m),parmmat(nextindex,0),parmmat(nextindex,1),TRUE,FALSE);
                    nextindex += 1;
                    
                }
            }
        }
    }
    

    double tempsum;
    double part1;
    double part2;
    double temprowsum;
    
    for(m=1; m<ngrid; m++){
        for(i=0; i<M; i++){
            for(j=0; j<M; j++){
                if(i!=j){
                    pp(i,j,m) = 0;
                    
                    for(k=0; k<M; k++){
                        if(k!=i){
                            
                            tempsum = 0;
                            for(s=0; s<m-1; s++){
                                part1 = pp(k,j,m-s-2);
                                
                                temprowsum = 0;
                                for(t=0; t<M; t++){
                                    temprowsum += cdfs(i,t,s);
                                }
                                    
                                part2 = densities(i,k,s)*cdfs(i,k,s)/temprowsum;
                                    
                                tempsum += part1 * part2;
                                
                            }
                            
                            pp(i,j,m) += h*tempsum;
                            
                            /*
                            double result=0;
                            for(int i=0; i<m-1; i++){
                                result += vec1(m-i-2) * vec2(i);
                            }
                            */
                            
                        }
                    }
                    
                }
            }
        }
        
        
        for(i=0; i<M; i++){
            temprowsum = 0;
            for(t = 0; t<M; t++) temprowsum += pp(i,t,m);
            if(temprowsum<1) pp(i,i,m) = 1 - temprowsum;
            else pp(i,i,m) = 0;
        }
        
    }
    
    return pp.slice(ngrid-1);
    
}









