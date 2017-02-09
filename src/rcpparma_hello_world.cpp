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


/*
arma::mat rcpparma_hello_world(int k) {
    arma::mat m1 = arma::eye<arma::mat>(k, k);
    arma::mat m2 = arma::eye<arma::mat>(k, k);
	                     
    return m1 + 3 * (m1 + m2);
}


// another simple example: outer product of a vector,
// returning a matrix
//
// [[Rcpp::export]]
arma::mat rcpparma_outerproduct(const arma::colvec & x) {
    arma::mat m = x * x.t();
    return m;
}

// and the inner product returns a scalar
//
// [[Rcpp::export]]
double rcpparma_innerproduct(const arma::colvec & x) {
    double v = arma::as_scalar(x.t() * x);
    return v;
}


//' use Rcpp::List to return both at the same time
//' @export
// [[Rcpp::export]]
Rcpp::List rcpparma_bothproducts(const arma::colvec & x) {
    arma::mat op = x * x.t();
    double    ip = arma::as_scalar(x.t() * x);
    return Rcpp::List::create(Rcpp::Named("outer")=op,
                              Rcpp::Named("inner")=ip);
}

*/



//using namespace Rcpp;
/*
arma::vec diagonal(arma::mat a){
    arma::vec result;
    arma::mat temp;
    //Rcpp::Rcout<< MAX(NAN,2) << std::endl; //this works
    //Rcpp::Rcout<< MAX(NAN,NAN) << std::endl; //this works
    temp = arma::inv(a.t()*a);
    result = a.diag();
    return(result);
}
*/

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
    if(y==0) result = p + (1-p) * exp(-theta);
    else result = (1-p) * R::dpois(y, theta, loga);
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

//shift the poisson to the right by shift
//[[Rcpp::export]]
arma::vec rshiftpois(int n, int theta, int shift){
    arma::vec result;
    result = Rcpp::rpois(n, theta) + shift; /*everything plus one to avoid zero for dwell time dist.*/
    return(result);
}

//shift the poisson to the right by shift
//[[Rcpp::export]]
double dshiftpois(int x, int theta, int shift, bool loga){
    double result;
    int temp = x - shift;
    result = R::dpois(temp, theta, loga); /*everything plus one to avoid zero for dwell time dist.*/
    return result;
}



//[[Rcpp::export]]
double dlogp(int x, double p, bool loga){
    double result;
    if(loga==FALSE) result = -exp(x*log(p))/(x*log(1-p));
    else result = log(-exp(x*log(p))/(x*log(1-p)));
    return result;
}


//[[Rcpp::export]]
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


//[[Rcpp::export]]
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
//[[Rcpp::export]]
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

//[[Rcpp::export]]
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


///////////////////////////////////////////
//[[Rcpp::export]]
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
            
        
      }//end section for zeroinflated poisson distribution
        
        // recover transition matrix
        
        a = hsmm_hmm (omega, dm, trunc);
        
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


