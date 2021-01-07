#include <Rcpp.h>
using namespace Rcpp;

double Cquantile(NumericVector x, double q) {
    NumericVector y = clone(x);
	NumericVector shell;
	double out;
    std::sort(y.begin(), y.end());
	shell = y[x.size()*(q - 0.000000001)];
	out = shell[0];
    return out;
}


double vecmax(NumericVector x) {
    // Rcpp supports STL-style iterators
    NumericVector::iterator it = std::max_element(x.begin(), x.end());
    // we want the value so dereference
    return *it;
}

int vecminInd(NumericVector x) {
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::min_element(x.begin(), x.end());
  // we want the value so dereference 
  return it - x.begin();
}


double minentropy(NumericVector x, IntegerVector y) {
    NumericVector x1 = clone(x);
	IntegerVector y1 = clone(y);

	double x20, x80, xin; 
	x20 = Cquantile(x1,0.20);
	x80 = Cquantile(x1,0.80);
	xin = (x80- x20)/19; 
 
	int xlen = x1.size();
 
	NumericVector cutpoints(20),ginilist(20); 

	double grp1sum,grp2sum, counter1,counter2;
	 
	for (int i = 0; i < 20; i++) {
		cutpoints[i] = x20+(i*xin);
	 
		grp1sum=0;
		grp2sum=0; 
		counter1=0;
		counter2=0;
		for (int j = 0; j < xlen; j++) {

			if(x1[j] > cutpoints[i]){
				grp1sum = grp1sum + y1[j]; 
				counter1 = counter1 + 1;
				}else{
				grp2sum = grp2sum + y1[j]; 	
				counter2 = counter2 + 1;	
			}
	}
	ginilist[i] =   (counter1/xlen)*(grp1sum/counter1)*(1-(grp1sum/counter1)) + (counter2/xlen)*(grp2sum/counter2)*(1-(grp2sum/counter2));
	}
	
	int tt;
	tt=vecminInd(ginilist);
	double output =  cutpoints[tt];
    return output;
}


// [[Rcpp::export]]
NumericMatrix getpairsOddXOut(NumericMatrix MatOut, NumericMatrix MatIn, IntegerVector Yout, int pch2){
  int nrows = MatOut.nrow();
  int ncol = MatOut.ncol();
  int nrowsIN = MatIn.nrow();

  NumericMatrix out(nrowsIN, pch2);
  
  NumericVector testrow(nrows);
  NumericVector testrow1b(nrows);

  NumericVector sortedtestrow(nrows);

  NumericVector testrow2(nrowsIN);

  double averagemean;
  
  int count =0;
  for (int i = 0; i < (ncol-1); i++) {
    for(int j = (i+1);  j < ncol ; j++){

      testrow = MatOut(_,i)/MatOut(_,j);
        for (int k = 0; k < nrows; k++) {
            if(MatOut(k,j)  == 0) testrow[k]=0;
			if((MatOut(k,i)  == 0) && (MatOut(k,j)  == 0)) testrow[k]= 1;
        }
		testrow1b = testrow;
        for (int k = 0; k < nrows; k++) {
            if( (MatOut(k,i)  != 0) && (MatOut(k,j)  == 0)) testrow[k]	= (vecmax(testrow1b));
        }
    
       averagemean = minentropy(testrow, Yout);
	  
	 for (int k = 0; k < nrowsIN; k++) {
        if((MatIn(k,i)  >  (averagemean * MatIn(k,j)))) testrow2[k]= 1;
        if((MatIn(k,i)  <= (averagemean * MatIn(k,j)))) testrow2[k]= 0;
        if((MatIn(k,i)==0) && (MatIn(k,j)==0)) testrow2[k]= rbinom(1, 1, 0.5)[0];
      }
      out(_,count) = testrow2;
      count= count+1;
    }
  }
  return out;
} 
  

  // [[Rcpp::export]]
  NumericMatrix getpairsOddXOutRandomOrder(NumericMatrix MatOut, NumericMatrix MatIn, IntegerVector Yout, int pch2){
    int nrows = MatOut.nrow();
    int ncol = MatOut.ncol();
    int nrowsIN = MatIn.nrow();
    
    NumericMatrix out(nrowsIN, pch2);
    
    NumericVector testrow(nrows);
    NumericVector testrow1b(nrows);
    
    NumericVector sortedtestrow(nrows);
    
    NumericVector testrow2(nrowsIN);
    
    double averagemean;

    NumericVector OutTemp1;
    NumericVector OutTemp2;
    
    NumericVector InTemp1;
    NumericVector InTemp2;
    
    
    int count =0;
    for (int i = 0; i < (ncol-1); i++) {
      for(int j = (i+1);  j < ncol ; j++){
        
        OutTemp1= MatOut(_,i);
        OutTemp2= MatOut(_,j);
        
        InTemp1= MatIn(_,i);
        InTemp2= MatIn(_,j);
      
        int movegrp = rbinom(1, 1, 0.5)[0];
        if(movegrp ==1){
          OutTemp1= MatOut(_,j);
          OutTemp2= MatOut(_,i);
          
          InTemp1= MatIn(_,j);
          InTemp2= MatIn(_,i);
        }
        
        testrow = OutTemp1/OutTemp2;
        for (int k = 0; k < nrows; k++) {
          if(OutTemp2[k] == 0) testrow[k]=0;
          if((OutTemp1[k]  == 0) && (OutTemp2[k] == 0)) testrow[k]= 1;
        }
        testrow1b = testrow;
        for (int k = 0; k < nrows; k++) {
          if( (OutTemp1[k] != 0) && (OutTemp2[k] == 0)) testrow[k]	= (vecmax(testrow1b));
        }
        
        averagemean = minentropy(testrow, Yout);
        
        for (int k = 0; k < nrowsIN; k++) {
          if((InTemp1[k]  >  (averagemean * InTemp2[k]))) testrow2[k]= 1;
          if((InTemp1[k]  <= (averagemean * InTemp2[k]))) testrow2[k]= 0;
          if((InTemp1[k]==0) && (MatIn(k,j)==0)) testrow2[k]= rbinom(1, 1, 0.5)[0];
        }
        out(_,count) = testrow2;
        count= count+1;
      }
    }
    return out;
  } 


  
  
  
  
  
  
  
// [[Rcpp::export]]
NumericMatrix RatioX(NumericMatrix MatIn, int pch2){
  int nrows = MatIn.nrow();
  int ncol = MatIn.ncol();
 
  NumericMatrix out(nrows, pch2);
  
  NumericVector testrow(nrows);
  NumericVector sortedtestrow(nrows);

  NumericVector testrow2(nrows);

  int count =0;
  for (int i = 0; i < (ncol-1); i++) {
    for(int j = (i+1);  j < ncol ; j++){

      testrow = MatIn(_,i)/MatIn(_,j);
        for (int k = 0; k < nrows; k++) {
            if(MatIn(k,j)  == 0) testrow[k]=0;
			if((MatIn(k,i)  == 0) && (MatIn(k,j)  == 0)) testrow[k]= 1;
        }

	     testrow2 = testrow;
        for (int k = 0; k < nrows; k++) {
          if( (MatIn(k,i)  != 0) && (MatIn(k,j)  == 0)) testrow[k] = (vecmax(testrow2));
        }
    
      out(_,count) = testrow;
      count= count+1;
    }
  }
  return out;
} 
  
class Comparator {
private:
    const Rcpp::NumericVector& ref;

    bool is_na(double x) const 
    {
        return Rcpp::traits::is_na<REALSXP>(x);    
    }

public:
    Comparator(const Rcpp::NumericVector& ref_)
        : ref(ref_)
    {}

    bool operator()(const int ilhs, const int irhs) const
    {
        double lhs = ref[ilhs], rhs = ref[irhs]; 
        if (is_na(lhs)) return false;
        if (is_na(rhs)) return true;
        return lhs < rhs;
    }
};


Rcpp::NumericVector avg_rank(Rcpp::NumericVector x)
{
    R_xlen_t sz = x.size();
    Rcpp::IntegerVector w = Rcpp::seq(0, sz - 1);
    std::sort(w.begin(), w.end(), Comparator(x));

    Rcpp::NumericVector r = Rcpp::no_init_vector(sz);
    for (R_xlen_t n, i = 0; i < sz; i += n) {
        n = 1;
        while (i + n < sz && x[w[i]] == x[w[i + n]]) ++n;
        for (R_xlen_t k = 0; k < n; k++) {
            r[w[i + k]] = i + (n + 1) / 2.;
        }
    }

    return r;
}
  


  
  
// [[Rcpp::export]]
NumericMatrix RankX(NumericMatrix MatIn, int pch2){
  int nrows = MatIn.nrow();
  int ncol = MatIn.ncol();
 
  NumericMatrix out(nrows, pch2);
  
  NumericVector testrow(nrows);
  NumericVector sortedtestrow(nrows);

  NumericVector testrow2(nrows);
  NumericVector testrow3(nrows);

  int count =0;
  for (int i = 0; i < (ncol-1); i++) {
    for(int j = (i+1);  j < ncol ; j++){

      testrow = MatIn(_,i)/MatIn(_,j);
        for (int k = 0; k < nrows; k++) {
            if(MatIn(k,j)  == 0) testrow[k]=0;
			if((MatIn(k,i)  == 0) && (MatIn(k,j)  == 0)) testrow[k]= 1;
        }

	     testrow2 = testrow;
        for (int k = 0; k < nrows; k++) {
          if( (MatIn(k,i)  != 0) && (MatIn(k,j)  == 0)) testrow[k] = (vecmax(testrow2));
        }
    
		testrow3 = avg_rank(testrow);
	
      out(_,count) = testrow3;
      count= count+1;
    }
  }
  return out;
} 
  

