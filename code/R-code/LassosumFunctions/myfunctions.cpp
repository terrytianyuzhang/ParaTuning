/**
 mylassosum
 functions.cpp
 Purpose: functions to perform mylassosum
 @author Peng Liu
 @version 0.1
 */
// [[Rcpp::interfaces(r, cpp)]]

#define ARMA_64BIT_WORD 1
#include <stdio.h>
#include <string>
#include <bitset>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

/**
 Opens a Plink binary files
 @s file name
 @BIT ifstream
 @return is plink file in major mode
 */

bool openPlinkBinaryFile(const std::string s, std::ifstream &BIT) {
  BIT.open(s.c_str(), std::ios::in | std::ios::binary);
  if (!BIT.is_open()) {
    throw "Cannot open the bed file";
  }

  // 2) else check for 0.99 SNP/Ind coding
  // 3) else print warning that file is too old
  char ch[1];
  BIT.read(ch, 1);
  std::bitset<8> b;
  b = ch[0];
  bool bfile_SNP_major = false;
  bool v1_bfile = true;
  // If v1.00 file format
  // Magic numbers for .bed file: 00110110 11011000 = v1.00 bed file
  // std::cerr << "check magic number" << std::endl;
  if ((b[2] && b[3] && b[5] && b[6]) && !(b[0] || b[1] || b[4] || b[7])) {
    // Next number
    BIT.read(ch, 1);
    b = ch[0];
    if ((b[0] && b[1] && b[3] && b[4]) && !(b[2] || b[5] || b[6] || b[7])) {
      // Read SNP/Ind major coding
      BIT.read(ch, 1);
      b = ch[0];
      if (b[0])
        bfile_SNP_major = true;
      else
        bfile_SNP_major = false;

      // if (bfile_SNP_major) std::cerr << "Detected that binary PED file is
      // v1.00 SNP-major mode" << std::endl;
      // else std::cerr << "Detected that binary PED file is v1.00
      // individual-major mode" << std::endl;

    } else
      v1_bfile = false;

  } else
    v1_bfile = false;
  // Reset file if < v1
  if (!v1_bfile) {
    Rcerr << "Warning, old BED file <v1.00 : will try to recover..."
          << std::endl;
    Rcerr << "  but you should --make-bed from PED )" << std::endl;
    BIT.close();
    BIT.clear();
    BIT.open(s.c_str(), std::ios::in | std::ios::binary);
    BIT.read(ch, 1);
    b = ch[0];
  }
  // If 0.99 file format
  if ((!v1_bfile) && (b[1] || b[2] || b[3] || b[4] || b[5] || b[6] || b[7])) {
    Rcerr << std::endl
          << " *** Possible problem: guessing that BED is < v0.99      *** "
          << std::endl;
    Rcerr << " *** High chance of data corruption, spurious results    *** "
          << std::endl;
    Rcerr
      << " *** Unless you are _sure_ this really is an old BED file *** "
      << std::endl;
    Rcerr << " *** you should recreate PED -> BED                      *** "
          << std::endl
          << std::endl;
    bfile_SNP_major = false;
    BIT.close();
    BIT.clear();
    BIT.open(s.c_str(), std::ios::in | std::ios::binary);
  } else if (!v1_bfile) {
    if (b[0])
      bfile_SNP_major = true;
    else
      bfile_SNP_major = false;
    Rcerr << "Binary PED file is v0.99" << std::endl;
    if (bfile_SNP_major)
      Rcerr << "Detected that binary PED file is in SNP-major mode"
            << std::endl;
      else
        Rcerr << "Detected that binary PED file is in individual-major mode"
              << std::endl;
  }
  return bfile_SNP_major;
}


// [[Rcpp::export]]
int myrepelnet(double lambda1, double lambda2, double gamma, arma::vec& diag1, arma::vec& diag2,
               arma::mat& X1, arma::mat& X2, arma::vec& r1, arma::vec& r2, //double thr,
               arma::vec& x, arma::vec& yhat1, arma::vec& yhat2, int trace, int maxiter,
               arma::Col<int>& startvec1, arma::Col<int>& endvec1,
               arma::Col<int>& startvec2, arma::Col<int>& endvec2,
               arma::vec& blocks1, arma::vec& blocks2)
{

  // Loop through all the SNPs
  int n1=X1.n_rows;
  int n2=X2.n_rows;
  int p=X1.n_cols;

  if(X1.n_cols != X2.n_cols) stop("X1.n_cols != X2.n_cols");
  if(r1.n_elem != p) stop("r1.n_elem != p");
  if(r2.n_elem != p) stop("r2.n_elem != p");
  if(x.n_elem != p) stop("x.n_elem != p");
  if(yhat1.n_elem != n1) stop("yhat1.n_elem != n1");
  if(yhat2.n_elem != n2) stop("yhat2.n_elem != n2");
  if(diag1.n_elem != p) stop("diag1.n_elem != p");
  if(diag2.n_elem != p) stop("diag2.n_elem != p");

  double disx, thr, del,t1,t2,t;
  int s1,s2;

  arma::vec yhat1touse;
  arma::vec yhat2touse;
  arma::vec yhat1del;
  arma::vec yhat2del;

  int conv=0; // initialize convergence indicator
  for(int k=0;k<maxiter ;k++) {
    disx=0.0;
    arma::vec x_old=x;
    for(int i=0;i < p; i++) {
      s1 = blocks1(i);
      s2 = blocks2(i);

      if(i==0) {
        arma::vec x1touse=x.subvec(startvec1(s1), endvec1(s1));
        yhat1touse=X1.cols(startvec1(s1), endvec1(s1)) * x1touse;
      } else if(s1 != blocks1(i-1)) {
        arma::vec x1touse=x.subvec(startvec1(s1), endvec1(s1));
        yhat1touse=X1.cols(startvec1(s1), endvec1(s1)) * x1touse;
      }

      if(i==0) {
        arma::vec x2touse=x.subvec(startvec2(s2), endvec2(s2));
        yhat2touse=X2.cols(startvec2(s2), endvec2(s2)) * x2touse;
      } else if(s2 != blocks2(i-1)){
        arma::vec x2touse=x.subvec(startvec2(s2), endvec2(s2));
        yhat2touse=X2.cols(startvec2(s2), endvec2(s2)) * x2touse;
      }

      double denom=gamma * diag1(i) + (1-gamma) * diag2(i) + lambda2; //double check this step does what it is supposed to do

      t1 = diag1(i) * x(i) + r1(i) - arma::dot(X1.col(i), yhat1touse);
      t2 = diag2(i) * x(i) + r2(i) - arma::dot(X2.col(i), yhat2touse);
      t = gamma* t1 + (1-gamma) * t2;
      // Think of it as r(j) - (dotproduct(X.col(j), yhat) - diag(j)*xj)
      if(std::abs(t)-lambda1 > 0.0) {
        x(i)=copysign(std::abs(t)-lambda1, t)/denom;
      } else{
        x(i)=0.0;
      }
      if(x(i)==x_old(i)) continue;
      del= x(i) - x_old(i);
      yhat1del= del* X1.col(i);
      yhat2del= del* X2.col(i);
      disx += del*del;

      //update
      yhat1touse += yhat1del;
      yhat2touse += yhat2del;
      yhat1 += yhat1del;
      yhat2 += yhat2del;
    }

    checkUserInterrupt(); //not sure what it does, just keep it for now
    //if(trace > 0) Rcout << "Iteration: " << k << "\n";

    //dlx = sqrt(accu(pow(x-x_old,2)));
    //dlx = max(abs(x-x_old)); // this step could be improved
    disx = sqrt(disx);
    thr=0.01*sqrt(accu(pow(x_old,2)));
    if(trace > 0) Rcout << "Iteration: " << k << "Distance: " << disx << " Threshold:"<< thr <<"\n";
    if(disx <= thr) {
      conv=1;
      break;
    }
  }

  //yhat1 = X1 * x;
  //yhat2 = X2 * x;

  return conv;
}


//' imports genotypeMatrix
//'
//' @param fileName location of bam file
//' @param N number of subjects
//' @param P number of positions
//' @param col_skip_pos which variants should we skip
//' @param col_skip which variants should we skip
//' @param keepbytes which bytes to keep
//' @param keepoffset what is the offset
//' @return an armadillo genotype matrix
//' @keywords internal
//'
// [[Rcpp::export]]
arma::mat genotypeMatrix(const std::string fileName, int N, int P,
                         arma::Col<int> col_skip_pos, arma::Col<int> col_skip,
                         arma::Col<int> keepbytes, arma::Col<int> keepoffset,
                         const int fillmissing) {

  std::ifstream bedFile;
  bool snpMajor = openPlinkBinaryFile(fileName+".bed", bedFile);

  if (!snpMajor)
    throw std::runtime_error("We currently have no plans of implementing the "
                               "individual-major mode. Please use the snp-major "
                               "format");

  int i = 0;
  int ii = 0;
  const bool colskip = (col_skip_pos.n_elem > 0);
  unsigned long long int Nbytes = ceil(N / 4.0);
  const bool selectrow = (keepbytes.n_elem > 0);
  int n, p, nskip;
  if (selectrow)
    n = keepbytes.n_elem;
  else
    n = N;

  if (colskip) {
    nskip = arma::accu(col_skip);
    p = P - nskip;
  }  else
    p = P;

  int j, jj, iii;

  arma::mat genotypes = arma::mat(n, p, arma::fill::zeros);
  std::bitset<8> b; // Initiate the bit array
  char ch[Nbytes];

  iii=0;
  while (i < P) {
     //Rcout << i << std::endl;
     Rcpp::checkUserInterrupt();
    if (colskip) {
      if (ii < col_skip.n_elem) {
        if (i == col_skip_pos[ii]) {
          bedFile.seekg(col_skip[ii] * Nbytes, bedFile.cur);
          i = i + col_skip[ii];
          ii++;
          continue;
        }
      }
    }

    bedFile.read(ch, Nbytes); // Read the information
    if (!bedFile)
      throw std::runtime_error(
          "Problem with the BED file...has the FAM/BIM file been changed?");

    j = 0;
    //Rcout << "i: " << i << "\n" << std::endl;
    //Rcout << "select row: " << selectrow << "\n" << std::endl;
    if (!selectrow) {
      for (jj = 0; jj < Nbytes; jj++) {
        b = ch[jj];

        int c = 0;
        while (c < 7 &&
               j < N) { // from the original PLINK: 7 because of 8 bits
          int first = b[c++];
          int second = b[c++];
          if (first == 0) {
            genotypes(j, iii) = (2 - second);
          }
          if(fillmissing == 0 && first == 1 && second == 0) genotypes(j, iii) =arma::datum::nan;
          j++;
        }
      }
    } else {
      for (jj = 0; jj < keepbytes.n_elem; jj++) {
        //Rcout << "j: " << j << "\n" << std::endl;
        //Rcout << "jj: " << jj << "\n" << std::endl;
        b = ch[keepbytes[jj]];
        //Rcout << "keepbytes[jj]:" << keepbytes[jj] << "\n" << std::endl;
        //Rcout << "b: " << b << "\n" << std::endl;

        int c = keepoffset[jj];
        //Rcout << "keepoffset: " << keepoffset << "\n" << std::endl;
        //Rcout << "keepoffset[jj]: " << keepoffset[jj] << "\n" << std::endl;
        //Rcout << "c: " << c << "\n" << std::endl;
        int first = b[c++];
        int second = b[c];
        if (first == 0) {
          genotypes(j, iii) = (2 - second);
        }
        if(fillmissing == 0 && first == 1 && second == 0) genotypes(j, iii) =arma::datum::nan;
        j++;
      }
    }
    i++;
    iii++;
  }
  return genotypes;
}


//' normalize genotype matrix
//'
//' @param genotypes a armadillo genotype matrix
//' @return standard deviation
//' @keywords internal
//'
// [[Rcpp::export]]
arma::vec normalize(arma::mat &genotypes)
{
  int k = genotypes.n_cols;
  int n = genotypes.n_rows;
  arma::vec sd(k);
  for (int i = 0; i < k; ++i) {
    double m = arma::mean(genotypes.col(i));
    arma::vec mm(n); mm.fill(m);
    sd(i) = arma::stddev(genotypes.col(i));
    // sd(i) = 1.0;
    //return its normalised version (ie. having unit 2-norm)
    genotypes.col(i) = arma::normalise(genotypes.col(i) - mm);
  }
  return sd;
}


//' Runs elnet with various parameters
//'
//' @param lambda1 a vector of lambdas (lambda2 is 0)
//' @param fileName the file name of the reference panel
//' @param r a vector of correlations
//' @param N number of subjects
//' @param P number of position in reference file
//' @param col_skip_posR which variants should we skip
//' @param col_skipR which variants should we skip
//' @param keepbytesR required to read the PLINK file
//' @param keepoffsetR required to read the PLINK file
//' @param thr threshold
//' @param x a numeric vector of beta coefficients
//' @param trace if >1 displays the current iteration
//' @param maxiter maximal number of iterations
//' @param Constant a constant to multiply the standardized genotype matrix
//' @return a list of results
//' @keywords internal
//'
// [[Rcpp::export]]
List myrunElnet(arma::vec& lambda, double shrink, double gamma,
                const std::string fileName1,const std::string fileName2,
                arma::vec& r1, arma::vec& r2, int N1, int N2, int P,
                arma::Col<int>& col_skip_pos1, arma::Col<int>& col_skip1,
                arma::Col<int>& keepbytes1, arma::Col<int>& keepoffset1,
                arma::Col<int>& col_skip_pos2, arma::Col<int>& col_skip2,
                arma::Col<int>& keepbytes2, arma::Col<int>& keepoffset2,
                arma::vec& x, int trace, int maxiter,
                arma::Col<int>& startvec1, arma::Col<int>& endvec1,
                arma::Col<int>& startvec2, arma::Col<int>& endvec2,
                arma::vec& blocks1, arma::vec& blocks2) {
  // a) read bed file
  // b) standardize genotype matrix
  // c) multiply by constatant factor
  // d) perfrom elnet

  Rcout << "ABC" << std::endl;

  int i,j;
  arma::mat X1 = genotypeMatrix(fileName1, N1, P, col_skip_pos1, col_skip1, keepbytes1,
                                keepoffset1, 1);
  arma::mat X2 = genotypeMatrix(fileName2, N2, P, col_skip_pos2, col_skip2, keepbytes2,
                                keepoffset2, 1);

  Rcout << "DEF" << std::endl;

  if (X1.n_cols != r1.n_elem) {
    throw std::runtime_error("Number of positions in reference file 1 is not "
                               "equal the number of regression coefficients");
  }

  if (X2.n_cols != r2.n_elem) {
    throw std::runtime_error("Number of positions in reference file 2 is not "
                               "equal the number of regression coefficients");
  }

  arma::vec sd1 = normalize(X1);
  arma::vec sd2 = normalize(X2);

  X1 *= sqrt(1.0 - shrink);
  X2 *= sqrt(1.0 - shrink);

  arma::Col<int> conv(lambda.n_elem);
  int len = r1.n_elem;

  arma::mat beta(len, lambda.n_elem);
  arma::mat pred1(X1.n_rows, lambda.n_elem); pred1.zeros();
  arma::mat pred2(X2.n_rows, lambda.n_elem); pred2.zeros();
  arma::vec out(lambda.n_elem);
  arma::vec loss(lambda.n_elem);
  //tianyu added these two to record training error for each population
  //these quantities will be used in parameter tuning step
  arma::vec trainerror1(lambda.n_elem);
  arma::vec trainerror2(lambda.n_elem);
  arma::vec diag1(r1.n_elem); diag1.fill(1.0 - shrink);
  arma::vec diag2(r2.n_elem); diag2.fill(1.0 - shrink);
  // Rcout << "HIJ" << std::endl;

  for(j=0; j < diag1.n_elem; j++) {
    if(sd1(j) == 0.0) diag1(j) = 0.0;
  }
  for(j=0; j < diag2.n_elem; j++) {
    if(sd2(j) == 0.0) diag2(j) = 0.0;
  }
  // Rcout << "LMN" << std::endl;

  arma::vec fbeta(lambda.n_elem);
  //arma::vec yhat1 = X1 * x;
  //arma::vec yhat2 = X2 * x;
  arma::vec yhat1(X1.n_rows);
  arma::vec yhat2(X2.n_rows);
  //arma::vec yhat(genotypes.n_rows);
  // yhat = genotypes * x;


  // Rcout << "Starting loop" << std::endl;

  for (i = 0; i < lambda.n_elem; ++i) {
    if (trace > 0)
      Rcout << "lambda: " << lambda(i) << "\n" << std::endl;
    // yhat1 and yhat2 are changed because we are passing by reference
    out(i) =
      myrepelnet(lambda(i), shrink, gamma, diag1, diag2, X1, X2, r1, r2, x, yhat1, yhat2, trace-1, maxiter,
               startvec1, endvec1, startvec2, endvec2, blocks1, blocks2);

    beta.col(i) = x;
    for(j=0; j < beta.n_rows; j++) {
      if(sd1(j) == 0.0 && sd2(j) == 0.0) beta(j,i)=beta(j,i) * shrink;
    }
    if (out(i) != 1) {
      throw std::runtime_error("Not converging.....");
    }
    pred1.col(i) = yhat1;
    pred2.col(i) = yhat2;

    float part1;
    part1 = arma::as_scalar(gamma * arma::sum(arma::pow(yhat1, 2)));
    float part2;
    part2 = arma::as_scalar(gamma * 2.0 * arma::sum(x % r1));

    Rcout << "This is first part:" << part1 << "\n";
    Rcout << "This is second part:" << part2 << "\n";

    loss(i) = arma::as_scalar(gamma * arma::sum(arma::pow(yhat1, 2)) -
      gamma * 2.0 * arma::sum(x % r1) + (1-gamma) * (arma::sum(arma::pow(yhat2, 2))) - (1-gamma) * 2.0 * arma::sum(x % r2));

    trainerror1(i) = arma::as_scalar(arma::sum(arma::pow(yhat1, 2)) - 2.0 * arma::sum(x % r1));
    trainerror2(i) = arma::as_scalar(arma::sum(arma::pow(yhat2, 2)) - 2.0 * arma::sum(x % r2));

    fbeta(i) =
      arma::as_scalar(loss(i) + 2.0 * arma::sum(arma::abs(x)) * lambda(i) +
      arma::sum(arma::pow(x, 2)) * shrink);
  }
  return List::create(Named("lambda") = lambda, //let's say there are 10 lambdas
                      Named("beta") = beta, //dimension: SNP# times 10
                      Named("conv") = out, //Tianyu guesses this is a convergence indicator, same length as lambda
                      Named("pred1") = pred1,
                      Named("pred2") = pred2,
                      Named("loss") = loss, //without the penalty on regression coefficients
                      Named("fbeta") = fbeta, //inlcude the penalty
                      Named("trainerror1") = trainerror1,
                      Named("trainerror1") = trainerror1,
                      Named("sd1")= sd1,
                      Named("sd2")= sd2);
}

//' Multiply genotypeMatrix by a matrix (sparse)
//'
//' @param fileName location of bam file
//' @param N number of subjects
//' @param P number of positions
//' @param input the matrix
//' @param col_skip_pos which variants should we skip
//' @param col_skip which variants should we skip
//' @param keepbytes which bytes to keep
//' @param keepoffset what is the offset
//' @return an armadillo genotype matrix
//' @keywords internal
//'
// [[Rcpp::export]]
arma::mat pgs(const std::string fileName, int N, int P, double shrink,
              arma::Col<int>& col_skip_pos, arma::Col<int>& col_skip,
              arma::Col<int>& keepbytes, arma::Col<int>& keepoffset,
              const arma::mat& x,
              const int trace) {
  arma::mat X = genotypeMatrix(fileName, N, P, col_skip_pos, col_skip, keepbytes,
                                keepoffset, 1);

  arma::vec sd = normalize(X);
  X *= sqrt(1.0 - shrink);
  arma::mat yhat = X * x;
  return yhat;
}

// [[Rcpp::export]]
arma::mat prod(const std::string fileName, int N, int P,
              arma::Col<int>& col_skip_pos, arma::Col<int>& col_skip,
              arma::Col<int>& keepbytes, arma::Col<int>& keepoffset,
              const arma::mat& x) {
  arma::mat X = genotypeMatrix(fileName, N, P, col_skip_pos, col_skip, keepbytes,
                               keepoffset, 1);

  //arma::vec sd = normalize(X);
  //X *= sqrt(1.0 - shrink);
  arma::mat yhat = X * x;
  return yhat;
}
