#include "RcppArmadillo.h"
#include "LogisModel.hpp"

using namespace std;
using namespace arma;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//


//LogisMVS2GVB: LogisModelVariableSelection2GroupViaVariationalBayes
ObjLogisMVS2GVB rcpparma_LogisMVS2GVB(mat & Xr, colvec & yr, Options* opts) {

	int n = Xr.n_rows;
	int p = Xr.n_cols;
	int dispF = opts->dispF;
	int dispEvery = opts->display_gap;
	int max_iter = opts->max_iter;
	double tol = opts->epsStopLogLik;

	int iter = 1;
	int i = 0;

	mat X(Xr.begin(), n, p, false);// pointer
	colvec y(yr.begin(), yr.size(), false);

	//vector definition
	rowvec Pr = arma::zeros(1, 4);

	colvec F = arma::zeros(1, 1);
	colvec Fold = arma::zeros(1, 1); // low bound
	Fold(0, 0) = -arma::datum::inf;
	colvec diff = arma::zeros(1, 1);
	diff(0, 0) = arma::datum::inf;

	colvec alpha_prior = arma::zeros(1, 1);
	alpha_prior(0, 0) = 0.1;
	colvec logalpha_prior = arma::ones(1, 1);

	colvec inv_Sigma2beta = arma::zeros(1, 1);
	colvec sigma_beta = arma::zeros(1, 1);
	colvec sigma_error = arma::zeros(1, 1);

	sigma_beta(0, 0) = 1;

	colvec mu = arma::zeros(p, 1);     // initial variational mean
	colvec s = arma::zeros(p, 1);
	double beta0 = 0;
	colvec u = arma::zeros(p, 1);
	colvec alpha = arma::ones(p, 1)*alpha_prior;
	rowvec low_bound = arma::zeros(1, max_iter);

	colvec xxsigma = arma::zeros(p, 1);
	colvec y_i = arma::zeros(n, 1);


	colvec term1 = arma::zeros(p, 1);
	colvec alphamu = arma::zeros(p, 1);
	colvec alphamu2 = arma::zeros(p, 1);
	colvec logpow = arma::zeros(p, 1);

	colvec xty = X.t()*y;                //p-by-1
	colvec diagxtx = arma::sum(X%X).t(); //p-by-1 vector   
	colvec y_tilde = X*(alpha%mu); //n-by-1 vector
	colvec sumX = (sum(X, 0)).t();

	colvec psi = arma::zeros(n, 1);
	double a = 0.25;
	colvec b = a*psi - 1.0 / (1.0 + exp(-psi));
	colvec c = 0.5*a*psi%psi - psi / (1.0 + exp(-psi)) + log(1.0 + exp(psi));

	while (iter<max_iter && fabs(diff(0, 0)) > tol)
	{
		logalpha_prior = log(alpha_prior / (1 - alpha_prior));

		inv_Sigma2beta(0, 0) = 1.0 / sigma_beta(0, 0);

		xxsigma = a*diagxtx + inv_Sigma2beta(0, 0);

		s = 1.0 / xxsigma;

		for (i = 0; i<p; i++)
		{
			y_i = y_tilde - alpha(i, 0)*mu(i, 0)*X(arma::span(0, n - 1), i);

			mu(i, 0) = arma::as_scalar((xty(i, 0) + X(arma::span(0, n - 1), i).t()*(b%y) - a*beta0*sumX(i) - a*X(arma::span(0, n - 1), i).t()*y_i) / xxsigma(i, 0));

			u(i, 0) = logalpha_prior(0, 0) + 0.5*log(s(i, 0)) - 0.5*log(sigma_beta(0, 0)) + 0.5*mu(i, 0)*mu(i, 0) / s(i, 0);

			alpha(i, 0) = arma::as_scalar(1 / (1 + exp(-u(i, 0))));

			y_tilde = y_i + alpha(i, 0)*mu(i, 0)*X(arma::span(0, n - 1), i);
		}

		// M step
		term1 = alpha % (s + mu%mu);
		alphamu = alpha%mu;

		sigma_beta(0, 0) = sum(term1) / sum(alpha);
		if (sigma_beta(0, 0)<1e-6)
			sigma_beta(0, 0) = 1e-6;

		alpha_prior = arma::sum(alpha) / p;
		
		beta0 = (sum((1 + b)%y) - a*sum(y_tilde)) / (n*a);

		psi = (y_tilde + beta0)%y;
		b = a*psi - 1 / (1 + exp(-psi));
		c = 0.5*a*psi%psi - psi / (1 + exp(-psi)) + log(1 + exp(psi));

		for (i = 0; i<p; i++)
		{
			logpow(i, 0) = log(pow((1 - alpha(i, 0)) / (1 - alpha_prior(0, 0)), 1 - alpha(i, 0)));
		}
		//the limit using of power

		low_bound(0, iter - 1) = arma::as_scalar( sum((y_tilde  + beta0)%((1 + b)%y)) - sum(c) 
			- 0.5*a*sum((y_tilde + beta0) % (y_tilde + beta0))
			- 0.5*a*sum((term1 - (alphamu) % (alphamu)) % diagxtx)
			- sum(alpha%log(alpha / alpha_prior(0, 0)))
			- sum(logpow)
			+ sum(0.5*alpha % (1 + log(1 / sigma_beta(0, 0)*s) - 1 / sigma_beta(0, 0)*(s + mu%mu))));

		// update
		F(0, 0) = low_bound(0, iter-1);
		diff(0, 0) = F(0, 0) - Fold(0, 0);

		// check monotone
		if (diff(0, 0)<0)
		{
			throw std::range_error("The lower bound is not monotone");
		}

		if (dispF == 1)
		{
			if (fmod(iter, dispEvery) == 0)
			{
				Pr(0, 0) = iter;
				Pr(0, 1) = F(0, 0);
				Pr(0, 2) = Fold(0, 0);
				Pr(0, 3) = diff(0, 0);
				Pr.print("***Iteration*******Fnew********Fold**********Diff***");
			}
		}

		Fold(0, 0) = F(0, 0);
		iter = iter + 1;

		if (iter == max_iter)
		{
			throw std::range_error("The iteration reaches maximal iteration");
		}

	}

	rowvec lowerbound = low_bound(0, arma::span(0, iter - 2));

	ObjLogisMVS2GVB obj;

	obj.vardist_gamma = alpha;
	obj.vardist_mu = mu;
	obj.vardist_sigma2beta = s;
	obj.sigma2beta = sigma_beta;
	obj.alpha = alpha_prior;
	obj.Lq = lowerbound;
	obj.beta0 = beta0;

	return obj;
}


// Logis4LowerBound function
double Logis4LowerBound(mat beta0, colvec y, colvec y2, mat logpow, colvec y_tilde1, colvec y_tilde2, 
	colvec diagxtx, colvec diagxtx2, mat Zalpha, colvec sigma_beta, mat mu, mat s, 
	double a1, colvec b1, colvec c1, double a2, colvec b2, colvec c2, int p)
{
	colvec mu1 = mu(arma::span(0, p - 1), 0);
	colvec mu2 = mu(arma::span(0, p - 1), 1);

	colvec s1 = s(arma::span(0, p - 1), 0);
	colvec s2 = s(arma::span(0, p - 1), 1);

	colvec Zalpha1 = Zalpha(arma::span(0, p - 1), 0);
	colvec Zalpha2 = Zalpha(arma::span(0, p - 1), 1);

	colvec beta1 = Zalpha1 % mu1;
	colvec beta2 = Zalpha2 % mu2;

	double X1_term1 = arma::as_scalar(sum((y_tilde1 + beta0(0,0)) % ((1 + b1) % y)) - sum(c1));
	double X2_term1 = arma::as_scalar(sum((y_tilde2 + beta0(0,1)) % ((1 + b2) % y2)) - sum(c2));
	double term2 = arma::as_scalar(-0.5*a1*sum(square(y_tilde1 + beta0(0, 0))) - 0.5*a2*sum(square(y_tilde2 + beta0(0,1))));
	double X1_term3 = arma::as_scalar(-0.5*a1*sum((Zalpha1 % (s1 + mu1 % mu1) - beta1 % beta1) % diagxtx));
	double X2_term3 = arma::as_scalar(-0.5*a2*sum((Zalpha2 % (s2 + mu2 % mu2) - beta2 % beta2) % diagxtx2));
	double term4 = arma::as_scalar(-sum(arma::sum(logpow, 0)));
	double X1_term5 = arma::as_scalar(0.5*sum(Zalpha1 % (1 + log(1 / sigma_beta(0, 0)*s1) - 1 / sigma_beta(0, 0)*(s1 + mu1 % mu1))));
	double X2_term5 = arma::as_scalar(0.5*sum(Zalpha2 % (1 + log(1 / sigma_beta(1, 0)*s2) - 1 / sigma_beta(1, 0)*(s2 + mu2 % mu2))));
	double lowerbound = X1_term1 + X2_term1 + term2 + X1_term3 + X2_term3 + term4 + X1_term5 + X2_term5;
	return lowerbound;
}


//LogisMVS4GVB: LogisModelVariableSelection4GroupViaVariationalBayes
ObjLogisMVS4GVB rcpparma_LogisMVS4GVB(mat & Xr, mat & Xr2, colvec & yr, colvec & yr2, Options* opts) {

	int n1 = Xr.n_rows, p = Xr.n_cols;
	int n2 = Xr2.n_rows;
	int dispF = opts->dispF;
	int dispEvery = opts->display_gap;
	int max_iter = opts->max_iter;
	double tol = opts->epsStopLogLik;
	int constraintalpha = opts->constraintalpha;
	int iter = 1;
	int i = 0;

	mat X(Xr.begin(), n1, p, false);// pointer
	mat X2(Xr2.begin(), n2, p, false);
	colvec y(yr.begin(), yr.size(), false);
	colvec y2(yr2.begin(), yr2.size(), false);

	//vector definition
	colvec F = arma::zeros(1, 1);
	colvec Fold = arma::zeros(1, 1); // low bound
	colvec diff = arma::zeros(1, 1);
	diff(0, 0) = arma::datum::inf;
	Fold(0, 0) = -arma::datum::inf;

	rowvec alpha_prior = arma::zeros(1, 4);
	alpha_prior(0, 0) = 0.79;
	alpha_prior(0, 1) = 0.1;
	alpha_prior(0, 2) = 0.1;
	alpha_prior(0, 3) = 0.01;

	rowvec logalpha_prior21 = arma::zeros(1, 3);
	rowvec logalpha_prior31 = arma::zeros(1, 3);
	rowvec logalpha_prior42 = arma::zeros(1, 3);
	colvec inv_Sigma2beta = arma::zeros(2, 1);
	colvec sigma_beta = arma::zeros(2, 1);
	sigma_beta(0, 0) = 1;
	sigma_beta(1, 0) = 1;

	mat Zalpha = arma::zeros(p, 2);
	mat mu = arma::zeros(p, 2);     // initial variational mean
	mat s = arma::zeros(p, 2);
	mat u = arma::zeros(p, 3);
	mat beta0 = arma::zeros(1, 2);

	mat alpha = arma::ones(p, 4);
	alpha = alpha%repmat(alpha_prior, p, 1);
	mat xxsigma = arma::zeros(p, 2);
	mat y_i1 = arma::zeros(n1, 1);
	mat y_i2 = arma::zeros(n2, 1);

	rowvec Pr = arma::zeros(1, 4);

	colvec tmp_alpha_prior = arma::zeros(2, 1);
	mat term1 = arma::zeros(p, 2);
	mat beta = arma::zeros(p, 2);
	mat logpow = arma::zeros(p, 4);

	rowvec low_bound = arma::zeros(1, max_iter);

	colvec xty1 = X.t()*y;                //p-by-1
	colvec xty2 = X2.t()*y2;                //p-by-1
	colvec diagxtx = arma::sum(X%X).t(); //p-by-1 vector
	colvec diagxtx2 = arma::sum(X2%X2).t(); //p-by-1 vector
	colvec sumX = (sum(X, 0)).t();
	colvec sumX2 = (sum(X2, 0)).t();

	Zalpha(arma::span(0, p - 1), 0) = alpha(arma::span(0, p - 1), 2) + alpha(arma::span(0, p - 1), 3);
	Zalpha(arma::span(0, p - 1), 1) = alpha(arma::span(0, p - 1), 1) + alpha(arma::span(0, p - 1), 3);

	colvec y_tilde1 = X*(Zalpha(arma::span(0, p - 1), 0) % mu(arma::span(0, p - 1), 0));
	colvec y_tilde2 = X2*(Zalpha(arma::span(0, p - 1), 1) % mu(arma::span(0, p - 1), 1));

	colvec psi1 = arma::zeros(n1,1);
	double a1 = 0.25;
	colvec b1 = a1*psi1 - 1 / (1 + exp(-psi1));
	colvec c1 = 0.5*a1*psi1%psi1 - psi1 / (1 + exp(-psi1)) + log(1 + exp(psi1));

	colvec psi2 = arma::zeros(n2, 1);
	double a2 = 0.25;
	colvec b2 = a2*psi2 - 1 / (1 + exp(-psi2));
	colvec c2 = 0.5*a2*psi2%psi2 - psi2 / (1 + exp(-psi2)) + log(1 + exp(psi2));

	while (iter<max_iter && fabs(diff(0, 0)) > tol)
	{
		logalpha_prior21(0, 0) = log(alpha_prior(0, 1) / (alpha_prior(0, 0)));
		logalpha_prior31(0, 0) = log(alpha_prior(0, 2) / (alpha_prior(0, 0)));
		logalpha_prior42(0, 0) = log(alpha_prior(0, 3) / (alpha_prior(0, 1)));

		inv_Sigma2beta = 1.0 / sigma_beta;

		xxsigma(arma::span(0, p - 1), 0) = a1*diagxtx + inv_Sigma2beta(0, 0);
		xxsigma(arma::span(0, p - 1), 1) = a2*diagxtx2 + inv_Sigma2beta(1, 0);

		s(arma::span(0, p - 1), 0) = 1.0 / xxsigma(arma::span(0, p - 1), 0);
		s(arma::span(0, p - 1), 1) = 1.0 / xxsigma(arma::span(0, p - 1), 1);

		for (i = 0; i<p; i++)
		{
			y_i1 = y_tilde1 - Zalpha(i, 0)*mu(i, 0)*X(arma::span(0, n1 - 1), i);
			y_i2 = y_tilde2 - Zalpha(i, 1)*mu(i, 1)*X2(arma::span(0, n2 - 1), i);

			mu(i, 0) = arma::as_scalar((xty1(i, 0) + X(arma::span(0, n1 - 1), i).t()*(b1%y) - a1*beta0(0,0)*sumX(i,0) - a1*X(arma::span(0, n1 - 1), i).t()*y_i1) / xxsigma(i, 0));
			mu(i, 1) = arma::as_scalar((xty2(i, 0) + X2(arma::span(0, n2 - 1), i).t()*(b2%y2) - a2*beta0(0,1)*sumX2(i,0) - a2*X2(arma::span(0, n2 - 1), i).t()*y_i2) / xxsigma(i, 1));

			u(i, 0) = logalpha_prior21(0, 0) + 0.5*log(s(i, 1)) - 0.5*log(sigma_beta(1, 0)) + 0.5*mu(i, 1)*mu(i, 1) / s(i, 1);
			u(i, 1) = logalpha_prior31(0, 0) + 0.5*log(s(i, 0)) - 0.5*log(sigma_beta(0, 0)) + 0.5*mu(i, 0)*mu(i, 0) / s(i, 0);
			u(i, 2) = logalpha_prior42(0, 0) + 0.5*log(s(i, 0)) - 0.5*log(sigma_beta(0, 0)) + 0.5*mu(i, 0)*mu(i, 0) / s(i, 0);

			alpha(i, 0) = 1 / (1 + exp(u(i, 0)) + exp(u(i, 1)) + exp(u(i, 0) + u(i, 2)));
			alpha(i, 1) = 1 / (1 + exp(-u(i, 0)) + exp(u(i, 1) - u(i, 0)) + exp(u(i, 2)));
			alpha(i, 2) = 1 / (1 + exp(u(i, 0) - u(i, 1)) + exp(-u(i, 1)) + exp(u(i, 0) + u(i, 2) - u(i, 1)));
			alpha(i, 3) = 1 / (1 + exp(-u(i, 2)) + exp(u(i, 1) - u(i, 0) - u(i, 2)) + exp(-u(i, 0) - u(i, 2)));

			Zalpha(i, 0) = alpha(i, 2) + alpha(i, 3);
			Zalpha(i, 1) = alpha(i, 1) + alpha(i, 3);

			y_tilde1 = y_i1 + Zalpha(i, 0)*mu(i, 0)*X(arma::span(0, n1 - 1), i);
			y_tilde2 = y_i2 + Zalpha(i, 1)*mu(i, 1)*X2(arma::span(0, n2 - 1), i);

		}
		// M step
		term1 = Zalpha % (s + mu%mu);  //p-by-2
		beta = Zalpha%mu;     //p-by-2

		sigma_beta(0, 0) = sum(term1(arma::span(0, p - 1), 0)) / sum(Zalpha(arma::span(0, p - 1), 0));
		if (sigma_beta(0, 0)<1e-6)
			sigma_beta(0, 0) = 1e-6;

		sigma_beta(1, 0) = sum(term1(arma::span(0, p - 1), 1)) / sum(Zalpha(arma::span(0, p - 1), 1));
		if (sigma_beta(1, 0)<1e-6)
			sigma_beta(1, 0) = 1e-6;

		alpha_prior = arma::sum(alpha, 0) / p;
		if (constraintalpha == 1)
		{
			tmp_alpha_prior(0) = alpha_prior(2) + alpha_prior(3);
			tmp_alpha_prior(1) = alpha_prior(1) + alpha_prior(3);
			alpha_prior(0) = prod(1 - tmp_alpha_prior);
			alpha_prior(1) = tmp_alpha_prior(1)*(1 - tmp_alpha_prior(0));
			alpha_prior(2) = tmp_alpha_prior(0)*(1 - tmp_alpha_prior(1));
			alpha_prior(3) = prod(tmp_alpha_prior);
		}



		beta0(0,0) = (sum((1 + b1) % y) - a1*sum(y_tilde1)) / (n1*a1);
		beta0(0,1) = (sum((1 + b2) % y2) - a2*sum(y_tilde2)) / (n2*a2);

		psi1 =(y_tilde1 + beta0(0,0)) % y;
		b1 = a1*psi1 - 1 / (1 + exp(-psi1));
		c1 = 0.5*a1*psi1%psi1 - psi1 / (1 + exp(-psi1)) + log(1 + exp(psi1));

		psi2 = (y_tilde2 + beta0(0,1)) % y2;
		b2 = a2*psi2 - 1 / (1 + exp(-psi2));
		c2 = 0.5*a2*psi2%psi2 - psi2 / (1 + exp(-psi2)) + log(1 + exp(psi2));

		for (i = 0; i<p; i++)
		{
			logpow(i, 0) = log(pow((alpha(i, 0)) / (alpha_prior(0, 0)), alpha(i, 0)));
			logpow(i, 1) = log(pow((alpha(i, 1)) / (alpha_prior(0, 1)), alpha(i, 1)));
			logpow(i, 2) = log(pow((alpha(i, 2)) / (alpha_prior(0, 2)), alpha(i, 2)));
			logpow(i, 3) = log(pow((alpha(i, 3)) / (alpha_prior(0, 3)), alpha(i, 3)));
		}

		//C language attention integer / integer = integer
		low_bound(0, iter - 1) = Logis4LowerBound(beta0, y, y2, logpow, y_tilde1, y_tilde2, 
			diagxtx, diagxtx2, Zalpha, sigma_beta, mu, s, a1, b1, c1, a2, b2, c2, p);

		// update
		F(0, 0) = low_bound(0, iter-1);
		diff = F(0, 0) - Fold(0, 0);

		// check monotone
		if (diff(0, 0)<0)
		{
			throw std::range_error("The lower bound is not monotone");
		}

		if (dispF == 1)
		{
			if (fmod(iter, dispEvery) == 0)
			{
				Pr(0, 0) = iter;
				Pr(0, 1) = F(0, 0);
				Pr(0, 2) = Fold(0, 0);
				Pr(0, 3) = diff(0, 0);
				Pr.print("***Iteration*******Fnew********Fold**********Diff***");
			}
		}

		Fold(0, 0) = F(0, 0);
		iter = iter + 1;

		if (iter == max_iter)
		{
			throw std::range_error("The iteration reaches maximal iteration");
		}

	}

	rowvec lowerbound = low_bound(0, arma::span(0, iter - 2));

	ObjLogisMVS4GVB obj;

	obj.vardist_gamma = Zalpha;
	obj.vardist_mu = mu;
	obj.vardist_sigma2beta = s;
	obj.sigma2beta = sigma_beta;
	obj.alpha = alpha_prior;
	obj.Lq = lowerbound;
	obj.beta0 = beta0;

	return obj;
}


//LogisWFMVS2GVB: LogisModelWithFixVariableSelection2GroupViaVariationalBayes
ObjLogisWFMVS2GVB rcpparma_LogisWFMVS2GVB(mat & Xr, mat & Zr, colvec & yr, Options* opts) {

	int n = Xr.n_rows;
	int p = Xr.n_cols;
	int q_ = Zr.n_cols;
	int q = q_ + 1;
	int dispF = opts->dispF;
	int dispEvery = opts->display_gap;
	int max_iter = opts->max_iter;
	double tol = opts->epsStopLogLik;
	int iter = 1;
	int i = 0;

	mat X(Xr.begin(), n, p, false);// pointer
	mat Z_(Zr.begin(), n, q_, false);// pointer
	mat Z = arma::ones(n, q); // add constant
	Z(span(0, n - 1), span(1,q_)) = Z_;
	colvec y(yr.begin(), yr.size(), false);

	//vector definition
	rowvec Pr = arma::zeros(1, 4);

	colvec F = arma::zeros(1, 1);
	colvec Fold = arma::zeros(1, 1); // low bound
	Fold(0, 0) = -arma::datum::inf;
	colvec diff = arma::zeros(1, 1);
	diff(0, 0) = arma::datum::inf;

	colvec alpha_prior = arma::zeros(1, 1);
	alpha_prior(0, 0) = 0.1;
	colvec logalpha_prior = arma::ones(1, 1);

	colvec inv_Sigma2beta = arma::zeros(1, 1);
	colvec sigma_beta = arma::zeros(1, 1);
	colvec sigma_error = arma::zeros(1, 1);

	sigma_beta(0, 0) = 1;

	colvec mu = arma::zeros(p, 1);     // initial variational mean
	colvec s = arma::zeros(p, 1);
	colvec u = arma::zeros(p, 1);
	colvec alpha = arma::ones(p, 1)*alpha_prior;
	rowvec low_bound = arma::zeros(1, max_iter);
	double beta0 = 0;

	colvec xxsigma = arma::zeros(p, 1);
	colvec y_i = arma::zeros(n, 1);

	colvec term1 = arma::zeros(p, 1);
	colvec beta = arma::zeros(p, 1);
	colvec alphamu2 = arma::zeros(p, 1);
	colvec logpow = arma::zeros(p, 1);

	colvec xty = X.t()*y;                //p-by-1
	colvec diagxtx = arma::sum(X%X).t(); //p-by-1 vector   
	colvec y_tilde = X*(alpha%mu);         //n-by-1 vector
	//mat xtx = X.t()*X;

	mat ZtZ = Z.t()*Z;
	mat Zty = Z.t()*y;
	mat Ztx = Z.t()*X;
	colvec sumX = (sum(X, 0)).t();

	colvec U = arma::zeros(q, 1);

	colvec psi = arma::zeros(n, 1);
	double a = 0.25;
	colvec b = a*psi - 1 / (1 + exp(-psi));
	colvec c = 0.5*a*psi%psi - psi / (1 + exp(-psi)) + log(1 + exp(psi));

	while (iter<max_iter && fabs(diff(0, 0)) > tol)
	{
		logalpha_prior = log(alpha_prior / (1 - alpha_prior));

		inv_Sigma2beta(0, 0) = 1.0 / sigma_beta(0, 0);

		xxsigma = a*diagxtx + inv_Sigma2beta(0, 0);

		s = 1.0 / xxsigma;

		for (i = 0; i<p; i++)
		{
			y_i = y_tilde - alpha(i, 0)*mu(i, 0)*X(arma::span(0, n - 1), i);

			mu(i, 0) = arma::as_scalar((xty(i, 0) + X(arma::span(0, n - 1), i).t()*(b%y) - a*beta0*sumX(i) - a*Ztx(arma::span(0, q - 1), i).t()*U - a*X(arma::span(0, n - 1), i).t()*y_i) / xxsigma(i, 0));

			u(i, 0) = logalpha_prior(0, 0) + 0.5*log(s(i, 0)) - 0.5*log(sigma_beta(0, 0)) + 0.5*mu(i, 0)*mu(i, 0) / s(i, 0);

			alpha(i, 0) = arma::as_scalar(1 / (1 + exp(-u(i, 0))));

			y_tilde = y_i + alpha(i, 0)*mu(i, 0)*X(arma::span(0, n - 1), i);
		}

		// M step
		term1 = alpha % (s + mu%mu);
		beta = alpha%mu;

		sigma_beta(0, 0) = sum(term1) / sum(alpha);
		if (sigma_beta(0, 0)<1e-6)
			sigma_beta(0, 0) = 1e-6;

		alpha_prior = arma::sum(alpha) / p;

		U = inv_sympd(ZtZ)*(1.0 / a*Z.t()*((1+b)%y) - Z.t()*y_tilde - beta0*(sum(Z, 0)).t());

		//beta0 = (sum((1 + b) % y) - a*sum(Z*U + y_tilde)) / (n*a);

		psi = (y_tilde + Z*U + beta0)%y;
		b = a*psi - 1 / (1 + exp(-psi));
		c = 0.5*a*psi%psi - psi / (1 + exp(-psi)) + log(1 + exp(psi));
		
		for (i = 0; i<p; i++)
		{
			logpow(i, 0) = log(pow((1 - alpha(i, 0)) / (1 - alpha_prior(0, 0)), 1 - alpha(i, 0)));
		}
		//the limit using of power

		low_bound(0, iter - 1) = arma::as_scalar(sum((y_tilde + Z*U + beta0) % ((1 + b) % y)) - sum(c)
			- a*beta0*sum(Z*U) - a*y_tilde.t()*Z*U - 0.5*a*(U.t()*ZtZ*U)
			- 0.5*a*sum((y_tilde + beta0) % (y_tilde + beta0))
			- 0.5*a*sum((term1 - (beta) % (beta)) % diagxtx)
			- sum(alpha%log(alpha / alpha_prior(0, 0)))
			- sum(logpow)
			+ sum(0.5*alpha % (1 + log(1 / sigma_beta(0, 0)*s) - 1 / sigma_beta(0, 0)*(s + mu%mu))));

		// update
		F(0, 0) = low_bound(0, iter - 1);
		diff(0, 0) = F(0, 0) - Fold(0, 0);
		
		// check monotone
		if (diff(0, 0)<0)
		{
			throw std::range_error("The lower bound is not monotone");
		}

		if (dispF == 1)
		{
			if (fmod(iter, dispEvery) == 0)
			{
				Pr(0, 0) = iter;
				Pr(0, 1) = F(0, 0);
				Pr(0, 2) = Fold(0, 0);
				Pr(0, 3) = diff(0, 0);
				Pr.print("***Iteration*******Fnew********Fold**********Diff***");
			}
		}

		Fold(0, 0) = F(0, 0);
		iter = iter + 1;

		if (iter == max_iter)
		{
			throw std::range_error("The iteration reaches maximal iteration");
		}

	}

	rowvec lowerbound = low_bound(0, arma::span(0, iter - 2));

	ObjLogisWFMVS2GVB obj;

	obj.vardist_gamma = alpha;
	obj.vardist_mu = mu;
	obj.vardist_sigma2beta = s;
	obj.sigma2beta = sigma_beta;
	obj.alpha = alpha_prior;
	obj.Lq = lowerbound;
	obj.u = U;
	//obj.beta0 = beta0;

	return obj;
}


// Logis4WFLowerBound function
double Logis4WFLowerBound(colvec beta0, mat Z, mat Z2, colvec psi1, colvec psi2, mat logpow, 
	colvec xty1, colvec xty2, colvec y_tilde1, colvec y_tilde2, 
	colvec diagxtx, colvec diagxtx2, int n1, int n2, mat Zalpha, colvec sigma_beta, mat mu, mat s, 
	double a1, colvec b1, colvec c1, double a2, colvec b2, colvec c2, int p, int q1,int q2, mat ZtZ, mat Zty, mat Ztx, 
	mat ZtZ2, mat Zty2, mat Ztx2, mat U)
{

	colvec mu1 = mu(arma::span(0, p - 1), 0);
	colvec mu2 = mu(arma::span(0, p - 1), 1);

	colvec s1 = s(arma::span(0, p - 1), 0);
	colvec s2 = s(arma::span(0, p - 1), 1);

	colvec Zalpha1 = Zalpha(arma::span(0, p - 1), 0);
	colvec Zalpha2 = Zalpha(arma::span(0, p - 1), 1);

	colvec beta1 = Zalpha1 % mu1;
	colvec beta2 = Zalpha2 % mu2;

	colvec U1 = U(arma::span(0, q1 - 1), 0);
	colvec U2 = U(arma::span(0, q2 - 1), 1);

	double X1_term1 = arma::as_scalar( sum( psi1 % (1 + b1) ) - sum(c1) );
	double X2_term1 = arma::as_scalar( sum( psi2 % (1 + b2) ) - sum(c2) );

	double X1_term2 = arma::as_scalar(-a1*beta0(0,0)*sum(Z*U1) - a1*y_tilde1.t()*Z*U1 - 0.5*a1*(U1.t()*ZtZ*U1) - 0.5*a1*sum(square(y_tilde1 + beta0(0,0)) ));
	double X2_term2 = arma::as_scalar(-a2*beta0(1,0)*sum(Z2*U2) - a2*y_tilde2.t()*Z2*U2 - 0.5*a2*(U2.t()*ZtZ2*U2) - 0.5*a2*sum(square(y_tilde2 + beta0(1,0)) ));
	
	double X1_term3 = arma::as_scalar(-0.5*a1*sum((Zalpha1 % (s1 + mu1 % mu1) - beta1 % beta1) % diagxtx));
	double X2_term3 = arma::as_scalar(-0.5*a2*sum((Zalpha2 % (s2 + mu2 % mu2) - beta2 % beta2) % diagxtx2));
	double term4 = arma::as_scalar(-sum(arma::sum(logpow, 0)));
	double X1_term5 = arma::as_scalar(0.5*sum(Zalpha1 % (1 + log(1 / sigma_beta(0, 0)*s1) - 1 / sigma_beta(0, 0)*(s1 + mu1 % mu1))));
	double X2_term5 = arma::as_scalar(0.5*sum(Zalpha2 % (1 + log(1 / sigma_beta(1, 0)*s2) - 1 / sigma_beta(1, 0)*(s2 + mu2 % mu2))));
	double lowerbound = X1_term1 + X2_term1 + X1_term2 + X2_term2 + X1_term3 + X2_term3 + term4 + X1_term5 + X2_term5;
	return lowerbound;
}


//LogisMVS4GVB: LogisModelVariableSelection4GroupViaVariationalBayes
ObjLogisWFMVS4GVB rcpparma_LogisWFMVS4GVB(mat & Xr, mat & Xr2, mat & Zr, mat & Zr2, colvec & yr, colvec & yr2, Options* opts) {

	int n1 = Xr.n_rows, p = Xr.n_cols;
	int n2 = Xr2.n_rows;
	int q1_ = Zr.n_cols;
	int q2_ = Zr2.n_cols;
	int q1 = q1_ + 1;
	int q2 = q2_ + 1;
	int dispF = opts->dispF;
	int dispEvery = opts->display_gap;
	int max_iter = opts->max_iter;
	double tol = opts->epsStopLogLik;
	int constraintalpha = opts->constraintalpha;
	int iter = 1;
	int i = 0;

	mat X(Xr.begin(), n1, p, false);// pointer
	mat X2(Xr2.begin(), n2, p, false);

	mat Z_(Zr.begin(), n1, q1_, false);// pointer
	mat Z = arma::ones(n1, q1);
	Z(span(0, n1 - 1), span(1, q1_)) = Z_;

	mat Z2_(Zr2.begin(), n2, q2_, false);
	mat Z2 = arma::ones(n2, q2);
	Z2(span(0, n2 - 1), span(1, q2_)) = Z2_;


	colvec y(yr.begin(), yr.size(), false);
	colvec y2(yr2.begin(), yr2.size(), false);

	//vector definition
	colvec F = arma::zeros(1, 1);
	colvec Fold = arma::zeros(1, 1); // low bound
	colvec diff = arma::zeros(1, 1);
	diff(0, 0) = arma::datum::inf;
	Fold(0, 0) = -arma::datum::inf;

	rowvec alpha_prior = arma::zeros(1, 4);
	alpha_prior(0, 0) = 0.79;
	alpha_prior(0, 1) = 0.1;
	alpha_prior(0, 2) = 0.1;
	alpha_prior(0, 3) = 0.01;

	rowvec logalpha_prior21 = arma::zeros(1, 3);
	rowvec logalpha_prior31 = arma::zeros(1, 3);
	rowvec logalpha_prior42 = arma::zeros(1, 3);
	colvec inv_Sigma2beta = arma::zeros(2, 1);
	colvec sigma_beta = arma::zeros(2, 1);
	sigma_beta(0, 0) = 1;
	sigma_beta(1, 0) = 1;

	mat Zalpha = arma::zeros(p, 2);
	mat mu = arma::zeros(p, 2);     // initial variational mean
	mat s = arma::zeros(p, 2);
	mat u = arma::zeros(p, 3);
	colvec beta0 = arma::zeros(2, 1);

	mat alpha = arma::ones(p, 4);
	alpha = alpha%repmat(alpha_prior, p, 1);
	mat xxsigma = arma::zeros(p, 2);
	mat y_i1 = arma::zeros(n1, 1);
	mat y_i2 = arma::zeros(n2, 1);

	rowvec Pr = arma::zeros(1, 4);

	colvec tmp_alpha_prior = arma::zeros(2, 1);
	mat term1 = arma::zeros(p, 2);
	mat beta = arma::zeros(p, 2);
	mat logpow = arma::zeros(p, 4);

	rowvec low_bound = arma::zeros(1, max_iter);

	colvec xty1 = X.t()*y;                //p-by-1
	colvec xty2 = X2.t()*y2;                //p-by-1
	colvec diagxtx = arma::sum(X%X).t(); //p-by-1 vector
	colvec diagxtx2 = arma::sum(X2%X2).t(); //p-by-1 vector

	Zalpha(arma::span(0, p - 1), 0) = alpha(arma::span(0, p - 1), 2) + alpha(arma::span(0, p - 1), 3);
	Zalpha(arma::span(0, p - 1), 1) = alpha(arma::span(0, p - 1), 1) + alpha(arma::span(0, p - 1), 3);

	colvec y_tilde1 = X*(Zalpha(arma::span(0, p - 1), 0) % mu(arma::span(0, p - 1), 0));
	colvec y_tilde2 = X2*(Zalpha(arma::span(0, p - 1), 1) % mu(arma::span(0, p - 1), 1));

	colvec psi1 = zeros(n1,1);
	double a1 = 0.25;
	colvec b1 = a1*psi1 - 1 / (1 + exp(-psi1));
	colvec c1 = 0.5*a1*psi1%psi1 - psi1 / (1 + exp(-psi1)) + log(1 + exp(psi1));

	colvec psi2 = zeros(n2, 1);
	double a2 = 0.25;
	colvec b2 = a2*psi2 - 1 / (1 + exp(-psi2));
	colvec c2 = 0.5*a2*psi2%psi2 - psi2 / (1 + exp(-psi2)) + log(1 + exp(psi2));

	mat ZtZ = Z.t()*Z;
	mat Zty = Z.t()*y;
	mat Ztx = Z.t()*X;

	int q = 0;
	if (q1 > q2)
	{
		q = q1;
	}
	else{
		q = q2;
	}
	mat U = arma::zeros(q, 2);

	mat ZtZ2 = Z2.t()*Z2;
	mat Zty2 = Z2.t()*y2;
	mat Ztx2 = Z2.t()*X2;
	colvec sumX = (sum(X, 0)).t();
	colvec sumX2 = (sum(X2, 0)).t();

	while (iter<max_iter && fabs(diff(0, 0)) > tol)
	{
		logalpha_prior21(0, 0) = log(alpha_prior(0, 1) / (alpha_prior(0, 0)));
		logalpha_prior31(0, 0) = log(alpha_prior(0, 2) / (alpha_prior(0, 0)));
		logalpha_prior42(0, 0) = log(alpha_prior(0, 3) / (alpha_prior(0, 1)));

		inv_Sigma2beta = 1.0 / sigma_beta;

		xxsigma(arma::span(0, p - 1), 0) = a1*diagxtx + inv_Sigma2beta(0, 0);
		xxsigma(arma::span(0, p - 1), 1) = a2*diagxtx2 + inv_Sigma2beta(1, 0);

		s(arma::span(0, p - 1), 0) = 1.0 / xxsigma(arma::span(0, p - 1), 0);
		s(arma::span(0, p - 1), 1) = 1.0 / xxsigma(arma::span(0, p - 1), 1);

		for (i = 0; i<p; i++)
		{
			y_i1 = y_tilde1 - Zalpha(i, 0)*mu(i, 0)*X(arma::span(0, n1 - 1), i);
			y_i2 = y_tilde2 - Zalpha(i, 1)*mu(i, 1)*X2(arma::span(0, n2 - 1), i);

			mu(i, 0) = arma::as_scalar((xty1(i, 0) + X(arma::span(0, n1 - 1), i).t()*(b1%y) - a1*beta0(0,0)*sumX(i)
				- a1*Ztx(arma::span(0, q1 - 1), i).t()*U(arma::span(0, q1 - 1), 0) - a1*X(arma::span(0, n1 - 1), i).t()*y_i1) / xxsigma(i, 0));
			mu(i, 1) = arma::as_scalar((xty2(i, 0) + X2(arma::span(0, n2 - 1), i).t()*(b2%y2) - a2*beta0(1,0)*sumX2(i)
				- a2*Ztx2(arma::span(0, q2 - 1), i).t()*U(arma::span(0, q2 - 1), 1) - a2*X2(arma::span(0, n2 - 1), i).t()*y_i2) / xxsigma(i, 1));

			u(i, 0) = logalpha_prior21(0, 0) + 0.5*log(s(i, 1)) - 0.5*log(sigma_beta(1, 0)) + 0.5*mu(i, 1)*mu(i, 1) / s(i, 1);
			u(i, 1) = logalpha_prior31(0, 0) + 0.5*log(s(i, 0)) - 0.5*log(sigma_beta(0, 0)) + 0.5*mu(i, 0)*mu(i, 0) / s(i, 0);
			u(i, 2) = logalpha_prior42(0, 0) + 0.5*log(s(i, 0)) - 0.5*log(sigma_beta(0, 0)) + 0.5*mu(i, 0)*mu(i, 0) / s(i, 0);

			alpha(i, 0) = 1 / (1 + exp(u(i, 0)) + exp(u(i, 1)) + exp(u(i, 0) + u(i, 2)));
			alpha(i, 1) = 1 / (1 + exp(-u(i, 0)) + exp(u(i, 1) - u(i, 0)) + exp(u(i, 2)));
			alpha(i, 2) = 1 / (1 + exp(u(i, 0) - u(i, 1)) + exp(-u(i, 1)) + exp(u(i, 0) + u(i, 2) - u(i, 1)));
			alpha(i, 3) = 1 / (1 + exp(-u(i, 2)) + exp(u(i, 1) - u(i, 0) - u(i, 2)) + exp(-u(i, 0) - u(i, 2)));

			Zalpha(i, 0) = alpha(i, 2) + alpha(i, 3);
			Zalpha(i, 1) = alpha(i, 1) + alpha(i, 3);

			y_tilde1 = y_i1 + Zalpha(i, 0)*mu(i, 0)*X(arma::span(0, n1 - 1), i);
			y_tilde2 = y_i2 + Zalpha(i, 1)*mu(i, 1)*X2(arma::span(0, n2 - 1), i);

		}
		// M step
		term1 = Zalpha % (s + mu%mu);  //p-by-2
		beta = Zalpha%mu;     //p-by-2

		sigma_beta(0, 0) = sum(term1(arma::span(0, p - 1), 0)) / sum(Zalpha(arma::span(0, p - 1), 0));
		if (sigma_beta(0, 0)<1e-6)
			sigma_beta(0, 0) = 1e-6;
		
		sigma_beta(1, 0) = sum(term1(arma::span(0, p - 1), 1)) / sum(Zalpha(arma::span(0, p - 1), 1));
		if (sigma_beta(1, 0)<1e-6)
			sigma_beta(1, 0) = 1e-6;
		
		alpha_prior = arma::sum(alpha, 0) / p;
		if (constraintalpha == 1)
		{
			tmp_alpha_prior(0) = alpha_prior(2) + alpha_prior(3);
			tmp_alpha_prior(1) = alpha_prior(1) + alpha_prior(3);
			alpha_prior(0) = prod(1 - tmp_alpha_prior);
			alpha_prior(1) = tmp_alpha_prior(1)*(1 - tmp_alpha_prior(0));
			alpha_prior(2) = tmp_alpha_prior(0)*(1 - tmp_alpha_prior(1));
			alpha_prior(3) = prod(tmp_alpha_prior);
		}

		U(arma::span(0, q1 - 1), 0) = inv_sympd(ZtZ)*(1.0 / a1*Z.t()*((1 + b1) % y) - Z.t()*y_tilde1 - beta0(0,0)*(sum(Z, 0)).t());
		U(arma::span(0, q2 - 1), 1) = inv_sympd(ZtZ2)*(1.0 / a2*Z2.t()*((1 + b2) % y2) - Z2.t()*y_tilde2 - beta0(1,0)*(sum(Z2, 0)).t());

		//beta0(0,0) = (sum((1 + b1) % y) - a1*sum(Z*U(arma::span(0, q - 1), 0) + y_tilde1)) / (n1*a1);
		//beta0(1,0) = (sum((1 + b2) % y2) - a2*sum(Z2*U(arma::span(0, q - 1), 1) + y_tilde2)) / (n2*a2);

		psi1 = (y_tilde1 + Z*U(arma::span(0, q1 - 1), 0) + beta0(0,0) )%y;
		b1 = a1*psi1 - 1 / (1 + exp(-psi1));
		c1 = 0.5*a1*psi1%psi1 - psi1 / (1 + exp(-psi1)) + log(1 + exp(psi1));
		
		psi2 = (y_tilde2 + Z2*U(arma::span(0, q2 - 1), 1) + beta0(1,0) ) % y2;
		b2 = a2*psi2 - 1 / (1 + exp(-psi2));
		c2 = 0.5*a2*psi2%psi2 - psi2 / (1 + exp(-psi2)) + log(1 + exp(psi2));

		
		for (i = 0; i<p; i++)
		{
			logpow(i, 0) = log(pow((alpha(i, 0)) / (alpha_prior(0, 0)), alpha(i, 0)));
			logpow(i, 1) = log(pow((alpha(i, 1)) / (alpha_prior(0, 1)), alpha(i, 1)));
			logpow(i, 2) = log(pow((alpha(i, 2)) / (alpha_prior(0, 2)), alpha(i, 2)));
			logpow(i, 3) = log(pow((alpha(i, 3)) / (alpha_prior(0, 3)), alpha(i, 3)));
		}

		//C language attention integer / integer = integer
		low_bound(0, iter - 1) = Logis4WFLowerBound(beta0, Z, Z2, psi1, psi2, logpow, 
			xty1, xty2, y_tilde1, y_tilde2, 
			diagxtx, diagxtx2, n1, n2, Zalpha, sigma_beta, mu, s, 
			a1, b1, c1, a2, b2, c2, p, q1,q2, ZtZ, Zty, Ztx, ZtZ2, Zty2, Ztx2, U);

		// update
		F(0, 0) = low_bound(0, iter - 1);
		diff = F(0, 0) - Fold(0, 0);

		// check monotone
		if (diff(0, 0)<0)
		{
			throw std::range_error("The lower bound is not monotone");
		}

		if (dispF == 1)
		{
			if (fmod(iter, dispEvery) == 0)
			{
				Pr(0, 0) = iter;
				Pr(0, 1) = F(0, 0);
				Pr(0, 2) = Fold(0, 0);
				Pr(0, 3) = diff(0, 0);
				Pr.print("***Iteration*******Fnew********Fold**********Diff***");
			}
		}

		Fold(0, 0) = F(0, 0);
		iter = iter + 1;

		if (iter == max_iter)
		{
			throw std::range_error("The iteration reaches maximal iteration");
		}

	}

	rowvec lowerbound = low_bound(0, arma::span(0, iter - 2));

	ObjLogisWFMVS4GVB obj;

	obj.vardist_gamma = Zalpha;
	obj.vardist_mu = mu;
	obj.vardist_sigma2beta = s;
	obj.sigma2beta = sigma_beta;
	obj.alpha = alpha_prior;
	obj.Lq = lowerbound;
	obj.u = U;
	return obj;
}
