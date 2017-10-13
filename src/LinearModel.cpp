#include "RcppArmadillo.h"
#include "LinearModel.hpp"
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

// C++ function
//LinearMVS4GVB: LinearModelVariableSelection4GroupViaVariationalBayes
ObjLinearMVS4GVB rcpparma_LinearMVS4GVB(mat & Xr, mat & Xr2, colvec & yr, colvec & yr2, Options* opts) {
	
	int n1 = Xr.n_rows, p = Xr.n_cols;
	int n2 = Xr2.n_rows;
	int dispF = opts->dispF;
	int dispEvery = opts->display_gap;
	int max_iter = opts->max_iter;
	double tol = opts->epsStopLogLik;
	int constraintalpha = opts->constraintalpha;
	int iter = 1;
	int i = 0;

	double normyytilde = 0;
	double normyytilde2 = 0;

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
	colvec sigma2e_Sigma2beta = arma::zeros(2, 1);
	colvec sigma_beta = arma::zeros(2, 1);
	colvec sigma_error = arma::zeros(2, 1);
	sigma_beta(0, 0) = var(y) / 2;
	sigma_error(0, 0) = var(y) / 2;
	sigma_beta(1, 0) = var(y2) / 2;
	sigma_error(1, 0) = var(y2) / 2;

	mat Zalpha = arma::zeros(p, 2);
	mat mu = arma::zeros(p, 2);     // initial variational mean
	mat s = arma::zeros(p, 2);
	mat u = arma::zeros(p, 3);
	mat alpha = arma::ones(p, 4);
	alpha = alpha%repmat(alpha_prior, p, 1);
	mat xxsigma = arma::zeros(p, 2);
	mat y_i1 = arma::zeros(n1, 1);
	mat y_i2 = arma::zeros(n2, 1);

	rowvec Pr = arma::zeros(1, 4);

	colvec tmp_alpha_prior = arma::zeros(2, 1);
	mat term1 = arma::zeros(p, 2);
	mat alphamu = arma::zeros(p, 2);
	mat alphamu2 = arma::zeros(p, 2);
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

	while (iter<max_iter && fabs(diff(0, 0)) > tol)
	{
		logalpha_prior21(0, 0) = log(alpha_prior(0, 1) / (alpha_prior(0, 0)));
		logalpha_prior31(0, 0) = log(alpha_prior(0, 2) / (alpha_prior(0, 0)));
		logalpha_prior42(0, 0) = log(alpha_prior(0, 3) / (alpha_prior(0, 1)));

		sigma2e_Sigma2beta = sigma_error / sigma_beta;

		xxsigma(arma::span(0, p - 1), 0) = diagxtx + sigma2e_Sigma2beta(0, 0);
		xxsigma(arma::span(0, p - 1), 1) = diagxtx2 + sigma2e_Sigma2beta(1, 0);

		s(arma::span(0, p - 1), 0) = sigma_error(0, 0) / xxsigma(arma::span(0, p - 1), 0);
		s(arma::span(0, p - 1), 1) = sigma_error(1, 0) / xxsigma(arma::span(0, p - 1), 1);

		for (i = 0; i<p; i++)
		{
			y_i1 = y_tilde1 - Zalpha(i, 0)*mu(i, 0)*X(arma::span(0, n1 - 1), i);
			y_i2 = y_tilde2 - Zalpha(i, 1)*mu(i, 1)*X2(arma::span(0, n2 - 1), i);

			mu(i, 0) = arma::as_scalar((xty1(i, 0) - X(arma::span(0, n1 - 1), i).t()*y_i1) / xxsigma(i, 0));
			mu(i, 1) = arma::as_scalar((xty2(i, 0) - X2(arma::span(0, n2 - 1), i).t()*y_i2) / xxsigma(i, 1));

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
		alphamu = Zalpha%mu;     //p-by-2
		alphamu2 = alphamu%alphamu;   //p-by-2

		normyytilde = arma::norm(y - y_tilde1, 2);
		normyytilde2 = arma::norm(y2 - y_tilde2, 2);

		sigma_error(0, 0) = (pow(normyytilde, 2) + arma::sum((term1(arma::span(0, p - 1), 0) - alphamu2(arma::span(0, p - 1), 0)) % diagxtx)) / n1;
		if (sigma_error(0, 0)<1e-6)
			sigma_error(0, 0) = 1e-6;

		sigma_error(1, 0) = (pow(normyytilde2, 2) + arma::sum((term1(arma::span(0, p - 1), 1) - alphamu2(arma::span(0, p - 1), 1)) % diagxtx2)) / n2;
		if (sigma_error(1, 0)<1e-6)
			sigma_error(1, 0) = 1e-6;

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

		for (i = 0; i<p; i++)
		{
			logpow(i, 0) = log(pow((alpha(i, 0)) / (alpha_prior(0, 0)), alpha(i, 0)));
			logpow(i, 1) = log(pow((alpha(i, 1)) / (alpha_prior(0, 1)), alpha(i, 1)));
			logpow(i, 2) = log(pow((alpha(i, 2)) / (alpha_prior(0, 2)), alpha(i, 2)));
			logpow(i, 3) = log(pow((alpha(i, 3)) / (alpha_prior(0, 3)), alpha(i, 3)));
		}

		//C language attention integer / integer = integer
		low_bound(0, iter-1) = -n1*0.5*log(2 * arma::datum::pi*sigma_error(0, 0)) - pow(normyytilde, 2) / sigma_error(0, 0) / 2
			- 1 / sigma_error(0, 0)*0.5*sum((term1(arma::span(0, p - 1), 0) - alphamu2(arma::span(0, p - 1), 0)) % diagxtx)
			- n2*0.5*log(2 * arma::datum::pi*sigma_error(1, 0)) - pow(normyytilde2, 2) / sigma_error(1, 0) / 2
			- 1 / sigma_error(1, 0)*0.5*sum((term1(arma::span(0, p - 1), 1) - alphamu2(arma::span(0, p - 1), 1)) % diagxtx2)
			- sum(arma::sum(logpow, 0))
			+ sum(0.5*Zalpha(arma::span(0, p - 1), 0) % (1 + log(1 / sigma_beta(0, 0)*s(arma::span(0, p - 1), 0)) - 1 / sigma_beta(0, 0)*(s(arma::span(0, p - 1), 0) + mu(arma::span(0, p - 1), 0) % mu(arma::span(0, p - 1), 0))))
			+ sum(0.5*Zalpha(arma::span(0, p - 1), 1) % (1 + log(1 / sigma_beta(1, 0)*s(arma::span(0, p - 1), 1)) - 1 / sigma_beta(1, 0)*(s(arma::span(0, p - 1), 1) + mu(arma::span(0, p - 1), 1) % mu(arma::span(0, p - 1), 1))));

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

		if (iter==max_iter)
		{
			throw std::range_error("The iteration reaches maximal iteration");
		}
	}

	rowvec lowerbound = low_bound(0, arma::span(0, iter - 2));

	ObjLinearMVS4GVB obj;

	obj.vardist_gamma = Zalpha;
	obj.vardist_mu = mu;
	obj.vardist_sigma2beta = s;
	obj.sigma2beta = sigma_beta;
	obj.alpha = alpha_prior;
	obj.sigma2e = sigma_error;
	obj.Lq = lowerbound;

	return obj;
}
 

//LinearMVS2GVB: LinearModelVariableSelection2GroupViaVariationalBayes
ObjLinearMVS2GVB rcpparma_LinearMVS2GVB(mat & Xr, colvec & yr, Options* opts) {

	int n = Xr.n_rows;
	int p = Xr.n_cols;
	int dispF = opts->dispF;
	int dispEvery = opts->display_gap;
	int max_iter = opts->max_iter;
	double tol = opts->epsStopLogLik;
	int iter = 1;
	int i = 0;
	double normyytilde = 0;

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

	colvec sigma2e_Sigma2beta = arma::zeros(1, 1);
	colvec sigma_beta = arma::zeros(1, 1);
	colvec sigma_error = arma::zeros(1, 1);

	sigma_beta(0, 0) = arma::var(y) / 2.0;
	sigma_error(0, 0) = arma::var(y) / 2.0;

	colvec mu = arma::zeros(p, 1);     // initial variational mean
	colvec s = arma::zeros(p, 1);
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
	colvec y_tilde = X*(alpha%mu);         //n-by-1 vector

	while (iter<max_iter && fabs(diff(0, 0)) > tol)
	{
		logalpha_prior = log(alpha_prior / (1 - alpha_prior));

		sigma2e_Sigma2beta = sigma_error / sigma_beta;

		xxsigma = diagxtx + sigma2e_Sigma2beta(0, 0);

		s = arma::repmat(sigma_error, p, 1) / xxsigma;

		for (i = 0; i<p; i++)
		{
			y_i = y_tilde - alpha(i, 0)*mu(i, 0)*X(arma::span(0, n - 1), i);

			mu(i, 0) = arma::as_scalar((xty(i, 0) - X(arma::span(0, n - 1), i).t()*y_i) / xxsigma(i, 0));

			u(i, 0) = logalpha_prior(0, 0) + 0.5*log(s(i, 0)) - 0.5*log(sigma_beta(0, 0)) + 0.5*mu(i, 0)*mu(i, 0) / s(i, 0);

			alpha(i, 0) = arma::as_scalar(1 / (1 + exp(-u(i, 0))));

			y_tilde = y_i + alpha(i, 0)*mu(i, 0)*X(arma::span(0, n - 1), i);

		}

		// M step
		term1 = alpha % (s + mu%mu);
		alphamu = alpha%mu;
		alphamu2 = alphamu%alphamu;
		normyytilde = arma::norm(y - y_tilde, 2);

		sigma_error = (pow(normyytilde, 2) + sum((term1 - alphamu2) % diagxtx)) / n;
		if (sigma_error(0, 0)<1e-6)
			sigma_error(0, 0) = 1e-6;

		sigma_beta(0, 0) = sum(term1) / sum(alpha);
		if (sigma_beta(0, 0)<1e-6)
			sigma_beta(0, 0) = 1e-6;

		alpha_prior = arma::sum(alpha) / p;

		for (i = 0; i<p; i++)
		{
			logpow(i, 0) = log(pow((1 - alpha(i, 0)) / (1 - alpha_prior(0, 0)), 1 - alpha(i, 0)));
		}
		//the limit using of power

		low_bound(0, iter-1) = -n*0.5*log(2 * arma::datum::pi*sigma_error(0, 0)) - normyytilde*normyytilde / sigma_error(0, 0) / 2
			- 1 / sigma_error(0, 0)*0.5*sum((term1 - alphamu2) % diagxtx)
			- sum(alpha%log(alpha / alpha_prior(0, 0)))
			- sum(logpow)
			+ sum(0.5*alpha % (1 + log(1 / sigma_beta(0, 0)*s) - 1 / sigma_beta(0, 0)*(s + mu%mu)));

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

	ObjLinearMVS2GVB obj;

	obj.vardist_gamma = alpha;
	obj.vardist_mu = mu;
	obj.vardist_sigma2beta = s;
	obj.sigma2beta = sigma_beta;
	obj.alpha = alpha_prior;
	obj.sigma2e = sigma_error;
	obj.Lq = lowerbound;

	return obj;
}
