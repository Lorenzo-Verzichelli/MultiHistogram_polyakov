#pragma once
#include <iostream>
#include <cmath>

class weighted_mean
{
private:
	double log_sumw_, log_sqrw_, log_mean_, log_sqrm_;
	unsigned int count_;

public:
	weighted_mean() : log_sumw_(1), log_sqrw_(1), log_mean_(1), log_sqrm_(1), count_(0) {}

	inline void reset(void);

	inline weighted_mean& add(double element, double weight = 1);
	inline weighted_mean& add_log_weig(double element, double log_weig);
	inline weighted_mean& add_log_elem(double log_elem, double weight = 1);
	inline weighted_mean& add_log_both(double log_elem, double log_weig);

	inline weighted_mean& remove(double element, double weight = 1);
	inline weighted_mean& remove_log_weig(double element, double log_weig);
	inline weighted_mean& remove_log_elem(double log_elem, double weight = 1);
	inline weighted_mean& remove_log_both(double log_elem, double log_weig);

	inline unsigned int count(void) const;

	inline double mean(void) const;
	inline double mean_sqr(void) const;
	inline double variance(bool unbiased = true) const;

	inline double sum_weight(void) const;
	inline double sum_weig_sqr(void) const;
	inline double sum_element(void) const;
	inline double sum_elem_sqr(void) const;
};

inline void weighted_mean::reset(void)
{
	log_sumw_ = 1;
	log_sqrw_ = 1;
	log_mean_ = 1;
	log_sqrm_ = 1;
	count_ = 0;
}

inline weighted_mean& weighted_mean::add(double element, double weight)
{
	return add_log_both(log(element), log(weight));
}

inline weighted_mean& weighted_mean::add_log_weig(double element, double log_weig)
{
	return add_log_both(log(element), log_weig);
}

inline weighted_mean& weighted_mean::add_log_elem(double log_elem, double weight)
{
	return add_log_both(log_elem, log(weight));
}

inline weighted_mean& weighted_mean::add_log_both(double log_elem, double log_weig)
{
	if (count_ == 0) {
		count_ = 1;
		log_sumw_ = log_weig;
		log_sqrw_ = 2 * log_weig;
		log_mean_ = log_elem + log_weig;
		log_sqrm_ = 2 * log_elem + log_weig;
		return *this;
	}

	double diff = log_weig - log_sumw_;
	if (diff < 0) log_sumw_ += log1p(exp(diff));
	else log_sumw_ = log_weig + log1p(exp(-diff));

	double log_ele2 = 2 * log_elem + log_weig;
	diff = log_ele2 - log_sqrm_;
	if (diff <= 0) log_sqrm_ += log1p(exp(diff));
	else log_sqrm_ = log_ele2 + log1p(exp(-diff));
	
	log_elem += log_weig;
	diff = log_elem - log_mean_;
	if (diff <= 0) log_mean_ += log1p(exp(diff));
	else log_mean_ = log_elem + log1p(exp(-diff));

	log_weig *= 2;
	diff = log_weig - log_sqrw_;
	if (diff <= 0) log_sqrw_ += log1p(exp(diff));
	else log_sqrw_ = log_weig + log1p(exp(-diff));

	count_++;
	return *this;
}

inline weighted_mean& weighted_mean::remove(double element, double weight)
{
	return remove_log_both(log(element), log(weight));
}

inline weighted_mean& weighted_mean::remove_log_weig(double element, double log_weig)
{
	return remove_log_both(log(element), log_weig);
}

inline weighted_mean& weighted_mean::remove_log_elem(double log_elem, double weight)
{
	return remove_log_both(log_elem, log(weight));
}

inline weighted_mean& weighted_mean::remove_log_both(double log_elem, double log_weig)
{
	double diff = log_weig - log_sumw_;
	if (diff <= 0) log_sumw_ += log1p(-exp(diff));
	else log_sumw_ = log_weig + log1p(-exp(-diff));

	double log_ele2 = 2 * log_elem + log_weig;
	diff = log_ele2 - log_sqrm_;
	if (diff <= 0) log_sqrm_ += log1p(-exp(diff));
	else log_sqrm_ = log_ele2 + log1p(-exp(-diff));

	log_elem += log_weig;
	diff = log_elem - log_mean_;
	if (diff <= 0) log_mean_ += log1p(-exp(diff));
	else log_mean_ = log_elem + log1p(-exp(-diff));

	log_weig *= 2;
	diff = log_weig - log_sqrw_;
	if (diff <= 0) log_sqrw_ += log1p(-exp(diff));
	else log_sqrw_ = log_weig + log1p(-exp(-diff));

	count_--;
	return *this;
}

inline unsigned int weighted_mean::count(void) const
{
	return count_;
}

inline double weighted_mean::mean(void) const
{
	return exp(log_mean_ - log_sumw_);
}

inline double weighted_mean::mean_sqr(void) const
{
	return exp(log_sqrm_ - log_sumw_);
}



inline double weighted_mean::variance(bool unbiased) const
{
	double diff = 2 * log_mean_ - log_sumw_ - log_sqrm_; //diff should be negative
	double log_var = log_sqrm_ + log1p(-exp(diff));

	if (!unbiased) {
		log_var -= log_sumw_;
		return exp(log_var);
	}

	diff = log_sqrw_ - 2 * log_sumw_; //negative
	
	log_var -= log_sumw_ + log1p(-exp(diff));

	return exp(log_var);
}

inline double weighted_mean::sum_weight(void) const
{
	return exp(log_sumw_);
}

inline double weighted_mean::sum_weig_sqr(void) const
{
	return exp(log_sqrw_);
}

inline double weighted_mean::sum_element(void) const
{
	return exp(log_mean_);
}

inline double weighted_mean::sum_elem_sqr(void) const
{
	return exp(log_sqrm_);
}

