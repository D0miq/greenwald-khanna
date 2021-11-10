#pragma once

#include "result.h"
#include <string>
#include <list>
#include <cmath>
#include <limits>
#include <vector>

template <typename T>
class Greenwald_khanna {
public:
	struct Tuple {
	public:
		Tuple(T v, long g, int delta) : v(v), g(g), delta(delta), min_bounds(v), max_bounds(v) {}
		T v;
		long g;
		int delta;
		T min_bounds;
		T max_bounds;
	};

	Greenwald_khanna(double epsilon) : m_epsilon(epsilon), m_one_devide_2e(1 / (2 * epsilon)), m_n(0) {}

	void insert(T value) {
		// Compress the data structure
		if (this->m_n > 0 && this->m_n % this->m_one_devide_2e == 0) {
			this->compress();
		}

		// Insert the value
		tuples_it it = this->find_insert_iterator(value);
		this->m_tuples.insert(it, Tuple(value, 1, this->compute_delta(it)));

		// Increment number of processed values
		this->m_n++;
	}

	Tuple query(long rank, long& low_bounds) const {
		long sum = 0;
		long range = this->m_epsilon * this->m_n;

		for (auto it = this->m_tuples.begin(); it != this->m_tuples.end(); it++) {
			sum += it->g;

			if (rank - sum <= range && sum + it->delta - rank <= range) { // TODO: Might be able to find easier condition
				low_bounds = sum;
				return (*it);
			}
		}
	}

	long n() const {
		return this->m_n;
	}

	double epsilon() const {
		return this->m_epsilon;
	}
private:
	typedef typename std::list<Tuple>::iterator tuples_it;

	void compress() {
		int i = this->m_tuples.size() - 2;
		auto it = this->m_tuples.end();
		std::advance(it, -2);

		long two_eps_n = 2 * this->m_epsilon * this->m_n;
		std::vector<long> bands = create_bands(two_eps_n);

		for (; it != this->m_tuples.begin(); it--) {
			auto next_it = std::next(it);

			if (bands[it->delta] <= bands[next_it->delta] && it->g + next_it->g + next_it->delta < two_eps_n) {
				next_it->g += it->g;

				if (next_it->min_bounds > it->min_bounds) {
					next_it->min_bounds = it->min_bounds;
				}

				if (next_it->max_bounds < it->max_bounds) {
					next_it->max_bounds = it->max_bounds;
				}

				this->m_tuples.erase(it);
			}

			i--;
		}
	}

	std::vector<long> create_bands(long two_eps_n) {
		long p = two_eps_n;
		long max_alpha = std::ceil(std::log2(two_eps_n));

		std::vector<long> bands(p + 1);
		bands[p] = 0;
		bands[0] = std::numeric_limits<long>::max();

		for (int alpha = 0; alpha < max_alpha; alpha++) {
			long two_to_alpha = 1 << alpha;
			long two_to_alpha_minus_one = 1 << (alpha - 1);

			long low_bound = p - two_to_alpha - (p % two_to_alpha);
			if (low_bound < 0) low_bound = 0;

			long high_bound = p - two_to_alpha_minus_one - (p % two_to_alpha_minus_one);

			for (long i = low_bound + 1; i <= high_bound; i++) {
				bands[i] = alpha;
			}
		}

		return bands;
	}

	tuples_it find_insert_iterator(T value) {

		// Find position where the value should be inserted

		//auto it = this->m_tuples.begin();
		//while (it != this->m_tuples.end() && it->v < value) {
		//	it++;
		//}

		//return it;

		// TODO - compare performance
		for (auto it = this->m_tuples.begin(); it != this->m_tuples.end(); it++) {
			if (it->v >= value) {
				return it;
			}
		}
	}

	int compute_delta(tuples_it it) {
		int delta = 0;
		if (it != this->m_tuples.begin() && it != this->m_tuples.end() // delta = 2*epsilon*n should be set only for items inside the list
			&& this->m_n > this->m_one_devide_2e) { // Modification from https://www.stevenengelhardt.com/2018/03/07/calculating-percentiles-on-streaming-data-part-2-notes-on-implementing-greenwald-khanna/ to maintain invariant 
			delta = floor(2 * this->m_epsilon * this->m_n) - 1;
		}

		return delta;
	}

	T min_value;
	T max_value;
	std::list<Tuple> m_tuples;
	int m_one_devide_2e;
	long m_n;
	double m_epsilon;
};