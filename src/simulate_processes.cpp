// simulate_processes.cpp

// [[Rcpp::depends(RcppParallel, RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

// Wrapped Cauchy density
inline double dwrappedcauchy(double theta, double mu, double rho) {
  return (1 - rho * rho) / (2 * M_PI * (1 + rho * rho - 2 * rho * cos(theta - mu)));
}

// Compute next probabilities for a single simulation
arma::vec next_probabilities(double mu, double rho, const arma::vec& a, const arma::vec& b, int m, double small_arc_threshold) {
  arma::vec probs(m, arma::fill::zeros);
  for (int i = 0; i < m; i++) {
    double arc_length = b[i] - a[i];
    if (arc_length >= small_arc_threshold) {
      int num_steps = 100;
      double step = arc_length / num_steps;
      double integral = 0.0;
      for (int j = 0; j <= num_steps; ++j) {
        double theta = a[i] + j * step;
        double weight = (j == 0 || j == num_steps) ? 0.5 : 1.0;
        integral += weight * dwrappedcauchy(theta, mu, rho);
      }
      probs[i] = integral * step;
    }
  }
  return probs / arma::accu(probs);
}

// Parallel worker
struct SimulationWorker : public Worker {
  const double rho;
  const int m;
  const int n;
  const double small_arc_threshold;
  const std::string mu_behavior;
  arma::mat& output;
  
  SimulationWorker(double rho, int m, int n, double small_arc_threshold, std::string mu_behavior, arma::mat& output)
    : rho(rho), m(m), n(n), small_arc_threshold(small_arc_threshold), mu_behavior(mu_behavior), output(output) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t process = begin; process < end; ++process) {
      arma::vec P_prev(m, arma::fill::ones);
      P_prev /= m;
      arma::vec a = arma::linspace(0, 2 * M_PI - 2 * M_PI / m, m);
      arma::vec b = a + 2 * M_PI / m;
      double mu = R::runif(0, 2 * M_PI);
      
      for (int iter = process + 1; iter <= n; ++iter) {
        if (mu_behavior == "random_each_iter") {
          mu = R::runif(0, 2 * M_PI);
        } else if (mu_behavior == "fixed_start_rotate_pi") {
          mu += M_PI;
        } else if (mu_behavior == "fixed_start_rotate_half_pi") {
          mu += M_PI / 2;
        }
        
        mu = fmod(mu, 2 * M_PI);
        if (mu < 0) mu += 2 * M_PI;
        
        arma::vec P_next = next_probabilities(mu, rho, a, b, m, small_arc_threshold);
        P_prev = P_next;
        a = arma::cumsum(arma::join_cols(arma::vec({0}), P_prev)) * 2 * M_PI;
        b = a.tail(m);
        a = a.head(m);
      }
      output.col(process) = P_prev;
    }
  }
};

// [[Rcpp::export]]
NumericMatrix simulate_processes_cpp(int m, double rho, int n, double small_arc_threshold, std::string mu_behavior) {
  arma::mat result(m, n);
  SimulationWorker worker(rho, m, n, small_arc_threshold, mu_behavior, result);
  parallelFor(0, n, worker);
  return wrap(result);
}

// [[Rcpp::export]]
NumericMatrix simulate_process_to_fixed_iteration(double rho, int m, int n_start, int n, double small_arc_threshold, std::string mu_behavior) {
  // n_start is ignored for now — simulates from n_start to n always
  return simulate_processes_cpp(m, rho, n, small_arc_threshold, mu_behavior);
}

