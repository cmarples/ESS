set(HEADER_TRNG
        bernoulli_dist.hpp
        beta_dist.hpp
        binomial_dist.hpp
        cauchy_dist.hpp
        chi_square_dist.hpp
        constants.hpp
        correlated_normal_dist.hpp
        cuda.hpp
        discrete_dist.hpp
        exponential_dist.hpp
        extreme_value_dist.hpp
        fast_discrete_dist.hpp
        gamma_dist.hpp
        generate_canonical.hpp
        geometric_dist.hpp
        hypergeometric_dist.hpp
        int_math.hpp
        int_types.hpp
        lagfib2plus.hpp
        lagfib2xor.hpp
        lagfib4plus.hpp
        lagfib4xor.hpp
        lcg64.hpp
        lcg64_count_shift.hpp
        lcg64_shift.hpp
        limits.hpp
        linear_algebra.hpp
        logistic_dist.hpp
        lognormal_dist.hpp
        math.hpp
        maxwell_dist.hpp
        minstd.hpp
        mrg2.hpp
        mrg3.hpp
        mrg3s.hpp
        mrg4.hpp
        mrg5.hpp
        mrg5s.hpp
        mrg_parameter.hpp
        mrg_status.hpp
        mt19937_64.hpp
        mt19937.hpp
        negative_binomial_dist.hpp
        normal_dist.hpp
        pareto_dist.hpp
        poisson_dist.hpp
        powerlaw_dist.hpp
        rayleigh_dist.hpp
        snedecor_f_dist.hpp
        special_functions.hpp
        student_t_dist.hpp
        tent_dist.hpp
        truncated_normal_dist.hpp
        twosided_exponential_dist.hpp
        uniform01_dist.hpp
        uniform_dist.hpp
        uniform_int_dist.hpp
        uniformxx.hpp
        utility.hpp
        weibull_dist.hpp
        xoshiro256plus.hpp
        yarn2.hpp
        yarn3.hpp
        yarn3s.hpp
        yarn4.hpp
        yarn5.hpp
        yarn5s.hpp
        zero_truncated_poisson_dist.hpp
        )

set(SOURCE_TRNG
        lcg64.cc
        lcg64_shift.cc
        lcg64_count_shift.cc
        minstd.cc
        mrg2.cc
        mrg3.cc
        mrg3s.cc
        mrg4.cc
        mrg5.cc
        mrg5s.cc
        mt19937_64.cc
        mt19937.cc
        xoshiro256plus.cc
        yarn2.cc
        yarn3.cc
        yarn3s.cc
        yarn4.cc
        yarn5.cc
        yarn5s.cc
        )

add_library(trng4 STATIC ${HEADER_TRNG} ${SOURCE_TRNG})
target_include_directories(trng4 PUBLIC ${PROJECT_SOURCE_DIR}/lib)