#include "src/domain/base_domain.cpp"
#include "src/utils.cpp"


bool test_grid_point(){
    auto grid_point = GridPoint(0.1, 0.1, std::complex<double>{10., 1.});
    bool all_passed = true;
    if (!is_close(grid_point.x, 0.1, 1e-12))
    {
        all_passed = false;
    }
    if (!is_close(grid_point.y, 0.1, 1e-12))
    {
        all_passed = false;
    }
    if (!is_close(grid_point.wave_function.real(), 10., 1e-12))
    {
        all_passed = false;
    }
    if (!is_close(grid_point.wave_function.imag(), 1., 1e-12))
    {
        all_passed = false;
    }
    return all_passed;
}

// bool test_base_domain_contructor(){
//     auto domain = BaseDomain(100, 100, 0, 10, 11);
//     bool all_passed = true;
//     if (!is_close(domain.get_t_start(), 0., 1e-12)){
//         all_passed = false;
//     }
//     if (!is_close(domain.get_t_end(), 10., 1e-12))
//     {
//         all_passed = false;
//     }
//     if (domain.get_num_times() != 11)
//     {
//         all_passed = false;
//     }
//     if (!is_close(domain.get_dt(), 1., 1e-12))
//     {
//         all_passed = false;
//     }
//     return all_passed;
// }
