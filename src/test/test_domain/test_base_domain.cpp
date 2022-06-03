#include "src/domain/base_domain.cpp"
#include "src/utils.cpp"


bool test_base_domain_contructor(){
    auto domain = BaseDomain(100, 100, 0, 10, 11);
    bool all_passed = true;
    if (!is_close(domain.get_t_start(), 0., 1e-12)){
        all_passed = false;
    }
    if (!is_close(domain.get_t_end(), 10., 1e-12))
    {
        all_passed = false;
    }
    if (domain.get_num_times() != 11)
    {
        all_passed = false;
    }
    if (!is_close(domain.get_dt(), 1., 1e-12))
    {
        all_passed = false;
    }
}
