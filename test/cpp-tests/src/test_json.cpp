#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"

using namespace ZonoOpt;

int main()
{
    const ZonoPtr Z = make_regular_zono_2D(3., 12);
    const std::string filename = "test_zono.json";
    to_json(*Z, filename);

    std::cout << *Z << std::endl;

    const ZonoPtr Z_read = from_json(filename);

    std::cout << *Z_read << std::endl;

    return 0;
}