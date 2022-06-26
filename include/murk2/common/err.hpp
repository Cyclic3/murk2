#pragma once

#include <stdexcept>
#include <string>

#define MURK2_UNIMPLEMENTED {throw std::logic_error{"The function " + std::string{__func__} + " at " __FILE__ ":" + std::to_string(__LINE__) + " is not implemented"};}
