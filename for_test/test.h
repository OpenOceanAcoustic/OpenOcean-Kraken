#ifndef OOK_FOR_TEST_TEST_H
#define OOK_FOR_TEST_TEST_H

#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace ook_test
{
    class TestRunner
    {
    public:
        void require(bool condition, const std::string &message)
        {
            ++checks_;
            if (!condition)
            {
                throw std::runtime_error(message);
            }
        }

        template <typename Fn>
        void requireThrows(Fn &&fn, const std::string &message)
        {
            ++checks_;
            try
            {
                fn();
            }
            catch (const std::exception &)
            {
                return;
            }
            throw std::runtime_error(message);
        }

        void requireNear(double actual, double expected, double tolerance, const std::string &message)
        {
            ++checks_;
            if (std::abs(actual - expected) > tolerance)
            {
                std::ostringstream oss;
                oss << message << " actual=" << actual << " expected=" << expected;
                throw std::runtime_error(oss.str());
            }
        }

        int checks() const
        {
            return checks_;
        }

    private:
        int checks_ = 0;
    };
}

#endif
