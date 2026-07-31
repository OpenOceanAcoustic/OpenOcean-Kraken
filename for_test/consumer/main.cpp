#include <OpenOceanKrakencKernelInterface.h>

#include <iostream>

int main()
{
    OpenOceanKrakenc::KernelInterface interface;
    if (interface.isClosed() || interface.getHardwareThreads() < 1)
    {
        return 1;
    }
    interface.close();
    if (!interface.isClosed())
    {
        return 2;
    }
    std::cout << "INSTALLED_CONSUMER_OK\n";
    return 0;
}
