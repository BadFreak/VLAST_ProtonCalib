#include "ProtonCalib.h"
#include <iostream>

int main(int argc, char **argv) {
    std::cout << "======= Starting : VLAST ProtonCalib ======= " << std::endl;

    // ProtonCalib *calib = new ProtonCalib();
    // calib->ReadTreeAndSetBranches(argv[1]);
    // calib->ProtonMIPStat(argv[2]);

    ProtonCalib *calib = new ProtonCalib();
    std::cout << "======= InitHist ======= " << std::endl;
    calib->InitHist(argv[1]);

    std::cout << "======= ProtonSelection file1 ======= " << std::endl;
    calib->ReadTreeAndSetBranches(argv[1]);
    calib->ProtonSelection(1);

    std::cout << "======= ProtonSelection file2 ======= " << std::endl;
    calib->ReadTreeAndSetBranches(argv[2]);
    calib->ProtonSelection(0.708087);

    std::cout << "======= CalibFileWrite ======= " << std::endl;
    calib->CalibFileWrite(argv[3]);

    std::cout << "======= Ending : VLAST ProtonCalib ======= " << std::endl;
    return 0;
}