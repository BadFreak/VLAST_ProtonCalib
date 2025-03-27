#include "ProtonCalib.h"
#include <iostream>

int main(int argc, char** argv) {
    std::cout << "======= Starting : VLAST ProtonCalib ======= " << std::endl;
	
	ProtonCalib *calib = new ProtonCalib();
    calib->ReadTreeAndSetBranches(argv[1]);
	calib->ProtonMIPStat(argv[2]);
	
	std::cout << "======= Ending : VLAST ProtonCalib ======= " << std::endl;
	return 0;
}