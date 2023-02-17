#include "parafoil.h"

using namespace std;

int main() {


    Parafoil prova("data.json");

    if (prova.boolControl == 1){
        prova.simulate_control();
    }
    else {
        prova.simulate();
    }

    return 0;
}

