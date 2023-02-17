#include "parafoil.h"

using namespace std;

int main() {


    Parafoil prova("data.json");

    prova.simulate_control();

    return 0;
}