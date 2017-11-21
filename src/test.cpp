#include "input_args.h"
int main(int argc, char **argv) {
    InputArgs inputArgs = InputArgs(); 
    inputArgs.read(argc, argv);
    cout << inputArgs.inputMode << endl;
    return 0;
}
