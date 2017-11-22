#include "input_args.h"
int main(int argc, char **argv) {
    InputArgs inputArgs = InputArgs(); 
    if (inputArgs.read(argc, argv)) {
        cout << inputArgs.inputMode << endl;
    } else {
        return 0;
    }
}
